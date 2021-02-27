# Calculate plant gps coordinates -- generally handy functions that are used within site-specific calculations (see estim_plantGPScoords_bySite.R)


##### geometry ########

# calculate missing diagnal corners of a rectangle given 2 points
find_missingDiag.corners <- function(x1, y1, x3, y3, flipdiag){

  # calculate the diagonal's midpoint
  midpoint.x <- (x1 + x3) / 2
  midpoint.y <- (y1 + y3) / 2

  #plot it
  xs <- c(x1, x3, midpoint.x)
  ys <- c(y1, y3, midpoint.y)
  point.name <- c("1","3","midpoint")
  df <- data.frame(point.name, xs, ys)
  ggplot(df, aes(x = xs, y = ys, label = point.name)) +
    geom_text() + coord_fixed()

  # calculate the radius of the circle
  radius <- sqrt(((x1-x3)/2)^2 + ((y1-y3)/2)^2)

  # find the two points that fall on the circle, aka the 2 other rectangle points
  if(flipdiag == T){
    x2.points <- c(midpoint.x - radius, midpoint.x + radius)
    #radius^2 = (x2-midpoint.x)^2 + (y2-midpoint.y)^2
    y2 <- sqrt(radius^2 - (x2.points-midpoint.x)^2) + midpoint.y
    new.point.rows <- data.frame(point.name = "new.point",
                                 xs = x2.points,
                                 ys = y2)
    df.out <- rbind(df, new.point.rows)
  }else{
    y2.points <- c(midpoint.y - radius, midpoint.y + radius)
    #radius^2 = (x2-midpoint.x)^2 + (y2-midpoint.y)^2
    x2 <- sqrt(radius^2 - (y2.points-midpoint.y)^2) + midpoint.x
    new.point.rows <- data.frame(point.name = "new.point",
                                 xs = x2,
                                 ys = y2.points)
    df.out <- rbind(df, new.point.rows)
  }

  #plot it again
  ggplot(df.out, aes(x = xs, y = ys, label = point.name)) +
    geom_text() + coord_fixed()
  df.out %>%
    filter(point.name == "new.point") %>%
    dplyr::rename('lat'='ys',
                  'lon'='xs') -> df.result

  return(df.result)

}

# calculate the slope and intercept of a line given 2 points
calc_line.eq <- function(x1,y1,x2,y2){
  slope <- (y2 - y1) / (x2 - x1)
  intercept <- y2 - (slope*x2)
  result <- data.frame(slope = slope, intercept = intercept)
  names(result) <- c("slope","intercept")
  
  return(result)
}

# calculate the intercept of a line given the slope and 1 point
calc_intercept_fromSlopePoint <- function(slope, x, y){
  intercept <- y - (slope * x)
  return(intercept)
}

# calculate 1 point given 2 intersecting lines
calc_point_fromLines <- function(slope1, intercept1, slope2, intercept2){
  new.x = (intercept2 - intercept1) / (slope1 - slope2)
  new.y = (slope1 * new.x) + intercept1
  result <- data.frame(new.x = unlist(new.x), new.y = unlist(new.y), row.names = NULL)
  return(result)
}

# calculate 1 point given an origin point, a line, and the distance along that line 
useLineDist_findNewXY <- function(dist, slope, intercept, x0, y0){
  
  x1 <- x0 + dist/sqrt(1+slope^2)
  y1 <- (slope* x1) + intercept
  names(x1) <- names(y1) <- NULL
  xy1 <- data.frame(x1=x1, y1=y1)
  # do I need to correct the new longitude based on latitude?
  #new_longitude = longitude + (dx / r_earth) * (180 / pi) / cos(latitude * pi/180)
  
  return(xy1)
}

# convert meters to degree (for gps coords)
load_m.per.deg <- function(){
  # number of meters per degree
  r_earth <- 6378000 #meters
  m_per_deg <- (2*pi/360) * r_earth
  return(m_per_deg)
}

# calculate gps coordinates of samples within a plot given the plot's (a) origin lat and lon, (b) x and y lines in terms of lat and lon, and (c) xy coordinates of each sample within the plot. 
# reverse arguments allow for changing the positive and negative direction of the plot's xy grid on top of the lat lon surface. 
calc_plotplant.coords <- function(origin, xax.line, yax.line, plant.samps, df, 
                                  reverse.xdist,
                                  reverse.yslope,
                                  reverse.ydist){
  
  #origin = origin.b1
  #xax.line = xax.b1 
  #yax.line = yax.b1 
  #plant.samps = samps.b1 
  #df = gps.df
  #reverse.xdist <- F # reverse.xdist = T if distance from origin along xax.line is in the negative direction
  #reverse.yslope <- T # reverse.yslope = T if the yax.line slope is negative
  #reverse.ydist <- T # reverse.y = T if distance from origin along yax.line is in the negative direction
  ####reverse.x <- F # reverse.x = T if the yax.line slope is positive???? get rid of this one.
  
  # -------------------------------------#
  # load values to convert meters to gps units
  m_per_deg <- load_m.per.deg()
  
  # -------------------------------------#
  # 1. Find the GPS coordinates that correspond to each plant sample along the xax.line from the plot origin
  if(reverse.xdist == T){rev.xax <- -1}else{rev.xax <- 1}
  xax.coords<- useLineDist_findNewXY(dist = plant.samps$Xcoord.m/m_per_deg * rev.xax, 
                                     slope = xax.line$slope, 
                                     intercept = xax.line$intercept, 
                                     x0 = origin$lon, 
                                     y0 = origin$lat)
  xax.coords %>%
    dplyr::rename('lat'='y1', 'lon'='x1') -> df.xax.coords
  df.xax.coords$point.name <- paste0("x", plant.samps$Samp)
  df.xax.coords$type <- "other"
  p1 <- ggplot(df, aes(x = lon, y = lat)) +
    geom_text(aes(label = point.name)) +
    geom_abline(slope = xax.line$slope, intercept = xax.line$intercept) +
    geom_abline(slope = yax.line$slope, intercept = yax.line$intercept) +
    geom_text(data = df.xax.coords, aes(label = point.name)) +
    coord_fixed()
  p1
  
  # -------------------------------------#
  # 2. Use the new gps coords along xax and add y distance in the yax direction
  # first, use the new coords to find the intercepts for the lines parallel to yax
  if(reverse.yslope == T){rev.pbs <- -1}else{rev.pbs <- 1}
  parallel.bs <- xax.coords$y1 + rev.pbs*(yax.line$slope * xax.coords$x1)
  parallel.bs
  p2 <- p1 +
    geom_abline(slope = yax.line$slope, intercept = parallel.bs, linetype = 2)
  p2
  
  # second, find new point along new line using y dist
  if(reverse.ydist == T){rev.y <- -1}else{rev.y <- 1}
  rev.x <- 1 ### get rid of this?????
  plant.coords<- useLineDist_findNewXY(dist = (plant.samps$Ycoord.m/m_per_deg * rev.y), 
                                         slope = yax.line$slope * rev.x, 
                                         intercept = parallel.bs, 
                                         x0 = xax.coords$x1, 
                                         y0 = xax.coords$y1)
  plant.coords %>%
    dplyr::rename('lat'='y1', 'lon'='x1') -> df.plant.coords
  df.plant.coords$point.name <- paste0("p", plant.samps$Samp)
  df.plant.coords$type <- "plant"
  p3 <- p2 + geom_text(data = df.plant.coords, aes(y = lat, x = lon, label = point.name))
  p3

  # -------------------------------------#
  # make a df that holds the plot corner coordinates, the plot x-axis coordinates, and the plant coordinates
  #rbind(df.plant.coords, df.xax.coords) -> plotplant.coords
  
  # save everything in a list
  plotplant.data <- list(plant.coords = df.plant.coords, 
                         p = p3)
  return(plotplant.data)
}

calc_plot_xylines <- function(df, point.name.chr, xax.slope, yax.slope){
  
  xax.b1.intercept <- calc_intercept_fromSlopePoint(slope = xax.slope, 
                                                    x = df[df$point.name == point.name.chr,"lon"], 
                                                    y = df[df$point.name == point.name.chr,"lat"])
  yax.b1.intercept <- calc_intercept_fromSlopePoint(slope = yax.slope, 
                                                    x = df[df$point.name == point.name.chr,"lon"], 
                                                    y = df[df$point.name == point.name.chr,"lat"])
  
  xax.line.b1 <- data.frame(slope = xax.slope,
                            intercept = as.numeric(xax.b1.intercept))
  yax.line.b1 <- data.frame(slope = yax.slope,
                            intercept = as.numeric(yax.b1.intercept))
  
  result <- list(xax = xax.line.b1, yax = yax.line.b1)
  
  return(result)
}


##### plotting ########

# plot the samples in lat/lon space to veryify the estimated locations 
check_plotplant.coords <- function(plotplant.data, plot.title){
  
  # rename the df
  df <- plotplant.data$plant.coords
  parallel.bs <- plotplant.data$parallel.bs
  yax.line <- plotplant.data$yax.line
  xax.line <- plotplant.data$xax.line
  
  # find the corner coordinates to draw the plot edges
  a <- df[df$point.name == "A",c("lon","lat")]
  b <- df[df$point.name == "B",c("lon","lat")]
  c <- df[df$point.name == "C",c("lon","lat")]
  d <- df[df$point.name == "D",c("lon","lat")]
  
  # draw the plot
  #plot.title <-"thing"
  p <- ggplot(df, aes(x = lon, y = lat, label = point.name, color = type)) +
    geom_label() +
    coord_fixed() +
    geom_segment(x = a$lon, y = a$lat, xend = b$lon, yend = b$lat) +
    geom_segment(x = b$lon, y = b$lat, xend = c$lon, yend = c$lat) +
    geom_segment(x = c$lon, y = c$lat, xend = d$lon, yend = d$lat) +
    geom_segment(x = d$lon, y = d$lat, xend = a$lon, yend = a$lat) +
    ggtitle(plot.title)
  p
  # add plant samples
  if(length(parallel.bs) == 4){
    p +
      geom_abline(slope = yax.line$slope, intercept = parallel.bs[1], color = "gray", linetype = 2) +
      geom_abline(slope = yax.line$slope, intercept = parallel.bs[2], color = "gray", linetype = 2) +
      geom_abline(slope = yax.line$slope, intercept = parallel.bs[3], color = "gray", linetype = 2) +
      geom_abline(slope = yax.line$slope, intercept = parallel.bs[4], color = "gray", linetype = 2) -> p
  }
  if(length(parallel.bs) == 8){
    p +
      geom_abline(slope = yax.line$slope, intercept = parallel.bs[1], color = "gray", linetype = 2) +
      geom_abline(slope = yax.line$slope, intercept = parallel.bs[2], color = "gray", linetype = 2) +
      geom_abline(slope = yax.line$slope, intercept = parallel.bs[3], color = "gray", linetype = 2) +
      geom_abline(slope = yax.line$slope, intercept = parallel.bs[4], color = "gray", linetype = 2) +
      geom_abline(slope = yax.line$slope, intercept = parallel.bs[5], color = "gray", linetype = 2) +
      geom_abline(slope = yax.line$slope, intercept = parallel.bs[6], color = "gray", linetype = 2) +
      geom_abline(slope = yax.line$slope, intercept = parallel.bs[7], color = "gray", linetype = 2) +
      geom_abline(slope = yax.line$slope, intercept = parallel.bs[8], color = "gray", linetype = 2) -> p
  }
  
  return(p)
}



