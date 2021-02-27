# Calculate plant gps coordinates -- site-specific functions

#site.path <- 'data/sites/DOE-NC-FIELD_SiteData_2019_06_19.xlsx'
#gps.coords <- load_gpscoords(site.path)
#plants <- load_plants(site.path)


# helpers
get.curr.gps <- function(gps.coords, plants, site.name){
  
  gps.coords %>%
    filter(Site == site.name) -> curr.site
  plants %>%
    filter(Site == site.name) %>%
    select(Samp, Xcoord.m, Ycoord.m) -> curr.samps
  
  result <- list(curr.site = curr.site, curr.samps = curr.samps)
  return(result)
}

make.gps.long <- function(curr.site){
  
  gps.cols <- colnames(curr.site)[grepl("gps",colnames(curr.site))]
  curr.site %>%
    select(multiplot.notes, gps.cols) %>%
    gather(key = "key", value = "value", -multiplot.notes) %>%
    mutate(lat.lon = ifelse(grepl("N", key), "lat","lon")) %>%
    separate(key, into = c("thing1","point.name","thing2"), sep=c(3,4)) %>%
    select(multiplot.notes, point.name, lat.lon, value) %>%
    spread(key = lat.lon, value = value) -> curr.site.long
  
  return(curr.site.long)
}

# combine sites
compile_sitesamp_gps <- function(gps.coords, plants){
  
  #load data from each site
  mhc.one <- mhc.one.gpscoords(gps.coords, plants)
  sfa.one <- sfa.one.gpscoords(gps.coords, plants)
  cgf.mon <- cgf.mon.gpscoords(gps.coords, plants)
  cgf.mxg <- cgf.mxg.gpscoords(gps.coords, plants)
  oto.mon <- oto.mon.gpscoords(gps.coords, plants)
  oto.mxt <- oto.mxt.gpscoords(gps.coords, plants)
  ccr.one <- ccr.one.gpscoords(gps.coords, plants)
  cre.mxg <- cre.mxg.gpscoords(gps.coords, plants)
  cre.mxt <- cre.mxt.gpscoords(gps.coords, plants)
  ucp.mxg <- ucp.mxg.gpscoords(gps.coords, plants)
  lco.mxt <- lco.mxt.gpscoords(gps.coords, plants)
  wbi.nrt <- wbi.nrt.gpscoords(gps.coords, plants)
  brf.one <- brf.one.gpscoords(gps.coords, plants)
  lwr.bho <- lwr.bho.gpscoords(gps.coords, plants)
  
  site.names <- unique(plants$Site)[-15] # exclude MAF-ONE-PRO (site w/only 3 samples)
  coord.list <- list(mhc.one$coords,
       sfa.one$coords,
       cgf.mon$coords,
       cgf.mxg$coords,
       oto.mon$coords, 
       oto.mxt$coords, 
       ccr.one$coords,
       cre.mxg$coords,
       cre.mxt$coords,
       ucp.mxg$coords,
       lco.mxt$coords,
       wbi.nrt$coords,
       brf.one$coords,
       lwr.bho$coords)
  names(coord.list) <- site.names
  coord.df <- list_to_df(coord.list)
  coord.df <- data.frame(coord.df, row.names = NULL)
  coord.df %>%
    dplyr::rename('Site'='source') %>%
    arrange(type) -> coord.df
  
  return(coord.df)

}


########## MHC-ONE-NCD #########################################################
mhc.one.gpscoords <- function(gps.coords, plants){
  
  #---------------------------------------------------------#
  # 1. Identify relevant data
  curr.gps <- get.curr.gps(gps.coords, plants, site.name = "MHC-ONE-NCD")
  gps.df <- make.gps.long(curr.gps$curr.site)
  # these point names correspond to the site corners, so rename these
  gps.df %>%
    select(point.name, lat, lon) %>%
    filter(!is.na(lat)) %>%
    mutate(type = "site.corner") %>%
    mutate(plot = NA) -> gps.df
  
  #---------------------------------------------------------#
  # 2. Identify missing site corners
  ggplot(gps.df, aes(y = lat, x = lon, label = point.name)) +
    geom_text() +
    coord_fixed()
  missing.corners <- find_missingDiag.corners(x1 = as.numeric(gps.df[gps.df$point.name == "D","lon"]), 
                                              y1 = as.numeric(gps.df[gps.df$point.name == "D","lat"]), 
                                              x3 = as.numeric(gps.df[gps.df$point.name == "B","lon"]), 
                                              y3 = as.numeric(gps.df[gps.df$point.name == "B","lat"]),
                                              flipdiag = T)
  tmp <- data.frame(point.name = c("A","C"),
             lat = missing.corners$lat,
             lon = missing.corners$lon,
             type = "site.corner",
             plot = NA)
  rbind(gps.df, tmp) -> gps.df
  ggplot(gps.df, aes(y = lat, x = lon, label = point.name)) +
    geom_text() +
    coord_fixed()
  
  #---------------------------------------------------------#
  # 3a. Find plant sample coords for plot 1 == southwest
  origin.b1 <- gps.df[gps.df$point.name == "D",]
  xax.b1 <- calc_line.eq(x1 = gps.df[gps.df$point.name == "D","lon"], 
                         y1 = gps.df[gps.df$point.name == "D","lat"],
                         x2 = gps.df[gps.df$point.name == "C","lon"],
                         y2 = gps.df[gps.df$point.name == "C","lat"])
  yax.b1 <- calc_line.eq(x1 = gps.df[gps.df$point.name == "D","lon"], 
                         y1 = gps.df[gps.df$point.name == "D","lat"],
                         x2 = gps.df[gps.df$point.name == "A","lon"],
                         y2 = gps.df[gps.df$point.name == "A","lat"])
  curr.gps$curr.samps %>%
    filter(Samp %in% 1:4) -> samps.b1
  coords.b1 <- calc_plotplant.coords(origin = origin.b1, 
                                     xax.line = xax.b1, 
                                     yax.line = yax.b1, 
                                     plant.samps = samps.b1, 
                                     df = gps.df,
                                     reverse.xdist = F,
                                     reverse.yslope = T,
                                     reverse.ydist = T)
  #coords.b1$p
  #coords.b1$plant.coords
  coords.b1$plant.coords %>%
    mutate(plot = "b1") -> coords.b1$plant.coords
  # add plot.corner coordinates
  corners.b1 <- data.frame(lon = origin.b1$lon, lat = origin.b1$lat,
                         point.name = "D.b1",
                         type = "plot.corner",
                         plot = "b1")
  coords.b1$plant.coords %>%
    rbind(corners.b1) -> data.b1
  
  #---------------------------------------------------------#
  # 3b. Find plant sample coords for plot 2 == northeast
  
  # i. find coordA.b2
  xax2.b2 <- calc_line.eq(x1 = gps.df[gps.df$point.name == "B","lon"],
                         y1 = gps.df[gps.df$point.name == "B","lat"],
                         x2 = gps.df[gps.df$point.name == "A","lon"],
                         y2 = gps.df[gps.df$point.name == "A","lat"])
  # move 11.2m from B along xax.b1
  # load values to convert meters to gps units
  m_per_deg <- load_m.per.deg()
  coordA.b2 <- useLineDist_findNewXY(dist=11.2/m_per_deg * -1,
                                   slope = xax2.b2$slope, intercept = xax2.b2$intercept,
                                   x0 = gps.df[gps.df$point.name == "B","lon"],
                                   y0 = gps.df[gps.df$point.name == "B","lat"])
  p.tmp <- coords.b1$p +
    annotate("text", 
             x = coordA.b2$x1, 
             y = coordA.b2$y1, 
             color = "orange", label = "coordA.b2")
  p.tmp
  
  # ii. find coordD.b2
  # move 6.6m from B along the y-axis 
  yax2.b2 <- calc_line.eq(x1 = gps.df[gps.df$point.name == "B","lon"],
                           y1 = gps.df[gps.df$point.name == "B","lat"],
                           x2 = gps.df[gps.df$point.name == "C","lon"],
                           y2 = gps.df[gps.df$point.name == "C","lat"])
  new.intercept <- coordA.b2$y1 - (yax2.b2$slope * coordA.b2$x1)
  coordD.b2 <- useLineDist_findNewXY(dist=6.6/m_per_deg,
                                  slope = yax2.b2$slope, intercept = new.intercept,
                                  x0 = coordA.b2$x1,
                                  y0 = coordA.b2$y1)
  p.tmp <- p.tmp +
    annotate("text", 
             x = coordD.b2$x1, 
             y = coordD.b2$y1, 
             color = "orange", label = "coordD.b2")
  p.tmp
  
  # iii. identify xax.b2 and yax.b2
  new.intercept.x <- coordD.b2$y1 - (xax2.b2$slope * coordD.b2$x1)
  p.tmp <- p.tmp +
    geom_abline(slope = xax2.b2$slope, intercept = as.numeric(new.intercept.x), color = "orange") +
    geom_abline(slope = yax2.b2$slope, intercept = as.numeric(new.intercept), color = "orange")
  p.tmp
  xax.b2 <- data.frame(slope = xax2.b2$slope, intercept = as.numeric(new.intercept.x))
  yax.b2 <- data.frame(slope = yax2.b2$slope, intercept = as.numeric(new.intercept))
  
  # iv. identify the x,y coords for each plant sampled in the plot
  curr.gps$curr.samps %>%
    filter(Samp %in% 5:8) -> samps.b2
  coords.b2 <- calc_plotplant.coords(origin = data.frame(lat = coordD.b2$y1, lon = coordD.b2$x1), 
                                     xax.line = xax.b2, 
                                     yax.line = yax.b2, 
                                     plant.samps = samps.b2, 
                                     df = gps.df,
                                     reverse.xdist = F,
                                     reverse.yslope = T,
                                     reverse.ydist = T)
  #coords.b2$p
  #coords.b2$plant.coords
  coords.b2$plant.coords %>%
    mutate(plot = "b2") -> coords.b2$plant.coords
  # add plot.corner coordinates
  corners.b2 <- data.frame(lon = c(coordA.b2$x1, 
                                   as.numeric(gps.df[gps.df$point.name == "B","lon"]), 
                                   coordD.b2$x1), 
                           lat = c(coordA.b2$y1, 
                                   as.numeric(gps.df[gps.df$point.name == "B","lat"]), 
                                   coordD.b2$y1),
                           point.name = c("A.b2","B.b2", "D.b2"),
                           type = "plot.corner",
                           plot = "b2")
  coords.b2$plant.coords %>%
    rbind(corners.b2) -> data.b2
  
  #---------------------------------------------------------#
  # 4. Pull everything together to export and plot
  gps.df %>%
    rbind(data.b1) %>%
    rbind(data.b2) -> coords
  
  coords %>%
    filter(type == "site.corner") %>%
    arrange(point.name) -> df.site.corners
  p <- ggplot(coords, aes(y = lat, x = lon, label = point.name, color = type)) +
    geom_text() +
    coord_fixed() +
    geom_polygon(data = df.site.corners, fill = NA)
  
  coord.list <- list(coords = coords, p = p)
  
  return(coord.list)
  
}
#mhc.one.gpscoords(gps.coords, plants)

########## SFA-ONE-PRO #########################################################
sfa.one.gpscoords <- function(gps.coords, plants){
  
  #---------------------------------------------------------#
  # 1. Identify relevant data
  curr.gps <- get.curr.gps(gps.coords, plants, site.name = "SFA-ONE-PRO")
  gps.df <- make.gps.long(curr.gps$curr.site)
  # these point names correspond to the site corners, so rename these
  gps.df %>%
    select(point.name, lat, lon) %>%
    mutate(type = "site.corner") %>%
    mutate(plot = NA) -> gps.df
  
  #---------------------------------------------------------#
  # 2. Identify missing site corners -- NA
  # ggplot(gps.df, aes(y = lat, x = lon, label = point.name)) +
  #   geom_text() +
  #   coord_fixed() # this isn't really a rectangle, but that's ok
  
  
  #---------------------------------------------------------#
  # 3. Find plant sample coords 
  origin <- gps.df[gps.df$point.name == "D",]
  xax <- calc_line.eq(x1 = gps.df[gps.df$point.name == "D","lon"], 
                         y1 = gps.df[gps.df$point.name == "D","lat"],
                         x2 = gps.df[gps.df$point.name == "C","lon"],
                         y2 = gps.df[gps.df$point.name == "C","lat"])
  yax <- calc_line.eq(x1 = gps.df[gps.df$point.name == "D","lon"], 
                         y1 = gps.df[gps.df$point.name == "D","lat"],
                         x2 = gps.df[gps.df$point.name == "A","lon"],
                         y2 = gps.df[gps.df$point.name == "A","lat"])
  curr.gps$curr.samps %>%
    filter(Samp %in% 1:8) -> samps
  coords <- calc_plotplant.coords(origin = origin, 
                                     xax.line = xax, 
                                     yax.line = yax, 
                                     plant.samps = samps, 
                                     df = gps.df,
                                     reverse.xdist = F,
                                     reverse.yslope = F,
                                     reverse.ydist = T)
  coords$plant.coords
  coords$plant.coords %>%
    mutate(plot = NA) -> coords$plant.coords
  # add plot.corner coordinates, if applicable -- NA
  # corners.b1 <- data.frame(lon = origin.b1$lon, lat = origin.b1$lat,
  #                          point.name = "D.b1",
  #                          type = "plot.corner",
  #                          plot = "b1")
  # coords.b1$plant.coords %>%
  #   rbind(corners.b1) -> data.b1
  coords$plant.coords -> data
  
  #---------------------------------------------------------#
  # 4. Pull everything together to export and plot
  gps.df %>%
    rbind(data) -> coords
  coords
  coords %>%
    filter(type == "site.corner") %>%
    arrange(point.name) -> df.site.corners
  p <- ggplot(coords, aes(y = lat, x = lon, label = point.name, color = type)) +
    geom_text() +
    coord_fixed() +
    geom_polygon(data = df.site.corners, fill = NA)
  p
  
  coord.list <- list(coords = coords, p = p)
  
  return(coord.list)
  
}
#sfa.one.gpscoords(gps.coords, plants)

########## CGF-MON-PRO #########################################################
cgf.mon.gpscoords <- function(gps.coords, plants){
  
  #---------------------------------------------------------#
  # 1. Identify relevant data
  curr.gps <- get.curr.gps(gps.coords, plants, site.name = "CGF-MON-PRO")
  gps.df <- make.gps.long(curr.gps$curr.site)
  # these point names correspond to the site corners, so rename these
  gps.df %>%
    select(point.name, lat, lon) %>%
    filter(!is.na(lat)) %>%
    mutate(type = "site.corner") %>%
    mutate(plot = NA) -> gps.df
  
  #---------------------------------------------------------#
  # 2. Identify missing site corners 
  ggplot(gps.df, aes(y = lat, x = lon, label = point.name)) +
    geom_text() +
    coord_fixed() # this isn't really a rectangle and missing B
  
  # a. Locate A* as a point along line perpendicular to D-C that crosses through D and has the same lat as A
  # first, find the xax
  xax.line <- calc_line.eq(x1 = gps.df[gps.df$point.name == "D","lon"], 
                           y1 = gps.df[gps.df$point.name == "D","lat"],
                           x2 = gps.df[gps.df$point.name == "C","lon"],
                           y2 = gps.df[gps.df$point.name == "C","lat"])
  # second, calc the yax
  # perpendicular to xax.line and through D
  slope <- -1 / xax.line$slope
  intercept <- calc_intercept_fromSlopePoint(slope = slope, 
                                             x = gps.df[gps.df$point.name == "D","lon"], 
                                             y = gps.df[gps.df$point.name == "D","lat"])
  yax.line <- data.frame(slope, intercept)
  names(yax.line)[2] <- "intercept"
  # third, solve for x.A
  x.A <- (gps.df[gps.df$point.name == "A","lat"] - yax.line$intercept) /  yax.line$slope
  gps.df[gps.df$point.name == "A","lon"] <- x.A
  
  # b. Find B*
  # first create a line parallel to D-A through C
  yax2intercept <- calc_intercept_fromSlopePoint(slope = yax.line$slope, 
                                                 x = gps.df[gps.df$point.name == "C","lon"], 
                                                 y = gps.df[gps.df$point.name == "C","lat"])
  # and create a line parallel to D-C through A
  xax2intercept <- calc_intercept_fromSlopePoint(slope = xax.line$slope, 
                                                 x = gps.df[gps.df$point.name == "A","lon"], 
                                                 y = gps.df[gps.df$point.name == "A","lat"])
  # find B at the intersection of those lines
  new.point <- calc_point_fromLines(slope1 = yax.line$slope, 
                                    intercept1 = yax2intercept, 
                                    slope2 = xax.line$slope, 
                                    intercept2 = xax2intercept)
  new.row <- data.frame(point.name = "B",
             lat = new.point$new.y,
             lon = new.point$new.x,
             type = "site.corner",
             plot = NA)
  gps.df <- rbind(gps.df, new.row)
  # ggplot(gps.df, aes(y = lat, x = lon, label = point.name)) +
  #   geom_text() +
  #   coord_fixed() 
  
  
  #---------------------------------------------------------#
  # 3. Find plant sample coords 
  origin <- gps.df[gps.df$point.name == "D",]
  xax <- calc_line.eq(x1 = gps.df[gps.df$point.name == "D","lon"], 
                      y1 = gps.df[gps.df$point.name == "D","lat"],
                      x2 = gps.df[gps.df$point.name == "C","lon"],
                      y2 = gps.df[gps.df$point.name == "C","lat"])
  yax <- calc_line.eq(x1 = gps.df[gps.df$point.name == "D","lon"], 
                      y1 = gps.df[gps.df$point.name == "D","lat"],
                      x2 = gps.df[gps.df$point.name == "A","lon"],
                      y2 = gps.df[gps.df$point.name == "A","lat"])
  curr.gps$curr.samps %>%
    filter(Samp %in% 1:8) -> samps
  coords <- calc_plotplant.coords(origin = origin, 
                                  xax.line = xax, 
                                  yax.line = yax, 
                                  plant.samps = samps, 
                                  df = gps.df,
                                  reverse.xdist = F,
                                  reverse.yslope = T,
                                  reverse.ydist = F)
  coords$p
  coords$plant.coords
  coords$plant.coords %>%
    mutate(plot = NA) -> coords$plant.coords
  # add plot.corner coordinates, if applicable -- NA
  # corners.b1 <- data.frame(lon = origin.b1$lon, lat = origin.b1$lat,
  #                          point.name = "D.b1",
  #                          type = "plot.corner",
  #                          plot = "b1")
  # coords.b1$plant.coords %>%
  #   rbind(corners.b1) -> data.b1
  coords$plant.coords -> data
  
  #---------------------------------------------------------#
  # 4. Pull everything together to export and plot
  gps.df %>%
    rbind(data) -> coords
  coords
  coords %>%
    filter(type == "site.corner") %>%
    arrange(point.name) -> df.site.corners
  p <- ggplot(coords, aes(y = lat, x = lon, label = point.name, color = type)) +
    geom_text() +
    coord_fixed() +
    geom_polygon(data = df.site.corners, fill = NA)
  p
  
  coord.list <- list(coords = coords, p = p)
  
  return(coord.list)
  
}
#cgf.mon.gpscoords(gps.coords, plants)

########## CGF-MXG-PRO #########################################################
cgf.mxg.gpscoords <- function(gps.coords, plants){
  
  #---------------------------------------------------------#
  # 1. Identify relevant data
  curr.gps <- get.curr.gps(gps.coords, plants, site.name = "CGF-MXG-PRO")
  gps.df <- make.gps.long(curr.gps$curr.site)
  # these point names correspond to the site corners, so rename these
  gps.df %>%
    select(point.name, lat, lon) %>%
    filter(!is.na(lat)) %>%
    mutate(type = "site.corner") %>%
    mutate(plot = NA) -> gps.df
  
  #---------------------------------------------------------#
  # 2. Identify missing site corners 
  ggplot(gps.df, aes(y = lat, x = lon, label = point.name)) +
    geom_text() +
    coord_fixed() # missing A and B
  
  # a. Find points A and B by calculating xax.line slope, taking the perpendicular, and going a distance of 10m
  xax.line <- calc_line.eq(x1 = gps.df[gps.df$point.name == "C","lon"],
                           y1 = gps.df[gps.df$point.name == "C","lat"],
                           x2 = gps.df[gps.df$point.name == "D","lon"],
                           y2 = gps.df[gps.df$point.name == "D","lat"])
  # perpendicular to xax.line and through D
  slope <- -1 / xax.line$slope
  intercept <- calc_intercept_fromSlopePoint(slope = slope, 
                                             x = gps.df[gps.df$point.name == "D","lon"], 
                                             y = gps.df[gps.df$point.name == "D","lat"])
  yax.line <- data.frame(slope, intercept)
  names(yax.line)[2] <- "intercept"
  # add point A 10 meters distance along yax.line
  # load values to convert meters to gps units
  m_per_deg <- load_m.per.deg()
  new.point <- useLineDist_findNewXY(dist = 10/m_per_deg, 
                                     slope = yax.line$slope, 
                                     intercept = yax.line$intercept, 
                                     x0 = gps.df[gps.df$point.name == "D","lon"], 
                                     y0 = gps.df[gps.df$point.name == "D","lat"])
  new.row <- data.frame(point.name = "A",
                        lat = new.point$y1, lon = new.point$x1, 
                        type = "site.corner",
                        plot = NA)
  gps.df <- rbind(gps.df, new.row)
  
  # b. add point B 10 meters distance along yax2.line
  intercept <- calc_intercept_fromSlopePoint(slope = yax.line$slope, 
                                             x = gps.df[gps.df$point.name == "C","lon"], 
                                             y = gps.df[gps.df$point.name == "C","lat"])
  new.point <- useLineDist_findNewXY(dist = 10/m_per_deg, 
                                     slope = yax.line$slope, 
                                     intercept = as.numeric(intercept), 
                                     x0 = gps.df[gps.df$point.name == "C","lon"], 
                                     y0 = gps.df[gps.df$point.name == "C","lat"])
  new.row <- data.frame(point.name = "B",
                        lat = new.point$y1, lon = new.point$x1, 
                        type = "site.corner",
                        plot = NA)
  gps.df <- rbind(gps.df, new.row)
  
  
  #---------------------------------------------------------#
  # 3. Find plant sample coords 
  origin <- gps.df[gps.df$point.name == "D",]
  xax <- calc_line.eq(x1 = gps.df[gps.df$point.name == "D","lon"], 
                      y1 = gps.df[gps.df$point.name == "D","lat"],
                      x2 = gps.df[gps.df$point.name == "C","lon"],
                      y2 = gps.df[gps.df$point.name == "C","lat"])
  yax <- calc_line.eq(x1 = gps.df[gps.df$point.name == "D","lon"], 
                      y1 = gps.df[gps.df$point.name == "D","lat"],
                      x2 = gps.df[gps.df$point.name == "A","lon"],
                      y2 = gps.df[gps.df$point.name == "A","lat"])
  curr.gps$curr.samps %>%
    filter(Samp %in% 1:8) -> samps
  coords <- calc_plotplant.coords(origin = origin, 
                                  xax.line = xax, 
                                  yax.line = yax, 
                                  plant.samps = samps, 
                                  df = gps.df,
                                  reverse.xdist = T,
                                  reverse.yslope = T,
                                  reverse.ydist = F)
  coords$p
  coords$plant.coords
  coords$plant.coords %>%
    mutate(plot = NA) -> coords$plant.coords
  # add plot.corner coordinates, if applicable -- NA
  # corners.b1 <- data.frame(lon = origin.b1$lon, lat = origin.b1$lat,
  #                          point.name = "D.b1",
  #                          type = "plot.corner",
  #                          plot = "b1")
  # coords.b1$plant.coords %>%
  #   rbind(corners.b1) -> data.b1
  coords$plant.coords -> data
  
  #---------------------------------------------------------#
  # 4. Pull everything together to export and plot
  gps.df %>%
    rbind(data) -> coords
  coords
  coords %>%
    filter(type == "site.corner") %>%
    arrange(point.name) -> df.site.corners
  p <- ggplot(coords, aes(y = lat, x = lon, label = point.name, color = type)) +
    geom_text() +
    coord_fixed() +
    geom_polygon(data = df.site.corners, fill = NA)
  p
  
  coord.list <- list(coords = coords, p = p)
  
  return(coord.list)
  
}
#cgf.mxg.gpscoords(gps.coords, plants)

########## OTO-MON-NCD #########################################################
oto.mon.gpscoords <- function(gps.coords, plants){
  
  #---------------------------------------------------------#
  # 1. Identify relevant data
  curr.gps <- get.curr.gps(gps.coords, plants, site.name = "OTO-MON-NCD")
  gps.df <- make.gps.long(curr.gps$curr.site)
  # these point names correspond to the site corners, so rename these
  gps.df %>%
    select(point.name, lat, lon) %>%
    filter(!is.na(lat)) %>%
    mutate(type = "site.corner") %>%
    mutate(plot = NA) %>%
    unique() -> gps.df
  
  #---------------------------------------------------------#
  # 2. Identify missing site corners 
  ggplot(gps.df, aes(y = lat, x = lon, label = point.name)) +
    geom_text() +
    coord_fixed() # missing A and B
  
  # a. points A and B are 6m away from D and C
  # find xax.line
  xax.line <- calc_line.eq(x1 = gps.df[gps.df$point.name == "C","lon"],
                           y1 = gps.df[gps.df$point.name == "C","lat"],
                           x2 = gps.df[gps.df$point.name == "D","lon"],
                           y2 = gps.df[gps.df$point.name == "D","lat"])
  # perpendicular to xax.line and through D
  slope <- -1 / xax.line$slope
  intercept <- calc_intercept_fromSlopePoint(slope = slope, 
                                             x = gps.df[gps.df$point.name == "D","lon"], 
                                             y = gps.df[gps.df$point.name == "D","lat"])
  yax.line <- data.frame(slope, intercept)
  names(yax.line)[2]<-"intercept"
  # add point A 6 meters distance along yax.line
  m_per_deg <- load_m.per.deg()
  new.point <- useLineDist_findNewXY(dist = 6/m_per_deg, 
                                     slope = yax.line$slope, 
                                     intercept = yax.line$intercept, 
                                     x0 = gps.df[gps.df$point.name == "D","lon"], 
                                     y0 = gps.df[gps.df$point.name == "D","lat"])
  new.row <- data.frame(point.name = "A",
                        lat = new.point$y1, 
                        lon = new.point$x1, 
                        type = "site.corner",
                        plot = NA)
  gps.df <- rbind(gps.df, new.row)
  
  # b. add point B 6 meters distance along yax2.line
  intercept <- calc_intercept_fromSlopePoint(slope = yax.line$slope, 
                                             x = gps.df[gps.df$point.name == "C","lon"], 
                                             y = gps.df[gps.df$point.name == "C","lat"])
  new.point <- useLineDist_findNewXY(dist = 6/m_per_deg, 
                                     slope = yax.line$slope, 
                                     intercept = intercept, 
                                     x0 = gps.df[gps.df$point.name == "C","lon"], 
                                     y0 = gps.df[gps.df$point.name == "C","lat"])
  new.row <- data.frame(point.name = "B",
                        lat = new.point$y1, 
                        lon = new.point$x1, 
                        type = "site.corner",
                        plot = NA)
  gps.df <- rbind(gps.df, new.row)
  # ggplot(gps.df, aes(y = lat, x = lon, label = point.name)) +
  #   geom_text() +
  #   coord_fixed() 
  
  #---------------------------------------------------------#
  # 3. Find plant sample coords 
  origin <- gps.df[gps.df$point.name == "D",]
  xax <- calc_line.eq(x1 = gps.df[gps.df$point.name == "D","lon"], 
                      y1 = gps.df[gps.df$point.name == "D","lat"],
                      x2 = gps.df[gps.df$point.name == "C","lon"],
                      y2 = gps.df[gps.df$point.name == "C","lat"])
  yax <- calc_line.eq(x1 = gps.df[gps.df$point.name == "D","lon"], 
                      y1 = gps.df[gps.df$point.name == "D","lat"],
                      x2 = gps.df[gps.df$point.name == "A","lon"],
                      y2 = gps.df[gps.df$point.name == "A","lat"])
  curr.gps$curr.samps %>%
    filter(Samp %in% 1:8) -> samps
  # the 3 sub plots are N fertilization after 12, 11, and 10 m along the string from D to C
  # Nfert1 origin is D
  samps[1:3, "Xcoord.m"]
  # Nfert2 origin in D + 12m
  samps[4:6, "Xcoord.m"] <- samps[4:6, "Xcoord.m"] + 12
  # Nfert2 origin in D + 12m + 11m
  samps[7:8, "Xcoord.m"] <- samps[7:8, "Xcoord.m"] + 12 + 11
  coords <- calc_plotplant.coords(origin = origin, 
                                  xax.line = xax, 
                                  yax.line = yax, 
                                  plant.samps = samps, 
                                  df = gps.df,
                                  reverse.xdist = F,
                                  reverse.yslope = T,
                                  reverse.ydist = F)
  coords$p
  coords$plant.coords
  coords$plant.coords %>%
    mutate(plot = NA) -> coords$plant.coords
  # add plot.corner coordinates, if applicable -- NA
  # corners.b1 <- data.frame(lon = origin.b1$lon, lat = origin.b1$lat,
  #                          point.name = "D.b1",
  #                          type = "plot.corner",
  #                          plot = "b1")
  # coords.b1$plant.coords %>%
  #   rbind(corners.b1) -> data.b1
  coords$plant.coords -> data
  
  # need to go back and adjust the full plot size
  # a. find C* by adding full distance from Nfert1 + Nfert2 + Nfert3
  full.dist <- 12 + 11 + 10
  new.point <- useLineDist_findNewXY(dist = full.dist/m_per_deg, 
                                     slope = xax.line$slope, 
                                     intercept = xax.line$intercept, 
                                     x0 = gps.df[gps.df$point.name == "D","lon"], 
                                     y0 = gps.df[gps.df$point.name == "D","lat"])
  gps.df[gps.df$point.name=="C",c("lat","lon")] <- c(new.point$y1, new.point$x1)
  
  # b. find B* by adding full distance from Nfert1 + Nfert2 + Nfert3
  new.intercept <- calc_intercept_fromSlopePoint(slope = xax.line$slope, 
                                                 x = gps.df[gps.df$point.name == "A","lon"], 
                                                 y = gps.df[gps.df$point.name == "A","lat"])
  new.point <- useLineDist_findNewXY(dist = full.dist/m_per_deg, 
                                     slope = xax.line$slope, 
                                     intercept = new.intercept, 
                                     x0 = gps.df[gps.df$point.name == "A","lon"], 
                                     y0 = gps.df[gps.df$point.name == "A","lat"])
  gps.df[gps.df$point.name=="B",c("lat","lon")] <- c(new.point$y1, new.point$x1)
  

  #---------------------------------------------------------#
  # 4. Pull everything together to export and plot
  gps.df %>%
    rbind(data) -> coords
  coords
  coords %>%
    filter(type == "site.corner") %>%
    arrange(point.name) -> df.site.corners
  p <- ggplot(coords, aes(y = lat, x = lon, label = point.name, color = type)) +
    geom_text() +
    coord_fixed() +
    geom_polygon(data = df.site.corners, fill = NA)
  p
  
  coord.list <- list(coords = coords, p = p)
  
  return(coord.list)
  
}
#oto.mon.gpscoords(gps.coords, plants)

########## OTO-MXT-NCD #########################################################
oto.mxt.gpscoords <- function(gps.coords, plants){
  
  #---------------------------------------------------------#
  # 1. Identify relevant data
  curr.gps <- get.curr.gps(gps.coords, plants, site.name = "OTO-MXT-NCD")
  gps.df <- make.gps.long(curr.gps$curr.site)
  # these point names correspond to the site corners, so rename these
  gps.df %>%
    select(point.name, lat, lon) %>%
    filter(!is.na(lat)) %>%
    mutate(type = "site.corner") %>%
    mutate(plot = NA) %>%
    unique() -> gps.df

  
  #---------------------------------------------------------#
  # 2. Identify missing site corners 
  ggplot(gps.df, aes(y = lat, x = lon, label = point.name)) +
    geom_text() +
    coord_fixed() 
  # these plot coordinates look very odd considering that the distance between D and A should be 6.5 meters and the distance between D and C should be 8.67m * 9
  # remake the plot using just point D (origin) and the yax.line from D to A
  
  # a. find D to A line
  yax.line <- calc_line.eq(x1 = gps.df[gps.df$point.name == "D","lon"],
                           y1 = gps.df[gps.df$point.name == "D","lat"],
                           x2 = gps.df[gps.df$point.name == "A","lon"],
                           y2 = gps.df[gps.df$point.name == "A","lat"])
  # calc new A -- 6.5m along yax.line from D
  m_per_deg <- load_m.per.deg()
  new.point <- useLineDist_findNewXY(dist = 6.5/m_per_deg, 
                                     slope = yax.line$slope, 
                                     intercept = yax.line$intercept, 
                                     x0 = gps.df[gps.df$point.name == "D","lon"], 
                                     y0 = gps.df[gps.df$point.name == "D","lat"])
  gps.df[gps.df$point.name=="A",c("lat","lon")] <- c(new.point$y1, new.point$x1)
  
  # b. find xax.line from D to C
  slope <- -1/yax.line$slope # 0
  intercept <- calc_intercept_fromSlopePoint(slope = slope, 
                                             x = gps.df[gps.df$point.name == "D","lon"], 
                                             y = gps.df[gps.df$point.name == "D","lat"])
  slope
  intercept # -Inf
  # yax is a horizontal line, so xax is just a vertical line
  # find C* by subtracting xdistance from D lat
  xdist <- (8.67 * 9)/m_per_deg
  C.lat <- gps.df[gps.df$point.name == "D","lat"] - xdist 
  C.lon <- gps.df[gps.df$point.name == "D","lon"]
  new.row <- data.frame(point.name = "C",
                        lat = C.lat, 
                        lon = C.lon, 
                        type = "site.corner",
                        plot = NA)
  gps.df <- rbind(gps.df, new.row)
  gps.df
  
  # c. find point B* using 6.5m away from C* along horizontal line, slope = 0
  B.lat <- C.lat
  B.lon <- C.lon + 6.5/m_per_deg
  gps.df[gps.df$point.name == "B", c("lat","lon")] <- c(B.lat, B.lon)
  
  
  #---------------------------------------------------------#
  # 3. Find plant sample coords 
  origin <- gps.df[gps.df$point.name == "D",]
  xax <- calc_line.eq(x1 = gps.df[gps.df$point.name == "D","lon"], 
                      y1 = gps.df[gps.df$point.name == "D","lat"],
                      x2 = gps.df[gps.df$point.name == "C","lon"],
                      y2 = gps.df[gps.df$point.name == "C","lat"])
  yax <- calc_line.eq(x1 = gps.df[gps.df$point.name == "D","lon"], 
                      y1 = gps.df[gps.df$point.name == "D","lat"],
                      x2 = gps.df[gps.df$point.name == "A","lon"],
                      y2 = gps.df[gps.df$point.name == "A","lat"])
  curr.gps$curr.samps %>%
    filter(Samp %in% 1:8) -> samps
  # coords <- calc_plotplant.coords(origin = origin, 
  #                                 xax.line = xax, 
  #                                 yax.line = yax, 
  #                                 plant.samps = samps, 
  #                                 df = gps.df,
  #                                 reverse.xdist = F,
  #                                 reverse.yslope = T,
  #                                 reverse.ydist = F)
  # this isn't working because one of the slopes is Inf
  X.dist <- samps$Xcoord.m/m_per_deg
  xax.lats <- origin$lat - X.dist
  Y.dist <- samps$Ycoord.m/m_per_deg
  yax.lons <- origin$lon + Y.dist
  coords.df <- data.frame(lon = yax.lons, 
                          lat = xax.lats, 
                          point.name = paste0("p",samps$Samp), 
                          type = "plant",
                          plot = NA)
  coords.df -> data
  
  #---------------------------------------------------------#
  # 4. Pull everything together to export and plot
  gps.df %>%
    rbind(data) -> coords
  coords
  coords %>%
    filter(type == "site.corner") %>%
    arrange(point.name) -> df.site.corners
  p <- ggplot(coords, aes(y = lat, x = lon, label = point.name, color = type)) +
    geom_text() +
    coord_fixed() +
    geom_polygon(data = df.site.corners, fill = NA)
  p
  
  coord.list <- list(coords = coords, p = p)
  
  return(coord.list)
  
}
#oto.mxt.gpscoords(gps.coords, plants)

########## CCR-ONE-NCD #########################################################
ccr.one.gpscoords <- function(gps.coords, plants){
  
  #---------------------------------------------------------#
  # 1. Identify relevant data
  curr.gps <- get.curr.gps(gps.coords, plants, site.name = "CCR-ONE-NCD")
  gps.df <- make.gps.long(curr.gps$curr.site)
  # these point names correspond to the site corners, so rename these
  gps.df %>%
    select(point.name, lat, lon) %>%
    filter(!is.na(lat)) %>%
    mutate(type = "site.corner") %>%
    mutate(plot = NA) %>%
    unique() -> gps.df
  
  #---------------------------------------------------------#
  # 2. Identify missing site corners 
  ggplot(gps.df, aes(y = lat, x = lon, label = point.name)) +
    geom_text() +
    coord_fixed() 
  # this looks like a trapazoid, but that's ok... just follow from D
  
  
  
  #---------------------------------------------------------#
  # 3. Find plant sample coords 
  origin <- gps.df[gps.df$point.name == "D",]
  xax <- calc_line.eq(x1 = gps.df[gps.df$point.name == "D","lon"], 
                      y1 = gps.df[gps.df$point.name == "D","lat"],
                      x2 = gps.df[gps.df$point.name == "C","lon"],
                      y2 = gps.df[gps.df$point.name == "C","lat"])
  yax <- calc_line.eq(x1 = gps.df[gps.df$point.name == "D","lon"], 
                      y1 = gps.df[gps.df$point.name == "D","lat"],
                      x2 = gps.df[gps.df$point.name == "A","lon"],
                      y2 = gps.df[gps.df$point.name == "A","lat"])
  curr.gps$curr.samps %>%
    filter(Samp %in% 1:8) -> samps
  coords <- calc_plotplant.coords(origin = origin,
                                  xax.line = xax,
                                  yax.line = yax,
                                  plant.samps = samps,
                                  df = gps.df,
                                  reverse.xdist = T,
                                  reverse.yslope = T,
                                  reverse.ydist = F)
  coords$p
  coords$plant.coords %>%
    mutate(plot = NA) -> data
  
  # need to go back and adjust B so the plants fit in the plot
  # use y slope and C point to find yax2 intercept
  yax2.intercept <- calc_intercept_fromSlopePoint(slope = yax$slope, 
                                                  x = gps.df[gps.df$point.name == "C","lon"], 
                                                  y = gps.df[gps.df$point.name == "C","lat"])
  # use x slope and A point to find xax2 intercept
  xax2.intercept <- calc_intercept_fromSlopePoint(slope = xax$slope, 
                                                  x = gps.df[gps.df$point.name == "A","lon"], 
                                                  y = gps.df[gps.df$point.name == "A","lat"])
  # find the intersection of xax2 and yax2
  new.point <- calc_point_fromLines(slope1 = yax$slope, 
                                    intercept1 = yax2.intercept, 
                                    slope2 = xax$slope, 
                                    intercept2 = xax2.intercept)
  gps.df[gps.df$point.name == "B",c("lon","lat")] <- c(new.point$new.x, new.point$new.y)
  
  
  #---------------------------------------------------------#
  # 4. Pull everything together to export and plot
  gps.df %>%
    rbind(data) -> coords
  coords
  coords %>%
    filter(type == "site.corner") %>%
    arrange(point.name) -> df.site.corners
  p <- ggplot(coords, aes(y = lat, x = lon, label = point.name, color = type)) +
    geom_text() +
    coord_fixed() +
    geom_polygon(data = df.site.corners, fill = NA)
  p
  
  coord.list <- list(coords = coords, p = p)
  
  return(coord.list)
  
}
#ccr.one.gpscoords(gps.coords, plants)

########## CRE-MXG-NCD #########################################################
cre.mxg.gpscoords <- function(gps.coords, plants){
  
  #---------------------------------------------------------#
  # 1. Identify relevant data
  curr.gps <- get.curr.gps(gps.coords, plants, site.name = "CRE-MXG-NCD")
  gps.df <- make.gps.long(curr.gps$curr.site)
  # these point names correspond to the site corners, so rename these
  gps.df %>%
    select(point.name, lat, lon) %>%
    filter(!is.na(lat)) %>%
    mutate(type = "site.corner") %>%
    mutate(plot = NA) %>%
    unique() -> gps.df
  
  #---------------------------------------------------------#
  # 2. Identify missing site corners 
  ggplot(gps.df, aes(y = lat, x = lon, label = point.name)) +
    geom_text() +
    coord_fixed() 
  # this is messed up looking, especially B. go with D's coordinates
  
  # a. find the xax.line (D to C)
  xax.line <- calc_line.eq(x1 = gps.df[gps.df$point.name == "D","lon"],
                           y1 = gps.df[gps.df$point.name == "D","lat"],
                           x2 = gps.df[gps.df$point.name == "C","lon"],
                           y2 = gps.df[gps.df$point.name == "C","lat"])
  
  # b. find the yax.line using the point D and the slope perpendicular to xax.line
  slope <- -1 / xax.line$slope
  intercept <- calc_intercept_fromSlopePoint(slope = slope, 
                                             x = gps.df[gps.df$point.name == "D","lon"], 
                                             y = gps.df[gps.df$point.name == "D","lat"])
  yax.line <- data.frame(slope, intercept)
  names(yax.line)[2] <- "intercept"
  
  # c. find A from yax.line and 77m distance
  m_per_deg <- load_m.per.deg()
  new.point <- useLineDist_findNewXY(dist = (77/m_per_deg)*-1, 
                                     slope = yax.line$slope, 
                                     intercept = yax.line$intercept, 
                                     x0 = gps.df[gps.df$point.name == "D","lon"], 
                                     y0 = gps.df[gps.df$point.name == "D","lat"])
  gps.df[gps.df$point.name == "A", c("lat","lon")] <- c(new.point$y1, new.point$x1)
  
  # d. find the yax2.line 
  yax2.intercept <- calc_intercept_fromSlopePoint(slope = slope, 
                                                  x = gps.df[gps.df$point.name == "C","lon"], 
                                                  y = gps.df[gps.df$point.name == "C","lat"])
  
  # find B from C and yax2 line
  new.point <- useLineDist_findNewXY(dist = (77/m_per_deg)*-1, 
                                     slope = yax.line$slope, 
                                     intercept = yax2.intercept, 
                                     x0 = gps.df[gps.df$point.name == "C","lon"], 
                                     y0 = gps.df[gps.df$point.name == "C","lat"])
  gps.df[gps.df$point.name == "B", c("lat","lon")] <- c(new.point$y1, new.point$x1)
  # ggplot(gps.df, aes(y = lat, x = lon, label = point.name)) +
  #   geom_text() +
  #   coord_fixed() 
  
  
  #---------------------------------------------------------#
  # 3. Find plant sample coords 
  origin <- gps.df[gps.df$point.name == "D",]
  xax <- calc_line.eq(x1 = gps.df[gps.df$point.name == "D","lon"], 
                      y1 = gps.df[gps.df$point.name == "D","lat"],
                      x2 = gps.df[gps.df$point.name == "C","lon"],
                      y2 = gps.df[gps.df$point.name == "C","lat"])
  yax <- calc_line.eq(x1 = gps.df[gps.df$point.name == "D","lon"], 
                      y1 = gps.df[gps.df$point.name == "D","lat"],
                      x2 = gps.df[gps.df$point.name == "A","lon"],
                      y2 = gps.df[gps.df$point.name == "A","lat"])
  curr.gps$curr.samps %>%
    filter(Samp %in% 1:8) -> samps
  coords <- calc_plotplant.coords(origin = origin,
                                  xax.line = xax,
                                  yax.line = yax,
                                  plant.samps = samps,
                                  df = gps.df,
                                  reverse.xdist = T,
                                  reverse.yslope = T,
                                  reverse.ydist = T)
  coords$p
  coords$plant.coords %>%
    mutate(plot = NA) -> data

  
  #---------------------------------------------------------#
  # 4. Pull everything together to export and plot
  gps.df %>%
    rbind(data) -> coords
  coords
  coords %>%
    filter(type == "site.corner") %>%
    arrange(point.name) -> df.site.corners
  p <- ggplot(coords, aes(y = lat, x = lon, label = point.name, color = type)) +
    geom_text() +
    coord_fixed() +
    geom_polygon(data = df.site.corners, fill = NA)
  p
  
  coord.list <- list(coords = coords, p = p)
  
  return(coord.list)
  
}
#cre.mxg.gpscoords(gps.coords, plants)

########## CRE-MXT-NCD #########################################################
cre.mxt.gpscoords <- function(gps.coords, plants){
  
  #---------------------------------------------------------#
  # 1. Identify relevant data
  curr.gps <- get.curr.gps(gps.coords, plants, site.name = "CRE-MXT-NCD")
  gps.df <- make.gps.long(curr.gps$curr.site)
  gps.df
  # these correspond to the origin of each plot at the site
  gps.df %>%
    separate(multiplot.notes, into = c("plot.num","thing"), sep = 1) %>%
    mutate(plot = paste0("b", plot.num)) %>%
    mutate(point.name = paste(plot, point.name, sep = ".")) %>%
    select(point.name, lat, lon, plot) %>%
    filter(!is.na(lat)) %>%
    mutate(type = "plot.corner") -> gps.df
    
  
  #---------------------------------------------------------#
  # 2. Identify missing site corners 
  p1 <- ggplot(gps.df, aes(y = lat, x = lon, label = point.name)) +
    geom_text() +
    coord_fixed() 
  p1
  
  # a. use the origins for plots 1 and 2 to create overall x axis line
  xax.line <- calc_line.eq(x1 = gps.df[gps.df$plot == "b1", "lon"],
                           y1 = gps.df[gps.df$plot == "b1", "lat"],
                           x2 = gps.df[gps.df$plot == "b2", "lon"],
                           y2 = gps.df[gps.df$plot == "b2", "lat"])
  names(xax.line) <- c("slope","intercept")
  #xax.line # this works for (1) overall site x axis, (2) block1 x axis, and (3) block2 xaxis 
  p2 <- p1 + geom_abline(slope = xax.line$slope, intercept = xax.line$intercept)
  p2
  
  # b. find the intercept for the xax.line specifically for block3 and block4
  xax.b3.intercept <- calc_intercept_fromSlopePoint(slope = xax.line$slope, 
                                                    x = gps.df[gps.df$plot == "b3", "lon"], 
                                                    y = gps.df[gps.df$plot == "b3", "lat"])
  xax.b4.intercept <- calc_intercept_fromSlopePoint(slope = xax.line$slope, 
                                                    x = gps.df[gps.df$plot == "b4", "lon"], 
                                                    y = gps.df[gps.df$plot == "b4", "lat"])
  
  # c. find the overall y axis line using the origin for block1 point and the slope perpendicular to the xax.line
  yax.slope <- -1/xax.line$slope
  yax.intercept <- calc_intercept_fromSlopePoint(slope = yax.slope, 
                                                 x = gps.df[gps.df$plot == "b1", "lon"], 
                                                 y = gps.df[gps.df$plot == "b1", "lat"])
  names(yax.intercept) <- "intercept"
  yax.line <- data.frame(slope = yax.slope, intercept = yax.intercept)
  p3 <- p2 + 
    geom_abline(slope = yax.slope, intercept = yax.intercept)
  p3
  
  # d. find the intercept for the yax.line specifically for block2, block3, and block4
  yax.b2.intercept <- calc_intercept_fromSlopePoint(slope = yax.line$slope, 
                                                    x = gps.df[gps.df$plot == "b2", "lon"], 
                                                    y = gps.df[gps.df$plot == "b2", "lat"])
  yax.b3.intercept <- calc_intercept_fromSlopePoint(slope = yax.line$slope, 
                                                    x = gps.df[gps.df$plot == "b3", "lon"], 
                                                    y = gps.df[gps.df$plot == "b3", "lat"])
  yax.b4.intercept <- calc_intercept_fromSlopePoint(slope = yax.line$slope, 
                                                    x = gps.df[gps.df$plot == "b4", "lon"], 
                                                    y = gps.df[gps.df$plot == "b4", "lat"])
  names(yax.b2.intercept) <- names(yax.b3.intercept) <- names(yax.b4.intercept) <- "intercept"
  p4 <- p3 + geom_abline(slope = yax.line$slope, intercept = yax.b2.intercept) + # y for b2
    geom_abline(slope = yax.line$slope, intercept = yax.b3.intercept) + # y for b3
    geom_abline(slope = yax.line$slope, intercept = yax.b4.intercept) # y for b4
  p4
  
  # summarize the site corners if possible
  new.row <- data.frame(point.name = "D", 
                        lat = as.numeric(gps.df[gps.df$point.name == "b1.D","lat"]),
                        lon = as.numeric(gps.df[gps.df$point.name == "b1.D","lon"]), 
                        plot = NA,
                        type = "site.corner")
  
  gps.df <- rbind(gps.df, new.row)
  
  
  #---------------------------------------------------------#
  # 3. Find plant sample coords 
  
  # block1
  xax.line.b1 <- data.frame(slope = xax.line$slope,
                            intercept = xax.line$intercept)
  yax.line.b1 <- data.frame(slope = yax.line$slope,
                            intercept = yax.line$intercept)
  curr.gps$curr.samps %>%
    filter(Samp %in% 1:2) -> samps.b1
  coords.b1 <- calc_plotplant.coords(origin = gps.df[gps.df$point.name == "b1.D",],
                                  xax.line = xax.line.b1,
                                  yax.line = yax.line.b1,
                                  plant.samps = samps.b1,
                                  df = gps.df,
                                  reverse.xdist = F,
                                  reverse.yslope = T,
                                  reverse.ydist = F)
  coords.b1$p
  coords.b1$plant.coords %>%
    mutate(plot = "b1") -> data.b1
 
  
  # block2
  xax.line.b2 <- data.frame(slope = xax.line$slope,
                            intercept = xax.line$intercept)
  yax.line.b2 <- data.frame(slope = yax.line$slope,
                            intercept = yax.b2.intercept)
  curr.gps$curr.samps %>%
    filter(Samp %in% 3:4) -> samps.b2
  coords.b2 <- calc_plotplant.coords(origin = gps.df[gps.df$point.name == "b2.D",],
                                     xax.line = xax.line.b2,
                                     yax.line = yax.line.b2,
                                     plant.samps = samps.b2,
                                     df = gps.df,
                                     reverse.xdist = F,
                                     reverse.yslope = T,
                                     reverse.ydist = F)
  coords.b2$p
  coords.b2$plant.coords %>%
    mutate(plot = "b2") -> data.b2
  
  
  # block3
  xax.line.b3 <- data.frame(slope = xax.line$slope,
                            intercept = as.numeric(xax.b3.intercept))
  yax.line.b3 <- data.frame(slope = yax.line$slope,
                            intercept = yax.b3.intercept)
  curr.gps$curr.samps %>%
    filter(Samp %in% 5:6) -> samps.b3
  coords.b3 <- calc_plotplant.coords(origin = gps.df[gps.df$point.name == "b3.D",],
                                     xax.line = xax.line.b3,
                                     yax.line = yax.line.b3,
                                     plant.samps = samps.b3,
                                     df = gps.df,
                                     reverse.xdist = F,
                                     reverse.yslope = T,
                                     reverse.ydist = F)
  coords.b3$p
  coords.b3$plant.coords %>%
    mutate(plot = "b3") -> data.b3
  
  #####
  # block4
  xax.line.b4 <- data.frame(slope = xax.line$slope,
                            intercept = as.numeric(xax.b4.intercept))
  yax.line.b4 <- data.frame(slope = yax.line$slope,
                            intercept = yax.b4.intercept)
  curr.gps$curr.samps %>%
    filter(Samp %in% 7:8) -> samps.b4
  coords.b4 <- calc_plotplant.coords(origin = gps.df[gps.df$point.name == "b4.D",],
                                     xax.line = xax.line.b4,
                                     yax.line = yax.line.b4,
                                     plant.samps = samps.b4,
                                     df = gps.df,
                                     reverse.xdist = F,
                                     reverse.yslope = T,
                                     reverse.ydist = F)
  coords.b4$p
  coords.b4$plant.coords %>%
    mutate(plot = "b4") -> data.b4

  
  # go back and estimate site corners
  # a. find C
  #first, find the point at the intersection of xax.line and yax.line.b3
  new.point <- calc_point_fromLines(slope1 = xax.line$slope, 
                                    intercept1 = xax.line$intercept, 
                                    slope2 = yax.line.b4$slope, 
                                    intercept2 = yax.line.b4$intercept)
  #second, calc 42m from new point along xax line
  m_per_deg <- load_m.per.deg()
  new.point <- useLineDist_findNewXY(dist = 42/m_per_deg, 
                                     slope = xax.line$slope, 
                                     intercept = xax.line$intercept, 
                                     x0 = new.point$new.x, 
                                     y0 = new.point$new.y)
  new.row <- data.frame(point.name = "C",
                        lat = new.point$y1,
                        lon = new.point$x1,
                        plot = NA,
                        type = "site.corner")
  gps.df <- rbind(gps.df, new.row)
  
  # b. find A for the whole site
  # first, find A for block3 by going 6.3m along yax.lineb3 from origin.b3
  new.point <- useLineDist_findNewXY(dist = 6.3/m_per_deg, 
                                     slope = yax.line.b3$slope, 
                                     intercept = yax.line.b3$intercept, 
                                     x0 = gps.df[gps.df$point.name == "b3.D","lon"], 
                                     y0 = gps.df[gps.df$point.name == "b3.D","lat"])
  # second, use A for block3 and xax slope to find xax2 intercept
  xax2.intercept <- calc_intercept_fromSlopePoint(slope = xax.line$slope, 
                                                  x = new.point$x1, 
                                                  y = new.point$y1)
  # third, find the intersection of xax2 line and yax line --> A
  new.point <- calc_point_fromLines(slope1 = yax.line$slope, 
                                    intercept1 = yax.line$intercept, 
                                    slope2 = xax.line$slope, 
                                    intercept2 = xax2.intercept)
  new.row <- data.frame(point.name = "A",
                        lat = new.point$new.y,
                        lon = new.point$new.x,
                        plot = NA,
                        type = "site.corner")
  gps.df <- rbind(gps.df, new.row)
  
  # c. find B for the whole site
  # first, find yax2 intercept using yax slope and point C
  yax2.intercept <- calc_intercept_fromSlopePoint(slope = yax.line$slope, 
                                                  x = gps.df[gps.df$point.name == "C","lon"], 
                                                  y = gps.df[gps.df$point.name == "C","lat"])
  # second, find the intersection of yax2 and xax2 lines --> B
  new.point <- calc_point_fromLines(slope1 = yax.line$slope, 
                                    intercept1 = yax2.intercept, 
                                    slope2 = xax.line$slope, 
                                    intercept2 = xax2.intercept)
  new.row <- data.frame(point.name = "B",
                        lat = new.point$new.y,
                        lon = new.point$new.x,
                        plot = NA,
                        type = "site.corner")
  gps.df <- rbind(gps.df, new.row)
  ggplot(gps.df, aes(y = lat, x = lon, label = point.name, color = type)) +
    geom_text() +
    coord_fixed()
  
  #---------------------------------------------------------#
  # 4. Pull everything together to export and plot
  gps.df %>%
    rbind(data.b1) %>%
    rbind(data.b2) %>%
    rbind(data.b3) %>%
    rbind(data.b4) -> coords
  coords %>%
    filter(type == "site.corner") %>%
    arrange(point.name) -> df.site.corners
  p <- ggplot(coords, aes(y = lat, x = lon, label = point.name, color = type)) +
    geom_text() +
    coord_fixed() +
    geom_polygon(data = df.site.corners, fill = NA)
  p
  
  coord.list <- list(coords = coords, p = p)
  
  return(coord.list)

}
#cre.mxt.gpscoords(gps.coords, plants)

########## UCP-MXG-NCD #########################################################
ucp.mxg.gpscoords <- function(gps.coords, plants){
  
  #---------------------------------------------------------#
  # 1. Identify relevant data
  curr.gps <- get.curr.gps(gps.coords, plants, site.name = "UCP-MXG-NCD")
  gps.df <- make.gps.long(curr.gps$curr.site)
  gps.df
  # these correspond to site corners
  gps.df %>%
    mutate(plot = NA) %>%
    select(point.name, lat, lon, plot) %>%
    filter(!is.na(lat)) %>%
    mutate(type = "site.corner") -> gps.df
  
  
  #---------------------------------------------------------#
  # 2. Identify missing site corners 
  p1 <- ggplot(gps.df, aes(y = lat, x = lon, label = point.name)) +
    geom_text() +
    coord_fixed() 
  p1
  
  # a. find B
  # first, define xax.line and yax.line
  xax.line <- calc_line.eq(x1 = gps.df[gps.df$point.name == "D","lon"],
                           y1 = gps.df[gps.df$point.name == "D","lat"],
                           x2 = gps.df[gps.df$point.name == "C","lon"],
                           y2 = gps.df[gps.df$point.name == "C","lat"])
  yax.line <- calc_line.eq(x1 = gps.df[gps.df$point.name == "D","lon"],
                           y1 = gps.df[gps.df$point.name == "D","lat"],
                           x2 = gps.df[gps.df$point.name == "A","lon"],
                           y2 = gps.df[gps.df$point.name == "A","lat"])
  # second, find yax2 intercept using C and yax slope
  yax2.intercept <- calc_intercept_fromSlopePoint(slope = yax.line$slope, 
                                                  x = gps.df[gps.df$point.name == "C","lon"], 
                                                  y = gps.df[gps.df$point.name == "C","lat"])
  # third, find B using yax2 line and 18m distance from C
  m_per_deg <- load_m.per.deg()
  new.point <- useLineDist_findNewXY(dist = (18/m_per_deg)*-1, 
                                     slope = yax.line$slope, 
                                     intercept = yax2.intercept, 
                                     x0 = gps.df[gps.df$point.name == "C","lon"], 
                                     y0 = gps.df[gps.df$point.name == "C","lat"])
  new.row <- data.frame(point.name = "B",
                        lat = new.point$y1,
                        lon = new.point$x1,
                        plot = NA,
                        type = "site.corner")
  gps.df <- rbind(gps.df, new.row)
  
  
  
  #---------------------------------------------------------#
  # 3. Find plant sample coords 
  origin <- gps.df[gps.df$point.name == "D",]
  xax <- calc_line.eq(x1 = gps.df[gps.df$point.name == "D","lon"], 
                      y1 = gps.df[gps.df$point.name == "D","lat"],
                      x2 = gps.df[gps.df$point.name == "C","lon"],
                      y2 = gps.df[gps.df$point.name == "C","lat"])
  yax <- calc_line.eq(x1 = gps.df[gps.df$point.name == "D","lon"], 
                      y1 = gps.df[gps.df$point.name == "D","lat"],
                      x2 = gps.df[gps.df$point.name == "A","lon"],
                      y2 = gps.df[gps.df$point.name == "A","lat"])
  curr.gps$curr.samps %>%
    filter(Samp %in% 1:8) -> samps
  coords <- calc_plotplant.coords(origin = origin,
                                  xax.line = xax,
                                  yax.line = yax,
                                  plant.samps = samps,
                                  df = gps.df,
                                  reverse.xdist = T,
                                  reverse.yslope = T,
                                  reverse.ydist = T)
  coords$p
  coords$plant.coords %>%
    mutate(plot = NA) -> data
  
  
  #---------------------------------------------------------#
  # 4. Pull everything together to export and plot
  gps.df %>%
    rbind(data) -> coords
  coords
  coords %>%
    filter(type == "site.corner") %>%
    arrange(point.name) -> df.site.corners
  p <- ggplot(coords, aes(y = lat, x = lon, label = point.name, color = type)) +
    geom_text() +
    coord_fixed() +
    geom_polygon(data = df.site.corners, fill = NA)
  p
  
  coord.list <- list(coords = coords, p = p)
  
  return(coord.list)
  
}
#ucp.mxg.gpscoords(gps.coords, plants)

########## LCO-MXT-COM #########################################################
lco.mxt.gpscoords <- function(gps.coords, plants){
  
  #---------------------------------------------------------#
  # 1. Identify relevant data
  curr.gps <- get.curr.gps(gps.coords, plants, site.name = "LCO-MXT-COM")
  gps.df <- make.gps.long(curr.gps$curr.site)
  gps.df
  # these correspond to site corners
  gps.df %>%
    mutate(plot = NA) %>%
    select(point.name, lat, lon, plot) %>%
    filter(!is.na(lat)) %>%
    mutate(type = "site.corner") -> gps.df
  
 
  #---------------------------------------------------------#
  # 2. Identify missing site corners 
  p1 <- ggplot(gps.df, aes(y = lat, x = lon, label = point.name)) +
    geom_text() +
    coord_fixed() 
  p1
  
  # define xax.line and yax.line
  xax.line <- calc_line.eq(x1 = gps.df[gps.df$point.name == "D","lon"],
                           y1 = gps.df[gps.df$point.name == "D","lat"],
                           x2 = gps.df[gps.df$point.name == "C","lon"],
                           y2 = gps.df[gps.df$point.name == "C","lat"])
  yax.line <- calc_line.eq(x1 = gps.df[gps.df$point.name == "D","lon"],
                           y1 = gps.df[gps.df$point.name == "D","lat"],
                           x2 = gps.df[gps.df$point.name == "A","lon"],
                           y2 = gps.df[gps.df$point.name == "A","lat"])
  # find yax2 intercept using C and yax slope
  yax2.intercept <- calc_intercept_fromSlopePoint(slope = yax.line$slope, 
                                                  x = gps.df[gps.df$point.name == "C","lon"], 
                                                  y = gps.df[gps.df$point.name == "C","lat"])
  # find xax2 intercept using A and xax slope
  xax2.intercept <- calc_intercept_fromSlopePoint(slope = xax.line$slope, 
                                                  x = gps.df[gps.df$point.name == "A","lon"], 
                                                  y = gps.df[gps.df$point.name == "A","lat"])
  # find B using the intersection of yax2 and xax2
  new.point <- calc_point_fromLines(slope1 = xax.line$slope, 
                                    intercept1 = xax2.intercept, 
                                    slope2 = yax.line$slope, 
                                    intercept2 = yax2.intercept)
  new.row <- data.frame(point.name = "B",
                        lat = new.point$new.y,
                        lon = new.point$new.x,
                        plot = NA,
                        type = "site.corner")
  gps.df <- rbind(gps.df, new.row)
  # ggplot(gps.df, aes(x = lon, y = lat, label = point.name)) +
  #   geom_text() +
  #   coord_fixed()
  
  
  #---------------------------------------------------------#
  # 3. Find plant sample coords 
  origin <- gps.df[gps.df$point.name == "D",]
  xax <- calc_line.eq(x1 = gps.df[gps.df$point.name == "D","lon"], 
                      y1 = gps.df[gps.df$point.name == "D","lat"],
                      x2 = gps.df[gps.df$point.name == "C","lon"],
                      y2 = gps.df[gps.df$point.name == "C","lat"])
  yax <- calc_line.eq(x1 = gps.df[gps.df$point.name == "D","lon"], 
                      y1 = gps.df[gps.df$point.name == "D","lat"],
                      x2 = gps.df[gps.df$point.name == "A","lon"],
                      y2 = gps.df[gps.df$point.name == "A","lat"])
  curr.gps$curr.samps %>%
    filter(Samp %in% 1:8) -> samps
  coords <- calc_plotplant.coords(origin = origin,
                                  xax.line = xax,
                                  yax.line = yax,
                                  plant.samps = samps,
                                  df = gps.df,
                                  reverse.xdist = T,
                                  reverse.yslope = T,
                                  reverse.ydist = F)
  coords$p
  coords$plant.coords %>%
    mutate(plot = NA) -> data
  
  
  #---------------------------------------------------------#
  # 4. Pull everything together to export and plot
  gps.df %>%
    rbind(data) -> coords
  coords
  coords %>%
    filter(type == "site.corner") %>%
    arrange(point.name) -> df.site.corners
  p <- ggplot(coords, aes(y = lat, x = lon, label = point.name, color = type)) +
    geom_text() +
    coord_fixed() +
    geom_polygon(data = df.site.corners, fill = NA)
  p
  
  coord.list <- list(coords = coords, p = p)
  
  return(coord.list)
  
}
#lco.mxt.gpscoords(gps.coords, plants)

########## WBI-NRT-NCS #########################################################
wbi.nrt.gpscoords <- function(gps.coords, plants){
  
  #---------------------------------------------------------#
  # 1. Identify relevant data
  curr.gps <- get.curr.gps(gps.coords, plants, site.name = "WBI-NRT-NCS")
  gps.df <- make.gps.long(curr.gps$curr.site)
  gps.df
  # id site corners and origins for each plot
  gps.df %>%
    mutate(plot = NA) %>%
    mutate(plot = ifelse(multiplot.notes == "Rep I", "b1", plot)) %>%
    mutate(plot = ifelse(multiplot.notes == "Rep II", "b2", plot)) %>%
    mutate(plot = ifelse(multiplot.notes == "Rep III", "b3", plot)) %>%
    filter(!is.na(lat)) %>%
    mutate(type = ifelse(is.na(plot), "site.corner", "plot.corner")) %>%
    mutate(point.name = ifelse(!is.na(plot), 
                               paste(plot, point.name, sep = "."), 
                               point.name)) %>%
    select(point.name, lat, lon, plot, type) -> gps.df
  
  #---------------------------------------------------------#
  # 2. Identify missing site corners 
  p1 <- ggplot(gps.df, aes(y = lat, x = lon, label = point.name)) +
    geom_text() +
    coord_fixed() 
  p1
  
  
  #---------------------------------------------------------#
  # 3. Find plant sample coords 
  
  # first, find overall x and y slopes
  xax.line <- calc_line.eq(x1 = gps.df[gps.df$point.name == "D","lon"],
                           y1 = gps.df[gps.df$point.name == "D","lat"],
                           x2 = gps.df[gps.df$point.name == "C","lon"],
                           y2 = gps.df[gps.df$point.name == "C","lat"])
  names(xax.line) <- c("slope","intercept")
  yax.line <- calc_line.eq(x1 = gps.df[gps.df$point.name == "D","lon"],
                           y1 = gps.df[gps.df$point.name == "D","lat"],
                           x2 = gps.df[gps.df$point.name == "A","lon"],
                           y2 = gps.df[gps.df$point.name == "A","lat"])
  names(yax.line) <- c("slope","intercept")
  p2 <- p1 +
    geom_abline(slope = xax.line$slope, intercept = xax.line$intercept) +
    geom_abline(slope = yax.line$slope, intercept = yax.line$intercept) 
  p2  
  
  
  ######
  # block1
  b1.xys <- calc_plot_xylines(df = gps.df, 
                              point.name.chr = "b1.D",
                              xax.slope = xax.line$slope,
                              yax.slope = yax.line$slope)
  curr.gps$curr.samps %>%
    filter(Samp %in% c(1,2,3)) -> samps.b1
  coords.b1 <- calc_plotplant.coords(origin = gps.df[gps.df$point.name == "b1.D",c("lon","lat")],
                                  xax.line = b1.xys$xax,
                                  yax.line = b1.xys$yax,
                                  plant.samps = samps.b1,
                                  df = gps.df,
                                  reverse.xdist = F,
                                  reverse.yslope = T,
                                  reverse.ydist = T)
  coords.b1$p
  coords.b1$p +
    xlim(c(-78.1023,-78.1020)) + ylim(c(34.767, 34.76725))
  coords.b1$plant.coords %>%
    mutate(plot = "b1") -> data.b1
  
  
  ######
  # block2
  b2.xys <- calc_plot_xylines(df = gps.df, 
                              point.name.chr = "b2.D",
                              xax.slope = xax.line$slope,
                              yax.slope = yax.line$slope)
  b2.xys
  curr.gps$curr.samps %>%
    filter(Samp %in% c(4,5,6)) -> samps.b2
  coords.b2 <- calc_plotplant.coords(origin = gps.df[gps.df$point.name == "b2.D",c("lon","lat")],
                                     xax.line = b2.xys$xax,
                                     yax.line = b2.xys$yax,
                                     plant.samps = samps.b2,
                                     df = gps.df,
                                     reverse.xdist = F,
                                     reverse.yslope = T,
                                     reverse.ydist = T)
  coords.b2$p
  coords.b2$p +
     xlim(c(-78.102,-78.1015)) + ylim(c(34.76725, 34.76775))
  coords.b2$plant.coords %>%
    mutate(plot = "b2") -> data.b2
  
  
  ######
  # block3
  b3.xys <- calc_plot_xylines(df = gps.df, 
                              point.name.chr = "b3.D",
                              xax.slope = xax.line$slope,
                              yax.slope = yax.line$slope)
  b3.xys
  curr.gps$curr.samps %>%
    filter(Samp %in% c(7,8)) -> samps.b3
  coords.b3 <- calc_plotplant.coords(origin = gps.df[gps.df$point.name == "b3.D",c("lon","lat")],
                                     xax.line = b3.xys$xax,
                                     yax.line = b3.xys$yax,
                                     plant.samps = samps.b3,
                                     df = gps.df,
                                     reverse.xdist = F,
                                     reverse.yslope = T,
                                     reverse.ydist = T)
  coords.b3$p
  # this doesn't look quite right since the plants have been plotted beyond the edge of the site boundary
  coords.b3$plant.coords %>%
    mutate(plot = "b3") -> data.b3
  
  
  #---------------------------------------------------------#
  # 4. Pull everything together to export and plot
  gps.df %>%
    rbind(data.b1) %>%
    rbind(data.b2) %>%
    rbind(data.b3) -> coords
  coords
  coords %>%
    filter(type == "site.corner") %>%
    arrange(point.name) -> df.site.corners
  p <- ggplot(coords, aes(y = lat, x = lon, label = point.name, color = type)) +
    geom_text() +
    coord_fixed() +
    geom_polygon(data = df.site.corners, fill = NA)
  p
  
  coord.list <- list(coords = coords, p = p)
  
  return(coord.list)
  
}
#wbi.nrt.gpscoords(gps.coords, plants)

########## BRF-ONE-COM #########################################################
brf.one.gpscoords <- function(gps.coords, plants){
  
  #---------------------------------------------------------#
  # 1. Identify relevant data
  curr.gps <- get.curr.gps(gps.coords, plants, site.name = "BRF-ONE-COM")
  gps.df <- make.gps.long(curr.gps$curr.site)
  gps.df
  # these are site corners
  gps.df %>%
    mutate(plot = NA) %>%
    filter(!is.na(lat)) %>%
    mutate(type = "site.corner") %>%
    select(point.name, lat, lon, plot, type) -> gps.df
  
  #---------------------------------------------------------#
  # 2. Identify missing site corners 
  p1 <- ggplot(gps.df, aes(y = lat, x = lon, label = point.name)) +
    geom_text() +
    coord_fixed() 
  p1
  
  # define xax.line and yax.line
  # y and x fall along gps grid, so find A by calc distance to A from D based on dist from C to B
  x.C <- gps.df[gps.df$point.name == "C","lon"]
  x.B <- gps.df[gps.df$point.name == "B","lon"]
  dist <- x.C - x.B
  x.D <- x.C
  x.A <- x.D - dist
  y.D <- gps.df[gps.df$point.name == "D","lat"]
  y.A <- y.D
  new.row <- data.frame(point.name = "A",
                        lat = y.A,
                        lon = x.A,
                        plot = NA,
                        type = "site.corner")
  gps.df <- rbind(gps.df, new.row)
  ggplot(gps.df, aes(x = lon, y = lat, label = point.name)) +
    geom_text() +
    coord_fixed()
  
  
  #---------------------------------------------------------#
  # 3. Find plant sample coords 
  curr.gps$curr.samps -> samps
  
  #calculate the gps coords for each plant sampled in the plot
  origin <- data.frame(lon = gps.df[gps.df$point.name == "D","lon"],
                       lat = gps.df[gps.df$point.name == "D","lat"])
  m_per_deg <- load_m.per.deg()
  X.dist <- samps$Xcoord.m/m_per_deg
  xax.lats <- origin$lat + X.dist
  Y.dist <- samps$Ycoord.m/m_per_deg
  yax.lons <- origin$lon - Y.dist
  data <- data.frame(point.name = paste0("p", samps$Samp),
                    lat = xax.lats, 
                    lon = yax.lons, 
                    plot = NA,
                    type = "plant")
  
  
  #---------------------------------------------------------#
  # 4. Pull everything together to export and plot
  gps.df %>%
    rbind(data) -> coords
  coords
  coords %>%
    filter(type == "site.corner") %>%
    arrange(point.name) -> df.site.corners
  p <- ggplot(coords, aes(y = lat, x = lon, label = point.name, color = type)) +
    geom_text() +
    coord_fixed() +
    geom_polygon(data = df.site.corners, fill = NA)
  p
  
  coord.list <- list(coords = coords, p = p)
  
  return(coord.list)
  
}
#brf.one.gpscoords(gps.coords, plants)

########## LWR-BHO-NCS #########################################################
lwr.bho.gpscoords <- function(gps.coords, plants){
  
  #---------------------------------------------------------#
  # 1. Identify relevant data
  curr.gps <- get.curr.gps(gps.coords, plants, site.name = "LWR-BHO-NCS")
  gps.df <- make.gps.long(curr.gps$curr.site)
  gps.df
  # these are site corners
  gps.df %>%
    mutate(plot = NA) %>%
    filter(!is.na(lat)) %>%
    mutate(type = "site.corner") %>%
    select(point.name, lat, lon, plot, type) -> gps.df
  
  
  #---------------------------------------------------------#
  # 2. Identify missing site corners 
  p1 <- ggplot(gps.df, aes(y = lat, x = lon, label = point.name)) +
    geom_text() +
    coord_fixed() 
  p1
  
  # define xax.line and yax.line
  xax.line <- calc_line.eq(x1 = gps.df[gps.df$point.name == "D","lon"],
                           y1 = gps.df[gps.df$point.name == "D","lat"],
                           x2 = gps.df[gps.df$point.name == "C","lon"],
                           y2 = gps.df[gps.df$point.name == "C","lat"])
  xax.line
  # yax.line <- calc_line.eq(x1 = gps.df[gps.df$point.name == "C","lon"],
  #                          y1 = gps.df[gps.df$point.name == "C","lat"],
  #                          x2 = gps.df[gps.df$point.name == "B","lon"],
  #                          y2 = gps.df[gps.df$point.name == "B","lat"])
  #yax.line # this doesn't match the shape in my notes. I'm going to force the yax slope to be perpendicular to the xax slope
  yax.slope <- -1/ xax.line$slope
  # find yax.line's intercept using the slope and point D 
  yax.intercept <- calc_intercept_fromSlopePoint(slope = yax.slope, 
                                                 x = gps.df[gps.df$point.name == "D","lon"], 
                                                 y = gps.df[gps.df$point.name == "D","lat"])
  # find A using yax line and distance from D
  m_per_deg <- load_m.per.deg()
  new.point <- useLineDist_findNewXY(dist = 27/m_per_deg, 
                                     slope = yax.slope, 
                                     intercept = yax.intercept, 
                                     x0 = gps.df[gps.df$point.name == "D","lon"], 
                                     y0 = gps.df[gps.df$point.name == "D","lat"])
  new.row <- data.frame(point.name = "A",
                        lat = new.point$y1,
                        lon = new.point$x1,
                        plot = NA,
                        type = "site.corner")
  
  gps.df <- rbind(gps.df, new.row)
  
  # update B using yax2 line and distance from C
  # first, find yax2 line intercept
  yax2.intercept <- calc_intercept_fromSlopePoint(slope = yax.slope, 
                                                  x = gps.df[gps.df$point.name == "C","lon"], 
                                                  y = gps.df[gps.df$point.name == "C","lat"])
  # second, find new B using yax2 line and distance from C
  new.point <- useLineDist_findNewXY(dist = 27/m_per_deg, 
                                     slope = yax.slope, 
                                     intercept = yax2.intercept, 
                                     x0 = gps.df[gps.df$point.name == "C","lon"], 
                                     y0 = gps.df[gps.df$point.name == "C","lat"])
  gps.df[gps.df$point.name == "B", c("lat","lon")] <- c(new.point$y1, new.point$x1)
  # ggplot(gps.df, aes(x = lon, y = lat, label = point.name)) +
  #   geom_text() +
  #   coord_fixed()
  

  #---------------------------------------------------------#
  # 3. Find plant sample coords 
  origin <- gps.df[gps.df$point.name == "D",]
  xax <- calc_line.eq(x1 = gps.df[gps.df$point.name == "D","lon"], 
                      y1 = gps.df[gps.df$point.name == "D","lat"],
                      x2 = gps.df[gps.df$point.name == "C","lon"],
                      y2 = gps.df[gps.df$point.name == "C","lat"])
  yax <- calc_line.eq(x1 = gps.df[gps.df$point.name == "D","lon"], 
                      y1 = gps.df[gps.df$point.name == "D","lat"],
                      x2 = gps.df[gps.df$point.name == "A","lon"],
                      y2 = gps.df[gps.df$point.name == "A","lat"])
  curr.gps$curr.samps -> samps
  coords <- calc_plotplant.coords(origin = origin,
                                  xax.line = xax,
                                  yax.line = yax,
                                  plant.samps = samps,
                                  df = gps.df,
                                  reverse.xdist = F,
                                  reverse.yslope = T,
                                  reverse.ydist = F)
  coords$p
  coords$plant.coords %>%
    mutate(plot = NA) -> data
  
  
  #---------------------------------------------------------#
  # 4. Pull everything together to export and plot
  gps.df %>%
    rbind(data) -> coords
  coords
  coords %>%
    filter(type == "site.corner") %>%
    arrange(point.name) -> df.site.corners
  p <- ggplot(coords, aes(y = lat, x = lon, label = point.name, color = type)) +
    geom_text() +
    coord_fixed() +
    geom_polygon(data = df.site.corners, fill = NA)
  p
  
  coord.list <- list(coords = coords, p = p)
  
  return(coord.list)
  
}
#lwr.bho.gpscoords(gps.coords, plants)




