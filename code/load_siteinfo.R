# load Site Info

# File: DOE-NC-FIELD_SiteInfo_2019_06_19.xlsx
# File: downloaded_climateData_2019_01_21.xlsx
# Copied from folder: HawkesLab/Projects/DOE LLNL/Data/NC Data/Site Info

#site.path <- 'data/sites/DOE-NC-FIELD_SiteData_2019_06_19.xlsx'

load_sites <- function(site.path){
  
  # read in the excel file
  excel_sheets(site.path)
  sites <- read_excel(path = site.path, sheet = "Sites", .name_repair = "universal")
  
  # check for typos
  unique(sites$Site)
  unique(sites$Short.Site)
  unique(sites$mono.mixed)
  unique(sites$stand.age.yrs)
  unique(sites$yr.estab)
  unique(sites$num.cultivars)
  unique(sites$cultivar)
  unique(sites$other.veg)
  unique(sites$pasture.yn)
  unique(sites$harvest.yn)
  unique(sites$mow.burn.yn)
  unique(sites$fert.yn)
  unique(sites$County)
  
  return(sites)
}

load_gpscoords <- function(site.path){
  
  # read in the excel file
  excel_sheets(site.path)
  gps.coords <- read_excel(path = site.path, sheet = "GPS.coords", .name_repair = "universal")
  
  # check for typos
  sites <- load_sites(site.path)
  unique(gps.coords$Site)[!unique(gps.coords$Site) %in% unique(sites$Site)] # need to change UPC to UCP
  gps.coords[gps.coords$Site == "UPC-MXG-NCD","Site"] <- "UCP-MXG-NCD"
  
  # switch cols to numeric
  gps.coords %>%
    mutate(length.m = as.numeric(length.m)) %>%
    mutate(width.m = as.numeric(width.m)) %>%
    mutate(gpsA.N = as.numeric(gpsA.N)) %>%
    mutate(gpsA.W = as.numeric(gpsA.W)) %>%
    mutate(gpsB.N = as.numeric(gpsB.N)) %>%
    mutate(gpsB.W = as.numeric(gpsB.W)) %>%
    mutate(gpsC.N = as.numeric(gpsC.N)) %>%
    mutate(gpsC.W = as.numeric(gpsC.W)) %>%
    mutate(gpsD.N = as.numeric(gpsD.N)) %>%
    mutate(gpsD.W = as.numeric(gpsD.W)) -> gps.coords
  
  return(gps.coords)
  
}

load_plants <- function(site.path){
  
  # read in the excel file
  excel_sheets(site.path)
  plants <- read_excel(path = site.path, sheet = "Plants", .name_repair = "universal")
  
  # check for typos
  sites <- load_sites(site.path)
  unique(plants$Site)[!unique(plants$Site) %in% unique(sites$Site)]
  unique(plants$Samp)
  
  # switch cols to numeric
  plants %>%
    mutate(Xcoord.m = as.numeric(Xcoord.m)) %>%
    mutate(Ycoord.m = as.numeric(Ycoord.m)) %>%
    mutate(gps.N = as.numeric(gps.N)) %>%
    mutate(gps.W = as.numeric(gps.W)) -> plants
  
  # update typo in plant plot coordinate -- SFA-ONE-PRO plant 8 should have Xcoord of 148, not 14.8
  plants[plants$Site == "SFA-ONE-PRO" & plants$Samp == "8","Xcoord.m"] <- 148
  
  return(plants)
  
}

calc_plotSize <- function(gps.coords){
  
  require(tidyverse)
  
  # reshape
  gps.coords %>%
    select(Site, 
           multiplot, multiplot.notes, 
           gpsA.N, gpsA.W,
           gpsB.N, gpsB.W,
           gpsC.N, gpsC.W,
           gpsD.N, gpsD.W,
           length.m, width.m) -> data
  
  # just deal with the sites that are 1 plot
  data %>%
    filter(multiplot == "n") %>%
    mutate(plotarea.m2 = length.m * width.m) -> data.mn
  
  # why is there missing data for one of these?
  data.mn %>%
    filter(is.na(plotarea.m2))
  # this is MAF-ONE-PRO, which isn't a fully-sampled site anyway
  
  # deal with sites with more than 1 plot
  data %>%
    filter(multiplot == "y") %>%
    select(Site, multiplot.notes, length.m, width.m) %>%
    filter(!is.na(length.m)) %>% # this gets rid of the row for WBI-NRT-NCS that was the outer experimental matrix dimensions for reference
    mutate(plotarea.m2 = length.m * width.m) -> data.my
  
  # take the average plot size for site with more than 1 plot
  data.my %>%
    group_by(Site) %>%
    summarize(mean = mean(plotarea.m2),
              n = length(plotarea.m2),
              se = sd(plotarea.m2)/sqrt(n)) -> summ.data.my
  
  # put plot size data back together
  data.mn %>%
    select(Site, plotarea.m2) %>%
    filter(!is.na(plotarea.m2)) %>% # gets rid of MAF, since it wasn't a fully-sampled plot
    mutate(numberOfplots = 1) %>%
    mutate(plotarea.m2.se = NA) -> tmp
  tmp1 <- data.frame(summ.data.my)
  tmp1 %>%
    dplyr::rename('plotarea.m2'='mean',
           'numberOfplots'='n',
           'plotarea.m2.se'='se') -> tmp1
  data.plotsize <- rbind(tmp, tmp1)

  return(data.plotsize)
  
}

calc_site_centroids <- function(gps.coords){
  
  require(tidyverse)
  require(geosphere)
  #source(file = 'code/helpers.R')
  
  # reshape
  num.rows <- dim(gps.coords)[1]
  gps.coords %>%
    select(Site, gpsA.N, gpsA.W,
           gpsB.N, gpsB.W,
           gpsC.N, gpsC.W,
           gpsD.N, gpsD.W) %>%
    mutate(subplotid = seq(1:num.rows)) %>% # add a unique number that identifies each plot
    gather(key = "corner.latlong", value = "coord", -c(Site, subplotid)) %>%
    separate(corner.latlong, into=c("corner","latlong")) %>%
    filter(!is.na(coord)) %>%
    spread(key = latlong, value = coord) %>%
    dplyr::rename(lat = "N",
           long = "W") -> data.gps
  
  #for each site, find the centroid
  SITE <- unique(data.gps$Site)
  centroid.list <- list()
  for(i in 1:length(SITE)){
    
    data.gps %>%
      filter(Site == SITE[i]) %>%
      select(long,lat) -> tmp
    tmp <- unique(tmp)
    tmp <- tmp[complete.cases(tmp),]
    
    if(dim(tmp)[1] > 2){
      cent <- centroid(tmp)
      centroid.list[[i]] <- data.frame(lon = cent[1], lat = cent[2])
    }
    
    if(dim(tmp)[1] == 2){
      centroid.list[[i]] <- data.frame(lon = sapply(tmp, mean)[1], lat = sapply(tmp, mean)[2], row.names = NULL)
    }
    
    if(dim(tmp)[1] == 1){
      centroid.list[[i]] <- data.frame(lon = NA, lat = NA)
    }
    
  }
  names(centroid.list) <- SITE
  centroid.df <- list_to_df(centroid.list)
  centroid.df %>%
    dplyr::rename('Site'='source') -> centroid.df
  df <- data.frame(centroid.df, row.names = NULL)
  df %>%
    select(Site, lat, lon) -> df
  
  return(df)

}




