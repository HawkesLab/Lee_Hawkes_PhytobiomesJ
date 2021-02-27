# Make the site matrix

# This script pulls in raw data from...
# data/sites/DOE-NC-FIELD_SiteData_2020_01_13.xlsx
# data/downloaded_climateData_2019_01_21.xlsx

thisfile <- "make_site_matrix.R"
# To create...
# data/data_intermediates/dataCleaningProducts/DOE-NC-FIELD_SiteData_Rcompiled.xlsx


# load libraries
library(tidyverse)
library(readxl); library(openxlsx) # openxlsx is needed to write the metadata file
library(gridExtra)
library(corrplot)
library(lubridate)

# load fxns
source('code/helpers.R')
sourceDir('code/')




# -------------------------#
#1. Read in the raw site-level data from this excel workbook: 'data/sites/DOE-NC-FIELD_SiteData_2019_06_19.xlsx'.
site.path <- 'data/sites/DOE-NC-FIELD_SiteData_2020_04_22.xlsx'
sites <- load_sites(site.path)

#sites$Short.Site
meta.notes
note4 <- data.frame(Date = "2020-01-13",
                    Initials = "ML",
                    Description = "Add short codes for site that were used to keep track of cultured fungal isolates",
                    R.file = thisfile,
                    stringsAsFactors = F)

# -------------------------#
#2. Load sheet = "GPS.coords"
gps.coords <- load_gpscoords(site.path)

# -------------------------#
#3. Use data in "GPS.coords" to do the following calculations:

# -------------------------#
#a. Calculate the size of the plot(s) that make up the site
site.plotsize <- calc_plotSize(gps.coords)

note1 <- data.frame(Date = "2020-01-13",
                    Initials = "ML",
                    Description = "Summarize the number of plots per site. Calculate the average plot area in m^2 per site using plot gps coordinates and plot length and width.",
                    R.file = thisfile,
                    stringsAsFactors = F)

def1 <- data.frame(Name = "numberOfplots",
                   Description = "Number of plots per site where plot is a contiguous stand of switchgrass", 
                   Unit = "NA")

def2 <- data.frame(Name = "plotarea.m2",
                   Description = "Average plot area in m2 per site", 
                   Unit = "m^2")

def3 <- data.frame(Name = "plotarea.m2.se",
                   Description = "Standard error plot area per site", 
                   Unit = "m^2")

# -------------------------#
#b. Calculate the site's centroid
site.centroids <- calc_site_centroids(gps.coords)

note2 <- data.frame(Date = "2020-01-13",
                    Initials = "ML",
                    Description = "Calculate the centroid (gps lat and lon) using site gps coordinates.",
                    R.file = thisfile,
                    stringsAsFactors = F)

def4 <- data.frame(Name = "lat",
                   Description = "Latitude", 
                   Unit = "decimal form")

def5 <- data.frame(Name = "lon",
                   Description = "Longitude", 
                   Unit = "decimal form")

# -------------------------#
#c. Load climate data based on site centroids. Data from excel workbook "downloaded_climateData_2019_01_21.xlsx", sheet = "Data"
clim.path <- 'data/sites/downloaded_climateData_2019_01_21.xlsx'
clim.data <- read_excel(path = clim.path, sheet = "Data", .name_repair = "universal")
clim.data %>%
  select(site, MAP.mm, MALT.C, MAT.C, MAHT.C) %>%
  dplyr::rename('Site'='site') -> clim.data

note3 <- data.frame(Date = "2020-01-13",
                    Initials = "ML",
                    Description = "Incorporate mean annual climate data from XXXXX for each site based on site centroids",
                    R.file = thisfile,
                    stringsAsFactors = F)

def6 <- data.frame(Name = "MAP.mm",
                   Description = "Mean annual precipitation in mm", 
                   Unit = "mm")

def7 <- data.frame(Name = "MALT.C",
                   Description = "Mean annual low temperature in degrees C", 
                   Unit = "C")

def8 <- data.frame(Name = "MAT.C",
                   Description = "Mean annual temperature in degrees C", 
                   Unit = "C")

def9 <- data.frame(Name = "MAHT.C",
                   Description = "Mean annual high temperature in degrees C", 
                   Unit = "C")



# -------------------------#
#4. Combine site-level information into one dataframe
sites %>%
  left_join(site.plotsize) %>%
  left_join(site.centroids) %>%
  left_join(clim.data) -> siteData

note5 <- data.frame(Date = "2020-01-13",
                    Initials = "ML",
                    Description = "Combine site-level data into one matrix",
                    R.file = thisfile,
                    stringsAsFactors = F)

# -------------------------#
#5. Cleaning

# -------------------------#
#a. Remove one site ("MAF-ONE-PRO") with only 3 plants sampled
siteData %>%
  filter(Site != "MAF-ONE-PRO") -> siteData

note6 <- data.frame(Date = "2020-01-13",
                    Initials = "ML",
                    Description = "Remove Site = MAF-ONE-PRO. This site only included three samples from a mixed grass stand",
                    R.file = thisfile,
                    stringsAsFactors = F)

# -------------------------#
#b. Cateorize a borderline site ("LCO-MXT-COM") into Ecoregion == Middle Atlantic Coastal Plain
selection <- siteData$Site == "LCO-MXT-COM"
siteData[selection, "Ecoregion"] <- "Middle Atlantic Coastal Plain"
siteData$Ecoregion <- factor(siteData$Ecoregion, levels = c("Blue Ridge","Piedmont",
                                      "Southeastern Plains","Middle Atlantic Coastal Plain"))

note7 <- data.frame(Date = "2020-01-13",
                    Initials = "ML",
                    Description = "Cateorize a borderline site LCO-MXT-COM into Ecoregion = Middle Atlantic Coastal Plain",
                    R.file = thisfile,
                    stringsAsFactors = F)

# -------------------------#
#c. Create a categorical variable for stand age in years 
siteData %>%
  mutate(stand.age.yrs.num = as.numeric(stand.age.yrs)) %>%
  mutate(stand.age.yrs.cat = ifelse(stand.age.yrs.num < 5, "lessThan5yrs", "5-10yrs")) %>%
  mutate(stand.age.yrs.cat = ifelse(stand.age.yrs.num > 10, "greaterThan10yrs", stand.age.yrs.cat)) -> siteData
selection <- siteData$stand.age.yrs == ">10"
siteData[selection, "stand.age.yrs.cat"] <- "greaterThan10yrs"
siteData$stand.age.yrs.cat <- factor(siteData$stand.age.yrs.cat, 
                                     levels = c("lessThan5yrs","5-10yrs","greaterThan10yrs"))

note8 <- data.frame(Date = "2020-01-13",
                    Initials = "ML",
                    Description = "Create a categorical variable for stand age with three levels: lessThan5yrs, 5-10yrs, and greatherThan10yrs",
                    R.file = thisfile,
                    stringsAsFactors = F)

def10 <- data.frame(Name = "stand.age.yrs.cat",
                   Description = "Categorical variable for stand age with three levels: lessThan5yrs, 5-10yrs, and greatherThan10yrs", 
                   Unit = "NA")

# -------------------------#
#d. Create a categorical variable for average plot area. Log-transform the value plotarea.m2 and divide into quantiles.
log.quants <- quantile(log(siteData$plotarea.m2))
exp(log.quants)
siteData %>%
  mutate(log.plotarea.m2 = log(plotarea.m2)) %>%
  mutate(plotarea.cat = ifelse(log.plotarea.m2 > log.quants[1], "small","very_small")) %>%
  mutate(plotarea.cat = ifelse(log.plotarea.m2 > log.quants[2], "medium", plotarea.cat)) %>%
  mutate(plotarea.cat = ifelse(log.plotarea.m2 > log.quants[3], "big", plotarea.cat)) %>%
  mutate(plotarea.cat = ifelse(log.plotarea.m2 > log.quants[4], "very_big", plotarea.cat)) -> siteData
siteData$plotarea.cat <- factor(siteData$plotarea.cat, 
                                levels = c("very_small","small","medium","big","very_big"))

note9 <- data.frame(Date = "2020-01-13",
                    Initials = "ML",
                    Description = "Create a categorical variable for average plot area per site with 5 levels based on log-transformed quantiles: very_small (0-25m2), small (25-504m2), medium (504-1487m2), big (1487-3417m2), and very_big (3417-33957m2)",
                    R.file = thisfile,
                    stringsAsFactors = F)

def11 <- data.frame(Name = "plotarea.cat",
                   Description = "Categorical variable for average plot area per site with 5 levels: very_small (0-25m2), small (25-504m2), medium (504-1487m2), big (1487-3417m2), and very_big (3417-33957m2)", 
                   Unit = "NA")

# -------------------------#
#e. Simplify two categorical variables (harvest and mow.burn) into one (harvest.mow.burn). In other words, just keep track of whether the site experienced some sort of biomass removal management -- harvesting, mowing, or burning.  This is coded as "yes" or "no"
# simplify harvest and mow.burn into one variable
siteData$harvest.mow.burn.yn <- c("unknown","yes","no","yes","yes","yes","no","no",
                                  "unknown","unknown","unknown","yes","yes","yes")

note10 <- data.frame(Date = "2020-01-13",
                    Initials = "ML",
                    Description = "Simplify two categorical variables (harvest and mow.burn) into one (harvest.mow.burn). In other words, just keep track of whether the site experienced some sort of biomass removal management -- harvesting, mowing, or burning",
                    R.file = thisfile,
                    stringsAsFactors = F)

def12 <- data.frame(Name = "harvest.mow.burn.yn",
                   Description = "Is the site regularly harvested, mowed, or burned? Yes/no", 
                   Unit = "NA")

# r-format cleaning...
# set up the levels for mono.mixed
siteData$mono.mixed <- factor(siteData$mono.mixed, levels = c("mono","mixed-grass","mixed-tree"))

# -------------------------#
#6. Only keep columns that might be reported in manuscript or used in analyses
siteData %>%
  select(Site, Site.name, Short.Site,
         sampling.day, sampling.month, sampling.year,
         Ecoregion,
         mono.mixed, 
         stand.age.yrs, stand.age.yrs.num, stand.age.yrs.cat,
         num.cultivars, cultivar, other.veg,
         pasture.yn, harvest.mow.burn.yn, fert.yn, mow.burn.notes, fert.notes,
         numberOfplots, 
         plotarea.m2, plotarea.m2.se, plotarea.cat,
         lat, lon,
         MAP.mm, MALT.C, MAT.C, MAHT.C,
         Site.address, County, Land.owner,
         Site.access.contact, Site.access.email, Site.access.phone) -> siteData

note11 <- data.frame(Date = "2020-01-13",
                    Initials = "ML",
                    Description = "Only keep columns that might be reported in manuscript or used in analyses",
                    R.file = thisfile,
                    stringsAsFactors = F)


# -------------------------#
#7. Annotate metadata from raw "DOE-NC-FIELD_SiteData_XXXX_XX_XX.xlsx" with new information from this script to create metadata in a "readme" for "DOE-NC-FIELD_SiteData_Rcompiled.csv"

#- transfer and update column descriptions
meta.cols <- read_excel(path = site.path, 
           sheet = "Metadata",
           range = anchored("B14", dim = c(65, 5)), col_names = T)

# pull column descriptions for those that match 
curr.meta.cols <- meta.cols[meta.cols$Name %in% colnames(siteData),c("Name","Description","Unit")]

# add new column info
colnames(siteData)[!colnames(siteData) %in% meta.cols$Name]
curr.meta.cols <- rbind(curr.meta.cols, def1, def2, def3, def4, def5, def6, def7, def8, def9, def10, def11, def12)

#- update the Notes meta info
meta.notes <- read_excel(path = site.path, 
           sheet = "Metadata",
           range = anchored("B80", dim = c(6, 3)), col_names = T)
meta.notes$R.file <- NA

# combine notes from this script
curr.meta.notes <- rbind(note1, note2, note3, note4, 
                         note5, note6, note7, note8, note9, 
                         note10, note11)
curr.meta.notes <- rbind(meta.notes, curr.meta.notes)
curr.meta.notes$Date <- as.character(as_date(curr.meta.notes$Date))

#- update the project and datasheet meta info
meta.project <- read_excel(path = site.path, 
           sheet = "Metadata",
           range = anchored("B3", dim = c(2, 3)), col_names = T)

meta.sheets <- read_excel(path = site.path, 
           sheet = "Metadata",
           range = anchored("B6", dim = c(7, 6)), col_names = T)

curr.meta.sheet <- data.frame(Name = "DOE-NC-FIELD_SiteData_Rcompiled.csv",
           Description = "Matrix that combines all environmental variables at the site-level and includes additional information about sites that will be useful to report in the manuscript. Raw data come from DOE-NC-FIELD_SiteData_2020_01_13.xlsx and downloaded_climateData_2019_01_21.xlsx. Data were combined using the R script make_site_matrix.R",
           `Data type` = "Not raw",
           `Hard copy location` = "NA",
           `Version info` = "Raw data come from DOE-NC-FIELD_SiteData_2020_01_13.xlsx and downloaded_climateData_2019_01_21.xlsx. Data were combined using the R script make_site_matrix.R",
           `Previous versions` = "NA")


# -------------------------#
#8. Combine data and metadata and export to xlsx file in nearby location

## Create Workbook object and add worksheets
wb <- createWorkbook()
addWorksheet(wb, "data")
writeData(wb, "data", siteData, startCol = 1, startRow = 1, rowNames = FALSE)

addWorksheet(wb, "metadata")
writeData(wb, "metadata", "Metadata", startCol = 1, startRow = 1, rowNames = FALSE)
# Project
writeData(wb, "metadata", "Project", startCol = 1, startRow = 3, rowNames = FALSE)
writeData(wb, "metadata", meta.project, startCol = 2, startRow = 3, rowNames = FALSE)
# Sheets
writeData(wb, "metadata", "Sheets", startCol = 1, startRow = 6, rowNames = FALSE)
writeData(wb, "metadata", curr.meta.sheet, startCol = 2, startRow = 6, rowNames = FALSE)
# Columns
writeData(wb, "metadata", "Columns", startCol = 1, startRow = 9, rowNames = FALSE)
writeData(wb, "metadata", curr.meta.cols, startCol = 2, startRow = 9, rowNames = FALSE)

9 + dim(curr.meta.cols)[1] + 1
# Notes
writeData(wb, "metadata", "Notes", startCol = 1, startRow = 45, rowNames = FALSE)
writeData(wb, "metadata", curr.meta.notes, startCol = 2, startRow = 45, rowNames = FALSE)

saveWorkbook(wb, "data_intermediates/dataCleaningProducts/DOE-NC-FIELD_SiteData_Rcompiled.xlsx", overwrite = TRUE)






