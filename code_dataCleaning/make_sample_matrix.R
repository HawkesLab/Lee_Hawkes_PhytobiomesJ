# Make the sample matrix

# This script pulls in raw data from...
# ...
thisfile <- "make_sample_matrix.R"
# To create...
# data/data_intermediates/dataCleaningProducts/DOE-NC-FIELD_SampleData_Rcompiled.xlsx

# load libraries
library(tidyverse)
library(readxl); library(openxlsx) # openxlsx is needed to write the metadata file
library(gridExtra)
library(corrplot)

# load fxns
source('code/helpers.R')
sourceDir('code/')




# -------------------------#
#1. Read in the raw biogeochemistry data from this excel workbook and combine sheets: 
#'data/sites/DOE-NC-FIELD_BGCData_2019_01_19.xlsx'.
bgc.path <- 'data/bgc/DOE-NC-FIELD_BGCData_2020_01_19.xlsx'
soilData <- load_soildata(bgc.path)
note1 <- data.frame(Date = "1/19/2020",
                    Initials = "ML",
                    Description = "Combine bgc data into one matrix",
                    R.file = thisfile)
note2 <- data.frame(Date = "1/19/2010",
                    Initials = "ML",
                    Description = "Renamed variables for coding: 
                    'watercontent'<-'gravimetric.soil.water.content..g.water...g.dry.soil.',
                    'gdry.gfresh'<-'moisture.multiplier..dry.fresh.',
                    'nh4'='NH4..ug.g.',
                    'no3'='NO3..ug.g.',
                    'mbc'='MBC..fum.unfum.',
                    'doc'='UNFUM.C..ug.g.',
                    'p.resin'='P..ug.g.',
                    'perc.sand'='..Sand',
                    'perc.clay'='..Clay....8hrs',
                    'perc.silt'='..Silt.....8hrs',
                    'usda.class'='USDA.Soil.Classification',
                    'ph = mean_pH'",
                    R.file = thisfile)
note_new <- data.frame(Date = "2/17/2020",
                    Initials = "ML",
                    Description = "Added EATS soil C and N data",
                    R.file = thisfile)

# -------------------------#
#2. Read in plant data from this excel workbook:
#'data/sites/DOE-NC-FIELD_BGCData_2019_01_19.xlsx'.

site.path <- 'data/sites/DOE-NC-FIELD_SiteData_2020_04_22.xlsx'
plants <- load_plants(site.path)
plants %>%
  mutate(SiteSamp = paste(Site, Samp, sep = "-")) %>%
  mutate(Samp = as.character(Samp)) %>%
  select(Site, Samp, SiteSamp, max.height.m, max.basalwidth.m, max.basallength.m) -> plants.c
note3 <- data.frame(Date = "1/20/2020",
                    Initials = "ML",
                    Description = "Correct typo in plant gps coordinate. SFA-ONE-PRO plant 8 should have Xcoord of 148, not 14.8",
                    R.file = thisfile)

# -------------------------#
#3. Calculate plant-specific gps coordinates
gps.coords <- load_gpscoords(site.path)
coord.df <- compile_sitesamp_gps(gps.coords, plants)
coord.df %>%
  filter(type == "plant") %>%
  separate(point.name, into = c("drop","Samp"), sep = 1) %>%
  mutate(SiteSamp = paste(Site, Samp, sep = "-")) %>%
  select(-c(drop, type)) %>%
  dplyr::rename('samp.lat'='lat',
                'samp.lon'='lon',
                'samp.plot' = 'plot')-> coord.df.c

# check on cluster of plants that are about the same distance apart
# head(coord.df.c)
# coord.df.c %>%
#   dplyr::select(samp.lon, samp.lat) -> samp.gps
# row.names(samp.gps) <- coord.df.c$SiteSamp
# hav.dist <- geodist(samp.gps, 
#                     paired = TRUE, 
#                     sequential = FALSE, pad = FALSE,
#                     measure = "haversine")
# colnames(hav.dist) <- row.names(samp.gps)
# row.names(hav.dist) <- row.names(samp.gps)
# hav.dist.df <- extract_uniquePairDists(hav.dist)
# hav.dist.df %>%
#   dplyr::rename('hav.dist.m'='dist') -> hav.dist.df
# View(hav.dist.df)


note4 <- data.frame(Date = "1/20/2020",
                    Initials = "ML",
                    Description = "Renamed variables for coding: 'samp.lat'='lat',
                'samp.lon'='lon',
                'samp.plot'='plot'",
                    R.file = thisfile)

# -------------------------#
#4. Combine sample-level BGC data, plant size, and plant gps into one dataframe

soilData %>%
  left_join(plants.c) %>%
  left_join(coord.df.c) -> sampleData

note5 <- data.frame(Date = "1/20/2020",
                    Initials = "ML",
                    Description = "Combine sample-level soil biogeochemistry data with plant size and plant gps coordinates into one dataframe",
                    R.file = thisfile)


# -------------------------#
#5. Cleaning

# -------------------------#
#a. Remove one site ("MAF-ONE-PRO") with only 3 plants sampled
sampleData %>%
  filter(Site != "MAF-ONE-PRO") -> sampleData

note6 <- data.frame(Date = "1/13/2020",
                    Initials = "ML",
                    Description = "Remove Site = MAF-ONE-PRO. This site only included three samples from a mixed grass stand",
                    R.file = thisfile)

# -------------------------#
#3. Only keep columns that might be reported in manuscript or used in analyses
sampleData %>%
  select(SiteSamp, Site, Samp,
         SOM, W.V, BS., Ac, CEC, ph,
         watercontent,
         P, K, Na, Ca, Cu, Mg, Mn, S, Zn, 
         nh4, no3, TIN,
         perc.C, perc.N,
         mbc, doc,
         p.resin,
         perc.sand, perc.clay, perc.silt, usda.class,
         max.height.m, max.basalwidth.m, max.basallength.m,
         samp.lat, samp.lon, samp.plot) -> sampleData

note7 <- data.frame(Date = "1/19/2020",
                     Initials = "ML",
                     Description = "Only keep columns that might be reported in manuscript or used in analyses",
                     R.file = thisfile)


# -------------------------#
#4. Annotate metadata 

# transfer and update column descriptions
meta.cols.bgc <- read_excel(path = bgc.path, 
                        sheet = "Metadata",
                        range = anchored("B23", dim = c(162, 5)), col_names = T)

meta.cols.site <- read_excel(path = site.path, 
                        sheet = "Metadata",
                        range = anchored("B14", dim = c(65, 5)), col_names = T)

meta.cols <- rbind(meta.cols.bgc, meta.cols.site)


# pull column descriptions for those that match 
curr.meta.cols <- unique(meta.cols[meta.cols$Name %in% colnames(sampleData),c("Name","Description","Unit")])
curr.meta.cols

# update column info for those that don't match
colnames(sampleData)[!colnames(sampleData) %in% meta.cols$Name]
# update re-named cols
meta.cols %>%
  filter(Name %in% c("HM%","W/V","BS%",
                     "mean_pH",
                     "gravimetric soil water content (g water / g dry soil)", 
                     "NH4 (ug/g)","NO3 (ug/g)",
                     "MBC (fum-unfum)","UNFUM C (ug/g)", 
                     "P (ug/g)", 
                     "% Sand", "% Clay -- 8hrs", "% Silt  -- 8hrs", "USDA Soil Classification")) -> fix.meta.cols
fix.meta.cols$Name <- c("SOM","W.V","BS.",
                        "ph",
                        "watercontent",
                        "nh4","no3",
                        "mbc","doc",
                        "p.resin",
                        "perc.sand","perc.clay","perc.silt","usda.class")

curr.meta.cols <- rbind(curr.meta.cols, fix.meta.cols[,c("Name","Description","Unit")])
curr.meta.cols

# add calculated variables
colnames(sampleData)[!colnames(sampleData) %in% curr.meta.cols$Name]
newvars.meta.cols <- data.frame(Name = c("samp.lat","samp.lon","samp.plot"),
           Description = c("Latitude of plant sample calculated based on sample XY coordinates within plot and plot corner GPS coordinates. 
                           See code/estim_plantGPScoords.R and estim_plantGPScoords_bySite.R",
                           
                           "Longitude of plant sample calculated based on sample XY coordinates within plot and plot corner GPS coordinates. 
                           See code/estim_plantGPScoords.R and estim_plantGPScoords_bySite.R",
                           
                           "Identification of which plot the sample was taken if there was >1 plot per site"),
           Unit = c("degree decimal coordinate", "degree decimal coordinate","NA"))

curr.meta.cols <- rbind(curr.meta.cols, newvars.meta.cols)
dim(curr.meta.cols); dim(sampleData)

#- update the Notes meta info
# combine notes from this script
curr.meta.notes <- rbind(note1, note2, note3, note4, note5, note6, note7)

#- update the project and datasheet meta info
meta.project <- read_excel(path = bgc.path, 
                           sheet = "Metadata",
                           range = anchored("B3", dim = c(2, 3)), col_names = T)

curr.meta.sheet <- data.frame(Name = "DOE-NC-FIELD_SampleData_Rcompiled.csv",
                              Description = "Matrix that combines all environmental variables at the sample-level and 
                              includes additional information that may useful to report in the manuscript. 
                              Raw data come from DOE-NC-FIELD_BGCData_2020_01_19.xlsx and DOE-NC-FIELD_SiteData_2020_01_13.xlsx. 
                              Data were cleaned using the R script make_sample_matrix.R",
                              `Data type` = "Not raw",
                              `Hard copy location` = "NA",
                              `Version info` = "Raw data come from DOE-NC-FIELD_BGCData_2020_01_19.xlsx and DOE-NC-FIELD_SiteData_2020_01_13.xlsx. 
                              Data were cleaned using the R script make_sample_matrix.R",
                              `Previous versions` = "NA")

####

# -------------------------#
#5. Combine data and metadata and export to xlsx file in nearby location

## Create Workbook object and add worksheets
wb <- createWorkbook()
addWorksheet(wb, "data")
writeData(wb, "data", sampleData, startCol = 1, startRow = 1, rowNames = FALSE)

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
# Notes
writeData(wb, "metadata", "Notes", startCol = 1, startRow = 47, rowNames = FALSE)
writeData(wb, "metadata", curr.meta.notes, startCol = 2, startRow = 47, rowNames = FALSE)

saveWorkbook(wb, "data_intermediates/dataCleaningProducts/DOE-NC-FIELD_SampleData_Rcompiled.xlsx", overwrite = TRUE)



