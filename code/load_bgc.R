# Load BGC data

# File: DOE-NC-FIELD_BGCData_2019_01_19.xlsx
# Copied from folder: HawkesLab/Projects/DOE LLNL/Data/NC Data/Biogeochem
#bgc.path <- 'data/bgc/DOE-NC-FIELD_BGCData_2020_01_19.xlsx'

load_soildata <- function(bgc.path){
  
  require(readxl)
  require(tidyverse)
  
  excel_sheets(bgc.path)
  
  ### ncda.cs data
  ncdacs <- read_excel(path = bgc.path, sheet = "NCDA&CS", .name_repair = "universal")
  ncdacs %>%
    dplyr::rename('SOM'='HM.') %>%
    select(-c(Initials, Site, Samp, ncda.Sample.ID, ncdda.notes)) -> ncdacs
  dim(ncdacs)
  
  ### remove this pH from ncda.cs and replace with measurement in lab
  ncdacs %>%
    select(-pH) -> ncdacs
  
  ### ph
  ph <- read_excel(path = bgc.path, sheet = "pH", .name_repair = "universal")
  ph %>%
    select(SiteSamp, mean_pH) %>%
    dplyr::rename('ph'='mean_pH')-> ph
  
  ### moisture
  moisture <- read_excel(path = bgc.path, sheet = "Moisture ", .name_repair = "universal")
  moisture %>%
    dplyr::rename('watercontent'='gravimetric.soil.water.content..g.water...g.dry.soil.') %>%
    dplyr::rename('gdry.gfresh'='moisture.multiplier..dry.fresh.') %>%
    select(SiteSamp, watercontent, gdry.gfresh) -> moisture
  dim(moisture)
  
  ### KCl
  kcl <- read_excel(path = bgc.path, sheet = "KCl", .name_repair = "universal")
  kcl %>%
    dplyr::rename('nh4'='NH4..ug.g.') %>%
    dplyr::rename('no3'='NO3..ug.g.') %>%
    filter(Site != "Blank") %>%
    mutate(nh4 = as.numeric(nh4)) %>%
    mutate(no3 = as.numeric(no3)) %>%
    mutate(TIN = as.numeric(TIN)) %>%
    select(SiteSamp, nh4, no3, TIN) -> kcl
  dim(kcl)
  
  ### K2SO4
  k2so4 <- read_excel(path = bgc.path, sheet = "K2SO4 ", .name_repair = "universal")
  k2so4 %>%
    dplyr::rename('mbc'='MBC..fum.unfum.') %>%
    dplyr::rename('doc'='UNFUM.C..ug.g.') %>%
    filter(Site != "Blank") %>%
    mutate(mbc = as.numeric(mbc)) %>%
    mutate(doc = as.numeric(doc)) %>%
    select(SiteSamp, mbc, doc) -> k2so4
  dim(k2so4)
  
  ### P resin
  presin <- read_excel(path = bgc.path, sheet = "P resin", .name_repair = "universal")
  presin %>%
    dplyr::rename('p.resin'='P..ug.g.') %>%
    filter(Site != "Blank") %>%
    filter(Site != "NA") %>%
    mutate(p.resin = as.numeric(p.resin)) %>%
    select(SiteSamp, p.resin) -> presin
  dim(presin)
  
  ### texture
  texture <- read_excel(path = bgc.path, sheet = "Texture", .name_repair = "universal")
  texture %>%
    dplyr::rename('perc.sand'='..Sand') %>%
    dplyr::rename('perc.clay'='..Clay....8hrs') %>%
    dplyr::rename('perc.silt'='..Silt.....8hrs') %>%
    dplyr::rename('usda.class'='USDA.Soil.Classification') %>%
    select(SiteSamp, perc.sand, perc.clay, perc.silt, usda.class) %>%
    mutate(perc.sand = as.numeric(perc.sand)) %>%
    mutate(perc.clay = as.numeric(perc.clay)) %>%
    mutate(perc.silt = as.numeric(perc.silt)) %>%
    filter(usda.class != "NA") -> texture
  dim(texture)
  
  ### Soil CN
  soilcn <- read_excel(path = bgc.path, sheet = "EATS", .name_repair = "universal")
  soilcn %>%
    select(SiteSamp, perc.C, perc.N) -> soilcn
  
  
  # combine sheets by SiteSamp
  ncdacs %>%
    full_join(ph) %>%
    full_join(moisture) %>%
    full_join(kcl) %>%
    full_join(k2so4) %>%
    full_join(presin) %>%
    full_join(texture) %>%
    full_join(soilcn) -> soilData
  
  return(soilData)

}

