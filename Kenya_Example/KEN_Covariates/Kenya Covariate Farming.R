
library(terra)
library(sf)
library(exactextractr)
library(rnaturalearth)
library(dplyr)

ken <- ne_countries(country = "Kenya", returnclass = "sf")
country.abbrev <- 'ken'

# functions -----
pop_weighted_mean <- function(values, coverage_fraction) {
  pop <- values$pop
  cov <- values$cov
  wmean <- sum(pop * cov * coverage_fraction, na.rm=TRUE) / 
    sum(pop * coverage_fraction, na.rm=TRUE)
  return(wmean)
}

pop_weighted_var <- function(values, coverage_fraction) {
  pop <- values$pop
  cov <- values$cov
  wvar <- var(pop * cov * coverage_fraction, na.rm=TRUE)
  return(wvar)
}

normalize <- function(vector){
  return((vector-mean(vector))/sd(vector))
}


# load geometry ------
setwd('/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/shapeFiles_gadm')

poly.adm2 <- st_read(dsn = paste0('gadm41_',country.abbrev,'_shp'), layer = paste0("gadm41_",country.abbrev,"_2"), options = "ENCODING=UTF-8")

poly.adm2$admin2 <- 1:nrow(poly.adm2)
poly.adm2$admin2.char <- paste0('admin2_',1:nrow(poly.adm2))
n_admin2 <- nrow(poly.adm2)

poly.adm1 <- poly.adm2 %>% group_by(NAME_1) %>% summarise(geometry = st_union(geometry))
poly.adm1$admin1 <- 1:nrow(poly.adm1)
poly.adm1$admin1.char <- paste0('admin1_',1:nrow(poly.adm1))
n_admin1 <- nrow(poly.adm1)

admin.key <- merge(as.data.frame(poly.adm2[,c('NAME_1','NAME_2','admin2','admin2.char')]),
                   as.data.frame(poly.adm1[,c('NAME_1','admin1','admin1.char')]),by='NAME_1')%>%
  dplyr::select(-c(geometry.x,geometry.y))

cmat_admin2 <- poly.adm2
cmat_admin1 <- poly.adm1

# Size of area ---------

sf_proj <- st_transform(poly.adm2, 6933)   # Cylindrical Equal Area (meters)
sf_proj$area_m2 <- st_area(sf_proj)
cmat_admin2$area <- as.numeric(sf_proj$area_m2) / 1e6

sf_proj <- st_transform(poly.adm1, 6933)   # Cylindrical Equal Area (meters)
sf_proj$area_m2 <- st_area(sf_proj)
cmat_admin1$area <- as.numeric(sf_proj$area_m2) / 1e6

cmat_admin2$area_log <- log(cmat_admin2$area)
cmat_admin1$area_log <- log(cmat_admin1$area)

# population density ------

pop_dens <- rast("/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/Kenya/Population/ken_u5_2022_1km.tif")

cmat_admin2$density <- exact_extract(pop_dens, poly.adm2, 'sum')/cmat_admin2$area
cmat_admin2$pop_var <- exact_extract(pop_dens, poly.adm2, 'variance')

cmat_admin1$density <- exact_extract(pop_dens, poly.adm1, 'sum')/cmat_admin1$area
cmat_admin1$pop_var <- exact_extract(pop_dens, poly.adm1, 'variance')

cmat_admin2$density_log <- log(cmat_admin2$density)
cmat_admin2$pop_var_log <- log(cmat_admin2$pop_var)
cmat_admin1$density_log <- log(cmat_admin1$density)
cmat_admin1$pop_var_log <- log(cmat_admin1$pop_var)

# nighttime lights -----
## NOAA VNL
setwd("/Users/alanamcgovern/Desktop/Research/KEN_Covariates")

global_health <- rast("VNL_v22_npp-j01_2022_global_vcmslcfg_c202303062300.median_masked.dat.tif")   # path to your VIIRS/VNL GeoTIFF

global_health <- resample(global_health, pop_dens)
r_stack <- c(pop_dens, global_health)
names(r_stack) <- c("pop", "cov")

cmat_admin2$nt_lights <- exact_extract(r_stack, poly.adm2, pop_weighted_mean)
cmat_admin2$nt_lights_var <- exact_extract(r_stack, poly.adm2, pop_weighted_var)

cmat_admin1$nt_lights <- exact_extract(r_stack, poly.adm1, pop_weighted_mean)
cmat_admin1$nt_lights_var <- exact_extract(r_stack, poly.adm1, pop_weighted_var)

cmat_admin2$nt_lights_log <- log(0.01 + cmat_admin2$nt_lights)
cmat_admin1$nt_lights_log <- log(0.01 + cmat_admin1$nt_lights)

# time to healthcare -----
## Malaria atlas project
global_travel <- rast("202001_Global_Motorized_Travel_Time_to_Healthcare_KEN.tiff")

global_travel <- resample(global_travel, pop_dens)
r_stack <- c(pop_dens, global_travel)
names(r_stack) <- c("pop", "cov")

cmat_admin2$tthc <- exact_extract(r_stack, poly.adm2, pop_weighted_mean)
cmat_admin2$tthc_var <- exact_extract(r_stack, poly.adm2, pop_weighted_var)

cmat_admin1$tthc <- exact_extract(r_stack, poly.adm1, pop_weighted_mean)
cmat_admin1$tthc_var <- exact_extract(r_stack, poly.adm1, pop_weighted_var)

cmat_admin2$tthc_log <- log(cmat_admin2$tthc)
cmat_admin1$tthc_log <- log(cmat_admin1$tthc)

# Precipitation (CHIRPS annual 2018)-------------
prec <- rast("chirps-v2.0.2022.tif")
prec_ken <- crop(prec, vect(ken))
prec_ken <- mask(prec_ken, vect(ken))
prec_ann <- sum(prec_ken, na.rm=TRUE)  # total annual mm

prec_ann <- resample(prec_ann, pop_dens)
r_stack <- c(pop_dens, prec_ann)
names(r_stack) <- c("pop", "cov")

cmat_admin2$precip <- exact_extract(r_stack, poly.adm2, pop_weighted_mean)
cmat_admin2$precip_var <- exact_extract(r_stack, poly.adm2, pop_weighted_var)

cmat_admin1$precip <- exact_extract(r_stack, poly.adm1, pop_weighted_mean)
cmat_admin1$precip_var <- exact_extract(r_stack, poly.adm1, pop_weighted_var)

# Temperature (WorldClim)--------------
# files <- list.files("wc2.1_2.5m_tavg/", pattern = "wc2.1_2.5m_tavg_.*\\.tif$", full.names = TRUE)
# tmean_monthly <- rast(files)
# 
# # Days per month (for non-leap year)
# days <- c(31,28,31,30,31,30,31,31,30,31,30,31)
# 
# # Compute weighted mean
# tmean_weighted <- sum(tmean_monthly * days) / sum(days)
# 
# # Save
# writeRaster(tmean_weighted, "wc2.1_2.5m_tavg/mean_annual_weighted.tif", overwrite = TRUE)

temp <- rast('wc2.1_2.5m_tavg/mean_annual_weighted.tif')
temp <- resample(temp, pop_dens)
r_stack <- c(pop_dens, temp)
names(r_stack) <- c("pop", "cov")

cmat_admin2$temp <- exact_extract(r_stack, poly.adm2, pop_weighted_mean)
cmat_admin2$temp_var <- exact_extract(r_stack, poly.adm2, pop_weighted_var)

cmat_admin1$temp <- exact_extract(r_stack, poly.adm1, pop_weighted_mean)
cmat_admin1$temp_var <- exact_extract(r_stack, poly.adm1, pop_weighted_var)

# Elevation (SRTM via WorldClim) -----------
elev <- rast('wc2.1_2.5m_elev.tif')

elev <- resample(elev, pop_dens)
r_stack <- c(pop_dens, elev)
names(r_stack) <- c("pop", "cov")

cmat_admin2$elev <- exact_extract(r_stack, poly.adm2, pop_weighted_mean)
cmat_admin2$elev_var <- exact_extract(r_stack, poly.adm2, pop_weighted_var)

cmat_admin1$elev <- exact_extract(r_stack, poly.adm1, pop_weighted_mean)
cmat_admin1$elev_var <- exact_extract(r_stack, poly.adm1, pop_weighted_var)

# Take log of variance covariates -----------
cmat_admin1$nt_lights_var_log <- log(cmat_admin1$nt_lights_var)
cmat_admin1$tthc_var_log <- log(cmat_admin1$tthc_var)
cmat_admin1$precip_var_log <- log(cmat_admin1$precip_var)
cmat_admin1$temp_var_log <- log(cmat_admin1$temp_var)
cmat_admin1$elev_var_log <- log(cmat_admin1$elev_var)

cmat_admin2$nt_lights_var_log <- log(cmat_admin2$nt_lights_var+0.1)
cmat_admin2$tthc_var_log <- log(cmat_admin2$tthc_var)
cmat_admin2$precip_var_log <- log(cmat_admin2$precip_var)
cmat_admin2$temp_var_log <- log(cmat_admin2$temp_var)
cmat_admin2$elev_var_log <- log(cmat_admin2$elev_var)

# Select covariates and normalize --------------

cmat_admin1 <- cmat_admin1 %>% dplyr::select(NAME_1, admin1, admin1.char, density_log,nt_lights_log,tthc_log,precip,temp,elev,
                              pop_var_log, nt_lights_var_log, tthc_var_log, precip_var_log, temp_var_log, elev_var_log, area_log)

cmat_admin2 <- cmat_admin2 %>% dplyr::select(NAME_2, admin2, admin2.char, density_log,nt_lights_log,tthc_log,precip,temp,elev,
                                             pop_var_log, nt_lights_var_log, tthc_var_log, precip_var_log, temp_var_log, elev_var_log, area_log)

cmat_admin1$density_log <- normalize(cmat_admin1$density_log)
cmat_admin1$nt_lights_log <- normalize(cmat_admin1$nt_lights_log)
cmat_admin1$tthc_log <- normalize(cmat_admin1$tthc_log)
cmat_admin1$precip <- normalize(cmat_admin1$precip)
cmat_admin1$temp <- normalize(cmat_admin1$temp)
cmat_admin1$elev <- normalize(cmat_admin1$elev)

cmat_admin1$area_log <- normalize(cmat_admin1$area_log)
cmat_admin1$pop_var_log <- normalize(cmat_admin1$pop_var_log)
cmat_admin1$nt_lights_var_log <- normalize(cmat_admin1$nt_lights_var_log)
cmat_admin1$tthc_var_log <- normalize(cmat_admin1$tthc_var_log)
cmat_admin1$precip_var_log <- normalize(cmat_admin1$precip_var_log)
cmat_admin1$temp_var_log <- normalize(cmat_admin1$temp_var_log)
cmat_admin1$elev_var_log <- normalize(cmat_admin1$elev_var_log)

cmat_admin2$density_log <- normalize(cmat_admin2$density_log)
cmat_admin2$nt_lights_log <- normalize(cmat_admin2$nt_lights_log)
cmat_admin2$tthc_log <- normalize(cmat_admin2$tthc_log)
cmat_admin2$precip <- normalize(cmat_admin2$precip)
cmat_admin2$temp <- normalize(cmat_admin2$temp)
cmat_admin2$elev <- normalize(cmat_admin2$elev)

cmat_admin2$area_log <- normalize(cmat_admin2$area_log)
cmat_admin2$pop_var_log <- normalize(cmat_admin2$pop_var_log)
cmat_admin2$nt_lights_var_log <- normalize(cmat_admin2$nt_lights_var_log)
cmat_admin2$tthc_var_log <- normalize(cmat_admin2$tthc_var_log)
cmat_admin2$precip_var_log <- normalize(cmat_admin2$precip_var_log)
cmat_admin2$temp_var_log <- normalize(cmat_admin2$temp_var_log)
cmat_admin2$elev_var_log <- normalize(cmat_admin2$elev_var_log)

# Save covariates ------
cmat_admin2 <- as.data.frame(cmat_admin2)
cmat_admin1 <- as.data.frame(cmat_admin1)

setwd("/Users/alanamcgovern/Desktop/Research/KEN_Covariates")
save(cmat_admin2,file='Kenya_admin2_covariates.rda')
save(cmat_admin1,file='Kenya_admin1_covariates.rda')
