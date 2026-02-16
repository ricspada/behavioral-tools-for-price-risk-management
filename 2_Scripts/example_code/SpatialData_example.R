# =========================================================================================
# install and load required packages
# =========================================================================================

## Install and load packages (not all are necessary, but all can be useful for further applications)
#packages <- c("sp","raster","rasterVis","maptools","mapview","maps","ncdf4","zoo","lubridate","rgdal","rgeos","plyr","dplyr","reshape2","biogeo","geosphere")
#
#package.check <- lapply(packages,FUN = function(x) {
#  if (!require(x, character.only = TRUE)) {
#    install.packages(x, dependencies = TRUE)
#    library(x, character.only = TRUE)
#  }
#}
#)

library(sf)
library(plyr)
library(sp)
library(raster)
library(ncdf4)
library(terra)
library(reshape2)

# =========================================================================================
# Abbreviations and explanation
# =========================================================================================

## variables (daily resolution)

# Rhires = precipitation in [mm]
# Tmax = maximum air temperature in [degree Celsius]
# Tmin = minimum air temperature in [degree Celsius]
# Tabs = mean air temperature in [degree Celsius]
# Srel = sunshine in [%] of daytime

## parameters to be defined below (that first has to be defined by author)

# start     = year of beginning in XXXX format (e.g. 2002)
# end       = year of ending in XXXX format (e.g. 2005)

# prepare extraction:
# all weather files must be in the following format: .nc (Network Common Data Form)
# all needed weather files (in .nc format) in a folder named "all_weatherdata"
# in all_weatherdata, gridded datasets for a specific variable must be in a subfolder (e.g. RhiresD/ for precipitation or TabsD/ for mean air temperature). D represents daily.
# know the PLZ-number of your municipalities of interest

# steps to do:
# - type in your information about the time serie (start and end) 
# - type in your information about the different variables you are interested in (format for ex.: c(RhiresD, TabsD))
# - set working directory as indicated below
# - load the files "PLZO_CSV_LV03.csv" and municipality shapefiles as indicated below
# - download all packages and load them
# - run the rest of the script (and make no changes unless corrections below)

# =========================================================================================
# Definition of variables
# =========================================================================================

# set start and end date of extracted time series in format YYYY e.g. 1990
start <- 2016 # format YYYY e.g. 1990
end <- 2021 # format YYYY e.g. 1990

# type in variables of interest
variables_of_interest <- c("TminD") # use the following format for daily values: c("RhiresD", "TmaxD", "TminD", "TabsD", "SrelD")

# set working directory to created folder "all_weatherdata".
setwd("W:/Riccardo/Extraction_weather_data/Code_Data/all_weatherdata")

# Enter your PLZ (zip-codes) of interest.
# Read PLZ (zip-codes) as a vector
loc <- read.csv("Ramseier_data_postcodes.csv")
PLZ.numbers <- unique(loc$pcode)

# Define the output directory
# Three outputs will be saved in the output directory
# i) list_weather_data.RData  : a list with structure [[weather variable]][PLZ,date]
# ii) all_weather_data.RData  : a data.frame with all weather variables (long form)
# iii) all_weather_data.csv   : a data.frame with all weather variables (long form)

output.dir <- "W:/Riccardo/Extraction_weather_data/Code_Data/results"

# =========================================================================================
# Load latest zip-code map
# =========================================================================================

# We will derive the coordinates of zip-code centers from this map
shapefile_PLZ <- read_sf(dsn = 'AMTOVZ_SHP_LV95/AMTOVZ_ZIP.shp')

# =========================================================================================
# run script and do not make any changes from here
# everything has been defined above - no need to adjust code below!
# =========================================================================================

# Visualization of selected municipalities in Switzerland
plot(shapefile_PLZ[which(shapefile_PLZ$ZIP4 %in% PLZ.numbers),], col="red")

# Get coordinates of each ZIP-code center 
PLZ.centroid <- st_centroid(shapefile_PLZ, byid = T)

# Get the municipalities of interest (dont joind the two data frames, only use PLZ.centroid)
PLZ.centroid <- PLZ.centroid[which(PLZ.centroid$ZIP4 %in% PLZ.numbers),]

# No duplicates (exclaves, etc.)
PLZ.centroid <- PLZ.centroid[!duplicated(PLZ.centroid$ZIP4),]

# Copying the centroid dataframe for later in the loop
PLZ.numbers.coordinates <- PLZ.centroid

# Convert it to spatial points
PLZ_spatial <- as(st_geometry(PLZ.centroid$geometry), "Spatial")

# Visual check on map
plot(PLZ_spatial, col = "blue")

#=========================================================================================
# Trying out file reading
#==========================================================================================
# Defining the crs
lv95_proj4 <- "+proj=somerc +lat_0=46.9524055555556 +lon_0=7.43958333333333 +k_0=1 +x_0=2600000 +y_0=1200000 +ellps=bessel +towgs84=674.374,15.056,405.346,0,0,0,0 +units=m +vunits=m +no_defs"

# Importing a weather data test file. Afterwards heck for errors regaring the crs
test <- raster('TminD/Test.nc', crs = lv95_proj4)
test

# Having a look at the meta data of the weather data file
test2 <- nc_open('TminD/Test.nc')

# Testing the extracting, check for errors
extracted_example <- extract(test, PLZ_spatial)

# =========================================================================================
# Extract weather data
# =========================================================================================

# Each element of the list contains the weather data
list_weather_data <- list() 

for (i in 1:length(variables_of_interest)) {
  
  # Get the weather data (which ones to download)
  if (variables_of_interest[i] == "RhiresD") {
    
    dir <- paste(variables_of_interest[i],"/", sep="")
    listoffiles <- list.files(path=dir ,pattern=".nc")
    
    years.data <- vector(length=length(listoffiles)) 
    for (f in 1:length(listoffiles)){
      years.data[f]<-substr(listoffiles[f],25,28)
    }
    
    files.to.import <- listoffiles[which(years.data==start):which(years.data==end)]
    
    # Download the datasets
    rawdata<- list()
    
    # coordinate reference system is for LV03
    for (j in 1:length(files.to.import)){
      rawdata[[j]]<-brick(paste(variables_of_interest[i],"/",files.to.import[j], sep=""))
      # Add coordinate reference system (LV03)
      crs(rawdata[[j]])<- lv95_proj4 
    }
    
    # Bring municipality center into the reference system of weather data
    locations <- PLZ_spatial
    #locations <- st_transform(PLZ_spatial, crs = crs(rawdata[[1]])) 
    
    # Get the weather data
    
    weatherdata <- extract(rawdata[[1]], locations)
    if(length(files.to.import) != 1){
      for (k in 2:length(files.to.import)){
        weatherdata.year <- extract(rawdata[[k]], locations)
        weatherdata <- cbind(weatherdata, weatherdata.year)
        rm(weatherdata.year)
      }
    } else {}
    
    years <- seq(start, end)
    days <- as.character(seq(as.Date(paste(years[1],1,1, sep = "-")), as.Date(paste(years[length(years)],12,31, sep = "-")), "day"))
    colnames(weatherdata) <- days
    rm(rawdata, dir, listoffiles, years.data, files.to.import)
    
    
  } else if (variables_of_interest[i] == "TmaxD"){
    dir <- paste(variables_of_interest[i],"/", sep="")
    listoffiles <- list.files(path=dir,pattern=".nc")
    
    years.data <- vector(length=length(listoffiles)) 
    for (f in 1:length(listoffiles)){
      years.data[f]<-substr(listoffiles[f],23,26)
    }
    
    files.to.import <- listoffiles[which(years.data==start):which(years.data==end)]
    
    # Download the datasets
    rawdata<- list()
    
    # coordinate reference system is for LV03
    for (j in 1:length(files.to.import)){
      rawdata[[j]]<-brick(paste(variables_of_interest[i],"/",files.to.import[j], sep=""))
      # Add coordinate reference system (LV03)
      crs(rawdata[[j]])<- lv95_proj4 
    }
    
    # Bring municipality center into the reference system of weather data
    locations <- PLZ_spatial
    #locations <- st_transform(PLZ_spatial, crs = crs(rawdata[[1]]))  
    
    # Get the weather data
    weatherdata <- extract(rawdata[[1]], locations)
    if(length(files.to.import) != 1){
      for (k in 2:length(files.to.import)){
        weatherdata.year <- extract(rawdata[[k]], locations)
        weatherdata <- cbind(weatherdata, weatherdata.year)
        rm(weatherdata.year)
      }
    } else {}
    
    years <- seq(start, end)
    days <- as.character(seq(as.Date(paste(years[1],1,1, sep = "-")), as.Date(paste(years[length(years)],12,31, sep = "-")), "day"))
    colnames(weatherdata) <- days
    rm(rawdata, dir, listoffiles, years.data, files.to.import)
    
    
  } else if (variables_of_interest[i] == "TminD") {
    dir <- paste(variables_of_interest[i],"/", sep="")
    listoffiles <- list.files(path=dir,pattern=".nc")
    
    years.data <- vector(length=length(listoffiles)) 
    for (f in 1:length(listoffiles)){
      years.data[f]<-substr(listoffiles[f],24,27)
    }
    
    files.to.import <- listoffiles[which(years.data==start):which(years.data==end)]
    
    # Download the datasets
    rawdata<- list()
    
    # coordinate reference system is for LV03
    for (j in 1:length(files.to.import)){
      rawdata[[j]]<-brick(paste(variables_of_interest[i],"/",files.to.import[j], sep=""))
      # Add coordinate reference system (LV03)
      crs(rawdata[[j]])<- lv95_proj4 
    }
    
    # Bring municipality center into the reference system of weather data
    locations <- PLZ_spatial
    #locations <- st_transform(PLZ_spatial, crs = crs(rawdata[[1]]))  
    
    # Get the weather data
    weatherdata <- extract(rawdata[[1]], locations)
    if(length(files.to.import) != 1){
      for (k in 2:length(files.to.import)){
        weatherdata.year <- extract(rawdata[[k]], locations)
        weatherdata <- cbind(weatherdata, weatherdata.year)
        rm(weatherdata.year)
      }
    } else {}
    
    years <- seq(start, end)
    days <- as.character(seq(as.Date(paste(years[1],1,1, sep = "-")), as.Date(paste(years[length(years)],12,31, sep = "-")), "day"))
    colnames(weatherdata) <- days
    rm(rawdata, dir, listoffiles, years.data, files.to.import)
    
  } else if (variables_of_interest[i] == "TabsD") {
    dir <- paste(variables_of_interest[i],"/", sep="")
    listoffiles <- list.files(path=dir,pattern=".nc")
    
    years.data <- vector(length=length(listoffiles)) 
    for (f in 1:length(listoffiles)){
      years.data[f]<-substr(listoffiles[f],23,26)
    }
    
    files.to.import <- listoffiles[which(years.data==start):which(years.data==end)]
    
    # Download the datasets
    rawdata<- list()
    
    # coordinate reference system is for LV03
    for (j in 1:length(files.to.import)){
      rawdata[[j]]<-brick(paste(variables_of_interest[i],"/",files.to.import[j], sep=""))
      # Add coordinate reference system (LV03)
      crs(rawdata[[j]])<- lv95_proj4 
    }
    
    # Bring municipality center into the reference system of weather data
    locations <- PLZ_spatial
    #locations <- st_transform(PLZ_spatial, crs = crs(rawdata[[1]])) 
    
    # Get the weather data
    weatherdata <- extract(rawdata[[1]], locations)
    if(length(files.to.import) != 1){
      for (k in 2:length(files.to.import)){
        weatherdata.year <- extract(rawdata[[k]], locations)
        weatherdata <- cbind(weatherdata, weatherdata.year)
        rm(weatherdata.year)
      }
    } else {}
    
    years <- seq(start, end)
    days <- as.character(seq(as.Date(paste(years[1],1,1, sep = "-")), as.Date(paste(years[length(years)],12,31, sep = "-")), "day"))
    colnames(weatherdata) <- days
    rm(rawdata, dir, listoffiles, years.data, files.to.import)
    
    
    # else it must be SrelD
  } else {
    dir <- paste(variables_of_interest[i],"/", sep="")
    listoffiles <- list.files(path=dir,pattern=".nc")
    
    years.data <- vector(length=length(listoffiles)) 
    for (f in 1:length(listoffiles)){
      years.data[f]<-substr(listoffiles[f],23,26)
    }
    
    files.to.import <- listoffiles[which(years.data==start):which(years.data==end)]
    
    # Download the datasets
    rawdata<- list()
    
    # coordinate reference system is for LV03
    for (j in 1:length(files.to.import)){
      rawdata[[j]]<-brick(paste(variables_of_interest[i],"/",files.to.import[j], sep=""))
      # Add coordinate reference system (LV03)
      crs(rawdata[[j]])<- lv95_proj4 
    }
    
    # Bring municipality center into the reference system of weather data
    locations <- PLZ_spatial
    #locations <- st_transform(PLZ_spatial, crs = crs(rawdata[[1]])) 
    
    # Get the weather data
    weatherdata <- extract(rawdata[[1]], locations)
    if(length(files.to.import) != 1){
      for (k in 2:length(files.to.import)){
        weatherdata.year <- extract(rawdata[[k]], locations)
        weatherdata <- cbind(weatherdata, weatherdata.year)
        rm(weatherdata.year)
      }
    } else {}
    
    years <- seq(start, end)
    days <- as.character(seq(as.Date(paste(years[1],1,1, sep = "-")), as.Date(paste(years[length(years)],12,31, sep = "-")), "day"))
    colnames(weatherdata) <- days
    rm(rawdata, dir, listoffiles, years.data, files.to.import)
  }
  
  # Save data from specific weather variable in list
  row.names(weatherdata) <- PLZ.numbers.coordinates$ZIP4
  list_weather_data[[i]] <- weatherdata
  rm(weatherdata)
  
}

# =========================================================================================
# Save weather data in the output directory
# ========================================================================================= 

# Transform each element of the list into long format (melt)
list_weather_data_melted <- lapply(list_weather_data, melt) 
for (i in 1:length(variables_of_interest)){
  colnames(list_weather_data_melted[[i]]) <- c("PLZ","date",variables_of_interest[i])
}  

# Merge all weather variables into one frame
weather_data_melted <- list_weather_data_melted[[1]]
if (length(variables_of_interest) !=1){
  for (i in 2:length(variables_of_interest)){
    weather_data_melted <- join(weather_data_melted,list_weather_data_melted[[i]], by=c("PLZ","date"))
  }
} else {}

# Save the weather data in the output directory  
save(list_weather_data, file=paste(output.dir,"/list_weather_data.RData",sep=""))
save(weather_data_melted, file=paste(output.dir,"/all_weather_data_melted.RData",sep=""))
write.csv(weather_data_melted, file=paste(output.dir,"/all_weather_data_melted.csv",sep=""))
