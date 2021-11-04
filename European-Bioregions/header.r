## ---------------------------
##
## Script name: header.r
##
## Purpose of script: This script contain functions and pre-processing code
## that are necessary to run the network-area code for the european bioregions.
##
## Authors: Dr Joao Braga
## Laboratory of Alpine Ecology
## Grenoble Alps University
## Grenoble, France
##
## and Dr Miguel Lurgi 
## Lecturer in Biosciences (Computational Ecology)
## Computational Ecology Lab - Department of Biosciences
## Swansea University, UK
## 
## Centre for Biodiversity Theory and Modelling
## Theoretical and Experimental Ecology Station, CNRS, France
##
## Date Created: November 2016
##
## Copyright (c) Joao Braga & Miguel Lurgi, 2016-2020
##
## ---------------------------
##
## Notes:
##
## This script is provided as supplementary material for the paper:
## Galiana, Lurgi, et al. (2020) The spatial scaling of ecological networks across the
## globe, Nature.
##
## The raw data used within this script was provided by 
## Dr Wilfried Thuiller
## Laboratory of Alpine Ecology
## Grenoble Alps University
## Grenoble, France
##
## ---------------------------

## This script reads a series of files from the home directory. Ensure these files are placed in the same
## directory from which this script is executed

# Libraries
source('../utils.r')
library(raster)

#################################################
# This section of the script was developed by Joao Braga
# Spatial networks 
# 10 KM
# 28/06/2016 Joao Braga
#################################################

# Metaweb with diet categories
load(file = "./raw-data/BARMdiet_binFUNDLINKS_50.RData")

#### This is the metaweb ######

#### Row and column names are species codes according to the species names dictionary
#### which is kept in 'SppID.txt'

head(BARMdiet.binary)   # Format: rows are predators ;  columns are preys VERY IMPORTANT!!!!

# the interaction matrix needs to be transposed to calculate matrics properly
BARMdiet.binary <- t(BARMdiet.binary)


#################################################
# Species codes and species names
# Function to Identify a spp by the code
whois <- function(SPPCODE = NULL, SPPNAME = NULL) {
  # Function to Identify a spp by the code
  
  if(is.null(SPPCODE) & is.null(SPPNAME)) stop("Must specify a species code or name(Genus_species)")
  if(!is.null(SPPCODE) & !is.null(SPPNAME)) stop("Must specify a species code or name(Genus_species)")
  
  SppID <- read.table(file = "SppID.txt", header = TRUE, stringsAsFactors = F)
  
  if(length(SPPCODE) > 1){
    SPPCODE <- paste0(SPPCODE, "$", collapse = "|")
  }
  if(length(SPPNAME) > 1){
    SPPNAME <- paste0(SPPNAME, "$", collapse = "|")
  }
  
  if(!is.null(SPPCODE))    who <- SppID[which(SppID$ID==SPPCODE),]$SPPname
  if(!is.null(SPPNAME))    who <- SppID[grepl(pattern = SPPNAME,x = SppID$SPPname),]
  
  return(who)
}



# Function to transform spp distribution (it can be network properties per pixel) into raster
fun.dbf2raster <- function(SPPPA, mask.dir = NULL){
  # SPPPA must be in this format - first colmun with CELL ID (PAGENALE) and second column with the value (to plot)
  if(is.null(mask.dir)) stop("Must specify the mask file directory")
  require(rgdal)
  require(raster)
  require(foreign)
  
  maskID <- read.dbf(list.files(path = mask.dir, full.names = TRUE, pattern = ".img.vat.dbf$"))
  maskk <- raster(x = list.files(path = mask.dir, full.names = TRUE, pattern = ".img$"))
  
  spp <- maskID
  spp$val <- NA
  spp$PageName <- as.character(spp$PageName)
  row.names(spp) <- spp$PageName
  
  SPPPA$PAGENAME <- as.character(SPPPA$PAGENAME)
  SPPPA[,2] <- as.numeric(as.character(SPPPA[,2]))
  row.names(SPPPA) <- SPPPA$PAGENAME
  
  cellID <- as.character(SPPPA$PAGENAME)
  if( nrow(spp[cellID,]) != nrow(SPPPA[cellID,])) stop("Cell IDs do not match")
  spp <- spp[cellID,] 
  spp$val <- SPPPA[,2]
  
  xx <- raster::values(maskk)
  
  if( length(xx[!is.na(xx)]) != nrow(spp)) stop("Mask size inadequate")
  xx[!is.na(xx)] <- spp$val    #[xx[!is.na(xx)]]
  
  raster::values(maskk) <- xx
  return(maskk)
}

whois(SPPCODE = "B122")
#      ID          SPPname
# 27 B122 Falco_rusticolus

whois(SPPNAME = "Falco")

#################################################
# Species presence/absence per pixel 
# e.g. at @10km
load(file = "./raw-data/MASTER.bin10000_allhab_tresh0.Rdata")
head(master) # rows are individual pixels; 1st columm (PAGENAME) is the unique name of the pixel, other columns are species

# Convert a dbf table to raster
# E.g. getting the spp "B122" distribution at @10k
sppDIST <- data.frame(PAGENAME = master$PAGENAME, SPP = master$M205) # 1st Column: PAGENAME, necessary to match the pixels; 2nd can be any value, in here is presence and absence, but it can be network properties
sppDIST$PAGENAME <- as.character(sppDIST$PAGENAME) # important: PAGENAME must be as character for this to work correctly
str(sppDIST)


sppRaster <- fun.dbf2raster(SPPPA = sppDIST, mask.dir = "./GIS-map/mask10k")
sppRaster

plot(sppRaster, main= whois(SPPCODE = "M205"))

#################################################
# Local webs (pixel level): sample of the metaweb for a given pixel species pool
pixel <- master[52000,]
head(pixel)
sum(pixel[-1])              # Note: first value is the pixel name
# 241 species in this pixel 

# Defining the species pool
sppSTRING <- names(pixel)[which( pixel==1, arr.ind=FALSE)]
head(sppSTRING)
length(sppSTRING)
# 241

# Sampling the metaweb
localWEB <- BARMdiet.binary[sppSTRING,sppSTRING] # Web adjacency matrix (predator x prey)
dim(localWEB)
# 241 241


############## Here ends the code provided by Joao Braga #####################

####### Here starts the data processing added by Miguel Lurgi  ##########
####### This is to see what the letters in the pagename mean

letters_in_codes <- c()

for(pn in sppDIST$PAGENAME){
  code <- gsub("[A-Z]", "", pn);
  
  if(! code %in% letters_in_codes){
    letters_in_codes <- append(letters_in_codes, code);
  }
}

letters_in_codes <- as.character(sort(as.numeric(letters_in_codes)))

for(pn in 1:(dim(sppDIST)[1])){
  pname <- sppDIST[pn,]$PAGENAME
  
  code <- gsub("[A-Z]", "", pname);
  
  sppDIST[pn,]$SPP <- which(letters_in_codes == code)
}

sppRaster <- fun.dbf2raster(SPPPA = sppDIST, mask.dir = "./GIS-map/mask10k/")
sppRaster

plot(sppRaster)


latitudes <- c()
longitudes <- c()

for(pn in master$PAGENAME){
  code_lat <- gsub("[0-9]", "", pn);
  code_lon <- gsub("[A-Z]", "", pn);
  
  if(! code_lat %in% latitudes)latitudes <- append(latitudes, code_lat);
  if(! code_lon %in% longitudes)longitudes <- append(longitudes, code_lon);
}

longitudes <- sort(as.numeric(longitudes))


###### the actual code for the runs starts here!!!!

#the extent of the raster
rows <- 589
cols <- 671

longitudes <- 1:cols
latitudes <- LETTERS

cur_lat <- length(latitudes)
done <- FALSE

for(f in LETTERS){
  for(s in LETTERS){
    cur_name <- paste(f,s,sep='')
    latitudes <- append(latitudes, cur_name)
    
    cur_lat <- cur_lat+1
    print(cur_lat)
    if(cur_lat == rows) done <- TRUE;
    
    if(done) break;
  }
  if(done) break;
}

####### These are species that are not on the interaction matrix (BARMdiet.binary) and hence not part of the food web
####### Some of this codes are also absent from the species names database
absent_species <- c("B511", "B512", "R250", "R242", "R215", "R216", "R89", "B164", "R249", "R237", "R238", "R71", "R129", "A64", "R119", "R82")



