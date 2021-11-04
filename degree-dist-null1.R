
## ---------------------------
##
## Script name: degree-dists-analysis.r
##
## Purpose of script: Similar to data-analysis.r, this script takes output tables from the data 
## aggregation procedures of all the datasets (generated using the corresponding data_processing_agg.r scripts)
## and analyses them to investigate the pattern of shapes in degree distributions across data sets and across
## spatial scales
##
## Author: Dr Miguel Lurgi and Dr Nuria Galiana
## Lecturer in Biosciences (Computational Ecology)
## Computational Ecology Lab - Department of Biosciences
## Swansea University, UK
## 
## and
##
## Postdoc at Centre for Biodiversity Theory and Modelling
## Theoretical and Experimental Ecology Station, CNRS, France
##
## Date Created: 6-12-2020
##
## Copyright (c) Miguel Lurgi, 2020
## Email: miguel.lurgi@swansea.ac.uk; galiana.nuria@gmail.com
##
## ---------------------------
##
## Notes:
##
## This script is provided as supplementary material for the paper:
## Galiana, Lurgi, et al. (2021) The spatial scaling of ecological networks across the
## globe.
##
## ---------------------------

## This script reads a series of files from the home directory. Ensure these files are placed in the same
## directory from which this script is executed

## Load required libraries
require(ggplot2)
require(RColorBrewer)
library(tibble)
library(stringr)

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

###### For each dataset we need to find out what is the model that is consistently selected across areas.
cur_fit_galpar <- read.csv('./Salix/fits-degree-galpar-null-model-1.csv')
cur_fit_galpar <- add_column(cur_fit_galpar[,-1], region = 'Galpar', .after = 'area')
modified_fit_table <- NULL
for(a in sort(unique(cur_fit_galpar$area))){
  temp_data <- eval(parse(text=paste0('subset(cur_fit_galpar, area == ', a, ')')))
  ranking <- table(temp_data$model)
  cur_out <- cbind(temp_data[which(temp_data$model == names(which(ranking == max(ranking)))[1]),], fraction=max(ranking)/sum(ranking))
  
  if(is.null(modified_fit_table)){
    modified_fit_table <- cur_out
  }else{
    modified_fit_table <- rbind(modified_fit_table, cur_out)
  }
}

write.csv(modified_fit_table, './Salix/fits-degree-dists-galpar-fraction-null1.csv')

cur_fit_salgal <- read.csv('./Salix/fits-degree-salgal-null-model-1.csv')
cur_fit_salgal <- add_column(cur_fit_salgal[,-1], region = 'Salgal', .after = 'area')
modified_fit_table <- NULL
for(a in sort(unique(cur_fit_salgal$area))){
  temp_data <- eval(parse(text=paste0('subset(cur_fit_salgal, area == ', a, ')')))
  ranking <- table(temp_data$model)
  cur_out <- cbind(temp_data[which(temp_data$model == names(which(ranking == max(ranking)))[1]),], fraction=max(ranking)/sum(ranking))
  
  if(is.null(modified_fit_table)){
    modified_fit_table <- cur_out
  }else{
    modified_fit_table <- rbind(modified_fit_table, cur_out)
  }
}

write.csv(modified_fit_table, './Salix/fits-degree-dists-salgal-fraction-null1.csv')


cur_fit_pyrenees <- read.csv('./Pyrenees/fits-degree-dists-pyrenees-null-model-1.csv')
cur_fit_pyrenees <- add_column(cur_fit_pyrenees[,-1], region = 'pyrenees', .after = 'area')
modified_fit_table <- NULL
for(a in sort(unique(cur_fit_pyrenees$area))){
  temp_data <- eval(parse(text=paste0('subset(cur_fit_pyrenees, area == ', a, ')')))
  ranking <- table(temp_data$model)
  cur_out <- cbind(temp_data[which(temp_data$model == names(which(ranking == max(ranking)))[1]),], fraction=max(ranking)/sum(ranking))
  
  if(is.null(modified_fit_table)){
    modified_fit_table <- cur_out
  }else{
    modified_fit_table <- rbind(modified_fit_table, cur_out)
  }
}

write.csv(modified_fit_table, './Pyrenees/fits-degree-dists-pyrenees-fraction-null1.csv')



### For europe the code is a bit more complex because the distributions are in separate files per area

european_bioregions <- c('mediterranean', 'anatolian', 'steppic', 'atlantic', 'pannonian', 'arctic', 'boreal', 'alpine', 'continental', 'blacksea')

for(bior in european_bioregions){
  cur_fit_data <- read.csv('./European-Bioregions/fits-degree-dists-europe-null-model-1.csv')[-1]
  if(bior=='blacksea'){
    cur_fit_data <- cur_fit_data[cur_fit_data$region=='BlackSea',]
  }else{
    cur_fit_data <- cur_fit_data[cur_fit_data$region==str_to_title(bior),]
  }
  modified_fit_table <- NULL
  for(a in sort(unique(cur_fit_data$area))){
    temp_data <- eval(parse(text=paste0('subset(cur_fit_data, area == ', a, ')')))
    ranking <- table(temp_data$model)
    cur_out <- cbind(temp_data[which(temp_data$model == names(which(ranking == max(ranking)))[1]),], fraction=max(ranking)/sum(ranking))
    
    if(is.null(modified_fit_table)){
      modified_fit_table <- cur_out
    }else{
      modified_fit_table <- rbind(modified_fit_table, cur_out)
    }
  }
  write.csv(modified_fit_table, paste0('./European-Bioregions/',bior,'/fits-degree-dists-',bior,'-fraction-null1.csv'))
}

#### now let's do the replicates #'Bristol'
datasets <- c('Chaco', 'Quercus', 'Nahuel', 'Bristol', 'Gottin-HP', 'Gottin-PP', 'Sanak', 'Garraf-HP', 'Garraf-PP', 'Garraf-PP-2', 'Olot', 'Montseny')
for(dset in datasets){
  folder <- dset
  if(dset %in% c('Gottin-HP', 'Gottin-PP')){
    folder <- 'Gottin'
  }
  if(dset %in% c('Garraf-HP', 'Garraf-PP', 'Garraf-PP-2', 'Olot', 'Montseny')){
    folder <- 'Garraf-Montseny-Olot'
  }
  
  cur_fit_data <- read.csv(paste0('./',folder,'/fits-degree-dists-',tolower(dset),'-null-model-1.csv'))[-1]
  if(dset == 'Chaco'){
    cur_fit_data$area <- rep(1:12, 100)
  }
  modified_fit_table <- NULL
  for(a in sort(unique(cur_fit_data$area))){
    temp_data <- eval(parse(text=paste0('subset(cur_fit_data, area == ', a, ')')))
    ranking <- table(temp_data$model)
    cur_out <- cbind(temp_data[which(temp_data$model == names(which(ranking == max(ranking)))[1]),], fraction=max(ranking)/sum(ranking))
    
    if(is.null(modified_fit_table)){
      modified_fit_table <- cur_out
    }else{
      modified_fit_table <- rbind(modified_fit_table, cur_out)
    }
  }
  write.csv(modified_fit_table, paste0('./',folder,'/fits-degree-dists-', tolower(dset),'-fraction-null1.csv'))
}


#### for mulder's dataset this is a bit different because of the 7 different categories
cur_all_data <- read.csv(paste0('./Soils/fits-degree-dists-soils-null-model-1.csv'))[-1]
modified_fit_table <- NULL
for(type in unique(cur_all_data$s_type)){
  cur_fit_data <- subset(cur_all_data, s_type == type)
  for(a in sort(unique(cur_fit_data$area))){
    temp_data <- eval(parse(text=paste0('subset(cur_fit_data, area == ', a, ')')))
    ranking <- table(temp_data$model)
    cur_out <- cbind(temp_data[which(temp_data$model == names(which(ranking == max(ranking)))[1]),], fraction=max(ranking)/sum(ranking))
    
    if(is.null(modified_fit_table)){
      modified_fit_table <- cur_out
    }else{
      modified_fit_table <- rbind(modified_fit_table, cur_out)
    }
  }
}

write.csv(modified_fit_table, paste0('./Soils/fits-degree-dists-soils-fraction-null1.csv'))


#### After that pre-processing, we read the models by fractions tables generated above and 
#### produce Figures 3e & 3f from the paper

### REPLICATES ###
########### For replicates

cur_fits <- read.csv('./Soils/fits-degree-dists-soils-fraction-null1.csv')
cur_fits <- cur_fits[c('s_type','area','model','a','b')]
cur_fits$s_type <- paste0('soils-',cur_fits$s_type)

datasets <- c('fits-degree-dists-chaco-fraction-null1.csv', 'fits-degree-dists-gottin-pp-fraction-null1.csv',
              'fits-degree-dists-gottin-hp-fraction-null1.csv', 'fits-degree-dists-quercus-fraction-null1.csv',
              'fits-degree-dists-montseny-fraction-null1.csv','fits-degree-dists-garraf-pp-2-fraction-null1.csv', 
              'fits-degree-dists-garraf-pp-fraction-null1.csv', 'fits-degree-dists-nahuel-fraction-null1.csv',
              'fits-degree-dists-garraf-hp-fraction-null1.csv', 'fits-degree-dists-olot-fraction-null1.csv',
              'fits-degree-dists-sanak-fraction-null1.csv', 'fits-degree-dists-bristol-fraction-null1.csv')

for(d in datasets){
  name <- strsplit(strsplit(d, 'fits-degree-dists-')[[1]][2], '.csv')[[1]]
  name <- strsplit(strsplit(name, '-fraction-null1')[[1]][1], '.csv')[[1]]
  
  if(name %in% c('chaco', 'quercus', 'nahuel', 'sanak', 'bristol')){
    file_loc <- paste0('./',firstup(name),'/',d)
  }else if(name %in% c("gottin-pp", "gottin-hp")){
    file_loc <- paste0('./Gottin/',d)
  }else{
    file_loc <- paste0('./Garraf-Montseny-Olot/',d)
  }
  cur_fits <- rbind(cur_fits, cbind(s_type=name, read.csv(file_loc)[c('area','model','a','b')]))
}

norm <- aggregate(cur_fits$area, by = list(cur_fits$s_type), max)
norm1 <- aggregate(cur_fits$area, by = list(cur_fits$s_type), min)
cur_fits$area_norm <- (cur_fits$area-norm1[match(cur_fits$s_type, norm1$Group.1),]$x)/(norm[match(cur_fits$s_type, norm$Group.1),]$x-norm1[match(cur_fits$s_type, norm1$Group.1),]$x)
cur_fits$model <- relevel(cur_fits$model, 'mod1')

cur_fits[which(cur_fits$s_type == 'soils-1'),]$s_type <- 'Soil 1'
cur_fits[which(cur_fits$s_type == 'soils-2'),]$s_type <- 'Soil 2'
cur_fits[which(cur_fits$s_type == 'soils-3'),]$s_type <- 'Soil 3'
cur_fits[which(cur_fits$s_type == 'soils-4'),]$s_type <- 'Soil 4'
cur_fits[which(cur_fits$s_type == 'soils-5'),]$s_type <- 'Soil 5'
cur_fits[which(cur_fits$s_type == 'soils-6'),]$s_type <- 'Soil 6'
cur_fits[which(cur_fits$s_type == 'soils-7'),]$s_type <- 'Soil 7'

cur_fits[which(cur_fits$s_type == 'sanak'),]$s_type <- 'Sanak'
cur_fits[which(cur_fits$s_type == 'nahuel'),]$s_type <- 'Nahuel'
cur_fits[which(cur_fits$s_type == 'quercus'),]$s_type <- 'Quercus'
cur_fits[which(cur_fits$s_type == 'olot'),]$s_type <- 'Olot'
cur_fits[which(cur_fits$s_type == 'montseny'),]$s_type <- 'Montseny'
cur_fits[which(cur_fits$s_type == 'bristol'),]$s_type <- 'Bristol'
cur_fits[which(cur_fits$s_type == 'gottin-pp'),]$s_type <- 'Gottin PP'
cur_fits[which(cur_fits$s_type == 'gottin-hp'),]$s_type <- 'Gottin HP'
cur_fits[which(cur_fits$s_type == 'garraf-pp'),]$s_type <- 'Garraf PP'
cur_fits[which(cur_fits$s_type == 'garraf-pp-2'),]$s_type <- 'Garraf PP2'
cur_fits[which(cur_fits$s_type == 'garraf-hp'),]$s_type <- 'Garraf HP'
cur_fits[which(cur_fits$s_type == 'chaco'),]$s_type <- 'Chaco'

cur_fits$s_type <- factor(as.character(cur_fits$s_type), levels=sort(as.character(unique(cur_fits$s_type)), decreasing=TRUE))


#### Run the code below to generate Figure 3e in the paper
png(filename = "degree_distributions_points.png", res= 300, width = 2300, height = 1600)
ggplot(cur_fits , aes(area_norm, s_type, color = factor(model))) +
  geom_point(size=3.5, shape=15) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"), axis.text=element_text(size=18),
        axis.title=element_text(size= 22 ,face="bold"),
        legend.title = element_blank(),
        legend.box.background = element_rect(colour = "black", size = 1.5),
        #legend.title = element_text(size=12, face="bold"),
        #legend.text = element_text(size = 10, face = "bold"),
        panel.background = element_rect(fill = "white", color = 'black', size = 1.5)) +
  scale_colour_brewer(palette = 'Dark2', name = "Model", labels = c("Truncated PL", "Exponential", "Power Law", 'Lognormal')) +
  ylab('') + xlab('Normalised area') +
  guides(color=guide_legend(override.aes=list(fill=NA))) 
dev.off()



#### BIOGEOGRAPHICAL ###

## We repeat the same procedure above for te Biogeographical data
cur_fit_boreal <- read.csv('./European-Bioregions/boreal/fits-degree-dists-boreal-fraction-null1.csv')[-1]
cur_fit_boreal <- cur_fit_boreal[match(unique(cur_fit_boreal$area), cur_fit_boreal$area),]
cur_fit_continental <- read.csv('./European-Bioregions/continental/fits-degree-dists-continental-fraction-null1.csv')[-1]
cur_fit_continental <- cur_fit_continental[match(unique(cur_fit_continental$area), cur_fit_continental$area),]
cur_fit_pannonian <- read.csv('./European-Bioregions/pannonian/fits-degree-dists-pannonian-fraction-null1.csv')[-1]
cur_fit_pannonian <- cur_fit_pannonian[match(unique(cur_fit_pannonian$area), cur_fit_pannonian$area),]
cur_fit_atlantic <- read.csv('./European-Bioregions/atlantic/fits-degree-dists-atlantic-fraction-null1.csv')[-1]
cur_fit_atlantic <- cur_fit_atlantic[match(unique(cur_fit_atlantic$area), cur_fit_atlantic$area),]
cur_fit_steppic <- read.csv('./European-Bioregions/steppic/fits-degree-dists-steppic-fraction-null1.csv')[-1]
cur_fit_steppic <- cur_fit_steppic[match(unique(cur_fit_steppic$area), cur_fit_steppic$area),]
cur_fit_blacksea <- read.csv('./European-Bioregions/blacksea/fits-degree-dists-blacksea-fraction-null1.csv')[-1]
cur_fit_blacksea <- cur_fit_blacksea[match(unique(cur_fit_blacksea$area), cur_fit_blacksea$area),]
cur_fit_mediterranean <- read.csv('./European-Bioregions/mediterranean/fits-degree-dists-mediterranean-fraction-null1.csv')[-1]
cur_fit_mediterranean <- cur_fit_mediterranean[match(unique(cur_fit_mediterranean$area), cur_fit_mediterranean$area),]
cur_fit_arctic <- read.csv('./European-Bioregions/arctic/fits-degree-dists-arctic-fraction-null1.csv')[-1]
cur_fit_arctic <- cur_fit_arctic[match(unique(cur_fit_arctic$area), cur_fit_arctic$area),]
cur_fit_anatolian <- read.csv('./European-Bioregions/anatolian/fits-degree-dists-anatolian-fraction-null1.csv')[-1]
cur_fit_anatolian <- cur_fit_anatolian[match(unique(cur_fit_anatolian$area), cur_fit_anatolian$area),]
cur_fit_alpine <- read.csv('./European-Bioregions/alpine/fits-degree-dists-alpine-fraction-null1.csv')[-1]
cur_fit_alpine <- cur_fit_alpine[match(unique(cur_fit_alpine$area), cur_fit_alpine$area),]
cur_fit_salgal <- read.csv('./Salix/fits-degree-dists-salgal-fraction-null1.csv')[-c(1,2)]
cur_fit_salgal <- cur_fit_salgal[match(unique(cur_fit_salgal$area), cur_fit_salgal$area),]
cur_fit_galpar <- read.csv('./Salix/fits-degree-dists-galpar-fraction-null1.csv')[-c(1,2)]
cur_fit_galpar <- cur_fit_galpar[match(unique(cur_fit_galpar$area), cur_fit_galpar$area),]
cur_fit_pyrenees <- read.csv('./Pyrenees/fits-degree-dists-pyrenees-fraction-null1.csv')[-c(1,2)]
cur_fit_pyrenees <- cur_fit_pyrenees[match(unique(cur_fit_pyrenees$area), cur_fit_pyrenees$area),]

cur_fits_bio <- rbind(cur_fit_galpar, cur_fit_salgal, cur_fit_pyrenees, cur_fit_boreal, cur_fit_continental, cur_fit_pannonian, cur_fit_alpine, cur_fit_anatolian, cur_fit_arctic, cur_fit_steppic, cur_fit_mediterranean, cur_fit_blacksea)

norm <- aggregate(cur_fits_bio$area, by = list(cur_fits_bio$region), max)
norm1 <- aggregate(cur_fits_bio$area, by = list(cur_fits_bio$region), min)

cur_fits_bio$area_norm <- (cur_fits_bio$area-norm1[match(cur_fits_bio$region, norm1$Group.1),]$x)/(norm[match(cur_fits_bio$region, norm$Group.1),]$x-norm1[match(cur_fits_bio$region, norm1$Group.1),]$x)
cur_fits_bio$region <- firstup(as.character(cur_fits_bio$region))
cur_fits_bio$region <- factor(cur_fits_bio$region, levels=sort(unique(cur_fits_bio$region), decreasing=TRUE))


#### Run the code below to generate Figure 3f in the paper.
png(filename = "degree_distributions_bars-bioregions.png", res= 300, width = 2300, height = 1600)
ggplot(cur_fits_bio , aes(area_norm, region, color = factor(model))) +
  geom_point(size=3) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"), axis.text=element_text(size=18),
        axis.title=element_text(size= 22 ,face="bold"),
        legend.title = element_blank(),
        legend.box.background = element_rect(colour = "black", size = 1.5),
        #legend.title = element_text(size=12, face="bold"),
        #legend.text = element_text(size = 10, face = "bold"),
        panel.background = element_rect(fill = "white", color = 'black', size = 1.5)) +
  scale_colour_brewer(palette = 'Dark2', name = "Model", labels = c("Truncated PL", "Exponential", "Power Law", 'Lognormal')) +
  ylab('') + xlab('Normalised area') +
  guides(color=guide_legend(override.aes=list(fill=NA))) 
dev.off()


### Code to reproduce figure S4 in supplementary material ###

ggplot(cur_fits[cur_fits$model=='mod3'& cur_fits$s_type==c('Nahuel','Chaco','Montseny','Gottin PP','Gottin HP','Olot','Quercus','Bristol','Garraf PP','Garraf HP'),] , aes(area_norm, a, color = factor(s_type))) +
  geom_smooth() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"), axis.text=element_text(size=10),
        axis.title=element_text(size= 18 ,face="bold"),
        legend.title = element_blank(),
        legend.box.background = element_rect(colour = "black", size = 1.5),
        #legend.title = element_text(size=12, face="bold"),
        #legend.text = element_text(size = 10, face = "bold"),
        panel.background = element_rect(fill = "white", color = 'black', size = 1.5)) +
  #scale_colour_brewer(palette = 'Dark2', name = "Model", labels = c("Truncated PL", "Exponential", "Power Law", 'Lognormal')) +
  ylab('Exponent') + xlab('Normalised area') + 
  facet_wrap(~s_type)+
  guides(color=guide_legend(override.aes=list(fill=NA))) 

ggplot(cur_fits_bio[cur_fits_bio$model=='mod1',] , aes(area_norm, a, color = factor(region))) +
  geom_smooth() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"), axis.text=element_text(size=10),
        axis.title=element_text(size= 18 ,face="bold"),
        legend.title = element_blank(),
        legend.box.background = element_rect(colour = "black", size = 1.5),
        #legend.title = element_text(size=12, face="bold"),
        #legend.text = element_text(size = 10, face = "bold"),
        panel.background = element_rect(fill = "white", color = 'black', size = 1.5)) +
  #scale_colour_brewer(palette = 'Dark2', name = "Model", labels = c("Truncated PL", "Exponential", "Power Law", 'Lognormal')) +
  ylab('Exponent') + xlab('Normalised area') + 
  facet_wrap(~region)+
  guides(color=guide_legend(override.aes=list(fill=NA))) 



### NULL MODEL ANALYSES ###

########### For replicates

cur_fits_null1_null1 <- read.csv('./Soils/fits-degree-dists-soils-null-model.csv')
cur_fits_null1 <- cur_fits_null1[c('s_type','area','model','a','b')]
cur_fits_null1$s_type <- paste0('soils-',cur_fits_null1$s_type)

datasets <- c('fits-degree-dists-chaco-fraction.csv', 'fits-degree-dists-gottin-pp-fraction.csv',
              'fits-degree-dists-gottin-hp-fraction.csv', 'fits-degree-dists-quercus-fraction.csv',
              'fits-degree-dists-montseny-fraction.csv','fits-degree-dists-garraf-pp-2-fraction.csv', 
              'fits-degree-dists-garraf-pp-fraction.csv', 'fits-degree-dists-nahuel-fraction.csv',
              'fits-degree-dists-garraf-hp-fraction.csv', 'fits-degree-dists-olot-fraction.csv',
              'fits-degree-dists-sanak-fraction.csv', 'fits-degree-dists-bristol-fraction.csv')

for(d in datasets){
  name <- strsplit(strsplit(d, 'fits-degree-dists-')[[1]][2], '.csv')[[1]]
  name <- strsplit(strsplit(name, '-fraction')[[1]][1], '.csv')[[1]]
  
  if(name %in% c('chaco', 'quercus', 'nahuel', 'sanak', 'bristol')){
    file_loc <- paste0('./',firstup(name),'/',d)
  }else if(name %in% c("gottin-pp", "gottin-hp")){
    file_loc <- paste0('./Gottin/',d)
  }else{
    file_loc <- paste0('./Garraf-Montseny-Olot/',d)
  }
  cur_fits_null1 <- rbind(cur_fits_null1, cbind(s_type=name, read.csv(file_loc)[c('area','model','a','b')]))
}

norm <- aggregate(cur_fits_null1$area, by = list(cur_fits_null1$s_type), max)
norm1 <- aggregate(cur_fits_null1$area, by = list(cur_fits_null1$s_type), min)
cur_fits_null1$area_norm <- (cur_fits_null1$area-norm1[match(cur_fits_null1$s_type, norm1$Group.1),]$x)/(norm[match(cur_fits_null1$s_type, norm$Group.1),]$x-norm1[match(cur_fits_null1$s_type, norm1$Group.1),]$x)
cur_fits_null1$model <- relevel(cur_fits_null1$model, 'mod1')

cur_fits_null1[which(cur_fits_null1$s_type == 'soils-1'),]$s_type <- 'Soil 1'
cur_fits_null1[which(cur_fits_null1$s_type == 'soils-2'),]$s_type <- 'Soil 2'
cur_fits_null1[which(cur_fits_null1$s_type == 'soils-3'),]$s_type <- 'Soil 3'
cur_fits_null1[which(cur_fits_null1$s_type == 'soils-4'),]$s_type <- 'Soil 4'
cur_fits_null1[which(cur_fits_null1$s_type == 'soils-5'),]$s_type <- 'Soil 5'
cur_fits_null1[which(cur_fits_null1$s_type == 'soils-6'),]$s_type <- 'Soil 6'
cur_fits_null1[which(cur_fits_null1$s_type == 'soils-7'),]$s_type <- 'Soil 7'

cur_fits_null1[which(cur_fits_null1$s_type == 'sanak'),]$s_type <- 'Sanak'
cur_fits_null1[which(cur_fits_null1$s_type == 'nahuel'),]$s_type <- 'Nahuel'
cur_fits_null1[which(cur_fits_null1$s_type == 'quercus'),]$s_type <- 'Quercus'
cur_fits_null1[which(cur_fits_null1$s_type == 'olot'),]$s_type <- 'Olot'
cur_fits_null1[which(cur_fits_null1$s_type == 'montseny'),]$s_type <- 'Montseny'
cur_fits_null1[which(cur_fits_null1$s_type == 'bristol'),]$s_type <- 'Bristol'
cur_fits_null1[which(cur_fits_null1$s_type == 'gottin-pp'),]$s_type <- 'Gottin PP'
cur_fits_null1[which(cur_fits_null1$s_type == 'gottin-hp'),]$s_type <- 'Gottin HP'
cur_fits_null1[which(cur_fits_null1$s_type == 'garraf-pp'),]$s_type <- 'Garraf PP'
cur_fits_null1[which(cur_fits_null1$s_type == 'garraf-pp-2'),]$s_type <- 'Garraf PP2'
cur_fits_null1[which(cur_fits_null1$s_type == 'garraf-hp'),]$s_type <- 'Garraf HP'
cur_fits_null1[which(cur_fits_null1$s_type == 'chaco'),]$s_type <- 'Chaco'

cur_fits_null1$s_type <- factor(as.character(cur_fits_null1$s_type), levels=sort(as.character(unique(cur_fits_null1$s_type)), decreasing=TRUE))


#### Run the code below to generate Figure 3e in the paper
png(filename = "degree_distributions_points.png", res= 300, width = 2300, height = 1600)
ggplot(cur_fits_null1 , aes(area_norm, s_type, color = factor(model))) +
  geom_point(size=.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"), axis.text=element_text(size=18),
        axis.title=element_text(size= 22 ,face="bold"),
        legend.title = element_blank(),
        legend.box.background = element_rect(colour = "black", size = 1.5),
        #legend.title = element_text(size=12, face="bold"),
        #legend.text = element_text(size = 10, face = "bold"),
        panel.background = element_rect(fill = "white", color = 'black', size = 1.5)) +
  scale_colour_brewer(palette = 'Dark2', name = "Model", labels = c("Truncated PL", "Exponential", "Power Law", 'Lognormal')) +
  ylab('') + xlab('Normalised area') +
  guides(color=guide_legend(override.aes=list(fill=NA))) 
dev.off()