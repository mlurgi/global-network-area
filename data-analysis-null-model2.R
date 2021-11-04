## ---------------------------
##
## Script name: data-analysis.r
##
## Purpose of script: This script takes the output tables from the null model 2 analyses for all datasets
## and analyses them in tandem. To be able to run this script all the 'data_processing_agg.r' scripts 
## embedded within each dataset folder must have been executed and the outputs thus obtained collated 
## in the same folder
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
require(dplyr)
require(purrr)
require(sars)

## Here we define a function used to extract normalised versions of the network properties
range_stats <- function(df){
  data.frame(areas = df$areas,
             species = (log10(df$species) - min(log10(df$species))),
             links = (log10(df$links) - min(log10(df$links))),
             links_per_sp = (log10(df$links_per_sp) - min(log10(df$links_per_sp))),
             indegree = (log10(df$indegree) - min(log10(df$indegree))),
             potential_indegree = (log10(df$potential_indegree) - max(log10(df$potential_indegree)))
  )
}

### LOADING DATA ###
### First we load the null model 2 output data from each dataset
### DATA FOR THE REPLICATES
output_garraf_hp_null2 <- read.table("./Garraf-Montseny-Olot/output-garraf-hp-null-model-2.csv", sep=",", header=TRUE)[-1]
output_garraf_pp_null2 <- read.table("./Garraf-Montseny-Olot/output-garraf-pp-null-model-2.csv", sep=",", header=TRUE)[-1]
output_garraf_pp2_null2 <- read.table("./Garraf-Montseny-Olot/output-garraf-pp-2-null-model-2.csv", sep=",", header=TRUE)[-1]
output_montseny_null2 <- read.table("./Garraf-Montseny-Olot/output-montseny-null-model-2.csv", sep=",", header=TRUE)[-1]
output_olot_null2 <- read.table("./Garraf-Montseny-Olot/output-olot-null-model-2.csv", sep=",", header=TRUE)[-1]

output_nahuel_null2 <- read.table("./Nahuel/output-nahuel-null-model-2.csv", sep=",", header=TRUE)[-1]

output_quercus_null2 <- read.table("./Quercus/output-quercus-null-model-2.csv", sep=",", header=TRUE)[-1]
colnames(output_quercus_null2)[4] <- "resources"
colnames(output_quercus_null2)[5] <- "consumers"

output_soils_null2 <- read.table("./Soils/output-soils-null-model-2.csv", sep=",", header=TRUE)[-1]
output_soils_null2$resources <- output_soils_null2$S_basal
output_soils_null2$consumers <- output_soils_null2$S_intermediate + output_soils_null2$S_top

output_soils1_null2 <- output_soils_null2[(output_soils_null2$site_type==1),]
output_soils2_null2 <- output_soils_null2[(output_soils_null2$site_type==2),]
output_soils3_null2 <- output_soils_null2[(output_soils_null2$site_type==3),]
output_soils4_null2 <- output_soils_null2[(output_soils_null2$site_type==4),]
output_soils5_null2 <- output_soils_null2[(output_soils_null2$site_type==5),]
output_soils6_null2 <- output_soils_null2[(output_soils_null2$site_type==6),]
output_soils7_null2 <- output_soils_null2[(output_soils_null2$site_type==7),]

output_gottin_hp_null2 <- read.table("./Gottin/output-gottin-hp-null-model-2.csv", sep=",", header=TRUE)[-1]
colnames(output_gottin_hp_null2)[5] <- "resources"
colnames(output_gottin_hp_null2)[6] <- "consumers"

output_gottin_pp_null2 <- read.table("./Gottin/output-gottin-pp-null-model-2.csv", sep=",", header=TRUE)[-1]
colnames(output_gottin_pp_null2)[5] <- "resources"
colnames(output_gottin_pp_null2)[6] <- "consumers"

output_chaco_null2 <- read.table("./Chaco/output-chaco-null-model-2.csv", sep=",", header=TRUE)[-1]
names(output_chaco_null2)[3] <- "areas"

output_sanak_null2 <- read.table("./Sanak/output-sanak-null-model-2.csv", sep=",", header=TRUE)[-1]

output_bristol_null2 <- read.table("./Bristol/output-bristol-null-model-2.csv", sep=",", header=TRUE)[-1]
names(output_bristol_null2)[2] <- "areas"

## Once we have read all datasets we standardise their format and merge them into a single data frame
all_datasets <- cbind(data.frame(dataset='Garraf HP'), output_garraf_hp_null2)
all_datasets <- rbind(all_datasets, cbind(data.frame(dataset='Garraf PP'), output_garraf_pp_null2))
all_datasets <- rbind(all_datasets, cbind(data.frame(dataset='Garraf PP2'), output_garraf_pp2_null2))
all_datasets <- rbind(all_datasets, cbind(data.frame(dataset='Montseny'), output_montseny_null2))
all_datasets <- rbind(all_datasets, cbind(data.frame(dataset='Olot'), output_olot_null2))

names(all_datasets)[5:6] <- c('resources', 'consumers')

all_datasets <- all_datasets[c('dataset', 'replicate', 'areas', 'species', 'links', 'links_per_sp', 'indegree', 'potential_indegree', 'resources', 'consumers')]

all_datasets <- rbind(all_datasets, cbind(data.frame(dataset='Nahuel'), output_nahuel_null2[c('replicate', 'areas', 'species', 'links', 'links_per_sp', 'indegree', 'potential_indegree', 'resources', 'consumers')]))
all_datasets <- rbind(all_datasets, cbind(data.frame(dataset='Quercus'), output_quercus_null2[c('replicate', 'areas', 'species', 'links', 'links_per_sp', 'indegree', 'potential_indegree', 'resources', 'consumers')]))

output_soils_null2$dataset <- paste0('Soil ', output_soils_null2$site_type)
all_datasets <- rbind(all_datasets, output_soils_null2[c('dataset', 'replicate', 'areas', 'species', 'links', 'links_per_sp', 'indegree', 'potential_indegree', 'resources', 'consumers')])

all_datasets <- rbind(all_datasets, cbind(data.frame(dataset='Gottin HP'), output_gottin_hp_null2[c('replicate', 'areas', 'species', 'links', 'links_per_sp', 'indegree', 'potential_indegree', 'resources', 'consumers')]))
all_datasets <- rbind(all_datasets, cbind(data.frame(dataset='Gottin PP'), output_gottin_pp_null2[c('replicate', 'areas', 'species', 'links', 'links_per_sp', 'indegree', 'potential_indegree', 'resources', 'consumers')]))

all_datasets <- rbind(all_datasets, cbind(data.frame(dataset='Chaco'), output_chaco_null2[c('replicate', 'areas', 'species', 'links', 'links_per_sp', 'indegree', 'potential_indegree', 'resources', 'consumers')]))
all_datasets <- rbind(all_datasets, cbind(data.frame(dataset='Bristol'), output_bristol_null2[c('replicate', 'areas', 'species', 'links', 'links_per_sp', 'indegree', 'potential_indegree', 'resources', 'consumers')]))
all_datasets <- rbind(all_datasets, cbind(data.frame(dataset='Sanak'), output_sanak_null2[c('replicate', 'areas', 'species', 'links', 'links_per_sp', 'indegree', 'potential_indegree', 'resources', 'consumers')]))

all_datasets$cr_ratio <- all_datasets$consumers/all_datasets$resources

#### In case you want to store the pre-processed data to have it all ready for the next time this script is run.
write.csv(all_datasets, file='merged_datasets_replicates_nullmodel2.csv')
#all_datasets_null2 <- read.table("./merged_datasets_replicates_nullmodel2.csv", sep=",", header=TRUE)

## Here we summarise and normalise the data
all_summaries_replicates <- all_datasets %>% 
  group_by(dataset, areas) %>% 
  dplyr::summarize(
    species = mean(species),
    links = mean(links),
    links_per_sp = mean(links_per_sp),
    indegree = mean(indegree),
    potential_indegree = mean(potential_indegree),
    cr_ratio = mean(cr_ratio)
  ) 

all_summaries_replicates_norm  <- all_summaries_replicates %>%
  split(.$dataset) %>%
  map_dfr(range_stats, .id = "dataset")

### DATA FROM BIOGEOGRAPHICAL DATASETS ### 
output_galpar <- read.table("./Salix/output-galpar-null-model-2.csv", sep=",", header=TRUE)[-1]
colnames(output_galpar)[4] <- "resources"
colnames(output_galpar)[5] <- "consumers"

output_salgal <- read.csv("./Salix/output-salgal-null-model-2.csv", sep=",", header=TRUE)[-1]
colnames(output_salgal)[4] <- "resources"
colnames(output_salgal)[5] <- "consumers"

output_pyrenees <- read.csv("./Pyrenees/output-pyrenees-null-model-2.csv", header=TRUE)[-1]
names(output_pyrenees)[2] <- "areas"


### Due to running constraints in terms of computational capacities we had to run different European
### Bioregions separately, so, there is a different output file for each.
### Replace the code below to read the output file 'output-european-bioregions.csv' once
### it has been generated accordingly

### Accordingly, the size of the outputs generated for each european bioregion is considerably large.
### So, an archive containing it would be too large to send over. For convenience hence, we provide
### one replicate of each of the european bioregion as a sample of the actual complete dataset.
### Those files are read in the following lines.

europe_datasets <- read.table("./European-Bioregions/output-europe-null-model-2.csv", sep=",", header=TRUE, stringsAsFactors = F)[-1]

names(europe_datasets)[c(1,3,7)] <- c('dataset', 'areas', 'indegree')
europe_datasets[is.na(europe_datasets)] <- 0

europe_datasets <- europe_datasets[c('dataset', 'replicate', 'areas', 'species', 'links', 'links_per_sp', 'indegree', 'resources', 'consumers')]

salgal <- cbind(data.frame(dataset='galpar'), output_galpar[c('rep', 'areas', 'species', 'links', 'links_per_sp', 'indegree', 'resources', 'consumers')])
galpar <- cbind(data.frame(dataset='salgal'), output_salgal[c('rep', 'areas', 'species', 'links', 'links_per_sp', 'indegree', 'resources', 'consumers')])

salgalpar <- rbind(salgal, galpar)
names(salgalpar)[2] <- 'replicate'

pyr <- cbind(data.frame(dataset='pyrenees'), output_pyrenees[c('replicate', 'areas', 'species', 'links', 'links_per_sp', 'indegrees', 'resources', 'consumers')])
names(pyr) <- names(europe_datasets)

#### Once we have extracted all the biogeographical datasets we merge them into a single data frame
biogeographical_datasets <- rbind(europe_datasets, salgalpar, pyr)

biogeographical_datasets$cr_ratio <- biogeographical_datasets$consumers/biogeographical_datasets$resources
#### In case you want to store the pre-processed data to have it all ready for the next time this script is run.
write.csv(biogeographical_datasets, file='merged_biogeographical-data_null2.csv')

#### Here we summarise and normalise the data

#### For the biogeographical regions this might take a while. Also, if not run for all replicates,
#### but using only the first replicate provided above instead for the European bioregions
#### would not yield the same results as show on the paper.

#### For convenience then, we have decided to provide a data table with this information ready to be loaded.
#### Run the following line instead of the two instructions below it to obtain the summarised and 
#### normalised data for the biogeographical datasets.
all_summaries_biogeo_norm <- read.csv('all_summaries_biogeographical.csv')

all_summaries_biogeo <- biogeographical_datasets %>% 
  group_by(dataset, areas) %>% 
  dplyr::summarize(
    species = mean(species),
    links = mean(links),
    links_per_sp = mean(links_per_sp),
    indegree = mean(indegree),
    potential_indegree = mean(potential_indegree),
    cr_ratio = mean(cr_ratio)
  ) 

all_summaries_biogeo_norm  <- all_summaries_biogeo %>%
  split(.$dataset) %>%
  map_dfr(range_stats, .id = "dataset")

##################################################################################
######## Once the data has been prepared we can proceed to fit the different models to 
######## the network-arearelationships (NARs)

#### the outcome of model fitting will be stored here
output_models_reps <- NULL

#### These are the properties that we want to analyse
properties <-  c('species', 'links', 'links_per_sp', 'indegree')

for(d in unique(all_summaries_replicates$dataset)){
  summarised_data <- subset(all_summaries_replicates, (dataset == d))
  for(p in properties){
    model_ranking <- sar_average(data=as.data.frame(summarised_data[ c('areas', p) ]), alpha_normtest = 0, alpha_homotest = 0)
    ##### This was changed to allow for multiple models to be selected and reported
    sum_ranking <- tryCatch({
      summary(model_ranking)
    }, warning = function(w) {
      NA
    }, error = function(e) {
      NA
    }, finally = {
      
    })
    
    if(is.na(sum_ranking)) next
    
    sum_ranking$Model_table <- sum_ranking$Model_table[order(sum_ranking$Model_table$AIC),]
    
    if(dim(sum_ranking$Model_table)[1] < 5){
      sel_models <- as.character(sum_ranking$Model_table$Model[1:dim(sum_ranking$Model_table)[1]])
    }else{
      sel_models <- as.character(sum_ranking$Model_table$Model[1:5])
    }
    idx <- 0
    
    if(!('power' %in% sel_models) & length(which(as.character(sum_ranking$Model_table$Model) == 'power')) != 0){
      idx_power <- which(as.character(sum_ranking$Model_table$Model) == 'power')
      sel_models <- append(sel_models, 'power')
    }
    
    for(cur_mod_name in sel_models){
      idx <- idx + 1
      if(idx > 5){
        idx <- idx_power
      }
      cur_model <- eval(parse(text=paste0('model_ranking$details$fits$', cur_mod_name)))
      
      mod_sum <- summary(cur_model)
      cur_model <- tryCatch({
        cbind(dataset=d, property=p, model=mod_sum$Model, ranking=idx, AIC=mod_sum$AIC, AkaikeWeight=subset(sum_ranking$Model_table, Model == cur_mod_name)$Weight, AICc=mod_sum$AICc, BIC=mod_sum$BIC, R2=mod_sum$R2, formula=as.character(mod_sum$formula), as.data.frame(mod_sum$Parameters[,1:4]))
      }, warning = function(w) {
        NA
      }, error = function(e) {
        NA
      }, finally = {
        
      })
      
      if(!is.na(cur_model)){
        cur_model$param <- mod_sum$parNames
        
        if(is.null(output_models_reps)){
          output_models_reps <- cur_model
        }else{
          output_models_reps <- rbind(output_models_reps, cur_model)
        }
      }
    }
  }
}

#### Here we store the results for the replicates datasets
write.csv(output_models_reps, file='fits-sars-replicates-nullmodel2.csv')


##### We repeat the same procedure above for the biogeographical datasets
#### the outcome of model fitting will be stored here
output_models_biogeo <- NULL

#### These are the properties that we want to analyse
properties <-  c('species', 'species', 'links', 'links_per_sp', 'indegree', 'potential_indegree')

for(d in unique(all_summaries_biogeo$dataset)){
  summarised_data <- subset(all_summaries_biogeo[complete.cases(all_summaries_biogeo),], (dataset == d))
  for(p in properties){
    model_ranking <- sar_average(data=as.data.frame(summarised_data[ c('areas', p) ]), alpha_normtest = 0, alpha_homotest = 0)
    ##### This was changed to allow for multiple models to be selected and reported
    sum_ranking <- tryCatch({
      summary(model_ranking)
    }, warning = function(w) {
      NA
    }, error = function(e) {
      NA
    }, finally = {
      
    })
    
    if(is.na(sum_ranking)) next
    
    sum_ranking$Model_table <- sum_ranking$Model_table[order(sum_ranking$Model_table$AIC),]
    
    if(dim(sum_ranking$Model_table)[1] < 5){
      sel_models <- as.character(sum_ranking$Model_table$Model[1:dim(sum_ranking$Model_table)[1]])
    }else{
      sel_models <- as.character(sum_ranking$Model_table$Model[1:5])
    }
    idx <- 0
    
    if(!('power' %in% sel_models) & length(which(as.character(sum_ranking$Model_table$Model) == 'power')) != 0){
      idx_power <- which(as.character(sum_ranking$Model_table$Model) == 'power')
      sel_models <- append(sel_models, 'power')
    }
    
    for(cur_mod_name in sel_models){
      idx <- idx + 1
      if(idx > 5){
        idx <- idx_power
      }
      cur_model <- eval(parse(text=paste0('model_ranking$details$fits$', cur_mod_name)))
      
      mod_sum <- summary(cur_model)
      cur_model <- tryCatch({
        cbind(dataset=d, property=p, model=mod_sum$Model, ranking=idx, AIC=mod_sum$AIC, AkaikeWeight=subset(sum_ranking$Model_table, Model == cur_mod_name)$Weight, AICc=mod_sum$AICc, BIC=mod_sum$BIC, R2=mod_sum$R2, formula=as.character(mod_sum$formula), as.data.frame(mod_sum$Parameters[,1:4]))
      }, warning = function(w) {
        NA
      }, error = function(e) {
        NA
      }, finally = {
        
      })
      
      if(!is.na(cur_model)){
        cur_model$param <- mod_sum$parNames
        
        if(is.null(output_models_biogeo)){
          output_models_biogeo <- cur_model
        }else{
          output_models_biogeo <- rbind(output_models_biogeo, cur_model)
        }
      }
    }
  }
}

#### Here we store the results for the biogeographical datasets
write.csv(output_models_biogeo, file='fits-sars-bioregions-nullmodel2.csv')

#### Lastly, we repeat the same code above but changed slightly to calculate the links-species relationships
#### We include this code only once. Repeat this for the all_summaries_biogeo
models_links_sp <- NULL
for(d in unique(all_summaries_replicates$dataset)){
  summarised_data <- subset(all_summaries_replicates, (dataset == d))
  model_ranking <- sar_average(data=as.data.frame(summarised_data[ c('species', 'links') ]), alpha_normtest = 0, alpha_homotest = 0)
  ##### This was changed to allow for multiple models to be selected and reported
  sum_ranking <- tryCatch({
    summary(model_ranking)
  }, warning = function(w) {
    NA
  }, error = function(e) {
    NA
  }, finally = {
    
  })
  
  if(is.na(sum_ranking)) next
  
  sum_ranking$Model_table <- sum_ranking$Model_table[order(sum_ranking$Model_table$AIC),]
  
  if(dim(sum_ranking$Model_table)[1] < 5){
    sel_models <- as.character(sum_ranking$Model_table$Model[1:dim(sum_ranking$Model_table)[1]])
  }else{
    sel_models <- as.character(sum_ranking$Model_table$Model[1:5])
  }
  idx <- 0
  
  if(!('power' %in% sel_models) & length(which(as.character(sum_ranking$Model_table$Model) == 'power')) != 0){
    idx_power <- which(as.character(sum_ranking$Model_table$Model) == 'power')
    sel_models <- append(sel_models, 'power')
  }
  
  for(cur_mod_name in sel_models){
    idx <- idx + 1
    if(idx > 5){
      idx <- idx_power
    }
    cur_model <- eval(parse(text=paste0('model_ranking$details$fits$', cur_mod_name)))
    plot(cur_model)
    mod_sum <- summary(cur_model)
    cur_model <- tryCatch({
      cbind(dataset=d, model=mod_sum$Model, ranking=idx, AIC=mod_sum$AIC, AkaikeWeight=subset(sum_ranking$Model_table, Model == cur_mod_name)$Weight, AICc=mod_sum$AICc, BIC=mod_sum$BIC, R2=mod_sum$R2, formula=as.character(mod_sum$formula), as.data.frame(mod_sum$Parameters[,1:4]))
    }, warning = function(w) {
      NA
    }, error = function(e) {
      NA
    }, finally = {
      
    })
    
    if(!is.na(cur_model)){
      cur_model$param <- mod_sum$parNames
      
      if(is.null(models_links_sp)){
        models_links_sp <- cur_model
      }else{
        models_links_sp <- rbind(models_links_sp, cur_model)
      }
    }
  }
}


# we store the output of the species-links
write.csv(models_links_sp, file = 'fits-species-links-relationships.csv')

#### Finally, use the code below to create the NAR figures in the paper
#### Change the dependent and indepedent variables to generate Figure 1 through 3 in the paper
#### Except Figure 3e and 3f.
require(ggplot2)
require(RColorBrewer)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

pd <- position_dodge(0.1)
ggplot(all_summaries_biogeo_norm[complete.cases(all_summaries_biogeo_norm) & all_summaries_biogeo_norm$dataset!='Macaronesia',], aes(x=log10(areas), y=(indegree), colour=dataset)) + 
  geom_line() +
  geom_point() + theme_bw() + ylab('log10(Indegree)') + xlab('log10(Area)')+
  scale_color_manual(values=getPalette(19))


pdf('links-sp-bioregions.pdf')
ggplot(all_summaries_replicates_norm, aes(x=species, y=links, colour=dataset)) + 
  geom_line(position=pd) +
  geom_point(position=pd) + theme_bw() +
  scale_color_manual(name='Dataset', values=getPalette(19)) + ylab('log10(links)')  + xlab('log10(species)')
dev.off()

### Plotting results for only sampled networks

all_summaries_sampled <- rbind(all_summaries_replicates_norm,all_summaries_biogeo_norm[all_summaries_biogeo_norm$dataset==c('Salgal','Galpar'),])
ggplot(all_summaries_replicates_norm, aes(x=log10(areas), y=(links), colour=dataset)) + 
  geom_line() +
  geom_point() + theme_bw() + ylab('log10(Links)') + xlab('log10(Area)')+
  scale_color_manual(values=getPalette(21))

