## ---------------------------
##
## Script name: data_processing_agg.r
##
## Purpose of script: This script constructs networks of ecological networks
## at different spatial scales according to the details of the relevant dataset.
## It then calculates a series of network properties which are in turn analysed 
## to build network-area relationships. These include fits of different functions
## to the degree distributions of the networks across spatial scales.
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
## The raw data used within this script was provided by 
## Professor Christian Mulder
## University of Catania
## Catania, Italy
##
## ---------------------------

## This script reads a series of files from the home directory. Ensure these files are placed in the same
## directory from which this script is executed

## Load required libraries
## This script relies on the libjvm java library to read Microsoft Excel files using the rJava
## and XLConnect R packages. Please make sure the line below spcifies the right location for that
## library on your local computer
# dyn.load('/Library/Java/JavaVirtualMachines/jdk-9.jdk/Contents/Home/lib/server/libjvm.dylib')

# dyn.load('/usr/lib/jvm/java-11-openjdk-amd64/lib/server/libjvm.so')

library(rJava)
library(XLConnect)
require(bipartite)
require(igraph)
require(AICcmodavg)
source('../utils.r')

## The section below reads the original raw data
metadata <- read.csv('metadata.csv', stringsAsFactors = F, header = T)
distances <- dist(metadata[c('Longitude.or.x.coordinate','Latitiude.or.y.coordinate')], method='euclidean')

excel <- loadWorkbook("./raw-data/interactions-1.xls")

# get sheet names
sheet_names <- getSheets(excel)
names(sheet_names) <- sheet_names

# put sheets into a list of data frames
sheet_list <- lapply(sheet_names, function(.sheet){readWorksheet(object=excel, .sheet)})
n.sites <- length(sheet_list)

##### because this data set is split into two different files we need to repeat
##### the above process for the other file
excel <- loadWorkbook("./raw-data/interactions-2.xls")

# get sheet names
sheet_names_2 <- getSheets(excel)
names(sheet_names_2) <- sheet_names_2

# put sheets into a list of data frames
sheet_list_2 <- lapply(sheet_names_2, function(.sheet){readWorksheet(object=excel, .sheet)})
n.sites <- n.sites + length(sheet_list_2)

sheet_list <- append(sheet_list, sheet_list_2)

sheet_names <- append(sheet_names, sheet_names_2)
sheet_names <- gsub('Monster', '', sheet_names)

#### Once the relevant raw data has been loaded we proceed with the network construction

#### first we define the categories (i.e. the different network types)
site_cats <- sort(unique(metadata$Site.type.code))
output <- NULL

#### use this flag to specify whether you want to run the null model(s) on this data
#### you will also need to comment / uncomment some lines of code below depending on which 
#### null model you wish to run
null_model <- FALSE


by_distance <- FALSE
n_reps <- 100
degree_dist_params <- NULL
#since we want to have one NAR per type of habitat we do this:
for(s_type in site_cats){
  cur_meta <- metadata[which(metadata$Site.type.code == s_type),]
  
  if(by_distance){
    replicates <- cur_meta$Site.ID
    #since we want to aggregate based on distances we calculate the distances among the current sites
    cur_dists_mat <- as.matrix(dist(cur_meta[c('Longitude.or.x.coordinate','Latitiude.or.y.coordinate')], method='euclidean'))
    
  }else{
    replicates <- 1:n_reps
  }
  
  first <- TRUE
  idxs <- match(cur_meta$Site.ID, sheet_names)
  
  ###### this check for NAs needs to be done because there are more networks
  ###### in the sheet list than sites on the metadata table. These networks
  ###### cannot be used for the analyses because they can't be matched to any
  ###### site type, therefore they are removed from the analyses.
  if(length(which(is.na(idxs))) != 0){
    idxs <- idxs[-which(is.na(idxs))]
  }
  for(i in idxs){
    n <- sheet_list[[i]]
    n <- n[c('Prooi', 'Predator')]
    n$Prooi <- as.character(n$Prooi)
    n$Predator <- as.character(n$Predator)
    
    local_net <- graph_from_edgelist(as.matrix(n), directed = T)
    local_net <- simplify(local_net)
    
    if(first){
      metaweb <- local_net
      first <- FALSE
    }else{
      metaweb <- igraph::union(metaweb, local_net)
    }
  }
  ##### and then, within each site type we repeat the network construction and analysis procedure
  ##### for several replicates
  for(r in replicates){
    print(r)
    ####### Here, we initialise a series of vector data structures to keep track of the 
    ####### network properties that are calculated over the networks at each locality
    species <- c()
    links <- c()
    connectances <- c()
    links_per_sp <- c()
    indegree <- c()
    outdegree <- c()
    sd_vul <- c()
    sd_gen <- c()
    basal <- c()
    top <- c()
    intermediate <- c()
    S_basal <- c()
    S_intermediate <- c()
    S_top <- c()
    omnivory <- c()
    mfcls <- c()
    modularity <- c()
    
    normalised_indegree <- c()
    normalised_outdegree <- c()
    
    potential_indegree <- c()
    
    ######## end of data structures
    
    ##### if sites are ordered by distance for aggregation, we do this
    if(by_distance){
      start_site <- r
      cur_dists <- cur_dists_mat[which(cur_meta$Site.ID == start_site),]
      names(cur_dists) <- cur_meta$Site.ID
      
      cur_dists <- sort(cur_dists)
      idxs <- names(cur_dists)
      
    }else{
      idxs <- sample(cur_meta$Site.ID)
    }
    
    idxs <- match(idxs, sheet_names)
    if(length(which(is.na(idxs))) != 0){
      idxs <- idxs[-which(is.na(idxs))]
    }
    
    first <- TRUE
    
    colors <- rev(rainbow(length(idxs), alpha=0.7))
    cc <- 1
    area_count <- 0
    for(i in idxs){
      ##### This section of the code reads the original raw data and builds a network
      ##### from it, storing it in the 'local_net' igraph structure
      area_count <- area_count + 1
      n <- sheet_list[[i]]
      n <- n[c('Prooi', 'Predator')]
      n$Prooi <- as.character(n$Prooi)
      n$Predator <- as.character(n$Predator)
      
      local_net <- graph_from_edgelist(as.matrix(n), directed = T)
      local_net <- simplify(local_net)
      
      
      ##### As areas get aggregated we store the regional (union) network in the
      ##### regional_net network structure
      if(first){
        regional_net_original <- local_net
        first <- FALSE
      }else{
        regional_net_original <- igraph::union(regional_net_original, local_net)
      }
      
      ##### This 'if' implements the two null models described in the paper.
      ##### Comment / uncomment each section appropriately
      if(null_model){
        ############### These lines are for the null models ######################
        random_net <- induced_subgraph(metaweb, sample(V(metaweb)$name, vcount(regional_net_original)))
        # random_net <- igraph::delete.vertices(random_net, names(which(igraph::degree(random_net) == 0)))
        
        #### for the second null model uncomment the line below (and comment the lines above)
        # random_net <- sample_gnm(vcount(regional_net_original), ecount(regional_net_original), directed=TRUE)
        # random_net <- igraph::delete.vertices(random_net, names(which(igraph::degree(random_net) == 0)))
        ##########################################################################
        
        regional_net <- random_net
      }else{
        regional_net <- regional_net_original
      }
      
      ##### Once we have the network corresponding to the current area we can calculate
      ##### a bunch of network properties. Those included in the paper and some more that
      ##### we also looked at
      S <- vcount(regional_net)
      species <- append(species, S)
      ls <- ecount(regional_net)
      links <- append(links, ls)
      connectances <- append(connectances, (2*ls/((S-1)*S)))
      links_per_sp <- append(links_per_sp, ls/S)
      
      n <- as_adjacency_matrix(regional_net)
      n <- as.matrix(n)
      
      indegree <- append(indegree, MeanGenerality(n))
      
      if(length(which(igraph::degree(regional_net, V(regional_net), mode='in') == 0)) != 0){
        normal_id <- mean(igraph::degree(regional_net, V(regional_net), mode='in')[-which(igraph::degree(regional_net, V(regional_net), mode='in') == 0)]/(ls/S))  
      }else{
        normal_id <- mean(igraph::degree(regional_net, V(regional_net), mode='in')/(ls/S))
      }
      
      normalised_indegree <- append(normalised_indegree, normal_id)
      predators <- names(which(igraph::degree(regional_net, mode='in') != 0))
      potential_indegree <- append(potential_indegree,  mean(igraph::degree(metaweb, predators, mode='in')))
      
      outdegree <- append(outdegree, MeanVulnerability(n))
      
      if(length(which(igraph::degree(regional_net, V(regional_net), mode='out') == 0)) != 0){
        normal_od <- mean(igraph::degree(regional_net, V(regional_net), mode='out')[-which(igraph::degree(regional_net, V(regional_net), mode='out') == 0)]/(ls/S))  
      }else{
        normal_od <- mean(igraph::degree(regional_net, V(regional_net), mode='out')/(ls/S))
      }
      
      normalised_outdegree <- append(normalised_outdegree, normal_od)
      sd_gen <- append(sd_gen, (SDGenerality(n)/(ls/S)))
      sd_vul <- append(sd_vul, (SDVulnerability(n)/(ls/S)))
      
      #####
      basal <- append(basal, FractionOfBasal(n))
      top <- append(top, FractionOfTop(n))
      intermediate <- append(intermediate, FractionOfIntermediate(n))
      S_basal <- append(S_basal, NumberOfBasal(n))
      S_top <- append(S_top, NumberOfTop(n))
      S_intermediate <- append(S_intermediate, NumberOfIntermediate(n))
      
      modularity <- append(modularity, max(walktrap.community(regional_net)$modularity))
      
      ########## this is for the degree distributions
      if(TRUE){
        degs <- igraph::degree(regional_net, mode='in')
        degs <- degs[-which(degs == 0)]
        
        #### a normalisation term
        norm_term <- 1 # ecount(regional_net)/vcount(regional_net)
        
        occur = as.vector(table(degs/norm_term))
        occur = occur/sum(occur)
        p = occur/sum(occur)
        y = rev(cumsum(rev(p)))
        x = as.numeric(names(table(degs/norm_term)))
        
        ##### comment out the following code to avoid drawing degree distributions
        # if(cc == 1){
        #   plot(x,y, log='xy', ylim=c(10^-3,1), xlim=c(10^0, 10^2), pch=16, cex=1.5, cex.lab=1.5, cex.axis=1.2, xlab='', ylab='log Pc(k)', col=colors[cc], tck=0, xaxt='n', yaxt='n')
        #   title(xlab="log k", line=2, cex.lab=1.5)
        #   box(lwd=1.5)
        #   x1 <- floor(log10(range(x)))
        #   # pow <- seq(x1[1], x1[2]+1)
        #   pow <- c(-1,0,1)
        #   
        #   ticksat <- as.vector(sapply(pow, function(p) (1:10)*10^p))[2:20]
        #   axis(1, 10^pow[2:3], tck=0.02, lwd=1.5, cex.axis=1.5, lwd.ticks=1.5, labels=c(expression(10^{0}), expression(10^{1})))
        #   axis(1, ticksat, labels=NA, tcl=-0.25, lwd=0, lwd.ticks=1, tck=0.01)
        #   
        #   
        #   y1 <- floor(log10(range(y)))
        #   # pow <- seq(y1[1], y1[2])
        #   pow <- c(-3,-2,-1,0)
        #   ticksat <- as.vector(sapply(pow, function(p) (1:10)*10^p))
        #   ticksat <- ticksat[3:30]
        #   axis(2, 10^pow[2:4], tck=0.02, lwd=1.5, cex.axis=1.5, lwd.ticks=1.5, labels=c(expression(10^{-2}), expression(10^{-1}), expression(10^{0})), las=1, mgp=c(1,.5,0))
        #   axis(2, ticksat, labels=NA, tcl=-0.25, lwd=0, lwd.ticks=1, tck=0.01)
        #   
        # }else{
        #   points(x,y,pch=16, cex.lab=1.5, cex.axis=1.5, col=colors[cc], cex=1.5)
        # }
        # 
        # cc <- cc + 1
        #### end of drawing code
        
        
        ###### To identify the shape of the degree distributions we fit four different models to it
        ###### the fitted models are: power law, exponential, truncated power law, and log-normal
        temp <- data.frame(x, y)
        failed <- tryCatch({
          mod1 <- tryCatch({
            ##### Truncated power law
            nls(y ~ ( (x^-a) *(exp(-x/b))), data = temp, start = list(a = .0001, b = 2), control=nls.control(maxiter = 1e3))
          }, error = function(e) {
            NA
          }, finally = {
          })
          mod2 <- tryCatch({
            ##### Exponential
            nls(y ~ (exp(-x/b)), data = temp, start = list(b = 2), control=nls.control(maxiter = 1000))
          }, error = function(e) {
            NA
          }, finally = {
          })
          mod3 <- tryCatch({
            ##### Power law
            nls(y ~ (x^-a), data = temp, start = list(a = .01), control=nls.control(maxiter = 1e3))
          }, error = function(e) {
            NA
          }, finally = {
          })
          mod4 <- tryCatch({
            ##### Log-normal
            nls(y ~ ( (1/ (x * b * sqrt(2*pi) )) * exp(- ( ((log(x) - a)^2) / (2*(b^2)) )) ), data = temp, start = list(a = .3, b = .3), control=nls.control(maxiter = 1000))
          }, error = function(e) {
            NA
          }, finally = {
          })
          FALSE
        }, warning = function(w) {
          FALSE
        }, error = function(e) {
          TRUE
        }, finally = {
          
        })
        
        ##### Once we have all model outputs we compare them using AIC and 
        ##### keep the result into the corresponding model fits data frame
        ##### depending on the type of model selected
        if(!failed){
          model_list <- list(mod1,mod2,mod3,mod4)
          names_list <- c('mod1','mod2','mod3','mod4')
          if(length(which(is.na(model_list))) != 0){
            names_list <- names_list[-which(is.na(model_list))]
            model_list <- model_list[-which(is.na(model_list))]  
          }
          names(model_list) <- names_list
          
          if(length(model_list) == 0) next
          
          if(length(model_list) == 1){
            model_name <- names_list[1]
            model <- summary(eval(as.symbol(model_name)))  
          }else{
            aic_comp <- aictab(model_list)
            model_name <- tolower(as.character(aic_comp$Modnames[1]))
            model <- summary(eval(as.symbol(model_name)))
          }
          
          if(model_name == 'mod1' | model_name == 'mod4'){
            cur_out <- data.frame(s_type, r, area=area_count, model=model_name, a=model$coefficients[1,1], a.std.err=model$coefficients[1,2], a.tval=model$coefficients[1,3], a.pval=model$coefficients[1,4], b=model$coefficients[2,1], b.std.err=model$coefficients[2,2], b.tval=model$coefficients[2,3], b.pval=model$coefficients[2,4])
          }else if(model_name == 'mod2'){
            cur_out <- data.frame(s_type, r, area=area_count, model=model_name, a=NA, a.std.err=NA, a.tval=NA, a.pval=NA, b=model$coefficients[1,1], b.std.err=model$coefficients[1,2], b.tval=model$coefficients[1,3], b.pval=model$coefficients[1,4])
          }else{
            cur_out <- data.frame(s_type, r, area=area_count, model=model_name, a=model$coefficients[1,1], a.std.err=model$coefficients[1,2], a.tval=model$coefficients[1,3], a.pval=model$coefficients[1,4], b=NA, b.std.err=NA, b.tval=NA, b.pval=NA)
          }
          
          if(is.null(degree_dist_params)){
            degree_dist_params <- cur_out
          }else{
            degree_dist_params <- rbind(degree_dist_params, cur_out)
          }
        }
      }
    }
    
    ##### Lastly, we store all the information for the network properties for the current replicate
    ##### into the data frame containing the output data for all the runs
    cur_out <- data.frame(site_type = rep(s_type,length(species)), replicate=rep(r,length(species)), areas=1:length(species), species, links, connectances, links_per_sp, indegree, outdegree, sd_gen, sd_vul, basal, top, intermediate, S_basal, S_top, S_intermediate, modularity, normalised_indegree, normalised_outdegree, potential_indegree) #, omnivory, mfcls) 
    if(is.null(output)){
      output <- cur_out
    }else{
      output <- rbind(output, cur_out)
    }
  }
}

##### After the runs are finished we store the output for the network properties across 
##### all areas (output) and the outcome of the analyses of the degree distributions (degree_dist_params)
##### into csv files that are readable for further analyses
write.csv(output, 'output-soils.csv')
write.csv(degree_dist_params, 'fits-degree-dists-soils.csv')
