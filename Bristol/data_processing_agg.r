## ---------------------------
##
## Script name: data_processing_agg.r
##
## Purpose of script: This script constructs networks of ecological networks
## at different spatial scales according to the details of the relevant dataset
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
## Dr Daniel Montoya
## Centre for Biodiversity Theory and Modelling
## Theoretical and Experimental Ecology Station, CNRS, France
##
## ---------------------------

## This script reads a series of files from the home directory. Ensure these files are placed in the same
## directory from which this script is executed

## Load required libraries
require(igraph)
require(bipartite)
require(AICcmodavg)
source('../utils.r')

## The section below reads the original raw data

#### Step 1: Read input data.
#### For this dataset networks are given as independent csv files with the adjacency
#### matrix inside
#### some metadata can be read
input_data <- read.csv('./metadata.csv', stringsAsFactors = F)
input_data_agg <- input_data[(input_data$Island.size > .5 & input_data$Island.size<=2),]
species_types <- read.csv('./species-types.csv', stringsAsFactors = F)
species_types$Basal <- gsub(' ', '.', species_types$Basal)
species_types$Intermediate <- gsub(' ', '.', species_types$Intermediate)
species_types$Top <- gsub(' ', '.', species_types$Top)

#### Step 2: We remove the islands with too few species (c('ST23' , 'PT6'))
#### because these do not represent real local communities
ids <- sort(unique(input_data$ISLAND.CODE))
removed_islands <- c('ST23' , 'PT6')
ids_agg <- sort(unique(input_data_agg$ISLAND.CODE))
ids_agg <- ids_agg[-which(ids_agg %in% removed_islands)]

#### Once the relevant raw data has been loaded we proceed with the network construction

#### Step 3: Since some of the metrics to be calculated require previous knowledge
#### of the metaweb, we calculate this network beforehand to be able to use it later on
first <- TRUE
for(net_id in ids_agg){
  #### here we read each local network
  cur_com <- read.csv(paste0('./raw-data/data',net_id,'.csv'), sep=';', stringsAsFactors = F)
  
  #### the names of rows and columns do not match... they need to match to be able to construct the network
  cur_com <- cur_com[-1]
  row.names(cur_com) <- names(cur_com)
  
  cur_com[is.na(cur_com)] <- 0
  n <- as.matrix(cur_com)
  
  remove <- which(rowSums(n)==0 & colSums(n)==0)
  if(length(remove) != 0){
    n <- n[-remove, -remove]
  }
  local_net <- graph_from_adjacency_matrix(n)
  
  #### for this particular dataset we need to remove some edges that should not appear twice
  #### first the ones to basal species
  for(b in intersect(V(local_net)$name, species_types$Basal)){
    to_remove <- incident(local_net, b , mode='in')
    local_net <- igraph::delete.edges(local_net, to_remove)
  }
  
  #### and the same thing for the top predators
  for(tp in intersect(V(local_net)$name, species_types$Top)){
    to_remove <- incident(local_net, tp , mode='out')
    local_net <- igraph::delete.edges(local_net, to_remove)
  }
  
  if(first){
    metaweb <- local_net 
    first <- FALSE
  }else{
    metaweb <- igraph::union(metaweb, local_net)
  }  
}

#### The algorithm below produces several different replicates of network aggregation of local
#### networks in an ever-increasing fashion to produce ecological networks at different spatial scales

bipartite <- FALSE
output <- NULL
degree_dist_params <- NULL
n.sites <- length(ids_agg)

null_model <- FALSE

#### This is to treat islands as replicates and proceed to aggregate them ###
replicates <- 100
for(r in 1:replicates){
  print(r)
  ####### Here, we initialise a series of vector data structures to keep track of the 
  ####### network properties that are calculated over the networks at each locality
  areas <- c()
  species <- c()
  links <- c()
  connectances <- c()
  links_per_sp <- c()
  indegree <- c()
  outdegree <- c()
  basal <- c() 
  top <- c() 
  intermediate <- c()
  S_basal <- c()
  S_intermediate <- c()
  S_top <- c()
  modularity <- c()
  sd_gen <- c()
  sd_vul <- c()
  normalised_indegree <- c()
  normalised_outdegree <- c()
  potential_indegree <- c()
  
  resources <- c()
  consumers <- c()
  
  ######## end of data structures
  
  ####### sites are randomly ordered for the random aggregation
  idxs <- sample(1:n.sites, n.sites)
  idxs <- ids_agg[idxs]
  first <- TRUE
  i <- 1
  colors <- rev(rainbow(length(idxs), alpha=0.7))
  cc <- 1
  counter <- 0
  for(net_id in idxs){
    #### Step 4: We create networks for each one of the local communities
    counter <- counter + 1  
    ##### This section of the code reads the original raw data and builds a network
    ##### from it, storing it in the 'local_net' igraph structure
    cur_com <- read.csv(paste0('./raw-data/data',net_id,'.csv'), sep=';', stringsAsFactors = F)
    
    #### the names of rows and columns do not match... they need to match to be able to construct the network
    cur_com <- cur_com[-1]
    row.names(cur_com) <- names(cur_com)
    
    cur_com[is.na(cur_com)] <- 0
    n <- as.matrix(cur_com)
    
    remove <- which(rowSums(n)==0 & colSums(n)==0)
    if(length(remove) != 0){
      n <- n[-remove, -remove]
    }
    local_net <- graph_from_adjacency_matrix(n)
    
    #### for this particular dataset we need to remove some edges that should not appear twice
    #### first the ones to basal species
    for(b in intersect(V(local_net)$name, species_types$Basal)){
      to_remove <- incident(local_net, b , mode='in')
      local_net <- igraph::delete.edges(local_net, to_remove)
    }
    
    #### and the same thing for the top predators
    for(tp in intersect(V(local_net)$name, species_types$Top)){
      to_remove <- incident(local_net, tp , mode='out')
      local_net <- igraph::delete.edges(local_net, to_remove)
    }
    
    ##### As areas get aggregated we store the regional (union) network in the
    ##### regional_net network structure
    if(first){
      regional_net_original <- local_net 
      first <- FALSE
    }else{
      regional_net_original <- igraph::union(regional_net_original, local_net)
    }  
    
    a <- input_data_agg[match(idxs[i], input_data_agg$ISLAND.CODE),]$Island.size
    if(i == 1){
      areas <- append(areas, a)
    }else{
      areas <- append(areas, a+areas[length(areas)])  
    }
    
    ##### This 'if' implements the two null models described in the paper.
    ##### Comment / uncomment each section appropriately
    if(null_model){
      ############### These lines are for the null models ######################
      random_net <- induced_subgraph(metaweb, sample(V(metaweb)$name, vcount(regional_net_original)))
      
      
      #### for the second null model uncomment the line below (and comment the lines above)
      # random_net <- sample_gnm(vcount(regional_net_original), ecount(regional_net_original), directed=TRUE)
      
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
    
    n <- as_adjacency_matrix(regional_net, sparse = F)
    
    indegree <- append(indegree, MeanGenerality(n))
    
    if(length(which(igraph::degree(regional_net, V(regional_net), mode='in') == 0)) != 0){
      normal_id <- mean(igraph::degree(regional_net, V(regional_net), mode='in')[-which(igraph::degree(regional_net, V(regional_net), mode='in') == 0)]/(ls/S))  
    }else{
      normal_id <- mean(igraph::degree(regional_net, V(regional_net), mode='in')/(ls/S))
    }
    
    normalised_indegree <- append(normalised_indegree, normal_id)
    
    outdegree <- append(outdegree, MeanVulnerability(n))
    
    if(length(which(igraph::degree(regional_net, V(regional_net), mode='in') == 0)) != 0){
      normal_outd <- mean(igraph::degree(regional_net, V(regional_net), mode='out')[-which(igraph::degree(regional_net, V(regional_net), mode='out') == 0)]/(ls/S))  
    }else{
      normal_outd <- mean(igraph::degree(regional_net, V(regional_net), mode='out')/(ls/S))
    }
    
    normalised_outdegree <- append(normalised_outdegree, normal_outd)
    
    predators <- names(which(igraph::degree(regional_net, mode='in')!=0))
    potential_indegree <- append(potential_indegree, mean(igraph::degree(metaweb, predators, mode='in')))
    
    sd_gen <- append(sd_gen, SDGenerality(n)/(ls/S))
    sd_vul <- append(sd_vul, SDVulnerability(n)/(ls/S))
    
    #####
    basal <- append(basal, FractionOfBasal(n))
    top <- append(top, FractionOfTop(n))
    intermediate <- append(intermediate, FractionOfIntermediate(n))
    S_basal <- append(S_basal, NumberOfBasal(n))
    S_top <- append(S_top, NumberOfTop(n))
    S_intermediate <- append(S_intermediate, NumberOfIntermediate(n))
    modularity <- append(modularity, max(walktrap.community(regional_net)$modularity))
    i <- i+1
    
    resources <- append(resources, length(which(igraph::degree(regional_net, mode='out') != 0)))
    consumers <- append(consumers, length(which(igraph::degree(regional_net, mode='in') != 0)))
    
    #### This is for the degree distributions
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
        cur_out <- data.frame(r, area=counter, area_original=areas[counter], model=model_name, a=model$coefficients[1,1], a.std.err=model$coefficients[1,2], a.tval=model$coefficients[1,3], a.pval=model$coefficients[1,4], b=model$coefficients[2,1], b.std.err=model$coefficients[2,2], b.tval=model$coefficients[2,3], b.pval=model$coefficients[2,4])
      }else if(model_name == 'mod2'){
        cur_out <- data.frame(r, area=counter, area_original=areas[counter], model=model_name, a=NA, a.std.err=NA, a.tval=NA, a.pval=NA, b=model$coefficients[1,1], b.std.err=model$coefficients[1,2], b.tval=model$coefficients[1,3], b.pval=model$coefficients[1,4])
      }else{
        cur_out <- data.frame(r, area=counter, area_original=areas[counter], model=model_name, a=model$coefficients[1,1], a.std.err=model$coefficients[1,2], a.tval=model$coefficients[1,3], a.pval=model$coefficients[1,4], b=NA, b.std.err=NA, b.tval=NA, b.pval=NA)
      }
      
      if(is.null(degree_dist_params)){
        degree_dist_params <- cur_out
      }else{
        degree_dist_params <- rbind(degree_dist_params, cur_out)
      }
    }
    
  }
  
  ##### Lastly, we store all the information for the network properties for the current replicate
  ##### into the data frame containing the output data for all the replicates
  cur_out <- data.frame(replicate=rep(r,length(species)), area=1:counter, area_original=areas, species, links, connectances, links_per_sp, indegree, outdegree, basal, top, intermediate, S_basal, S_intermediate, S_top, sd_gen, sd_vul, modularity, normalised_indegree, normalised_outdegree, consumers, resources, potential_indegree) 
  if(is.null(output)){
    output <- cur_out
  }else{
    output <- rbind(output, cur_out)
  }
  
}


##### After the runs are finished we store the output for the network properties across 
##### all areas (output) and the outcome of the analyses of the degree distributions (degree_dist_params)
##### into csv files that are readable for further analyses
write.csv(output, file = 'output-bristol.csv')
write.csv(degree_dist_params, file = 'fits-degree-dists-bristol.csv')
