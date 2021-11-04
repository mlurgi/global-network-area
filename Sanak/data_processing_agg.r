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
## Dr Spencer Wood
## Outdoor Recreation & Data Lab
## University of Washington, USA
##
## ---------------------------

## This script reads a series of files from the home directory. Ensure these files are placed in the same
## directory from which this script is executed

## Load required libraries
require(igraph)
require(AICcmodavg)
source('../utils.r')

## The section below reads the original raw data

#### Step 1: Read input data.
#### In this example case the interactions are reported as adjacency lists per area
red.frag <- read.table("./raw-data/WoodEtal_Append1_v2.csv",sep=",",header=T,na.strings="", stringsAsFactors = F)
input_data <- red.frag

#### Step 2: From this table we create networks for each one of the local communities
location <- sort(unique(input_data[input_data$WebScale=="Q",]$WebID))

n.sites <- length(location)

#### Once the relevant raw data has been loaded we proceed with the network construction

## First we construct the metaweb
first <- TRUE
for(l in location){
  cur_com <- input_data[which(input_data$WebID == l),]
  
  local_net <- make_empty_graph()
  for(int in 1:dim(cur_com)[1]){
    prey <- cur_com[int,]$PreyName
    if(! prey %in% V(local_net)$name){
      local_net <- local_net + vertices(prey)
    }
    pred <- cur_com[int,]$PredName
    if(! pred %in% V(local_net)$name){
      local_net <- local_net + vertices(pred)
    }
    local_net[prey, pred] <- 1
  }
  
  if(first){
    metaweb <- local_net 
    first <- FALSE
  }else{
    metaweb <- igraph::union(metaweb, local_net)
  } 
}

###### metaweb construction ends here! ######


#### Data structures to store the outcomes of the replicates
output <- NULL
degree_dist_params <- NULL
replicates <- 100
null_model <- FALSE

#### The algorithm below produces several different replicates of network aggregation of local
#### networks in an ever-increasing fashion to produce ecological networks at different spatial scales
for(r in 1:replicates){
  print(r)
  ####### Here, we initialise a series of vector data structures to keep track of the 
  ####### network properties that are calculated over the networks at each locality
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
  consumers <- c()
  resources <- c()
  potential_indegree <- c()
  #### end of data structures
  
  #### random order of sites for random aggregation
  idxs <- sample(1:n.sites, n.sites)
  idxs <- location[idxs]
  first <- TRUE
  colors <- rev(rainbow(length(idxs), alpha=0.7))
  cc <- 1
  counter <- 0
  for(l in idxs){
    ##### This section of the code reads the original raw data and builds a network
    ##### from it, storing it in the 'local_net' igraph structure
    counter <- counter + 1
    cur_com <- input_data[which(input_data$WebID == l),]
    
    local_net <- make_empty_graph()
    for(int in 1:dim(cur_com)[1]){
      prey <- cur_com[int,]$PreyName
      if(! prey %in% V(local_net)$name){
        local_net <- local_net + vertices(prey)
      }
      pred <- cur_com[int,]$PredName
      if(! pred %in% V(local_net)$name){
        local_net <- local_net + vertices(pred)
      }
      local_net[prey, pred] <- 1
    }
    
    ##### As areas get aggregated we store the regional (union) network in the
    ##### whole_net network structure
    if(first){
      whole_net_original <- local_net 
      first <- FALSE
    }else{
      whole_net_original <- igraph::union(whole_net_original, local_net)
    }
    
    
    ##### This 'if' implements the two null models described in the paper.
    ##### Comment / uncomment each section appropriately  
    if(null_model){
      ############### These lines are for the null models ######################
      random_net <- induced_subgraph(metaweb, sample(V(metaweb)$name, vcount(whole_net_original)))
      # random_net <- igraph::delete.vertices(random_net, names(which(igraph::degree(random_net) == 0)))

      #### for the second null model uncomment the line below (and comment the lines above)
      # random_net <- sample_gnm(vcount(whole_net_original), ecount(whole_net_original), directed=TRUE)
      # random_net <- igraph::delete.vertices(random_net, names(which(igraph::degree(random_net) == 0)))
      ##########################################################################
      
      whole_net <- random_net
    }else{
      whole_net <- whole_net_original
    }
    
    #### Step 4: Calculate the network features for the current community
    
    ##### Once we have the network corresponding to the current area we can calculate
    ##### a bunch of network properties. Those included in the paper and some more that
    ##### we also looked at
    webID <- unique(cur_com$WebID)
    S <- vcount(whole_net)
    species <- append(species, S)
    ls <- ecount(whole_net)
    links <- append(links, ls)
    C <- (2*ls/((S-1)*S))
    connectances <- append(connectances, C)
    links_per_sp <- append(links_per_sp, ls/S)
    
    n <- as_adjacency_matrix(whole_net,sparse=F)
    
    indegree <- append(indegree, MeanGenerality(n))
    
    if(length(which(igraph::degree(whole_net, V(whole_net), mode='in') == 0)) != 0){
      normal_id <- mean(igraph::degree(whole_net, V(whole_net), mode='in')[-which(igraph::degree(whole_net, V(whole_net), mode='in') == 0)]/(ls/S))  
    }else{
      normal_id <- mean(igraph::degree(whole_net, V(whole_net), mode='in')/(ls/S))
    }
    
    normalised_indegree <- append(normalised_indegree, normal_id)
    
    consumers <- append(consumers, length(which(igraph::degree(whole_net, mode='in') != 0)))
    resources <- append(resources, length(which(igraph::degree(whole_net, mode='out') != 0)))
    
    outdegree <- append(outdegree, MeanVulnerability(n))
    
    if(length(which(igraph::degree(whole_net, V(whole_net), mode='out') == 0)) != 0){
      normal_outd <- mean(igraph::degree(whole_net, V(whole_net), mode='out')[-which(igraph::degree(whole_net, V(whole_net), mode='out') == 0)]/(ls/S))  
    }else{
      normal_outd <- mean(igraph::degree(whole_net, V(whole_net), mode='out')/(ls/S))
    }
    
    normalised_outdegree <- append(normalised_outdegree, normal_outd)
    
    sd_gen <- append(sd_gen, (SDGenerality(n)/(ls/S)))
    sd_vul <- append(sd_vul, (SDVulnerability(n)/(ls/S)))
    modularity <- append(modularity, max(walktrap.community(whole_net)$modularity))
    
    #####
    basal <- append(basal, FractionOfBasal(n))
    top <- append(top, FractionOfTop(n))
    intermediate <- append(intermediate, FractionOfIntermediate(n))
    S_basal <- NumberOfBasal(n)
    S_top <- NumberOfTop(n)
    S_intermediate <- NumberOfIntermediate(n)
    
    predators <- names(which(igraph::degree(whole_net, mode='in') != 0))
    potential_indegree <- append(potential_indegree,  mean(igraph::degree(metaweb, predators, mode='in')))
    
    #### This is for the degree distributions
    degs <- igraph::degree(whole_net, mode='in')
    degs <- degs[-which(degs == 0)]
    
    #### a normalisation term
    norm_term <- 1 # ecount(regional_net)/vcount(regional_net)
    
    occur = as.vector(table(degs/norm_term))
    occur = occur/sum(occur)
    p = occur/sum(occur)
    y = rev(cumsum(rev(p)))
    x = as.numeric(names(table(degs/norm_term)))
    
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
        cur_out <- data.frame(r, area=counter, model=model_name, a=model$coefficients[1,1], a.std.err=model$coefficients[1,2], a.tval=model$coefficients[1,3], a.pval=model$coefficients[1,4], b=model$coefficients[2,1], b.std.err=model$coefficients[2,2], b.tval=model$coefficients[2,3], b.pval=model$coefficients[2,4])
      }else if(model_name == 'mod2'){
        cur_out <- data.frame(r, area=counter, model=model_name, a=NA, a.std.err=NA, a.tval=NA, a.pval=NA, b=model$coefficients[1,1], b.std.err=model$coefficients[1,2], b.tval=model$coefficients[1,3], b.pval=model$coefficients[1,4])
      }else{
        cur_out <- data.frame(r, area=counter, model=model_name, a=model$coefficients[1,1], a.std.err=model$coefficients[1,2], a.tval=model$coefficients[1,3], a.pval=model$coefficients[1,4], b=NA, b.std.err=NA, b.tval=NA, b.pval=NA)
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
  cur_out <- data.frame(replicate=rep(r,length(species)), areas=1:length(species), species, links, connectances, links_per_sp, indegree, outdegree, top, basal, intermediate, S_top, S_basal, S_intermediate, modularity, sd_gen, sd_vul, normalised_indegree, normalised_outdegree, consumers, resources, potential_indegree) #, indegree_normalized, outdegree, outdegree_normalized, sd_gen, sd_vul, basal, top, intermediate, S_basal, S_top, S_intermediate, omnivory, mfcls)
  
  if(is.null(output)){
    output <- cur_out
  }else{
    output <- rbind(output, cur_out)
  }
  
}

##### After the runs are finished we store the output for the network properties across 
##### all areas (output) and the outcome of the analyses of the degree distributions (degree_dist_params)
##### into csv files that are readable for further analyses
write.csv(output, file='output-sanak.csv')
write.csv(degree_dist_params, file='fits-degree-dists-sanak.csv')



