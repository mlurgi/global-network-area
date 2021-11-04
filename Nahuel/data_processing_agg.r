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
## Dr Diego Vazquez
## Argentine Institute for Dryland Research
## Mendoza, Argentina
##
## ---------------------------

## This script reads a series of files from the home directory. Ensure these files are placed in the same
## directory from which this script is executed

## Load required libraries
require(igraph)
require(bipartite)
require(AICcmodavg)
source("../utils.r")

## The section below reads the original raw data
vazag = t(read.csv("./raw-data/vaz_ag_matr_f.txt", sep="\t", header=FALSE))
vazcl = t(read.csv("./raw-data/vaz_cl_matr_f.txt", sep="\t", header=FALSE))
vazll = t(read.csv("./raw-data/vaz_ll_matr_f.txt", sep="\t", header=FALSE))
vazmh = t(read.csv("./raw-data/vaz_mh_matr_f.txt", sep="\t", header=FALSE))
vazmnh = t(read.csv("./raw-data/vaz_mnh_matr_f.txt", sep="\t", header=FALSE))
vazqh = t(read.csv("./raw-data/vaz_qh_matr_f.txt", sep="\t", header=FALSE))
vazqnh = t(read.csv("./raw-data/vaz_qnh_matr_f.txt", sep="\t", header=FALSE))
vazs = t(read.csv("./raw-data/vaz_s_matr_f.txt", sep="\t", header=FALSE))

vaz = list(vazag,vazcl,vazll,vazmh,vazmnh,vazqh,vazqnh,vazs)

n.sites=length(vaz)
areas=1:n.sites #Define areas

#### Ww build the metaweb
first <- TRUE
for (i in areas){
  n <- vaz[[i]]
  colnames(n) <- paste0('Pol',1:dim(vaz[[1]])[2])
  rownames(n) <- paste0('Plant',1:dim(vaz[[1]])[1])
  
  r_remove <- which(rowSums(n)==0)
  if(length(r_remove) != 0){
    n <- n[-r_remove,]
  }
  c_remove <- which(colSums(n)==0)
  if(length(c_remove) != 0){
    n <- n[, -c_remove]
  }
  
  local_net <- graph_from_incidence_matrix(n, directed = T, mode='out')
  
  if(first){
    metaweb <- local_net 
    first <- FALSE
  }else{
    metaweb <- igraph::union(metaweb, local_net)
    V(metaweb)$type <- V(metaweb)$type_1
    V(metaweb)$type[which(!is.na(V(metaweb)$type_2))] <- V(metaweb)$type_2[which(! is.na(V(metaweb)$type_2))]
  } 
}

################## metaweb construction ends here!!!!! ##################

#### Once the relevant raw data has been loaded we proceed with the network construction
#### Data structures to store outputs
output <- NULL
replicates <- 100
degree_dist_params <- NULL

null_model <- TRUE

for(r in 1:replicates){
  ####### Here, we initialise a series of vector data structures to keep track of the 
  ####### network properties that are calculated over the networks at each locality
  counter <- 0
  species <- c()
  consumers <- c()
  resources <- c()
  links <- c()
  connectances <- c()
  links_per_sp <- c()
  indegree <- c()
  outdegree <- c()
  modularity <- c()
  normalised_indegree <- c()
  normalised_outdegree <- c()
  potential_indegree <- c()
  sd_gen <- c()
  sd_vul <- c()
  ######## end of data structures
  
  ##### Random order of sites for random aggregation
  f_indx <- sample(1:n.sites, n.sites)
  first <- TRUE
  
  colors <- rev(rainbow(length(areas), alpha=0.7))
  cc <- 1
  
  for (i in f_indx){
    #### Step 4: We create the interactions network for each local community
    
    ##### This section of the code reads the original raw data and builds a network
    ##### from it, storing it in the 'local_net' igraph structure
    counter <- counter + 1
    n <- vaz[[i]]
    colnames(n) <- paste0('Pol',1:dim(vaz[[1]])[2])
    rownames(n) <- paste0('Plant',1:dim(vaz[[1]])[1])
  
    r_remove <- which(rowSums(n)==0)
    if(length(r_remove) != 0){
      n <- n[-r_remove,]
    }
    c_remove <- which(colSums(n)==0)
    if(length(c_remove) != 0){
      n <- n[, -c_remove]
    }
    
    local_net <- graph_from_incidence_matrix(n, directed = T, mode='out')
    
    ##### As areas get aggregated we store the regional (union) network in the
    ##### regional_net_original network structure
    if(first){
      regional_net_original <- local_net 
      first <- FALSE
    }else{
      regional_net_original <- igraph::union(regional_net_original, local_net)
      V(regional_net_original)$type <- V(regional_net_original)$type_1
      V(regional_net_original)$type[which(!is.na(V(regional_net_original)$type_2))] <- V(regional_net_original)$type_2[which(! is.na(V(regional_net_original)$type_2))]
    } 
    
    ##### This 'if' implements the two null models described in the paper.
    ##### Comment / uncomment each section appropriately  
    if(null_model){
      ############### These lines are for the null models ######################
      # random_net <- induced_subgraph(metaweb, sample(V(metaweb)$name, vcount(regional_net_original)))
      # 
      # while(ecount(random_net) < 1){
      #   random_net <- induced_subgraph(metaweb, sample(V(metaweb)$name, vcount(regional_net_original)))
      # }
      # 
      # random_net_bip <- graph.bipartite(bipartite.mapping(random_net)$type, as.vector(t(igraph::get.edges(random_net, 1:length(E(random_net))))), directed=T)
      # V(random_net_bip)$name <- V(random_net)$name
      # 
      #### for the second null model uncomment the line below (and comment the lines above)
      failed <- TRUE
      while(failed){
        failed <- tryCatch({
          resource_species <- sample(vcount(regional_net_original)-1,1)
          consumer_species <- vcount(regional_net_original) - resource_species
          random_net_bip <- sample_bipartite( resource_species, consumer_species,  type='gnm', m=ecount(regional_net_original), directed=TRUE)
          FALSE
        }, warning = function(w) {
          FALSE
        }, error = function(e) {
          TRUE
        }, finally = {

        })
      }
       ##########################################################################
      regional_net <- random_net_bip
    }else{
      regional_net <- regional_net_original
    }
    
    ##### Once we have the network corresponding to the current area we can calculate
    ##### a bunch of network properties. Those included in the paper and some more that
    ##### we also looked at
    n <- as_incidence_matrix(regional_net)
    S <- vcount(regional_net)
    Sc <- length(which(V(regional_net)$type == TRUE))
    Sr <- length(which(V(regional_net)$type == FALSE))
    
    species <- append(species, S)
    resources <- append(resources, Sr)
    consumers <- append(consumers, Sc)
    ls <- ecount(regional_net)
    links <- append(links, ls)
    connec <- ls / (Sr*Sc)
    connectances <- append(connectances, connec)
    links_per_sp <- append(links_per_sp, ls/S)
    
    indeg <- ls/Sc
    indegree <- append(indegree, indeg)
    outdeg <- ls/Sr
    outdegree <- append(outdegree, outdeg)
    
    
    normalised_indegree <- append(normalised_indegree, indeg/(ls/S))
    normalised_outdegree <- append(normalised_outdegree, outdeg/(ls/S))
    
    pollinators <- names(which(igraph::degree(regional_net, mode='in') != 0))
    potential_indegree <- append(potential_indegree,  mean(igraph::degree(metaweb, pollinators, mode='in')))
    
    sd_gen <- append(sd_gen, sd((specieslevel(n, 'degree', 'higher')/(ls/S))$degree))
    sd_vul <- append(sd_vul, sd((specieslevel(n, 'degree', 'lower')/(ls/S))$degree))
    
    cur_mod <- tryCatch({
      computeModules(n, forceLPA = T)@likelihood
    }, warning = function(w) {
      NA
    }, error = function(e) {
      NA
    }, finally = {
      
    })
    
    modularity <- append(modularity, cur_mod)
    
    #### This is for the degree distributions
    degs <- igraph::degree(regional_net, mode='in')
    degs <- degs[-which(degs == 0)]
    
    ##### normalisation term
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
  cur_out <- data.frame(replicate=r, areas=1:8, species, consumers, resources, links, connectances, links_per_sp, indegree, outdegree, modularity, normalised_indegree, normalised_outdegree, potential_indegree, sd_gen, sd_vul) 
  
  if(is.null(output)){
    output <- cur_out
  }else{
    output <- rbind(output, cur_out)
  }
}


##### After the runs are finished we store the output for the network properties across 
##### all areas (output) and the outcome of the analyses of the degree distributions (degree_dist_params)
##### into csv files that are readable for further analyses
write.csv(output, file='output-nahuel.csv')
write.csv(degree_dist_params, file='fits-degree-dists-nahuel.csv')






