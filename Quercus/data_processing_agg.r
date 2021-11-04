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
## Authors: Dr Miguel Lurgi and Nuria Galiana

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
## Professor Tomas Roslin
## Swedish University of Agricultural Sciences
## Sweden
##
## ---------------------------

## This script reads a series of files from the home directory. Ensure these files are placed in the same
## directory from which this script is executed

## Load required libraries
## This script relies on the libjvm java library to read Microsoft Excel files using the rJava
## and XLConnect R packages. Please make sure the line below spcifies the right location for that
## library on your local computer
# dyn.load('/Library/Java/JavaVirtualMachines/jdk-9.jdk/Contents/Home/lib/server/libjvm.dylib')

dyn.load('/usr/lib/jvm/java-11-openjdk-amd64/lib/server/libjvm.so')

library(rJava)
library(XLConnect)
require(igraph)
require(bipartite)
require(AICcmodavg)

## The section below reads the original raw data
metadata <- read.csv('metadata.csv',  stringsAsFactors = F)
metadata <- metadata[1:22,1:4]
by_lat <- FALSE

metaweb <- make_empty_graph(directed=TRUE)
for(t in metadata$Tree){
  excel <- loadWorkbook(paste0("./raw-data/Web",t," 2006 Kaartinen.xlsx"))
  # get sheet names
  sheet_names <- getSheets(excel)
  names(sheet_names) <- sheet_names
  sheet_list <- lapply(sheet_names, function(.sheet){readWorksheet(object=excel, .sheet)})
  
  int_mat <- sheet_list$Sheet1
  
  row.names(int_mat) <- int_mat$Col1
  int_mat <- int_mat[-1]
  
  #### remove the column corresponding to unparasitised individuals
  if(length(which(colnames(int_mat) %in% c('Unpar', 'unpar'))) > 0){
    int_mat <- int_mat[-which(colnames(int_mat) %in% c('Unpar', 'unpar'))]
  }
  
  local_net <- graph_from_incidence_matrix(int_mat, directed = T, mode='out')
  metaweb <- igraph::union(metaweb, local_net)
}

#### Once the relevant raw data has been loaded we proceed with the network construction
#### Data structures to keep record of the network metrics across replicates
output <- NULL
replicates <- 100
degree_dist_params <- NULL
null_model <- FALSE

#### The algorithm below produces several different replicates of network aggregation of local
#### networks in an ever-increasing fashion to produce ecological networks at different spatial scales
for(r in 1:replicates){
  print(r)
  if(by_lat==TRUE){
    metadata <- metadata[order(metadata$North),]
  }else{
    n.sites <- length(metadata$North)
    idxs <- sample(1:n.sites, n.sites)
    metadata <- metadata[order(metadata$North[idxs]),]
  }
  
  ####### Here, we initialise a series of vector data structures to keep track of the 
  ####### network properties that are calculated over the networks at each locality
  whole_net <- make_empty_graph(directed=TRUE)
  species <- c()
  hosts <- c()
  parasites <- c()
  links <- c()
  connectances <- c()
  links_per_sp <- c()
  indegree <- c()
  outdegree <- c()
  modularity <- c()
  normalised_indegree <- c()
  sd_gen <- c()
  potential_indegree <- c()
  #### end of data structures
  
  counter <- 0
  colors <- rev(rainbow(length(metadata$Tree), alpha=0.7))
  cc <- 1
  for(t in metadata$Tree){
    ##### This section of the code reads the original raw data and builds a network
    ##### from it, storing it in the 'local_net' igraph structure
    counter <- counter + 1
    excel <- loadWorkbook(paste0("./raw-data/Web",t," 2006 Kaartinen.xlsx"))
    # get sheet names
    sheet_names <- getSheets(excel)
    names(sheet_names) <- sheet_names
    sheet_list <- lapply(sheet_names, function(.sheet){readWorksheet(object=excel, .sheet)})
    
    int_mat <- sheet_list$Sheet1
    
    row.names(int_mat) <- int_mat$Col1
    int_mat <- int_mat[-1]
    
    #### remove the columns corresponding to unparasitised individuals
    if(length(which(colnames(int_mat) %in% c('Unpar', 'unpar'))) > 0){
      int_mat <- int_mat[-which(colnames(int_mat) %in% c('Unpar', 'unpar'))]
    }
    
    local_net <- graph_from_incidence_matrix(int_mat, directed = T, mode='out')
    
    ##### As areas get aggregated we store the regional (union) network in the
    ##### whole_net and g_bipart_original network structures
    whole_net <- igraph::union(whole_net, local_net)
    whole_net <- igraph::delete.vertices(whole_net, names(which(igraph::degree(whole_net) == 0)))
    
    g_bipart_original <- graph.bipartite(bipartite.mapping(whole_net)$type, as.vector(t(igraph::get.edges(whole_net, 1:length(E(whole_net))))), directed=T)
    V(g_bipart_original)$name <- V(whole_net)$name
    
    
    ##### This 'if' implements the two null models described in the paper.
    ##### Comment / uncomment each section appropriately  
    if(null_model){
      ############### These lines are for the null models ######################
      random_net <- induced_subgraph(metaweb, sample(V(metaweb)$name, vcount(g_bipart_original)))
      # # random_net <- igraph::delete.vertices(random_net, names(which(igraph::degree(random_net) == 0)))
      # 
      while(ecount(random_net) < 1){
        random_net <- induced_subgraph(metaweb, sample(V(metaweb)$name, vcount(g_bipart_original)))
      #   # random_net <- igraph::delete.vertices(random_net, names(which(igraph::degree(random_net) == 0)))
      }
      # 
      random_net_bip <- graph.bipartite(bipartite.mapping(random_net)$type, as.vector(t(igraph::get.edges(random_net, 1:length(E(random_net))))), directed=T)
      V(random_net_bip)$name <- V(random_net)$name
      
      #### for the second null model uncomment the line below (and comment the lines above)
      # failed <- TRUE
      # while(failed){
      #   failed <- tryCatch({
      #     resource_species <- sample(vcount(g_bipart_original)-1,1)
      #     consumer_species <- vcount(g_bipart_original) - resource_species
      #     random_net_bip <- sample_bipartite( resource_species, consumer_species,  type='gnm', m=ecount(g_bipart_original), directed=TRUE)
      #     FALSE
      #   }, warning = function(w) {
      #     FALSE
      #   }, error = function(e) {
      #     TRUE
      #   }, finally = {
      # 
      #   })
      # }
      ##########################################################################
      
      g_bipart <- random_net_bip
    }else{
      g_bipart <- g_bipart_original
    }
    
    ##### Once we have the network corresponding to the current area we can calculate
    ##### a bunch of network properties. Those included in the paper and some more that
    ##### we also looked at
    int_mat <- as_incidence_matrix(g_bipart)
    S <- vcount(g_bipart)
    Sc <- length(which(V(g_bipart)$type == TRUE))
    Sr <- length(which(V(g_bipart)$type == FALSE))
    
    species <- append(species, S)
    hosts <- append(hosts, Sr)
    parasites <- append(parasites, Sc)
    
    ls <- ecount(g_bipart)
    links <- append(links, ls)
    connec <- ls / (Sr*Sc)
    connectances <- append(connectances, connec)
    links_per_sp <- append(links_per_sp, ls/S)
    
    indeg <- ls/Sc
    indegree <- append(indegree, indeg)
    outdeg <- ls/Sr
    outdegree <- append(outdegree, outdeg)
    
    normalised_indegree <- append(normalised_indegree, indeg/(ls/S))
    # normalised_outdegree <- append(normalised_outdegree, outdeg/(ls/S))
    
    pollinators <- names(which(igraph::degree(g_bipart, mode='in') != 0))
    potential_indegree <- append(potential_indegree,  mean(igraph::degree(metaweb, pollinators, mode='in')))
    
    sd_gen <- append(sd_gen, sd((specieslevel(int_mat, 'degree', 'higher')/(ls/S))$degree))
    
    cur_mod <- tryCatch({
      computeModules(int_mat, forceLPA = T)@likelihood
    }, warning = function(w) {
      NA
    }, error = function(e) {
      NA
    }, finally = {
      
    })
    
    modularity <- append(modularity, cur_mod)
    
    #### This is for the degree distributions
    if(TRUE){
      degs <- igraph::degree(g_bipart, mode='in')
      degs <- degs[-which(degs == 0)]

      #### normalisation term
      norm_term <- 1 # ecount(regional_net)/vcount(regional_net)

      occur = as.vector(table(degs/norm_term))
      occur = occur/sum(occur)
      p = occur/sum(occur)
      y = rev(cumsum(rev(p)))
      x = as.numeric(names(table(degs/norm_term)))

      ##### comment out the following code to avoid drawing degree distributions
      # if(cc == 1){
      #   plot(x,y, log='xy', ylim=c(10^-2.5,1), xlim=c(10^-0.9, 10^1), pch=16, cex=1.5, cex.lab=1.5, cex.axis=1.2, xlab='', ylab='log Pc(k)', col=colors[cc], tck=0, xaxt='n', yaxt='n')
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
          nls(y ~ ( (x^-a) *(exp(-x/b))), data = temp, start = list(a = 1, b = 2), control=nls.control(maxiter = 1e3))
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
  }

  ##### Lastly, we store all the information for the network properties for the current replicate
  ##### into the data frame containing the output data for all the replicates  
  cur_out <- data.frame(replicate=rep(r,length(species)), areas=1:length(species), species, hosts, parasites, links, connectances, links_per_sp, indegree, outdegree, modularity, normalised_indegree, sd_gen, potential_indegree)  #, indegree_normalized, outdegree, outdegree_normalized, sd_gen, sd_vul, basal, top, intermediate, S_basal, S_top, S_intermediate, omnivory, mfcls)
  
  if(is.null(output)){
    output <- cur_out
  }else{
    output <- rbind(output, cur_out)
  }  
} 

##### After the runs are finished we store the output for the network properties across 
##### all areas (output) and the outcome of the analyses of the degree distributions (degree_dist_params)
##### into csv files that are readable for further analyses
write.csv(output, "output-quercus.csv")
write.csv(degree_dist_params, "fits-degree-dists-quercus.csv")


require(ggplot2)
# raw_data <- read.csv('output-quercus.csv', stringsAsFactors = FALSE)
# output <- read.csv('output-quercus-null-model-1.csv', stringsAsFactors = FALSE)
ggplot(output, aes(areas, links)) + geom_point() +
  geom_point(data=raw_data, mapping=aes(areas, links, color='red'))


