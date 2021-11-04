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
## Professor Tomas Roslin
## Swedish University of Agricultural Sciences
## Sweden
##
## Sections of the code below for reading raw data files were provided by Prof Roslin
## ---------------------------

## This script reads a series of files from the home directory. Ensure these files are placed in the same
## directory from which this script is executed

## Load required libraries

if(!require(magrittr)){ 
  install.packages("magrittr")
}

require(magrittr)
require(reshape2)
require(igraph)
require(bipartite)
require(AICcmodavg)
require(viridis)

## The section below reads the original raw data
####### THE FOLLOWING CODE WAS PROVIDED BY TOMAS ROSLIN TO READ DATA FILES #######
unlink("./raw-data/csv", recursive = TRUE) 
unlink("./raw-data/rdata", recursive = TRUE)

setwd('./raw-data')
source("../lib/format4R.r")
get_formatData("Salix_webs.csv")
setwd('..')

df_site <- readRDS("./raw-data/rdata/df_site.rds") 
str(df_site, strict.width="cut")

df_interact <- readRDS("./raw-data/rdata/df_interact.rds") 
str(df_interact, strict.width="cut")

site_interact <- merge(df_site, df_interact, by="REARING_NUMBER") 
head(site_interact)

#### Here we create the metawebs for the two networks
mweb_salgal <- df_interact[,c("RSAL","RGALLER")] %>% unique 
igr_salgal <- data.frame(
  from = mweb_salgal$RSAL,
  to = mweb_salgal$RGAL
) %>% graph_from_data_frame(directed=TRUE)
#

id <- df_interact$RPAR!="none"
mweb_galpar <- df_interact[id,c("RPAR","RGALLER")] %>% unique 
igr_salpar <- data.frame(
  from = mweb_galpar$RGAL,
  to = mweb_galpar$RPAR
) %>% graph_from_data_frame(directed=TRUE)

####### HERE ENDS THE CODE PROVIDED BY TOMAS ROSLIN #######


#### Once the relevant raw data has been loaded we proceed with the network construction
####### FROM NOW ON THE CODE THAT WE DID FOR NARS #######

###### this was added to be able to aggregate by gradient (north-south or south-north)
locations <- unique(site_interact[c('SITE','NDECDEG')])
north_south <- TRUE
locations <- locations[order(locations$NDECDEG, decreasing = north_south),]

random <- FALSE

## CHANGE THE VALUE FOR THIS VARIABLE TO SWITCH BETWEEN NETWORK TYPES
cur_dataset <- 'galpar'    #### Acceptable values: salgal  ,  galpar

if(cur_dataset == 'salgal'){
  metaweb <- graph_from_edgelist(as.matrix(mweb_salgal))
  metaweb_bip <- graph.bipartite(bipartite.mapping(metaweb)$type, as.vector(t(igraph::get.edges(metaweb, 1:length(E(metaweb))))), directed=T)
  V(metaweb_bip)$name <- V(metaweb)$name
  metaweb <- metaweb_bip
}else{
  metaweb <- graph_from_edgelist(as.matrix(mweb_galpar))
  metaweb_bip <- graph.bipartite(bipartite.mapping(metaweb)$type, as.vector(t(igraph::get.edges(metaweb, 1:length(E(metaweb))))), directed=T)
  V(metaweb_bip)$name <- V(metaweb)$name
  metaweb <- metaweb_bip
}

######## This code is to do the aggregation of the locations expanding from a randomly chosen
######## centre which is different for every replicate
replicates <- 1
output <- NULL
degree_dist_params <- NULL


null_model <- FALSE


for(rep in 1:replicates){
  locations <- locations[order(locations$NDECDEG, decreasing = north_south),]
  print(rep)
  if(random){
    locations <- locations[sample(dim(locations)[1]),]
  }else{
    locations$dist <- abs(locations$NDECDEG - sample(locations$NDECDEG, 1))
    locations <- locations[order(locations$dist),]
  }
  
  
  whole_net <- make_empty_graph(directed=T)
  
  ####### Here, we initialise a series of vector data structures to keep track of the 
  ####### network properties that are calculated over the networks at each locality
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
  normalised_outdegree <- c()
  
  sd_gen <- c()
  sd_vul <- c()
  potential_indegree <- c()
  
  ######## end of data structures
  
  counter <- 0
  # colors <- rev(rainbow(length(locations$SITE), alpha=0.7))
  colors <- plasma(length(locations$SITE))
  cc <- 1

  for(s in locations$SITE){
    ##### This section of the code reads the original raw data and builds a network
    ##### from it, storing it in the 'local_net' igraph structure
    print(counter)
    counter <- counter + 1
    cur_comm <- df_interact[which(site_interact$SITE == s),]
    
    if(cur_dataset == 'salgal'){
      bip_salgal <- cur_comm[,c("RSAL","RGALLER")] %>% table 
      int_mat <- as.matrix(bip_salgal)
    }
    else if(cur_dataset == 'galpar'){
      bip_galpar <- cur_comm[,c("RGALLER","RPAR")] %>% table
      int_mat <- as.matrix(bip_galpar)
    }
    
    new_col_names <- names(which(colSums(int_mat) != 0))
    new_row_names <- names(which(rowSums(int_mat) != 0))
    
    ##### this had to be done because there were networks with only one interaction
    if( length(new_row_names) == 1 ) {
      int_mat <- matrix((int_mat[which(rowSums(int_mat) != 0),which(colSums(int_mat) != 0)]), nrow=length(new_row_names), ncol=length(new_col_names))
      colnames(int_mat) <- new_col_names
      rownames(int_mat) <- new_row_names
    }else if(length(new_col_names) == 1) {
      int_mat <- matrix((int_mat[which(rowSums(int_mat) != 0),which(colSums(int_mat) != 0)]), nrow=length(new_row_names), ncol=length(new_col_names))
      colnames(int_mat) <- new_col_names
      rownames(int_mat) <- new_row_names
    }else{
      int_mat <- int_mat[which(rowSums(int_mat) != 0),which(colSums(int_mat) != 0)]
    }
    
    local_net <- graph_from_incidence_matrix(as.matrix(int_mat), directed = T, mode='out')
    
    if('none' %in% V(local_net)$name){
      local_net <- igraph::delete.vertices(local_net, 'none')
      
      if(length(which(igraph::degree(local_net) == 0)) != 0){
        local_net <- igraph::delete.vertices(local_net, names(which(igraph::degree(local_net) == 0)))
      }
    }
    
    ##### As areas get aggregated we store the regional (union) network in the whole_net
    ##### (and its bipartite representation g_bipart) network structure
    if(vcount(local_net) == 0 | ecount(local_net) == 0) next
    whole_net <- igraph::union(whole_net, local_net)
    
    g_bipart_original <- graph.bipartite(bipartite.mapping(whole_net)$type, as.vector(t(igraph::get.edges(whole_net, 1:length(E(whole_net))))), directed=T)
    V(g_bipart_original)$name <- V(whole_net)$name
    
    
    ##### This 'if' implements the two null models described in the paper.
    ##### Comment / uncomment each section appropriately
    if(null_model){
      ############### These lines are for the null models ######################
      # random_net <- induced_subgraph(metaweb, sample(V(metaweb)$name, vcount(g_bipart_original)))
      # 
      # while(ecount(random_net) < 1){
      #   random_net <- induced_subgraph(metaweb, sample(V(metaweb)$name, vcount(g_bipart_original)))
      # }
      # 
      # random_net_bip <- graph.bipartite(bipartite.mapping(random_net)$type, as.vector(t(igraph::get.edges(random_net, 1:length(E(random_net))))), directed=T)
      # V(random_net_bip)$name <- V(random_net)$name
      
      #### for the second null model uncomment the line below (and comment the lines above)
      failed <- TRUE
      while(failed){
        failed <- tryCatch({
          resource_species <- sample(vcount(g_bipart_original)-1,1)
          consumer_species <- vcount(g_bipart_original) - resource_species
          random_net_bip <- sample_bipartite( resource_species, consumer_species,  type='gnm', m=ecount(g_bipart_original), directed=TRUE)
          FALSE
        }, warning = function(w) {
          FALSE
        }, error = function(e) {
          TRUE
        }, finally = {

        })
      }
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
    normalised_outdegree <- append(normalised_outdegree, outdeg/(ls/S))
    
    pollinators <- names(which(igraph::degree(whole_net, mode='in') != 0))
    
    if(cur_dataset == 'salgal'){
      potential_indegree <- append(potential_indegree,  mean(igraph::degree(igr_salgal, pollinators, mode='in')))
    }else if(cur_dataset == 'galpar'){
      potential_indegree <- append(potential_indegree,  mean(igraph::degree(igr_salpar, pollinators, mode='in')))
    }
    
    sd_gen <- append(sd_gen, sd((specieslevel(int_mat, 'degree', 'higher')/(ls/S))$degree))
    sd_vul <- append(sd_vul, sd((specieslevel(int_mat, 'degree', 'lower')/(ls/S))$degree))
    
    cur_mod <- tryCatch({
      computeModules(int_mat, forceLPA = T)@likelihood
    }, warning = function(w) {
      NA
    }, error = function(e) {
      NA
    }, finally = {
      
    })
    
    modularity <- append(modularity, cur_mod)
    
    if(counter == 1 | counter%%37 == 0){
      degs <- igraph::degree(g_bipart, mode='in')
      degs <- degs[-which(degs == 0)]
      
      ###################### this is for the degree distributions
      norm_term <- ecount(g_bipart)/vcount(g_bipart)
      
      occur = as.vector(table(degs/norm_term))
      occur = occur/sum(occur)
      p = occur/sum(occur)
      y = rev(cumsum(rev(p)))
      x = as.numeric(names(table(degs/norm_term)))
      
      ##### comment out the following code to avoid drawing degree distributions
      if(counter == 1){
        plot(x,y, log='xy', ylim=c(10^-2.5,1), xlim=c(10^-.7, 10^1.7), pch=16, cex=1.5, cex.lab=1.7, cex.axis=1.2, xlab='', ylab=expression(bold('Pc(k)')), col=colors[cc], tck=0, xaxt='n', yaxt='n', line=3)
        title(xlab=expression(bold("links (k)")), line=2, cex.lab=1.7, font=2)
        box(lwd=2.5)
        x1 <- floor(log10(range(x)))
        # pow <- seq(x1[1], x1[2]+1)
        pow <- c(-1,0,1,2)
        
        ticksat <- as.vector(sapply(pow, function(p) (1:10)*10^p))
        axis(1, 10^pow[2:4], tck=0.02, lwd=1.5, cex.axis=1.8, lwd.ticks=1.5, labels=c(expression(10^{0}), expression(10^{1}), expression(10^{2})), col.axis='#4d4d4d')
        axis(1, ticksat, labels=NA, tcl=-0.25, lwd=0, lwd.ticks=1, tck=0.01, font=2)
        
        
        y1 <- floor(log10(range(y)))
        # pow <- seq(y1[1], y1[2])
        pow <- c(-3,-2,-1,0)
        ticksat <- as.vector(sapply(pow, function(p) (1:10)*10^p))
        ticksat <- ticksat[3:30]
        axis(2, 10^pow[2:4], tck=0.02, lwd=1.5, cex.axis=1.8, lwd.ticks=1.5, labels=c(expression(10^{-2}), expression(10^{-1}), expression(10^{0})), las=1, mgp=c(1,.5,0), col.axis='#4d4d4d')
        axis(2, ticksat, labels=NA, tcl=-0.25, lwd=0, lwd.ticks=1, tck=0.01)
  
      }else{
        points(x,y,pch=16, cex.lab=1.5, cex.axis=1.5, col=colors[cc], cex=1.5)
      }
  
      cc <- cc + 37 
    }
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
        cur_out <- data.frame(rep, area=counter, model=model_name, a=model$coefficients[1,1], a.std.err=model$coefficients[1,2], a.tval=model$coefficients[1,3], a.pval=model$coefficients[1,4], b=model$coefficients[2,1], b.std.err=model$coefficients[2,2], b.tval=model$coefficients[2,3], b.pval=model$coefficients[2,4])
      }else if(model_name == 'mod2'){
        cur_out <- data.frame(rep, area=counter, model=model_name, a=NA, a.std.err=NA, a.tval=NA, a.pval=NA, b=model$coefficients[1,1], b.std.err=model$coefficients[1,2], b.tval=model$coefficients[1,3], b.pval=model$coefficients[1,4])
      }else{
        cur_out <- data.frame(rep, area=counter, model=model_name, a=model$coefficients[1,1], a.std.err=model$coefficients[1,2], a.tval=model$coefficients[1,3], a.pval=model$coefficients[1,4], b=NA, b.std.err=NA, b.tval=NA, b.pval=NA)
      }
      
      if(is.null(degree_dist_params)){
        degree_dist_params <- cur_out
      }else{
        degree_dist_params <- rbind(degree_dist_params, cur_out)
      }
    }
    ######## end of degree distributions model fit and selection
  }
  
  
  ##### Lastly, we store all the information for the network properties for the current replicate
  ##### into the data frame containing the output data for all the runs
  cur_out <- data.frame(rep, areas=1:length(species), species, hosts, parasites, links, connectances, links_per_sp, indegree, outdegree, normalised_indegree,
                        potential_indegree, sd_gen, sd_vul, modularity)
  if(is.null(output)){
    output <- cur_out
  }else{
    output <- rbind(output, cur_out)
  }

}

##### After the runs are finished we store the output for the network properties across 
##### all areas (output) and the outcome of the analyses of the degree distributions (degree_dist_params)
##### into csv files that are readable for further analyses
write.csv(output, file = paste0('output-', cur_dataset, '.csv'))
write.csv(degree_dist_params, file = paste0('fits-degree-', cur_dataset, '.csv'))


