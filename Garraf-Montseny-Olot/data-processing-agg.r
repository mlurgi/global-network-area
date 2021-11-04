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
## Dr Jordi Bosch
## Centre for Ecological Research and Forestry Applications
## Autonomous University of Barcelona, Catalonia
##
## ---------------------------

## This script reads a series of files from the home directory. Ensure these files are placed in the same
## directory from which this script is executed

## Load required libraries
## This script relies on the libjvm java library to read Microsoft Excel files using the rJava
## and XLConnect R packages. Please make sure the line below spcifies the right location for that
## library on your local computer
# dyn.load('/Library/Java/JavaVirtualMachines/jdk-9.jdk/Contents/Home/lib/server/libjvm.dylib')

library(rJava)
library(XLConnect)
require(bipartite)
require(igraph)
require(AICcmodavg)
require(viridis)


## This script was designed to read both the Garraf, Montseny, and Olot datasets
## as such, before running the script, one must specify the specific dataset to be read
## this is done by assigning to the following global variable the name of the desired 
## dataset, amongst those commented out next to the line

cur_dataset <- 'garraf-pp'   ## possible values: garraf-hp , garraf-pp , garraf-pp-2, montseny, olot

## The section below reads the original raw data
excel <- loadWorkbook(paste0("./raw-data/",cur_dataset,".xlsx"))

# get sheet names
sheet_names <- getSheets(excel)
names(sheet_names) <- sheet_names

# put sheets into a list of data frames
sheet_list <- lapply(sheet_names, function(.sheet){readWorksheet(object=excel, .sheet)})
# We remove the last sheet from the original file because it contains the metaweb 
# (hence the -1 in the line below). We also remove the first sheet (see index 2 below)

#### Garraf PP2 does not have meta info nor metaweb, so the first and last sheets in the excel file need
#### to be considered

n.sites <- length(sheet_list) - 1
first_index <- 2

if(cur_dataset == 'garraf-pp-2'){
  n.sites <- length(sheet_list)
  first_index <- 1
}

print(n.sites)

##### because there are some metrics (e.g. potential indegree) that require previous knowledge of the metaweb
##### we first go through all the datasets to create the metaweb
first <- TRUE
for(i in first_index:n.sites){
  n <- sheet_list[[i]]
  # row.names(n) <- n$Col1
  row.names(n) <- n[,1]
  n <- n[-1]
  n[n!=0] <- 1
  
  if(cur_dataset == 'garraf-pp-2'){
    n <- t(n)
  }
  
  if(cur_dataset == 'montseny'){
    rownames(n) <- paste0('Plant-', rownames(n))
    colnames(n) <- paste0('Pol-', colnames(n))
  }
  local_net <- graph_from_incidence_matrix(n, directed = T, mode='out')
  local_net <- igraph::delete.vertices(local_net, which(igraph::degree(local_net, V(local_net)) == 0))
  
  if(first){
    metaweb <- local_net 
    first <- FALSE
  }else{
    metaweb <- igraph::union(metaweb, local_net)
    V(metaweb)$type <- V(metaweb)$type_1
    V(metaweb)$type[which(!is.na(V(metaweb)$type_2))] <- V(metaweb)$type_2[which(! is.na(V(metaweb)$type_2))]
  }
}

#### Once the relevant raw data has been loaded and the metaweb has been constructed, 
#### we proceed with the network construction

#### These structures are used to keep the results of each replicate
replicates <- 100
output <- NULL
degree_dist_params <- NULL

null_model <- FALSE

for(r in 1:replicates){
  print(r)
  
  ####### Here, we initialise a series of vector data structures to keep track of the 
  ####### network properties that are calculated over the networks at each locality
  species <- c()
  species_tl1 <- c()
  species_tl2 <- c()
  links <- c()
  connectances <- c()
  links_per_sp <- c()
  indegree <- c()
  outdegree <- c()
  normalised_indegree <- c()
  normalised_outdegree <- c()
  modularity <- c()
  
  potential_indegree <- c()
  sd_gen <- c()
  
  ######## end of data structures
  
  output_for_plot <- NULL
  idxs <- sample(first_index:n.sites, length(first_index:n.sites))
  first <- TRUE
  
  print(idxs)
  
  ##### These are the colours for the degree distributions (see code below)
  # colors <- rainbow(length(idxs), alpha=0.7)
  colors <- plasma(length(idxs))
  cc <- 1
  counter <- 0
  
  for(i in idxs){
    ##### This section of the code reads the original raw data and builds a network
    ##### from it, storoing it in the 'local_net' igraph structure
    counter <- counter + 1
    n <- sheet_list[[i]]
    
    row.names(n) <- n[,1]
    n <- n[-1]
    n[n!=0] <- 1
    
    if(cur_dataset == 'garraf-pp-2'){
      n <- t(n)
    }
    
    if(cur_dataset == 'montseny'){
      rownames(n) <- paste0('Plant-', rownames(n))
      colnames(n) <- paste0('Pol-', colnames(n))
    }

    local_net <- graph_from_incidence_matrix(n, directed = T, mode='out')
    local_net <- igraph::delete.vertices(local_net, which(igraph::degree(local_net, V(local_net)) == 0))
    
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
      # # random_net <- igraph::delete.vertices(random_net, names(which(igraph::degree(random_net) == 0)))
      # # 
      # while(ecount(random_net) < 1){
      #   random_net <- induced_subgraph(metaweb, sample(V(metaweb)$name, vcount(regional_net_original)))
      #   # random_net <- igraph::delete.vertices(random_net, names(which(igraph::degree(random_net) == 0)))
      # }
      #  
      # random_net_bip <- graph.bipartite(bipartite.mapping(random_net)$type, as.vector(t(igraph::get.edges(random_net, 1:length(E(random_net))))), directed=T)
      # V(random_net_bip)$name <- V(random_net)$name
      
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
    
    ##### This is 'if' is for avoiding drawing the degree distributions of all areas
    ##### Only the first area and then areas which are multiples of 3 are drawn.
    if(counter == 1 | counter%%3 == 0){
    # if(TRUE){
      degs <- igraph::degree(regional_net, mode='in')
      degs <- degs[-which(degs == 0)]
  
      ###################### this is for extracting, drawing and analysing the degree distributions ############################
      norm_term <- ecount(regional_net)/vcount(regional_net)
  
      occur = as.vector(table(degs/norm_term))
      occur = occur/sum(occur)
      p = occur/sum(occur)
      y = rev(cumsum(rev(p)))
      x = as.numeric(names(table(degs/norm_term)))
  
      if(i == idxs[1]){
        plot(x,y, log='xy', ylim=c(10^-2.5,1), xlim=c(10^-0.3, 10^1), pch=16, cex=1.5, cex.lab=1.7, cex.axis=1.2, xlab='', ylab=expression(bold('Pc(k)')), col=colors[cc], tck=0, xaxt='n', yaxt='n', line=3)
        title(xlab=expression(bold("links (k)")), line=2, cex.lab=1.7, font=2)
        box(lwd=2.5)
        x1 <- floor(log10(range(x)))
        # pow <- seq(x1[1], x1[2]+1)
        pow <- c(-1,0,1)

        ticksat <- as.vector(sapply(pow, function(p) (1:10)*10^p))[2:20]
        axis(1, 10^pow[2:3], tck=0.02, lwd=1.5, cex.axis=1.8, lwd.ticks=1.5, labels=c(expression(10^{0}), expression(10^{1})), col.axis='#4d4d4d')
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
      cc <- cc + 3
      #### end of drawing code
      
      temp <- data.frame(x, y)
      
      if(is.null(output_for_plot)){
        output_for_plot <- cbind(temp, data.frame(area=counter))
      }else{
        output_for_plot <- rbind(output_for_plot, cbind(temp, data.frame(area=counter)))
      }
      
      failed <- tryCatch({
        mod1 <- nls(y ~ ( (x^-a) *(exp(-x/b))), data = temp, start = list(a = .001, b = 20), control=nls.control(maxiter = 1000))
        mod2 <- nls(y ~ (exp(-x/b)), data = temp, start = list(b = 1))
        mod3 <- nls(y ~ (x^-a), data = temp, start = list(a = 1))
        mod4 <- nls(y ~ ( (1/ (x * b * sqrt(2*pi) )) * exp(- ( ((log(x) - a)^2) / (2*(b^2)) )) ), data = temp, start = list(a = .3, b = .3), control=nls.control(maxiter = 1000))
        FALSE
      }, warning = function(w) {
        FALSE
      }, error = function(e) {
        TRUE
      }, finally = {
        
      })
      
      if(!failed){
        aic_comp <- aictab(list(mod1,mod2,mod3, mod4))
        model_name <- row.names(aic_comp)[1]
        model_name <- tolower(as.character(aic_comp$Modnames[1]))
        model <- summary(eval(as.symbol(model_name)))
      
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
    
    ##### Once we have the network corresponding to the current area we can calculate
    ##### a bunch of network properties. Those included in the paper and some more that
    ##### we also looked at
    n <- get.incidence(regional_net)
    S <- vcount(regional_net)
    Sc <- length(which(V(regional_net)$type == TRUE))
    Sr <- length(which(V(regional_net)$type == FALSE))

    potential_indegree <- append(potential_indegree, mean(igraph::degree(metaweb, V(regional_net)[which(V(regional_net)$type == TRUE)]$name, mode='in')))
    
    species <- append(species, S)
    species_tl1 <- append(species_tl1, Sr)
    species_tl2 <- append(species_tl2, Sc)
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
    
    sd_gen <- append(sd_gen, sd((specieslevel(n, 'degree', 'higher')/(ls/S))$degree))
    
    cur_mod <- tryCatch({
      computeModules(n, forceLPA = T)@likelihood
    }, warning = function(w) {
      NA
    }, error = function(e) {
      NA
    }, finally = {

    })

    modularity <- append(modularity, cur_mod)
  }
  
  ##### Lastly, we store all the information for the network properties for the current replicate
  ##### into the data frame containing the output data for all the runs
  cur_out <- data.frame(replicate=rep(r,length(species)), areas=1:length(species), species, species_tl1, species_tl2, links, connectances, links_per_sp, indegree, outdegree, modularity, normalised_indegree, normalised_outdegree, potential_indegree, sd_gen) 
  if(is.null(output)){
    output <- cur_out
  }else{
    output <- rbind(output, cur_out)
  }
  
}

##### After the runs are finished we store the output for the network properties across 
##### all areas (output) and the outcome of the analyses of the degree distributions (degree_dist_params)
##### into csv files that are readable for further analyses
write.csv(output, paste0("output-", cur_dataset,"-null-model-2.csv"))
write.csv(degree_dist_params, paste0("fits-degree-dists-",cur_dataset,"-null-model-2.csv"))
