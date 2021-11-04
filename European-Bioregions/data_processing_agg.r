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
## For the specific case of this script the aggregation procedure includes a
## clustering algorithm developed by the author that merges cells on a raster map
## in an increasing, spiral-like, shape giving rise to patches of contiguous area.
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
## Dr Wilfried Thuiller
## Laboratory of Alpine Ecology
## Grenoble Alps University
## Grenoble, France
##
## ---------------------------

## This script reads a series of files from the home directory. Ensure these files are placed in the same
## directory from which this script is executed

## Load required libraries and required functions for the code below
require(raster)
require(rgdal)
require(igraph)
require(cheddar)
require(AICcmodavg)

## The script sourced below reads the original raw data
source('header.r')

#################################################
# This script is for running the network area 
# analysis per bioregion
#################################################

##### this is for building the network-area by bioregion

##### first the properties of the whole network at the bioregion level
##### next we get the network-area relationships per bioregion
europeRaster <- raster(x="./GIS-map/mask10k/reference_grid_10km.img")
cells_info <- read.dbf('./GIS-map/mask10k/reference_grid_10km.img.vat.dbf')

shape <- readOGR(dsn = "./GIS-map/BiogeoRegions2016_shapefile", layer = "BiogeoRegions2016")

region <- data.frame(PAGENAME = master$PAGENAME, SPP = 0)
regions_to_remove <- c('outside')

##### This code is for running several iterations in parallel
# library(foreach)
# library(doParallel)
# cl <- makeCluster(2)
# registerDoParallel(cl)

start_time <- proc.time()
agg_types <- c('spiral') #c('random', 'spiral')
output_bioregions_types <- NULL


start_time <- proc.time()

for(agg_type in agg_types){
  print(agg_type)
  replicates <- 100
  
  cur_dat <- NULL
  
  #### let's make it a bit fun and parallelise! 
  # cur_dat <- foreach(r = 1:replicates, .combine=rbind) %dopar%{
  for(r in 1:replicates){
    print(r)
    #### If you are paralellising, uncomment these library imports
    # require(igraph)
    # require(cheddar)
    # require(raster)
    # output <- NULL
    degree_dist_params <- NULL
    #### and here we extract the cell values
    for(i in shape@data$PK_UID){
      if(shape@data[which(shape@data$PK_UID == i),]$short_name %in% regions_to_remove) next;
      
      cur_pol <- shape[which(shape@data$PK_UID == i),]
      region_name <- as.character(shape@data[which(shape@data$PK_UID == i),]$code)
      cur_ids <- extract(europeRaster, cur_pol)[[1]]
      
      ##### here we fetch the cell codes from the database that matches the cell names on the raster
      cur_codes <- as.character(cells_info[which(cells_info$Value %in% cur_ids),]$PageName)
      
      ##### here we obtain the metaweb for the bioregion using the information of the cells
      meta_comm <- master[which(master$PAGENAME %in% cur_codes),]
      
      if(length(which(is.na(rowSums(meta_comm[-1])))) > 0){
        meta_comm <- meta_comm[-which(is.na(rowSums(meta_comm[-1]))),]  
      }
      
      meta_comm <- colSums(meta_comm[-1])
      meta_comm[meta_comm > 1] <- 1
      sum(meta_comm)
      # Defining the species pool at the bioregion level
      species_pool <- names(meta_comm)[which(meta_comm == 1)]
      species_pool <- setdiff(species_pool, absent_species)
      
      M <- BARMdiet.binary[species_pool,species_pool] # Web adjacency matrix (prey x predator)
      
      metaweb <- graph_from_adjacency_matrix(M)
      ######################### here ends the metaweb ##############################
      
      
      ##### this part of the code sets up different variables to keep track of the spatial aggregation
      ##### of cells in the raster map according to the aggregation heuristic to be followed: random or spiral
      
      if(agg_type == 'random') cur_codes <- sample(cur_codes)
      if(agg_type == 'random') j <- 1
      if(agg_type == 'spiral'){
        ### choose a random cell from the bioregion
        cur_cellid <- sample(cur_codes, 1)
        
        y_coord <- match(gsub("[0-9]", "", cur_cellid), latitudes);
        x_coord <- match(gsub("[A-Z]", "", cur_cellid), longitudes);
        
        x_index <- x_coord - 1;
        y_index <- y_coord - 1;
        x_count <- 1
        y_count <- 0
        x_next <- 1
        y_next <- 1
        x_offset <- -1
        y_offset <- 1
      }
      
      cells_to_traverse <- length(cur_codes) #rows*cols
      traversed_cells <- 1
      visited_cells <- c()
      traversed_set <- c()
      
      ################## end of aggregation data structures setup ##################
      
      ####### Here, we initialise a series of vector data structures to keep track of the 
      ####### network properties that are calculated over the networks at each locality
      species <- c()
      connectances <- c()
      links <- c()
      links_per_sp <- c()
      indegrees <- c()
      outdegrees <- c()
      basals <- c()
      tops <- c()
      intermediates <- c()
      sd_gen <- c()
      sd_vul <- c()
      S_tops <- c()
      S_intermediates <- c()
      S_basals <- c()
      modularities <- c()
      normalised_indegree <- c()
      normalised_outdegree <- c()
      consumers <- c()
      resources <- c()
      potential_indegree <- c()
      
      areas <- c()
      
      ######## end of data structures 
      
      print(region_name)
      print(paste0('cells to traverse = ', cells_to_traverse))
      
      sppSTRING_prev <- NULL
      ##### To ensure we visit all the cells in the current bioregion, we loop through
      ##### cells in sequence until the maximum number of 'cells to traverse' is reached
      while(length(visited_cells) < cells_to_traverse){
        if(agg_type == 'random'){
          cur_cellid <- cur_codes[j]
          j <- j + 1
        }else if(agg_type == 'spiral') {
          cur_cellid <- paste0(latitudes[y_coord],longitudes[x_coord])
        }
        
        traversed_set <- append(traversed_set, cur_cellid)
        skip <- FALSE
        
        ###### if the current cell has already been visited, we skip it
        if((! cur_cellid %in% cur_codes) | (cur_cellid %in% visited_cells)) skip <- TRUE
        if(! skip){
          if(is.na(sum(master[which(master$PAGENAME == cur_cellid),-1]))){
            skip <- TRUE
            cells_to_traverse <- cells_to_traverse - 1
          }
        }
        
        ###### if the current cell was not skipped we build the network for the current aggregation
        ###### and calculate all network properties
        if(! skip){
          visited_cells <- append(visited_cells, cur_cellid)
          
          #### This line extracts from the main database the species presence / absence information for
          #### all the cells corresponding to the current aggregation
          cur_comm <- master[which(master$PAGENAME %in% visited_cells),]
          
          cur_comm <- colSums(cur_comm[-1])
          cur_comm[cur_comm > 1] <- 1
          sum(cur_comm)
          # Defining the species pool
          sppSTRING <- names(cur_comm)[which(cur_comm == 1)]
          sppSTRING <- setdiff(sppSTRING, absent_species)
          
          ##### if the number of species is smaller than one we can't calculate any properties
          if(length(sppSTRING) <= 1){
            S <- length(sppSTRING)
            ls <- NA
            C <- NA
            l.s <- NA
            gen <- NA
            vul <- NA
            f_basal <- NA
            f_top <- NA
            f_int <- NA
            sd_g <- NA
            sd_v <- NA
            s_top <- NA
            s_int <- NA
            s_basal <- NA
            modularity <- NA
            normal_id <- NA
            normal_od <- NA
            
            cons <- NA
            resrcs <- NA
            
            pot_indeg <- NA
            
            sppSTRING_prev <- sppSTRING
          }else if(is.null(sppSTRING_prev) | length(setdiff(sppSTRING, sppSTRING_prev)) != 0){
            ##### to obtain the current network we sample the metaweb using the set of species
            ##### found in the current aggregation of cells
            M <- BARMdiet.binary[sppSTRING,sppSTRING] # Web adjacency matrix (predator x prey)
            
            ##### Once we have the network corresponding to the current area we can calculate
            ##### a bunch of network properties. Those included in the paper and some more that
            ##### we also looked at
            S <- length(sppSTRING)
            ls <- sum(M)
            C <- (2*ls/((S-1)*S))
            l.s <- ls/S
            gen <- Generality(M)
            vul <- Vulnerability(M)
            f_basal <- FractionOfBasal(M)
            f_top <- FractionOfTop(M)
            f_int <- FractionOfIntermediate(M)
            sd_g <- SDGenerality(M)/l.s
            sd_v <- SDVulnerability(M)/l.s
            s_top <- NumberOfTop(M)
            s_int <- NumberOfIntermediate(M)
            s_basal <- NumberOfBasal(M)
            
            ##### to calculate modularity we use the Walktrap modularity algorithm
            cur_net <- graph_from_adjacency_matrix(M)
            modularity <- tryCatch({
              max(walktrap.community(cur_net)$modularity)
            }, warning = function(w) {
              NA
            }, error = function(e) {
              NA
            }, finally = {
            })
            
            ##### normalisation of indegrees and outdegrees is done using the average links per species as normalisation term
            normal_id <- mean(degree(cur_net, V(cur_net), mode='in')[-which(degree(cur_net, V(cur_net), mode='in') == 0)]/l.s)
            normal_od <- mean(degree(cur_net, V(cur_net), mode='out')[-which(degree(cur_net, V(cur_net), mode='out') == 0)]/l.s)
            
            cons <- length(which(degree(cur_net, V(cur_net), mode='in') != 0))
            resrcs <- length(which(degree(cur_net, V(cur_net), mode='out') != 0))
            
            pot_indeg <- mean(degree(metaweb, names(which(degree(cur_net, mode='in') != 0)), mode='in'))
            sppSTRING_prev <- sppSTRING
            
            ########## this is for the degree distributions
            local_net <- graph_from_adjacency_matrix(M)
            degs <- igraph::degree(local_net, mode='all')
            if(length(which(degs == 0)) != 0){
              degs <- degs[-which(degs == 0)]
            }
            
            #### a normalisation term
            norm_term <- 1 #ecount(local_net)/vcount(local_net)
            
            occur = as.vector(table(degs/norm_term))
            occur = occur/sum(occur)
            p = occur/sum(occur)
            y = rev(cumsum(rev(p)))
            x = as.numeric(names(table(degs/norm_term)))
            temp <- data.frame(x,y)
            
            ###### To identify the shape of the degree distributions we fit four different models to it
            ###### the fitted models are: power law, exponential, truncated power law, and log-normal
            failed <- tryCatch({
              mod1 <- tryCatch({
                ##### Truncated power law
                nls(y ~ ( (x^-a) *(exp(-x/b))), data = temp, start = list(a = .1, b = 20), control=nls.control(maxiter = 1e3))
              }, error = function(e) {
                NA
              }, finally = {
              })
              mod2 <- tryCatch({
                ##### Exponential
                nls(y ~ (exp(-x/b)), data = temp, start = list(b = 35), control=nls.control(maxiter = 1000))
              }, error = function(e) {
                NA
              }, finally = {
              })
              mod3 <- tryCatch({
                ##### Power law
                nls(y ~ (x^-a), data = temp, start = list(a = .1), control=nls.control(maxiter = 1e3))
              }, error = function(e) {
                NA
              }, finally = {
              })
              mod4 <- tryCatch({
                ##### Log-normal
                nls(y ~ ( (1/ (x * b * sqrt(2*pi) )) * exp(- ( ((log(x) - a)^2) / (2*(b^2)) )) ), data = temp, start = list(a = 1.5, b = .3), control=nls.control(maxiter = 1000))
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
            if(! failed){
              model_list <- list(mod1,mod2,mod3,mod4)
              names_list <- c('mod1','mod2','mod3','mod4')
              if(length(which(is.na(model_list))) != 0){
                names_list <- names_list[-which(is.na(model_list))]
                model_list <- model_list[-which(is.na(model_list))]
              }
              names(model_list) <- names_list
              
              if(length(model_list) != 0){
                
                if(length(model_list) == 1){
                  model_name <- names_list[1]
                  model <- summary(eval(as.symbol(model_name)))
                }else{
                  aic_comp <- aictab(model_list)
                  model_name <- tolower(as.character(aic_comp$Modnames[1]))
                  model <- summary(eval(as.symbol(model_name)))
                }
                
                if(model_name == 'mod1' | model_name == 'mod4'){
                  cur_out <- data.frame(area=length(visited_cells), region=region_name,
                                        model=model_name, a=model$coefficients[1,1],
                                        a.std.err=model$coefficients[1,2], a.tval=model$coefficients[1,3],
                                        a.pval=model$coefficients[1,4], b=model$coefficients[2,1],
                                        b.std.err=model$coefficients[2,2], b.tval=model$coefficients[2,3],
                                        b.pval=model$coefficients[2,4])
                }else if(model_name == 'mod2'){
                  cur_out <- data.frame(area=length(visited_cells), region=region_name,
                                        model=model_name, a=NA, a.std.err=NA, a.tval=NA, a.pval=NA,
                                        b=model$coefficients[1,1], b.std.err=model$coefficients[1,2],
                                        b.tval=model$coefficients[1,3], b.pval=model$coefficients[1,4])
                }else{
                  cur_out <- data.frame(area=length(visited_cells), region=region_name,
                                        model=model_name, a=model$coefficients[1,1],
                                        a.std.err=model$coefficients[1,2], a.tval=model$coefficients[1,3],
                                        a.pval=model$coefficients[1,4], b=NA, b.std.err=NA, b.tval=NA,
                                        b.pval=NA)
                }
                
                if(is.null(degree_dist_params)){
                  degree_dist_params <- cur_out
                }else{
                  degree_dist_params <- rbind(degree_dist_params, cur_out)
                }
                
                ######## end of degree distributions model fit and selection
              }
            }
            
          }
          
          if(is.null(sppSTRING_prev))
            sppSTRING_prev <- sppSTRING
          
          ###### here we append the current network properties to the data structures (vectors)
          ###### that keep track of the changes in properties across areas
          species <- append(species, S)
          links <- append(links, ls)
          connectances <- append(connectances, C)
          links_per_sp <- append(links_per_sp, l.s)
          indegrees <- append(indegrees, gen)
          outdegrees <- append(outdegrees, vul)
          basals <- append(basals, f_basal)
          tops <- append(tops, f_top)
          intermediates <- append(intermediates, f_int)
          sd_gen <- append(sd_gen, sd_g)
          sd_vul <- append(sd_vul, sd_v)
          S_tops <- append(S_tops, s_top)
          S_intermediates <- append(S_intermediates, s_int)
          S_basals <- append(S_basals, s_basal)
          modularities <- append(modularities, modularity)
          normalised_indegree <- append(normalised_indegree, normal_id)
          normalised_outdegree <- append(normalised_outdegree, normal_od)
          
          consumers <- append(consumers, cons)
          resources <- append(resources, resrcs)
          
          potential_indegree <- append(potential_indegree, pot_indeg)
          
          areas <- append(areas, length(visited_cells))
        }
        
        ##### If we are performing spiral aggregation of cells across the map, the
        ##### following code performs a series of tricks to ensure the next cell
        ##### that we are going to look at is the correct one on the spiral order
        if(agg_type == 'spiral'){
          if(x_count > 0){
            x_index <- (x_index+x_offset) %% cols
            x_coord <- x_index + 1;
            x_count <- x_count - 1
            
            if(x_count == 0){
              x_next <- x_next + 1
              if(x_offset == 1) x_offset <- -1
              else if(x_offset == -1) x_offset <- 1
              y_count <- y_next
              
            }
            
          }else if(y_count > 0){
            y_index <- (y_index + y_offset) %% rows
            y_coord <- y_index + 1;
            y_count <- y_count - 1
            
            if(y_count == 0){
              y_next <- y_next + 1
              if(y_offset == 1) y_offset <- -1
              else if(y_offset == -1) y_offset <- 1
              x_count <- x_next
            }
          }
        }
        
        ###### end of spiral algorithm
        
        ###### we keep track of how many cells we have visited so far (for the while loop above)
        traversed_cells <- traversed_cells + 1
        print(length(visited_cells))
        
      }
      
      ##### Lastly, we store all the information for the network properties for the current replicate
      ##### into the data frame containing the output data for all the runs
      cur_out <- data.frame(region=rep(region_name, length(species)), area=1:length(species), species, links, links_per_sp, connectances,
                            indegrees, outdegrees, basals, tops, intermediates, sd_gen, sd_vul,
                            S_tops, S_intermediates, S_basals, modularities, normalised_indegree, normalised_outdegree, 
                            resources, consumers, potential_indegree)
      
      if(is.null(output)) output <- cur_out
      else output <- rbind(output, cur_out)
      
      #### write the output of the degree distribution fits to disk
      write.csv(degree_dist_params, file = paste0('fits-degree-dist-spiral-europe-rep-',r,'.csv')) 
    }
    
    #### write the output of the network properties for the current replicate to disk
    write.csv(output, file = paste0('output-bioregions-spiral-new-rep-',r,'.csv'))
    
    # cbind(data.frame(replicate=rep(r,dim(output)[1])),output)
    
    if(is.null(cur_dat))
      cur_dat <- cbind(data.frame(replicate=rep(r,dim(output)[1])),output)
    else
      cur_dat <- rbind(cur_dat, cbind(data.frame(replicate=rep(r,dim(output)[1])),output))
  }
  
  cur_dat <- cbind(data.frame(agg_type = rep(agg_type,dim(cur_dat)[1])), cur_dat)
  
  if(is.null(output_bioregions_types)) output_bioregions_types <- cur_dat
  else output_bioregions_types <- rbind(output_bioregions_types, cur_dat)
  
}

stop_time <- proc.time()
print(stop_time - start_time)
stopCluster(cl)

##### After the runs are finished we store the output for the network properties across 
##### all areas (output_bioregions_types) and the outcome of the analyses of the degree distributions (degree_dist_params)
##### into csv files that are readable for further analyses
write.csv(output_bioregions_types, file = 'output-european-bioregions.csv')
write.csv(degree_dist_params, file = 'fits-degree-dists-european-bioregions.csv')


