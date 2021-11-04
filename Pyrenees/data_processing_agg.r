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
## Dr Bernat Claramunt Lopez
## Centre for Ecological Research and Forestry Applications
## Autonomous University of Barcelona
## Barcelona, Catalonia
##
## ---------------------------

## This script reads a series of files from the home directory. Ensure these files are placed in the same
## directory from which this script is executed

## Load required libraries and required functions for the code below
require(raster)
require(rgdal)
require(igraph)
require(stringr)
require(AICcmodavg)
source('../utils.r')

#### This function makes the translation between the notation of cell labels in 
#### the original layer to the raster naming convention used in the code
TransformCoords <- function(x){
  cur_let <- gsub("[[:digit:]]","", x)
  cur_num <- as.numeric(str_extract(x, "[0-9]+"))
  
  if(cur_num <= 10){
    first_lett <- 'C'
  }else{
    first_lett <- 'D'
  }
  if(cur_let %in% LETTERS[1:5]){
    sec_lett <- 'H'
  }else{
    sec_lett <- 'G'
  }
 
  actual_num <- (15-(match(cur_let, LETTERS)))%%10
  actual_num <- actual_num + ((cur_num-1)%%10)*10
  actual_num <- formatC(actual_num, width=2, flag=0)
  
  return(paste0(first_lett, sec_lett, actual_num))
}

## The section below reads the original raw data
#### First we look at the map to get an order for the cells to be aggregated by
shape <- readOGR(dsn = "./GIS-map", layer='quadPIRINEU')
r <- raster(ncol=17, nrow=10)
extent(r) <- extent(shape)
rp <- raster::rasterize(shape, r, shape$ATRIBUT)

order_in_raster <- getValues(rp)
order_in_raster <- order_in_raster[-which(is.na(order_in_raster))]
order_in_shape <- as.numeric(shape$ATRIBUT)

cols <- 17
rows <- 10

latitudes <- LETTERS[1:rows]
longitudes <- 1:cols

mapped_cells <- apply(expand.grid(latitudes, as.character(longitudes)), 1, paste0, collapse='')

#### this is the order of the cells in the map (from north-west to south-east)
north_south <- FALSE
if(north_south){
  cells_order <- (shape$ATRIBUT[match(order_in_raster, order_in_shape)])  
}else{      #### from south to north
  cells_order <- rev(shape$ATRIBUT[match(order_in_raster, order_in_shape)])  
}

#### Now that we know the order of the cells we traverse them and get the networks from
cells_order <- as.character(cells_order)

metaweb <- NULL
for(cell in (cells_order)){
  if(cell == 'CH04') next
  local_net <- read.graph(paste0('./raw-data/networks-per-quad/network-',cell,'.graphml.xml'), format = 'graphml')
  V(local_net)$name <- V(local_net)$id
  local_net <- igraph::delete.vertices(local_net, 'Invertebrates')
  print(cell)
  if(is.null(metaweb)){
    metaweb <- local_net
  }else{
    metaweb <- igraph::union(metaweb, local_net)
  }
}

metaweb <- igraph::delete.vertices(metaweb, names(which(igraph::degree(metaweb) == 0)))

#### Once the relevant raw data has been loaded we proceed with the network construction
#### Data structures to store the results of the replicates
output <- NULL
agg_type <- 'spiral'
degree_dist_params <- NULL
replicates <- 100
null_model <- FALSE

#### The algorithm below produces several different replicates of network aggregation of local
#### networks in an ever-increasing fashion to produce ecological networks at different spatial scales
for(replicate in 1:replicates){
  print(replicate)
  ####### Here, we initialise a series of vector data structures to keep track of the 
  ####### network properties that are calculated over the networks at each locality
  species <- c()
  connectances <- c()
  links <- c()
  links_per_sp <- c()
  indegrees <- c()
  outdegrees <- c()
  maxsims <- c()
  basals <- c()
  tops <- c()
  intermediates <- c()
  sd_gen <- c()
  sd_vul <- c()
  overlaps <- c()
  S_tops <- c()
  S_intermediates <- c()
  S_basals <- c()
  modularities <- c()
  normalised_indegree <- c()
  normalised_outdegree <- c()
  consumers <- c()
  resources <- c()
  potential_indegree <- c()
  ######## end of data structures
  
  output_for_plot <- NULL
  cells_to_traverse <- length(cells_order) - 1
  visited_cells <- c()
  regional_net_original <- NULL
  colors <- rev(rainbow(94, alpha=0.7))
  cc <- 1
  counter <- 0
  
  
  if(agg_type == 'random'){
    cur_codes <- sample(mapped_cells)
    j <- 1
  }else if(agg_type == 'spiral'){
    ### choose a random cell from the map
    cur_cellid <- sample(mapped_cells, 1)
    
    ### CH02 and CH72 have disconnected networks, so we avoid starting on that one
    while(!TransformCoords(cur_cellid) %in% cells_order | TransformCoords(cur_cellid) %in% c('CH02', 'CH72')){
      cur_cellid <- sample(mapped_cells, 1) 
    }
    
    #### This part of the code sets some variables needed for the spiral aggregation
    #### algorithm to keep track of the cells to visit
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
  
  while(length(visited_cells) < cells_to_traverse){
    if(agg_type == 'random'){
      cur_cellid <- cur_codes[j]
      j <- j + 1
    }else if(agg_type == 'spiral') {
      cur_cellid <- paste0(latitudes[y_coord],longitudes[x_coord])
    }
    
    cur_cellid <- TransformCoords(cur_cellid)
    if(cur_cellid %in% cells_order & cur_cellid != 'CH04' & !cur_cellid %in% visited_cells){
      visited_cells <- append(visited_cells, cur_cellid)
      
      # This is to plot the visited cells to track the progress
      # colors_shape <- rep('black', 94)
      # colors_shape[which(as.character(shape$ATRIBUT) %in% visited_cells)] <- 'red'
      # plot(shape, col=colors_shape)
      # 
      
      ##### This section of the code reads the original raw data and builds a network
      ##### from it, storing it in the 'local_net' igraph structure
      local_net <- read.graph(paste0('./raw-data/networks-per-quad/network-',cur_cellid,'.graphml.xml'), format = 'graphml')
      V(local_net)$name <- V(local_net)$id
      local_net <- igraph::delete.vertices(local_net, 'Invertebrates')
      
      if(length(which(igraph::degree(local_net) == 0)) != 0){
        local_net <- igraph::delete.vertices(local_net, names(which(igraph::degree(local_net) == 0)))
      }
      
      ##### As areas get aggregated we store the regional (union) network in the
      ##### regional_net_original network structure
      if(is.null(regional_net_original)){
        regional_net_original <- local_net
      }else{
        regional_net_original <- igraph::union(regional_net_original, local_net)
      }
      
      ##### This 'if' implements the two null models described in the paper.
      ##### Comment / uncomment each section appropriately  
      if(null_model){
        ############### These lines are for the null models ######################
        random_net <- induced_subgraph(metaweb, sample(V(metaweb)$name, vcount(regional_net_original)))
        # random_net <- igraph::delete.vertices(random_net, names(which(igraph::degree(random_net) == 0)))

        while(ecount(random_net) < 1){
          random_net <- induced_subgraph(metaweb, sample(V(metaweb)$name, vcount(regional_net_original)))
          # random_net <- igraph::delete.vertices(random_net, names(which(igraph::degree(random_net) == 0)))
        }
        
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
      M <- get.adjacency(regional_net)
      S <- vcount(regional_net)
      ls <- ecount(regional_net)
      C <- (2*ls/((S-1)*S))
      l.s <- ls/S
      gen <- Generality(M)
      vul <- Vulnerability(M)
      max_sim <- 0 #MaximumSimilarity(M)
      f_basal <- FractionOfBasal(M)
      f_top <- FractionOfTop(M)
      f_int <- FractionOfIntermediate(M)
      sd_g <- SDGenerality(M)/(l.s)
      sd_v <- SDVulnerability(M)/(l.s)
      pred_overlap <- CalculatePredatorOverlap(M)
      s_top <- NumberOfTop(M)
      s_int <- NumberOfIntermediate(M)
      s_basal <- NumberOfBasal(M)
      
      predators <- names(which(igraph::degree(regional_net, mode='in') != 0))
      pot_indeg <- mean(igraph::degree(metaweb, predators, mode='in'))

      modularity <- tryCatch({
        max(walktrap.community(regional_net)$modularity)
      }, warning = function(w) {
        NA
      }, error = function(e) {
        NA
      }, finally = {
      })
      
      normal_id <- mean(igraph::degree(regional_net, V(regional_net), mode='in')[-which(igraph::degree(regional_net, V(regional_net), mode='in') == 0)]/l.s)
      normal_od <- mean(igraph::degree(regional_net, V(regional_net), mode='out')[-which(igraph::degree(regional_net, V(regional_net), mode='out') == 0)]/l.s)
      
      cons <- length(which(igraph::degree(regional_net, mode='in') != 0))
      rsrcs <- length(which(igraph::degree(regional_net, mode='out') != 0))
      
      #### and we insert the calculated values into the corresponding vectors
      species <- append(species, S)
      links <- append(links, ls)
      connectances <- append(connectances, C)
      links_per_sp <- append(links_per_sp, l.s)
      indegrees <- append(indegrees, gen)
      outdegrees <- append(outdegrees, vul)
      maxsims <- append(maxsims, max_sim)
      basals <- append(basals, f_basal)
      tops <- append(tops, f_top)
      intermediates <- append(intermediates, f_int)
      sd_gen <- append(sd_gen, sd_g)
      sd_vul <- append(sd_vul, sd_v)
      overlaps <- append(overlaps, pred_overlap)
      S_tops <- append(S_tops, s_top)
      S_intermediates <- append(S_intermediates, s_int)
      S_basals <- append(S_basals, s_basal)
      modularities <- append(modularities, modularity)
      normalised_indegree <- append(normalised_indegree, normal_id)
      normalised_outdegree <- append(normalised_outdegree, normal_od)
      resources <- append(resources, rsrcs)
      consumers <- append(consumers, cons)
      
      potential_indegree <- append(potential_indegree, pot_indeg)
      degs <- igraph::degree(regional_net, mode='all')
      
      #### This is for the degree distributions
      norm_term <- ecount(regional_net)/vcount(regional_net)
      occur = as.vector(table(degs/norm_term))
      occur = occur/sum(occur)
      p = occur/sum(occur)
      y = rev(cumsum(rev(p)))
      x = as.numeric(names(table(degs/norm_term)))
      
      ##### comment out the following code to avoid drawing degree distributions
      # if(cc == 1){
      #   plot(x,y, log='xy', ylim=c(10^-3,1), xlim=c(10^-2, 10^2), pch=16, cex=1.5, cex.lab=1.5, cex.axis=1.2, xlab='', ylab='log Pc(k)', col=colors[cc], tck=0, xaxt='n', yaxt='n')
      #   title(xlab='log k', line=2, cex.lab=1.5)
      # 
      #   box(lwd=1.5)
      #   x1 <- floor(log10(range(x)+1))
      #   pow <- seq(x1[1], x1[2])
      # 
      #   ticksat <- as.vector(sapply(pow, function(p) (1:10)*10^p))
      #   axis(1, 10^pow, tck=0.02, lwd=1.5, lwd.ticks=1.5, cex.axis=1.5, labels=c(expression(10^{0}), expression(10^{1})))
      #   # axis(1, 10^pow, tck=0.02, lwd=1.5, lwd.ticks=1.5, labels=c(expression(10^{0}), expression(10^{1}), expression(10^{2})))
      #   axis(1, ticksat, labels=NA, tcl=-0.25, lwd=0, lwd.ticks=1, tck=0.01)
      # 
      # 
      #   y1 <- floor(log10(range(y)))
      #   pow <- seq(y1[1], y1[2])
      #   ticksat <- as.vector(sapply(pow, function(p) (1:10)*10^p))
      #   axis(2, 10^pow, tck=0.02, lwd=1.5, lwd.ticks=1.5, cex.axis=1.5, labels=c(expression(10^{-2}), expression(10^{-1}), expression(10^{0})), las=1, mgp=c(1,.5,0))
      #   axis(2, ticksat, labels=NA, tcl=-0.25, lwd=0, lwd.ticks=1, tck=0.01)
      # 
      # }else{
      #   points(x,y,pch=16, cex.lab=1.5, cex.axis=1.5, col=colors[cc], cex=1.5)
      # }

      cc <- cc + 1
      counter <- counter + 1
      
      #### end of drawing code
      
      ###### To identify the shape of the degree distributions we fit four different models to it
      ###### the fitted models are: power law, exponential, truncated power law, and log-normal
      temp <- data.frame(x, y)
      if(is.null(output_for_plot)){
        output_for_plot <- cbind(temp, data.frame(area=counter))
      }else{
        output_for_plot <- rbind(output_for_plot, cbind(temp, data.frame(area=counter)))
      }
      
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
          cur_out <- data.frame(replicate, area=length(visited_cells), model=model_name, a=model$coefficients[1,1], a.std.err=model$coefficients[1,2], a.tval=model$coefficients[1,3], a.pval=model$coefficients[1,4], b=model$coefficients[2,1], b.std.err=model$coefficients[2,2], b.tval=model$coefficients[2,3], b.pval=model$coefficients[2,4])
        }else if(model_name == 'mod2'){
          cur_out <- data.frame(replicate, area=length(visited_cells), model=model_name, a=NA, a.std.err=NA, a.tval=NA, a.pval=NA, b=model$coefficients[1,1], b.std.err=model$coefficients[1,2], b.tval=model$coefficients[1,3], b.pval=model$coefficients[1,4])
        }else{
          cur_out <- data.frame(replicate, area=length(visited_cells), model=model_name, a=model$coefficients[1,1], a.std.err=model$coefficients[1,2], a.tval=model$coefficients[1,3], a.pval=model$coefficients[1,4], b=NA, b.std.err=NA, b.tval=NA, b.pval=NA)
        }
        
        if(is.null(degree_dist_params)){
          degree_dist_params <- cur_out
        }else{
          degree_dist_params <- rbind(degree_dist_params, cur_out)
        }
      }
    }
    
    ##### This part implements the algorithm to figure out the next cell in the order
    ##### for the spatial aggregation
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
  }
  
  ##### Lastly, we store all the information for the network properties for the current replicate
  ##### into the data frame containing the output data for all the replicates  
  cur_out <- data.frame(replicate, area=1:length(species), species, links, connectances, links_per_sp, indegrees, outdegrees, maxsims, basals, tops,
                        intermediates, sd_gen, sd_vul, overlaps, S_tops, S_intermediates, S_basals,
                        modularities, normalised_indegree, normalised_outdegree, resources, consumers, potential_indegree)
  
  if(is.null(output)){
    output <- cur_out
  }else{
    output <- rbind(output, cur_out)
  }
  
}

##### After the runs are finished we store the output for the network properties across 
##### all areas (output) and the outcome of the analyses of the degree distributions (degree_dist_params)
##### into csv files that are readable for further analyses
write.csv(output, file='output-pyrenees.csv')
write.csv(degree_dist_params, file='fits-degree-dists.csv')

# 
# require(ggplot2)
# raw_data <- read.csv('output-pyrenees.csv', stringsAsFactors = FALSE)
# output <- read.csv('output-pyrenees-null-model-1.csv', stringsAsFactors = FALSE)
# ggplot(output, aes(area, links)) + geom_point() +
#   geom_point(data=raw_data, mapping=aes(area, links, color='red'))
# 


#### with this code we can plot the degree distributions

# require(ggplot2)
# require(RColorBrewer)
# getPalette = colorRampPalette(brewer.pal(9, "Set1"))
# 
# pd <- position_dodge(1.2)
# 
# png(filename = "degree_distributions_pyrenees.png", res= 300, width = 2000, height = 1600)
# ggplot(subset(output_for_plot, area %in% seq(3,93,18)), aes(x, y, colour=as.factor(area))) + 
#   geom_jitter(size=3, width = 0, height = 0.05) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         axis.line = element_line(colour = "black"), axis.text=element_text(size=22),
#         axis.title=element_text(size= 22 ,face="bold"),
#         legend.title = element_blank(),
#         legend.box.background = element_rect(colour = "black", size = 1.5),
#         #legend.title = element_text(size=12, face="bold"),
#         #legend.text = element_text(size = 10, face = "bold"),
#         panel.background = element_rect(fill = "white", color = 'black', size = 1.5)) +
#   scale_color_manual(name='Dataset', values=getPalette(10)) + 
#   ylab('Pc(k)')  + xlab('links (k)') + scale_x_log10(breaks=c(1,10), labels=c(expression(10^{0}), expression(10^{1}))) + scale_y_log10() +
#   annotation_logticks(mid=unit(0.1, "cm")) +
#   guides(color=guide_legend(override.aes=list(fill=NA))) 
# dev.off()
# 


##################### HERE ENDS THE SPIRAL AGGREGATION CODE ###########################
