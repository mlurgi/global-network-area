## ---------------------------
##
## Script name: utils.r
##
## Purpose of script: This script implements a series of functions that are used to quantify
## food web properties
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
## Date Created: September-2013
##
## Copyright (c) Miguel Lurgi, 2013-2020
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

## Load required libraries
library(hash)
library(igraph)
library(cheddar)


Generality <- function(M){
  return(sum(colSums(M))/sum((colSums(M)!=0)));
}

Vulnerability <- function(M){
  return(sum(rowSums(M))/sum((rowSums(M)!=0)));
}

NormalisedGenerality <- function(M){
  return((sum(colSums(M))/sum(colSums(M)!=0))/(sum(M)/dim(M)[1]))
}

NormalisedVulnerability <- function(M){
  return((sum(rowSums(M))/sum(rowSums(M)!=0))/(sum(M)/dim(M)[1]))
}

InDegree <- TrophicGenerality <- NumberOfResources <- function(M){
  return(colSums(M));
}

OutDegree <- TrophicVulnerability <- NumberOfCosumers <- function(M){
  return(rowSums(M));
}

Degree <- function(M){
  return(InDegree(M)+OutDegree(M));
}

MeanGenerality <- function(M){
  return(mean(colSums(M)));
}

MeanVulnerability <- function(M){
  return(mean(rowSums(M)));
}

SDGenerality <- function(M){
  return(sd(colSums(M)));
}

SDVulnerability <- function(M){
  return(sd(rowSums(M)));
}

FractionOfBasal <- function(M){
  M_temp <- M;
  diag(M_temp) <- 0;
  
  b_sps <- sum(which(InDegree(M_temp) == 0) %in% which(OutDegree(M_temp) >= 1));
  
  return(b_sps / dim(M)[1]);
}

NumberOfBasal <- function(M){
  M_temp <- M;
  diag(M_temp) <- 0;
  
  b_sps <- sum(which(InDegree(M_temp) == 0) %in% which(OutDegree(M_temp) >= 1));
  
  return(b_sps);
}

FractionOfTop <- function(M){
  M_temp <- M;
  diag(M_temp) <- 0;
  
  t_sps <- sum(which(InDegree(M_temp) >= 1) %in% which(OutDegree(M_temp) == 0));
  
  return(t_sps / dim(M)[1]);
}

NumberOfTop <- function(M){
  M_temp <- M;
  diag(M_temp) <- 0;
  
  t_sps <- sum(which(InDegree(M_temp) >= 1) %in% which(OutDegree(M_temp) == 0));
  
  return(t_sps);
}

FractionOfIntermediate <- function(M){
  M_temp <- M;
  diag(M_temp) <- 0;
  
  i_sps <- sum(which(InDegree(M_temp) >= 1) %in% which(OutDegree(M_temp) >= 1));
  
  return(i_sps / dim(M)[1]);
}

NumberOfIntermediate <- function(M){
  M_temp <- M;
  diag(M_temp) <- 0;
  
  i_sps <- sum(which(InDegree(M_temp) >= 1) %in% which(OutDegree(M_temp) >= 1));
  
  return(i_sps);
}

FractionOfCannibalism <- function(M){
  return(sum(diag(M)) / dim(M)[1]);
}

MaximumSimilarity <- function(M){
  S <- dim(M)[1];
  similarity <- 0;
  
  for(i in 1:S){
    for(j in i:S){
      if(i == j) next;
      similarity <- similarity + ((sum(M[,i] & M[,j]) + sum(M[i,] & M[j,])) / (sum(M[,i] | M[,j]) + sum(M[i,] | M[j,])) );
    }
  }
  
  return(similarity/S);
}

Omnivory <- function(M){
  #Use PredationMatrixToLinks() to create a Cheddar community from a predation
  # matrix
  M <- t(M)
  node <- 1:dim(M)[1];
  for(n in 1:length(node)){
    node[n] <- paste(node[n],'-');
  }
  
  pm <- matrix(M, ncol=dim(M)[2], dimnames=list(node, node), byrow=TRUE);
  community <- Community(nodes=data.frame(node=node), trophic.links=PredationMatrixToLinks(pm), properties=list(title='Community'));
  community <- RemoveCannibalisticLinks(community, title='community');
  
  #community is a cheddar community
  Fractionomnivory <- FractionOmnivorous(community)
  return(Fractionomnivory)
}

MeanFoodChainLength <- function(M){
  # Use PredationMatrixToLinks() to create a Cheddar community from a predation
  # matrix
  M <- t(M)
  node <- 1:dim(M)[1];
  for(n in 1:length(node)){
    node[n] <- paste(node[n],'-');
  }
  
  pm <- matrix(M, ncol=dim(M)[2], dimnames=list(node, node), byrow=TRUE);
  community <- Community(nodes=data.frame(node=node), trophic.links=PredationMatrixToLinks(pm), properties=list(title='Community'));
  community <- RemoveCannibalisticLinks(community, title='community');
  chain.stats <- TrophicChainsStats(community)
  ch_lens <- (chain.stats$chain.lengths + 1)
  
  return(sum(ch_lens)/length(ch_lens));
}


##### This function calculates the connectivity of the 
##### predator-predator undirected network in which a link
##### is shared between predators that share at least one prey
CalculatePredatorOverlap <- function(M){
  cols <- dim(M)[2];
  M1 <- matrix(0, cols, cols);
  for(r in 1:dim(M)[1]){
    links <- which(M[r,] == 1);
    for(l in links){
      M1[l, links] <- 1;
    }
  }
  
  M1[lower.tri(M1)] <- 0;
  diag(M1) <- 0;
  
  preds <- which((colSums(M1) + rowSums(M1)) != 0)
  M1 <- M1[preds, preds]
  cols <- dim(M1)[2];
  
  return((sum(M1)) / ((cols) * ((cols-1)/2)));
}



AveragePPMR <- function(FW){
  d <- dim(FW$M)[1];
  index <- which(FW$M == 1);
  size <- length(index);
  
  ratio <- 0;
  for(i in index){
    predator <- i%/%d + 1;
    prey <- i%%d;
    if(length(FW$BS[predator]) == 0 || length(FW$BS[prey]) == 0){
      size <- size - 1;
      next;
    }
    ratio <- ratio + log10(((FW$BS[predator])**0.75)/FW$BS[prey]);
  }
  return(ratio/size);
}

FindAllPaths <- function(graph_paths, start, end, path=c()){
  if(!(start %in% V(graph_paths))){
    return(list());
  }
  
  path <- append(path, start);
  
  if(start == end){
    return(list((path)));
  }
  
  paths <- list();
  for(node in neighbors(graph_paths, start, 'out')){  #we get the successors (predators) of the node 'start'
    if(!(node %in% path)){
      newpaths <- FindAllPaths(graph_paths, node, end, path);
      for(newpath in newpaths){
        paths[[length(paths)+1]] <- newpath;
      }
    }
  }
  
  return(paths);
}

FoodChainStats <- function(M){
  M <- t(M)
  node <- paste('S', 1:ncol(M));
  colnames(M) <- rownames(M) <- node;
  
  community <- Community(nodes=data.frame(node=node), trophic.links=PredationMatrixToLinks(M), properties=list(title='Community'));
  community <- RemoveCannibalisticLinks(community, title='community')
  
  cheddar.stats <- TrophicChainsStats(community);
  
  stats <- list();
  
  stats$overall_mean_length <- mean(cheddar.stats$chain.lengths);
  stats$overall_number_of_paths <- length(cheddar.stats$chain.lengths);
  
  stats$path_length_sd <- sd(cheddar.stats$chain.lengths);
  
  return(stats);
  
}

FoodChainStatsOld <- function(M){
  M_temp <- M;
  diag(M_temp) <- 0;
  graph_paths <- graph.adjacency(M_temp);
  
  stats <- list();
  
  stats$overall_mean_length <- 0;
  stats$overall_number_of_paths <- 0;
  
  stats$path_length_variance <- 0.0;
  stats$path_length_sd <- 0.0;
  stats$number_of_paths <- hash(V(graph_paths), 0);
  stats$trophic_positions <- hash(V(graph_paths), 0.0);
  
  if(!is.dag(graph_paths)){
    print('the graph is not a DAG');
    #return(stats);
  }
  
  top_order_nodes <- topological.sort(graph_paths, mode='out');
  
  ##if there are cycles, the topological sort won't include the nodes involved in cycles...
  ##so, we check what these nodes are:
  missing_nodes <- setdiff(V(graph_paths), top_order_nodes);
  ## and we append those at the end of the nodes array...
  
  if(length(missing_nodes) > 0){
    top_order_nodes <- append(top_order_nodes, missing_nodes);
  }
  
  paths <- hash();
  
  basals <- intersect(which(InDegree(M_temp) == 0) , which(OutDegree(M_temp) >= 1));
  
  for(n in basals){
    .set(stats$trophic_positions, n, 1.0);
    .set(stats$number_of_paths, n, 0);
    .set(paths, n, list(c(n)));
  }
  
  for(n in top_order_nodes){
    n_hash <- as.character(n);
    
    if(n %in% basals){
      next;
    }
    
    .set(paths, n, list());
    predecs <- neighbors(graph_paths, n, 'in');
    
    for(predecessor in predecs){
      predecessor_hash <- as.character(predecessor);
      
      if( length(paths[[predecessor_hash]]) == 0 ){
        
        print(paste('there are no paths for', predecessor));
        
        all_pred_paths <- list();
        for(b in basals){
          pred_paths <- FindAllPaths(graph_paths, b, predecessor);
          for(p in pred_paths){
            all_pred_paths[[length(all_pred_paths)+1]] <- p;
          }
        }
        .set(paths, predecessor, all_pred_paths);  
      }
      
      pred_paths <- paths[[predecessor_hash]];
      
      for(p in pred_paths){
        #with this we avoid counting cycles...
        if(n %in% p){
          next;
        }
        p_new <- append(p, n);
        paths[[n_hash]][[length(paths[[n_hash]])+1]] <- p_new;  
      }   
    }
    
    mean_length <- 0;
    .set(stats$number_of_paths, n, length(paths[[n_hash]]));
    
    if(stats$number_of_paths[[n_hash]] == 0){
      .set(stats$trophic_positions, n, 0);
      next;
    }
    
    for(path in paths[[n_hash]]){
      mean_length <- mean_length + length(path) - 1;
    }
    
    stats$overall_mean_length <- stats$overall_mean_length + mean_length;
    stats$overall_number_of_paths <- stats$overall_number_of_paths + stats$number_of_paths[[n_hash]];
    .set(stats$trophic_positions, n, (mean_length/stats$number_of_paths[[n_hash]]) + 1);
    
  }
  
  if(stats$overall_number_of_paths == 0){
    stats$overall_mean_length <- 0;
  }else{
    stats$overall_mean_length <- stats$overall_mean_length/stats$overall_number_of_paths;
  }
  
  for(sp in keys(paths)){
    if(as.integer(sp) %in% basals){
      next;
    }
    
    for(path in paths[[sp]]){
      stats$path_length_variance <- stats$path_length_variance + (((length(path)-1) - stats$overall_mean_length) ** 2);
    }   
  }
  
  if(stats$overall_number_of_paths == 0){
    stats$path_length_variance <- 0.0;
    stats$path_length_sd <- 0.0;
  }else{
    stats$path_length_variance <- stats$path_length_variance / stats$overall_number_of_paths;
    stats$path_length_sd <- sqrt(stats$path_length_variance);
  }
  
  print('done!');
  return(stats);
}


CalculateFoodWebStats <- function(FW){
  g <- graph.adjacency(FW$M);
  
  FW$stats$S <- FW$S;
  FW$stats$C <- sum(FW$M)/((dim(FW$M)[1])**2);
  
  FW$stats$LxS <- sum(FW$M)/FW$S;
  FW$stats$MeanPPMR <- AveragePPMR(FW);
  
  FW$stats$degree <- Degree(FW$M);
  FW$stats$in_deg <- InDegree(FW$M);
  FW$stats$out_deg <- OutDegree(FW$M);
  
  FW$stats$sd_gen <- SDGenerality(FW$M);
  FW$stats$sd_vul <- SDVulnerability(FW$M);
  
  FW$stats$B <- FractionOfBasal(FW$M);
  FW$stats$I <- FractionOfIntermediate(FW$M);
  FW$stats$T <- FractionOfTop(FW$M);
  
  FW$stats$Ca <- FractionOfCannibalism(FW$M);
  
  if(igraph::is.connected(g, mode='weak')){
    comm <- spinglass.community(g);
    FW$stats$M <- modularity(comm);
    FW$stats$c_member <- membership(comm);  
  }else{
    FW$stats$M <- -1.0;
  }
  
  FW$stats$MxSim <- MaximumSimilarity(FW$M);
  #T-B T-I I-I I-B
  
  FW$stats$dd <- degree.distribution(g, mode="total", cumulative=TRUE);
  FW$stats$in_dd <- degree.distribution(g, mode="in", cumulative=TRUE);
  FW$stats$out_dd <- degree.distribution(g, mode="out", cumulative=TRUE);
  return(FW);
}

CalculateCommunityStats <- function(FW){
  FW$comm_stats$total_biomass <- sum(FW$Final);
  
  sums_trs <- rep(0,4);
  shannon <- 0;
  
  for(s in 1:FW$S){
    fraction <- (FW$final[s]/FW$comm_stats$total_biomass);
    shannon <- shannon + (fraction * log(fraction));
    sums_trs[FW$TrLevels[s] + 1] <- sums_trs[FW$TrLevels[s] + 1] + FW$final[s];
  }
  
  FW$comm_stats$shannon <- -shannon;
  return(FW);
}


NormalizeMatrix <- function(M){
  colsum_M <- colSums(M);
  colsum_M[colsum_M==0] <- 1;
  return(t(t(M)/colsum_M));
  
}


#calculate the trophic positions using the prey-averaged trophic level metric from:
# Richard J. Williams and Neo D. Martinez (2004). American Naturalist.
# Limits to trophic levels and omnivory in complex food webs: theory and data.
  
TrophicPositions <- function(FW){
  S <- dim(FW$M)[1];
  M <- NormalizeMatrix(FW$M);
  if(det(diag(S) - t(M)) != 0){
    FW$TP <- solve(diag(S)-t(M), rep(1,S));
    
  }else{
    tmp <- diag(S);
    for(i in 1:9){
      tmp <- tmp %*% t(M) + diag(S);
    }
    FW$TP <- tmp %*% rep(1,S);
  }
  FW$W <- M;
  FW <- TrophicLevels(FW);
  return(FW);
}


#we can also obtain the trophic level of each species
#e.g. basal, herbivore, omnivore, top predator
TrophicLevels <- function(FW){
  S <- dim(FW$M)[1];
  TrLevels <- rep(-1, S)
  
  #we find the neighbours of each node in the network
  M_temp <- FW$M;
  diag(M_temp) <- 0;
  
  #all the species with trophic position 1 are basal species, and therefore
  #belong to the trophic level 0
  TrLevels[FW$TP<=1.0] <- 0;
  producers <- which(TrLevels == 0);
  
  for(i in 1:S){
    if(TrLevels[i] == -1){
      herb <- TRUE;
      top <- TRUE;
      
      #for each species we verify two things:
      #first: if any of its prey does not belong to trophic level 0,
      #then it is not a herbivore
      if(sum(TrLevels[which(M_temp[,i] != 0)]) != 0) herb <- FALSE;
      
      #second: if it has any predators, then it is not a top predator
      if(sum(M_temp[i,]) > 0) top <- FALSE;
      
      #after we've found out whether it is a herbivore or a top
      #predator, or none of those, assign the value accordingly
      if(herb){
        TrLevels[i] = 1;
      }else if(top){
        TrLevels[i] = 3;
      }else{
        TrLevels[i] = 2;
      }
      
    }
  }
  
  FW$TrLevels <- TrLevels;
  return(FW);
  
}
