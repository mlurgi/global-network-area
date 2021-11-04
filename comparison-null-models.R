#### This code is to compare the outputs from null model 1 and 2 with the outputs of the original data.


library(stringr)
library(dplyr)
require(purrr)
require(sars)
library(tidyr)

## Here we define a function used to extract normalised versions of the network properties
range_stats <- function(df){
  data.frame(areas = df$areas,
             species = (log10(df$species) - min(log10(df$species))),
             links = (log10(df$links) - min(log10(df$links))),
             links_per_sp = (log10(df$links_per_sp) - min(log10(df$links_per_sp))),
             indegree = (log10(df$indegree) - min(log10(df$indegree))),
             cr_ratio = (log10(df$cr_ratio) - max(log10(df$cr_ratio)))
  )
}

substract_max <- function(df){
  data.frame(areas = df[df$areas<max(df$areas),]$areas,
             species = df[df$areas<max(df$areas),]$species,
             links = df[df$areas<max(df$areas),]$links,
             links_per_sp = df[df$areas<max(df$areas),]$links_per_sp,
             indegree = df[df$areas<max(df$areas),]$indegree,
             cr_ratio = df[df$areas<max(df$areas),]$cr_ratio
  )
}




# We load the already processed data for the null models and the original data

all_datasets_original <- read.table("./merged_datasets_replicates.csv", sep=",", header=TRUE)
all_datasets_original$comparison <- 'original'
all_datasets_null1 <- read.table("./merged_datasets_replicates_nullmodel1.csv", sep=",", header=TRUE)
all_datasets_null1$comparison <- 'null1'
all_datasets_null2 <- read.table("./merged_datasets_replicates_nullmodel2.csv", sep=",", header=TRUE)
all_datasets_null2$comparison <- 'null2'

all_datasets <- rbind(all_datasets_original,all_datasets_null1,all_datasets_null2)

## Here we summarise and normalise the data
all_summaries_replicates <- all_datasets %>% 
  group_by(dataset, areas, comparison) %>% 
  dplyr::summarize(
    species = mean(species),
    links = mean(links),
    links_per_sp = mean(links_per_sp),
    indegree = mean(indegree),
    potential_indegree = mean(potential_indegree),
    cr_ratio = mean(cr_ratio)
  ) 

all_summaries_replicates_norm  <- all_summaries_replicates %>%
  split(list(.$dataset, .$comparison)) %>%
  map_dfr(range_stats, .id = c("dataset"))%>% 
  separate(dataset, into = c("dataset","comparison"), 
           sep = "\\.")

all_summaries_replicates_no_mean  <- all_datasets %>%
  split(list(.$dataset, .$comparison)) %>%
  map_dfr(substract_max, .id = c("dataset"))%>%
  separate(dataset, into = c("dataset","comparison"), 
           sep = "\\.")


# Calculating the 95% confidence interval

all_summaries_replicates_CI <- all_datasets_original %>% 
  group_by(dataset, areas) %>% 
  summarise(mean.sp = mean(species, na.rm = TRUE),
            sd.sp = sd(species, na.rm = TRUE),
            mean.lk = mean(links, na.rm = TRUE),
            sd.lk = sd(links, na.rm = TRUE),
            mean.lks = mean(links_per_sp, na.rm = TRUE),
            sd.lks = sd(links_per_sp, na.rm = TRUE),
            mean.in = mean(indegree, na.rm = TRUE),
            sd.in = sd(indegree, na.rm = TRUE),
            mean.ratio = mean(cr_ratio, na.rm = TRUE),
            sd.ratio = sd(cr_ratio, na.rm = TRUE)) %>%
  mutate(lower.ci.sp = mean.sp - 2*sd.sp,
         upper.ci.sp = mean.sp + 2*sd.sp,
         lower.ci.lk = mean.lk - 2*sd.lk,
         upper.ci.lk = mean.lk + 2*sd.lk,
         lower.ci.lks = mean.lks - 2*sd.lks,
         upper.ci.lks = mean.lks + 2*sd.lks,
         lower.ci.in = mean.in - 2*sd.in,
         upper.ci.in = mean.in + 2*sd.in,
         lower.ci.ratio = mean.ratio - 2*sd.ratio,
         upper.ci.ratio = mean.ratio + 2*sd.ratio)



### Checking how many of the null model observations fall within the 95% confident interval

output_CI_null1 <- NULL
output_CI_null2 <- NULL

for(d in unique(all_datasets$dataset)){
  tmp_data <- all_datasets[all_datasets$dataset==d,]
  for(a in unique(tmp_data$areas)){
    tmp_CI <- as.data.frame(all_summaries_replicates_CI[all_summaries_replicates_CI$areas==a & all_summaries_replicates_CI$dataset==d,])
    tmp_df_null1 <- tmp_data[tmp_data$areas==a & tmp_data$comparison=='null1',]
    tmp_df_null2 <- tmp_data[tmp_data$areas==a & tmp_data$comparison=='null2',]
    
    species_CI_null1 <- dim(tmp_df_null1[tmp_df_null1$species>=tmp_CI$lower.ci.sp & tmp_df_null1$species<=tmp_CI$upper.ci.sp,])[1]/length(tmp_df_null1$species)
    links_CI_null1 <- dim(tmp_df_null1[tmp_df_null1$links>=tmp_CI$lower.ci.lk & tmp_df_null1$links<=tmp_CI$upper.ci.lk,])[1]/length(tmp_df_null1$links)
    links_sp_CI_null1 <- dim(tmp_df_null1[tmp_df_null1$links_per_sp>=tmp_CI$lower.ci.lks & tmp_df_null1$links_per_sp<=tmp_CI$upper.ci.lks,])[1]/length(tmp_df_null1$links_per_sp)
    indegree_CI_null1 <- dim(tmp_df_null1[tmp_df_null1$indegree>=tmp_CI$lower.ci.in & tmp_df_null1$indegree<=tmp_CI$upper.ci.in,])[1]/length(tmp_df_null1$indegree)
    ratio_CI_null1 <- dim(tmp_df_null1[tmp_df_null1$cr_ratio>=tmp_CI$lower.ci.ratio & tmp_df_null1$cr_ratio<=tmp_CI$upper.ci.ratio,])[1]/length(tmp_df_null1$cr_ratio)
    
    species_CI_null2 <- dim(tmp_df_null2[tmp_df_null2$species>=tmp_CI$lower.ci.sp & tmp_df_null2$species<=tmp_CI$upper.ci.sp,])[1]/length(tmp_df_null2$species)
    links_CI_null2 <- dim(tmp_df_null2[tmp_df_null2$links>=tmp_CI$lower.ci.lk & tmp_df_null2$links<=tmp_CI$upper.ci.lk,])[1]/length(tmp_df_null2$links)
    links_sp_CI_null2 <- dim(tmp_df_null2[tmp_df_null2$links_per_sp>=tmp_CI$lower.ci.lks & tmp_df_null2$links_per_sp<=tmp_CI$upper.ci.lks,])[1]/length(tmp_df_null2$links_per_sp)
    indegree_CI_null2 <- dim(tmp_df_null2[tmp_df_null2$indegree>=tmp_CI$lower.ci.in & tmp_df_null2$indegree<=tmp_CI$upper.ci.in,])[1]/length(tmp_df_null2$indegree)
    ratio_CI_null2 <- dim(tmp_df_null2[tmp_df_null2$cr_ratio>=tmp_CI$lower.ci.ratio & tmp_df_null2$cr_ratio<=tmp_CI$upper.ci.ratio,])[1]/length(tmp_df_null2$cr_ratio)
    
    cur_output_null1<- data.frame(dataset=d, area=a, species=species_CI_null1, links=links_CI_null1, links_per_sp=links_sp_CI_null1, indegree=indegree_CI_null1, ratio=ratio_CI_null1)
    cur_output_null2<- data.frame(dataset=d, area=a, species=species_CI_null2, links=links_CI_null2, links_per_sp=links_sp_CI_null2,indegree=indegree_CI_null2, ratio=ratio_CI_null2)
    
    if(is.null(output_CI_null1)){
      output_CI_null1 <- cur_output_null1
    }else{
      output_CI_null1 <- rbind(output_CI_null1, cur_output_null1)
    }  
    if(is.null(output_CI_null2)){
      output_CI_null2 <- cur_output_null2
    }else{
      output_CI_null2 <- rbind(output_CI_null2, cur_output_null2)
    }  
  }
}


output_CI_mean_null1 <- output_CI_null1 %>% 
  group_by(dataset) %>% 
  dplyr::summarize(
    species = mean(species),
    links = mean(links),
    links_per_sp = mean(links_per_sp),
    indegree = mean(indegree),
    cr_ratio = mean(ratio)
  ) 


output_CI_mean_null2 <- output_CI_null2 %>% 
  group_by(dataset) %>% 
  dplyr::summarize(
    species = mean(species),
    links = mean(links),
    links_per_sp = mean(links_per_sp),
    indegree = mean(indegree),
    cr_ratio = mean(ratio)
  ) 

write.csv(as.data.frame(output_CI_mean_null1), file='confidence_intervals-replicates-nullmodel1.csv')
write.csv(as.data.frame(output_CI_mean_null2), file='confidence_intervals-replicates-nullmodel2.csv')


### Plotting the comparison between null models and the original data 
properties <-  c('species', 'links', 'links_per_sp', 'indegree', 'cr_ratio')

# Normalised plots
for(d in unique(all_summaries_replicates_norm$dataset)){
  tmp_data <- all_summaries_replicates_norm[all_summaries_replicates_norm$dataset==d,]
  for(p in properties){
    temp_plot <- ggplot(tmp_data, aes(x=areas, y=eval(parse(text=paste0('tmp_data$', p))), colour=comparison)) + 
      geom_line()+
      geom_point() + ggtitle(paste0(d,'-',p)) + ylab(p)
    ggsave(temp_plot, file=paste0("plot_", d, p,".png"), width = 14, height = 10, units = "cm")
  }
}

## For the non-normalised plots
for(d in unique(all_summaries_replicates_no_mean$dataset)){
  tmp_data <- all_summaries_replicates_no_mean[all_summaries_replicates_no_mean$dataset==d,]
  for(p in properties){
    temp_plot <- ggplot(tmp_data, aes(x=areas, y=eval(parse(text=paste0('tmp_data$', p))), colour=comparison)) + 
      geom_smooth() + geom_point()+ ggtitle(paste0(d,'-',p)) + ylab(p)
    ggsave(temp_plot, file=paste0("plot_", d, p,"no_norm_smooth.png"), width = 14, height = 10, units = "cm")
  }
}


### Links-species scaling ###
properties <-  c('links')

# Normalised plots
for(d in unique(all_summaries_replicates_norm$dataset)){
  tmp_data <- all_summaries_replicates_norm[all_summaries_replicates_norm$dataset==d,]
  for(p in properties){
    temp_plot <- ggplot(tmp_data, aes(x=species, y=eval(parse(text=paste0('tmp_data$', p))), colour=comparison)) + 
      geom_line()+
      geom_abline(slope = 1) + geom_abline(slope = 2)+
      geom_point() + ggtitle(paste0(d,'-',p)) + ylab(p)
    ggsave(temp_plot, file=paste0("plot_l.s.scaling", d,".png"), width = 14, height = 10, units = "cm")
  }
}

## For the non-normalised plots
for(d in unique(all_summaries_replicates_no_mean$dataset)){
  tmp_data <- all_summaries_replicates_no_mean[all_summaries_replicates_no_mean$dataset==d,]
  for(p in properties){
    temp_plot <- ggplot(tmp_data, aes(x=log10(species), y=log10(eval(parse(text=paste0('tmp_data$', p)))), colour=comparison)) + 
      geom_smooth(method = 'lm') + geom_point()+
      #geom_abline(slope = 1) + geom_abline(slope = 2)+
      ggtitle(paste0(d,'-',p)) + ylab(p)
    ggsave(temp_plot, file=paste0("plot_l.s.scaling", d,"no_norm_smooth.png"), width = 14, height = 10, units = "cm")
  }
}


#### BIOGEOGRAPHICAL ####


# We load the already processed data for the null models and the original data

all_datasets_bio_original <- read.table("./merged_biogeographical-data.csv", sep=",", header=TRUE)
all_datasets_bio_original$comparison <- 'original'
all_datasets_bio_null1 <- read.table("./merged_biogeographical-data_null1.csv", sep=",", header=TRUE)
all_datasets_bio_null1$comparison <- 'null1'
all_datasets_bio_null2 <- read.table("./merged_biogeographical-data_null2.csv", sep=",", header=TRUE)
all_datasets_bio_null2$comparison <- 'null2'

all_datasets_bio <- rbind(all_datasets_bio_original,all_datasets_bio_null1,all_datasets_bio_null2)
all_datasets_bio <- all_datasets_bio[complete.cases(all_datasets_bio),]



## Here we summarise and normalise the data
all_summaries_bio <- all_datasets_bio %>% 
  group_by(dataset, areas, comparison) %>% 
  dplyr::summarize(
    species = mean(species),
    links = mean(links),
    links_per_sp = mean(links_per_sp),
    indegree = mean(indegree),
    cr_ratio = mean(cr_ratio)
  ) 

all_summaries_bio_norm  <- all_summaries_bio %>%
  split(list(.$dataset, .$comparison)) %>%
  map_dfr(range_stats, .id = c("dataset"))%>% 
  separate(dataset, into = c("dataset","comparison"), 
           sep = "\\.")

all_summaries_bio_no_mean  <- all_datasets_bio %>%
  split(list(.$dataset, .$comparison)) %>%
  map_dfr(substract_max, .id = c("dataset"))%>%
  separate(dataset, into = c("dataset","comparison"), 
           sep = "\\.")




# Calculating the 95% confidence interval
all_summaries_bio_CI <- all_datasets_bio_original %>% 
  group_by(dataset, areas) %>% 
  summarise(mean.sp = mean(species, na.rm = TRUE),
            sd.sp = sd(species, na.rm = TRUE),
            mean.lk = mean(links, na.rm = TRUE),
            sd.lk = sd(links, na.rm = TRUE),
            mean.lks = mean(links_per_sp, na.rm = TRUE),
            sd.lks = sd(links_per_sp, na.rm = TRUE),
            mean.in = mean(indegree, na.rm = TRUE),
            sd.in = sd(indegree, na.rm = TRUE),
            mean.ratio = mean(cr_ratio, na.rm = TRUE),
            sd.ratio = sd(cr_ratio, na.rm = TRUE)) %>%
  mutate(lower.ci.sp = mean.sp - 2*sd.sp,
         upper.ci.sp = mean.sp + 2*sd.sp,
         lower.ci.lk = mean.lk - 2*sd.lk,
         upper.ci.lk = mean.lk + 2*sd.lk,
         lower.ci.lks = mean.lks - 2*sd.lks,
         upper.ci.lks = mean.lks + 2*sd.lks,
         lower.ci.in = mean.in - 2*sd.in,
         upper.ci.in = mean.in + 2*sd.in,
         lower.ci.ratio = mean.ratio - 2*sd.ratio,
         upper.ci.ratio = mean.ratio + 2*sd.ratio)



### Checking how many of the null model observations fall within the 95% confident interval

output_CI_null1 <- NULL
output_CI_null2 <- NULL

for(d in unique(all_datasets_bio$dataset)){
  tmp_data <- all_datasets_bio[all_datasets_bio$dataset==d,]
  for(a in unique(tmp_data[tmp_data$areas<=200,]$areas)){
    tmp_CI <- as.data.frame(all_summaries_bio_CI[all_summaries_bio_CI$areas==a & all_summaries_bio_CI$dataset==d,])
    tmp_df_null1 <- tmp_data[tmp_data$areas==a & tmp_data$comparison=='null1',]
    tmp_df_null2 <- tmp_data[tmp_data$areas==a & tmp_data$comparison=='null2',]
    
    species_CI_null1 <- dim(tmp_df_null1[tmp_df_null1$species>=tmp_CI$lower.ci.sp & tmp_df_null1$species<=tmp_CI$upper.ci.sp,])[1]/length(tmp_df_null1$species)
    links_CI_null1 <- dim(tmp_df_null1[tmp_df_null1$links>=tmp_CI$lower.ci.lk & tmp_df_null1$links<=tmp_CI$upper.ci.lk,])[1]/length(tmp_df_null1$links)
    links_sp_CI_null1 <- dim(tmp_df_null1[tmp_df_null1$links_per_sp>=tmp_CI$lower.ci.lks & tmp_df_null1$links_per_sp<=tmp_CI$upper.ci.lks,])[1]/length(tmp_df_null1$links_per_sp)
    indegree_CI_null1 <- dim(tmp_df_null1[tmp_df_null1$indegree>=tmp_CI$lower.ci.in & tmp_df_null1$indegree<=tmp_CI$upper.ci.in,])[1]/length(tmp_df_null1$indegree)
    ratio_CI_null1 <- dim(tmp_df_null1[tmp_df_null1$cr_ratio>=tmp_CI$lower.ci.ratio & tmp_df_null1$cr_ratio<=tmp_CI$upper.ci.ratio,])[1]/length(tmp_df_null1$cr_ratio)
    
    species_CI_null2 <- dim(tmp_df_null2[tmp_df_null2$species>=tmp_CI$lower.ci.sp & tmp_df_null2$species<=tmp_CI$upper.ci.sp,])[1]/length(tmp_df_null2$species)
    links_CI_null2 <- dim(tmp_df_null2[tmp_df_null2$links>=tmp_CI$lower.ci.lk & tmp_df_null2$links<=tmp_CI$upper.ci.lk,])[1]/length(tmp_df_null2$links)
    links_sp_CI_null2 <- dim(tmp_df_null2[tmp_df_null2$links_per_sp>=tmp_CI$lower.ci.lks & tmp_df_null2$links_per_sp<=tmp_CI$upper.ci.lks,])[1]/length(tmp_df_null2$links_per_sp)
    indegree_CI_null2 <- dim(tmp_df_null2[tmp_df_null2$indegree>=tmp_CI$lower.ci.in & tmp_df_null2$indegree<=tmp_CI$upper.ci.in,])[1]/length(tmp_df_null2$indegree)
    ratio_CI_null2 <- dim(tmp_df_null2[tmp_df_null2$cr_ratio>=tmp_CI$lower.ci.ratio & tmp_df_null2$cr_ratio<=tmp_CI$upper.ci.ratio,])[1]/length(tmp_df_null2$cr_ratio)
    
    cur_output_null1<- data.frame(dataset=d, area=a, species=species_CI_null1, links=links_CI_null1, links_per_sp=links_sp_CI_null1, indegree=indegree_CI_null1, ratio=ratio_CI_null1)
    cur_output_null2<- data.frame(dataset=d, area=a, species=species_CI_null2, links=links_CI_null2, links_per_sp=links_sp_CI_null2,indegree=indegree_CI_null2, ratio=ratio_CI_null2)
    
    if(is.null(output_CI_null1)){
      output_CI_null1 <- cur_output_null1
    }else{
      output_CI_null1 <- rbind(output_CI_null1, cur_output_null1)
    }  
    if(is.null(output_CI_null2)){
      output_CI_null2 <- cur_output_null2
    }else{
      output_CI_null2 <- rbind(output_CI_null2, cur_output_null2)
    }  
  }
}


output_CI_mean_null1 <- output_CI_null1 %>% 
  group_by(dataset) %>% 
  dplyr::summarize(
    species = mean(species),
    links = mean(links),
    links_per_sp = mean(links_per_sp),
    indegree = mean(indegree),
    cr_ratio = mean(ratio)
  ) 


output_CI_mean_null2 <- output_CI_null2 %>% 
  group_by(dataset) %>% 
  dplyr::summarize(
    species = mean(species),
    links = mean(links),
    links_per_sp = mean(links_per_sp),
    indegree = mean(indegree),
    cr_ratio = mean(ratio)
  ) 

write.csv(as.data.frame(output_CI_mean_null1[output_CI_mean_null1$dataset!='Macaronesia',]), file='confidence_intervals-bio-nullmodel1.csv')
write.csv(as.data.frame(output_CI_mean_null2[output_CI_mean_null1$dataset!='Macaronesia',]), file='confidence_intervals-bio-nullmodel2.csv')



### Plotting the comparison between null models and the original data 

properties <-  c('species','cr_ratio', 'links', 'links_per_sp', 'indegree')

# Normalised plots
for(d in unique(all_summaries_bio_norm$dataset)){
  tmp_data <- all_summaries_bio_norm[all_summaries_bio_norm$dataset==d,]
  for(p in properties){
    temp_plot <- ggplot(tmp_data, aes(x=areas, y=eval(parse(text=paste0('tmp_data$', p))), colour=comparison)) + 
      geom_line()+
      geom_point() + ggtitle(paste0(d,'-',p)) + ylab(p) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
         axis.line = element_line(colour = "black"), axis.text=element_text(size=18),
         axis.title=element_text(size= 22 ,face="bold"),
         legend.title = element_blank(),
         legend.box.background = element_rect(colour = "black", size = 1.5),
         #legend.title = element_text(size=12, face="bold"),
         #legend.text = element_text(size = 10, face = "bold"),
         panel.background = element_rect(fill = "white", color = 'black', size = 1.5))
    ggsave(temp_plot, file=paste0("plot_", d, p,".png"), width = 14, height = 10, units = "cm")
  }
}

## For the non-normalised plots
datasets <- c('Quercus')
for(d in datasets){
  tmp_data <- all_summaries_replicates_no_mean[all_summaries_replicates_no_mean$dataset==d,]
  for(p in properties){
    temp_plot <- ggplot(tmp_data, aes(x=areas, y=eval(parse(text=paste0('tmp_data$', p))), colour=comparison)) + 
      geom_smooth() + ggtitle(paste0(d,'-',p)) + ylab(p) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black"), axis.text=element_text(size=14),
            axis.title=element_text(size= 17 ,face="bold"),
            legend.title = element_blank(),
            legend.box.background = element_rect(colour = "black", size = 1.5),
            legend.position = 'none',
            #legend.title = element_text(size=12, face="bold"),
            #legend.text = element_text(size = 10, face = "bold"),
            panel.background = element_rect(fill = "white", color = 'black', size = 1.5))
    ggsave(temp_plot, file=paste0("plot_", d, p,"no_norm_smooth.png"), width = 10, height = 10, units = "cm")
  }
}


### Links-species scaling ###
properties <-  c('links')

# Normalised plots
for(d in datasets){
  tmp_data <- all_summaries_replicates_norm[all_summaries_replicates_norm$dataset==d,]
  for(p in properties){
    temp_plot <- ggplot(tmp_data, aes(x=species, y=eval(parse(text=paste0('tmp_data$', p))), colour=comparison)) + 
      geom_line()+
      geom_abline(slope = 1) + geom_abline(slope = 2)+
      geom_point() + ylab(p)+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black"), axis.text=element_text(size=14),
            axis.title=element_text(size= 17 ,face="bold"),
            legend.title = element_blank(),
            legend.box.background = element_rect(colour = "black", size = 1.5),
            #legend.position = 'none',
            #legend.title = element_text(size=12, face="bold"),
            #legend.text = element_text(size = 10, face = "bold"),
            panel.background = element_rect(fill = "white", color = 'black', size = 1.5))
    ggsave(temp_plot, file=paste0("plot_l.s.scaling", d,".png"), width = 12, height = 10, units = "cm")
  }
}

## For the non-normalised plots
for(d in unique(all_summaries_bio_no_mean$dataset)){
  tmp_data <- all_summaries_bio_no_mean[all_summaries_bio_no_mean$dataset==d,]
  for(p in properties){
    temp_plot <- ggplot(tmp_data, aes(x=log10(species), y=log10(eval(parse(text=paste0('tmp_data$', p)))), colour=comparison)) + 
      geom_smooth(method = 'lm') + geom_point()+
      #geom_abline(slope = 1) + geom_abline(slope = 2)+
      ggtitle(paste0(d,'-',p)) + ylab(p)
    ggsave(temp_plot, file=paste0("plot_l.s.scaling", d,"no_norm_smooth.png"), width = 14, height = 10, units = "cm")
  }
}

