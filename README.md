# global-network-area

## Dr Miguel Lurgi - miguel.lurgi@swansea.ac.uk
Lecturer in Biosciences. Computational Ecology Lab, Department of Biosciences, Swansea University, UK
## 
## Dr Nuria Galiana - galiana.nuria@gmail.com
Centre for Biodiversity Theory and Modelling. Theoretical and Experimental Ecology Station, CNRS, France

Date Created: 25-05-2020 - Copyright (c) Miguel Lurgi, 2020

This project contains the data and source code used to produce the results
presented in our paper entitled:

'The spatial scaling of ecological networks across the globe' (Galiana, Lurgi et al. (2021))

In this README file we describe the contents of the archive (i.e. files and folders contained here) 
and provide a brief description of how the source code scripts should be executed to obtain the figures
and other results presented in the paper.

DATA FILES:
The archive is organised into a series of folders each containing the raw data and source code appropriate
for analysing that data set. Each data set is described in the Supplementary Material of the paper and the 
folder have been named according to the names of each data set as presented in the Supplementary Material.
Thus, for example, the folder title 'Sanak' contains the raw data for the spatial networks in the Sanak
dataset as well as the script used to process those data in order to construct networks of ecological
interactions at different spatial scales and quantify a suite of network properties at each of these scales.

Within each folder named after the corresponding dataset, a subfolder named 'raw-data' contains the actual
data as supplied by the different authors (as specified in the Supplementary Material of the paper). For
some datasets in which also a metadata table has been supplied, those metadata are outside of the raw-data
directory. In other instances, such as for example, the European Bioregions and Pyrenees datasets, further extra
files needed to analyse the data are also included. In these cases the GIS layers for maps onto which the
data is projected to obtain the desired spatial scales.

SOURCE CODE FILES:
1.- 'data_processing_agg.r': Within each dataset folder, a script called 'data_processing_agg.r' is provided. 
This script contains the code necessary to read and process the raw data for each dataset and construct 
ecological networks at different spatial scales. 
Each of these networks are then analysed using a collection of network metrics that are collated
into a final output file. The script is designed to run a spatial aggregation (i.e. spatial scaling) several
times (= 100 replicates). Aside from the network structure output, which is stored in output-XXXX.csv the 
script also fits several functions to the networks degree distributions and stores the outcome of those 
analyses in another output file: fits-degree-dists-XXXX.csv; where XXXX is the name of the corresponding dataset.

2.- 'data-analysis.r': This script, located within the parent folder, reads the output files from each of the 
datasets and produces the figures and performs the statistical analyses presented in the paper. It draws network-area
relationships and fits a battery of scaling functions to them in order to identify the scaling parameters
of each network property. Similar to 'data-analysis.r' we have 'data-analysis-null-model1.r' and 'data-analysis-null-model1.r' 
that read the output files for each null model from each dataset and performs the statistical analyses and the figures 
presented in the paper. The comparison between the null models and the original data can be done using 'comparison-null-models.r'.  
Finally, 'degree-dists-analysis.r', 'degree-dist-null1.r' and 'degree-dist-null2.r', analyse the output files of the degree distributions from each dataset for the original data and the null models, respectively.

3.- 'utils.r': An additional script containing a series of food web property measurement functions is provided
to aid with the calculations of network properties across datasets. This script is sourced from the 
'data_processing_agg.r' files above.

NOTE: notice that through the scripts the 'regional' domain category is called 'replicates'. 


HOW TO RUN:
1.- From within each dataset folder independently, run each one of the 'data_processing_agg.r' scripts 
located inside. Follow specific instructions (supplied as comments within the source code) when necessary.

2.- Executing the 'data_processing_agg.r' script should create two output files that are placed within that
folder: output-XXXX.csv and fits-degree-dists-XXXX.csv; where XXXX is replaced by the name of the 
corresponding dataset. These are the pre-processed data that are analysed afterwards to produce the final 
outputs. Make sure those files are created after the source code script has been executed.

3.- Go to the main (i.e. parent) folder and execute the script: 'data-analysis.r' (or any of the similar scripts, such as 
'data-analysis-null-model1.r'). In there you can see where each of the figures / results presented in the paper are generated.


Additionally, we provide .csv files with the outputs of the already processed data in case you want to directly produce 
the figures or perform the analyses without having to generate the network-area relationships. So 'merged_biogeographical-data.csv' 
and 'merged_datasets_replicates.csv' (and the respective null model .csv files) contain for all data sets the value of each network property at each spatial extent. Finally, 'merged_biogeographical-data-random.csv' contains the values of each network property at each spatial extent for the biogeographical data, where the different spatial units were randomly aggregated.   

I hope you enjoy the code!

If you have any questions / run into any issues, please contact us at miguel.lurgi@swansea.ac.uk or galiana.nuria@gmail.com







