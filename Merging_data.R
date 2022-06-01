#Read in normal cetacean array 
#Read in disarticualted array (see 'Mirroring_disarticulated.R')
#Make sure dims are the same

#Bind the two data sets 

#Bind the 3 datasets 
final_dataset=abind::abind(manual_skull, disarticulated_ptsarray, along = 3)

#What order is this in? 
View(dimnames(final_dataset)[[3]])

#change the order to alphabetical for all 
final_dataset=final_dataset[,,sort(dimnames(final_dataset)[[3]])]


#Read in the species data 

species_data_ALL=read.csv("E:/Smithsonian Postdoc/Year 1/Mandible scans/full mandibles asymmetry test/species_data_ALL.csv", header=TRUE, sep=",")

#pull the phylo names from the species data 
Full_names=species_data_ALL$binomial

#makes the species phylo.names the names of the array 
dimnames(final_dataset)[[3]]<-Full_names

