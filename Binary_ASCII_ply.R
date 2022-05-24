
#script for converting binary to ASCII .ply files
#rm(list=ls())

install.packages("Rvcg")
install.packages("rgl")
library(Rvcg)
library(rgl)

#set working directory with ascii or binary ply files
#put the .ply files you want to change (for reading into R) in the .ply binary folder (not the .ply folder) 
setwd("X:xxxx/xxxx/ply binary")

#set output folder
#The extra forward slash at the end tells R that that is the output folder
outputfolder<-("X:xxxx/xxxx/ply/")

meshlist<-dir(pattern='.ply',recursive=F)

for (i in 1:length(meshlist)){
  x<-vcgImport(file=meshlist[i]) #import meshes
  
  vcgPlyWrite(x,filename=paste(outputfolder,meshlist[i],sep=""),binary=FALSE)#export
}

#Next chunk shouldn't be needed but just incase...
#Remove ply from the filename if needed 
filelist <- dir(pattern='.ply', recursive=F)
gsub("ply", "", filelist)

outputfolder<-("X:xxxx/xxxx/ply/")

#To visualise
#Set the directory to the one with the .ply in it 
checkLM(newpts,path="./ply/",suffix=".ply",pt.size=1,render="s",alpha=1)


#========================================================================================

#For taking a snapshot for publications 
#install.packages('Morpho')
library(Morpho)
library(rgl)

#don't forget .ply
Lissodelphis=ply2mesh(file="X:xxxxx/ply/Lissodelphis borealis USNM 550188.ply")

shade3d(Lissodelphis, col="white") #"white" or bone1 (don't use the "" for the latter) can use bone2 or bone3 also
rgl.snapshot(filename = "X:xxxxx/ply/Lissodelphis borealis USNM 550188.png") #the object name and then underscore followed by the new file name 


