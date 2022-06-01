
Library(geomorph)
Library(rgl)

#========================================#
#      1. READ IN THE MANUAL LMS         #  #Or read in our data set 'manual_skull_LMs.R' if you are not using your own data and skip this part
#========================================#


###=========== LOADING DATA SET 1: WHOLE LANDMARKED SKULL 
#Read in landmarks manually placed on the whole skull 

ntaxa <- 74 ## number of specimens (extant only) - NB can also put this in the code below (x,y,z) 
#data set .pts from Checkpoint

ptslist<-dir(pattern='.pts',recursive=T)
ptsarray<-array(dim=c(32,3,74)) #dim=c(number of landmarks and semilandmarks, number of dimensions, number of specimens)
for(i in 1:length(ptslist))
{
  ptsarray[,,i]<-as.matrix(read.table(file=ptslist[i],skip=2,header=F,sep="",row.names=1))
}

#Need the .plys for this 
#[3] stays the same 
dimnames(ptsarray)[3]<-list(
  substr(dir("./pts",pattern=".pts"),1,(nchar(dir("./pts",pattern=".pts"))-4)))
arraylm<-ptsarray #this is your array


##### MISSING LANDMARKS #########

arraylm[which(arraylm==9999)] <- NA
arraylm <- estimate.missing(arraylm,method="TPS")

#careful not to have missing data 
#check here and also check all of the LMs are numbered correctly 
text3d(arraylm[,,22], text=1:dim(arraylm)[1])

#This would be to change the names according to the species (as above with .ply)
#dimnames(arraylm)[3]=species

#let's call arraylm 'manual skull' to differentiate from the mirrored skull LMs
manual_skull <- arraylm

#############################
#                           #
#   PROCRUSTES THE DATA     #  Or read in our data set 'mirror_skull_LMs.R' if you are not using your own data and skip this part 
#                           #
#############################

manual_skull=gpagen(manual_skull) #Remove non-shape aspects 
manual_coords=manual_skull$coords #Subset out the coords 
size=manual_skull$Csize #centroid size
#PlotTangentSpace if you want to see a quick morphospace of the skulls 

#Proc_full <- manual_coords

PCA <- geomorph::gm.prcomp(manual_coords)
plot(PCA, axis1 = 1, axis2 = 2) #to plot 
summary(PCA) # take a look at the PCs (now called comp)


#to look at specimen position 
text(PCA$x[,1],PCA$x[,2], pos = 4)


#-- Plot with 3D data (example with PC1max)
plotRefToTarget( PCA$shapes$shapes.comp1$min, PCA$shapes$shapes.comp1$min,
                 method = "points",axes = F) 

spheres3d(PCA$shapes$shapes.comp1$min, radius= .0008,color = "grey")
spheres3d(PCA$shapes$shapes.comp1$max, radius= .001,color = "grey")
spheres3d(PCA$shapes$shapes.comp2$min, radius=.01,color = "grey")
spheres3d(PCA$shapes$shapes.comp2$max, radius=.01,color = "red")

text3d(final_procrusted_ARCHS[frontal,,1], text=frontal)


#Plot the LMs on the mesh 
#DO NOT USE PROCRUSTED DATA 

atarfa=ply2mesh(file="E:/Smithsonian Postdoc/Year 1/Mandible scans/full mandibles asymmetry test/ASYMM DATA/PTS full mand/ply/Delphinapterus leucas NHM 1933.10.13.4.ply")
shade3d(atarfa,col='white')
spheres3d(manual_coords[,,1], radius = 4, color = 'darkgreen')
#spheres3d(final_dataset[c(19:24, 85:90, 369:458, 1359:1448),,48], radius = 4, color = 'pink')



#Manual_skull_AB is the Procrusted coords dataset of manually placed landmarks
#AB denotes manually LMed skull, AC denotes computer mirrored skull

############################################################
#                                                          #
#     2. Load second dataset - half LM skull to mirror     #  Or read in our data set 'mirror_skull_LMs.R' if you are not using your own data and skip to line 114
#                                                          #                          
############################################################


#Manually cut previous data set (to computer mirror)

mirror_coords <- manual_coords[c(1:18),,]
spheres3d(mirror_coords[,,1], radius = 0.004, color = 'darkgreen')
text3d(mirror_coords[,,1], text = 1:18)
#spheres3d(final_dataset[c(19:24, 85:90, 369:458, 1359:1448),,48], radiu
#Check 



#MIRROR THESE LANDMARKS over the central line plane of the skull 

########## SYMMETRISATION TO IMPROVE THE SHAPE ANALYSES #########################


#Make a midline - from non procrusted data 

midline<-as.integer(c(9, 10, 12, 13)) # LM that are on the midline + parasphenoid curve points + NO patch point
#got length(midline)= 9 points on the midline

left.lm <- c(1:8,11,14:18)
#exclude midline points. Last number = last number of newpts 

lengmatrice=dim(mirror_coords)[1]*2-length(midline)#-length(nasalfrontal) #should be the length with the both sides, 1 is the column and 2 
#just means that we are duplicating the data to be on both sides of the skull 

Matrice=array(NA,dim = c(lengmatrice,3,74)) #3 is the dimensions (x, y, z), 2 is specimen number 
Matrice[1:dim(mirror_coords)[1],,]=mirror_coords

#left.lm <- c(1:37,39,41:47,50,52,53,57:60,62:66)
#left.lm <- c(2,3,5:18,21:37,39,41:47,50,52,53,57:60,62:66)
#exclude midline points. Last number = last number of newpts 

#Check left.lm and midline [left.lm,,x] = species number
spheres3d(mirror_coords[left.lm,,1],radius=0.004) #left LMs
spheres3d(mirror_coords[midline,,1],radius=0.004,col='red') #midline

right.lm <- c(19:32) #left.lm +1:lenmatrice

bilat.landmarks <- cbind(left.lm, right.lm) #one column is the rows of the right side LM and the other column the rows of the left side

MirroredAC=mirrorfill(A=Matrice,  l1=midline, l2=bilat.landmarks) # the NA rows are now filled so that the numbers are the same on both
#sides of the skull 
MirroredAC
#deformGrid3d(MirroredAC[67:123,,2], Matrice[,,2], ngrid=0) #This shows you the new mirroed landmarks 

#These visualisations are done before Procrusted data 
#This shows the original landmarks

atarfa=ply2mesh(file="E:/Smithsonian Postdoc/Year 1/Mandible scans/full mandibles asymmetry test/ASYMM DATA/PTS full mand/ply/Delphinapterus leucas NHM 1933.10.13.4.ply")
shade3d(atarfa,col='white')
spheres3d(MirroredAC[,,22],col= 'red', radius=0.004)
spheres3d(manual_coords[,,22], col = 'green' ,radius=0.004)
#check dimensions

#############################
#                           #
#   PROCRUSTES THE DATA     #
#                           #
#############################

MirroredAC=gpagen(MirroredAC) #Remove non-shape aspects 
Proc_mirrored=MirroredAC$coords #Subset out the coords 


##########################
#                        #
#    Asymmetry test      #
#                        #
##########################


#==========

if(!require(devtools)) install.packages("devtools")
library(devtools)
install_github("TGuillerme/landvR")
library(landvR)

############# LANDVR ###################
#Check the values
MirroredAC[1,,2] == manual_coords[1,,2]
# FALSE FALSE FALSE
## This is due to rounding, in fact they have the same 9 digits - round them 
round(MirroredAC[1,,2], digits = 9) == round(manual_coords[1,,2], digits = 9)
# TRUE TRUE TRUE

#I’ve updated landvR to version 0.3 where the coordinates.difference function now have a tolerance optional argument.
#You can use the following to get the 0 difference results:

differences_between_lms <- coordinates.difference(coordinates = MirroredAC[,,1],
                                                  reference = manual_coords[,,1],
                                                  type = "spherical",
                                                  rounding = 9)

#Remove errornous missing landmarks (these should be zero because they are static)
#differences_between_lms[[1]][1:18, 1:3] <- c(0.000000, 0.000000, 0.000000)


#Ellen's own colour function 
colfunc <- colorRampPalette(c("red", "yellow", "white"))
colfunc(10)
plot(rep(1,10),col=colfunc(10),pch=19,cex=3)

get.col.spectrum <- landvR::procrustes.var.plot(manual_coords[,,22], MirroredAC[,,22], col.val = differences_between_lms[[1]][,1], col = colfunc)

test=differences_between_lms[[1]][,1] #this is a test for specimen 1 to look at the differences between lms 
test

##### LOOKING AT AN AVERAGE SPECIMEN ######
N=32 #number of landmarks 
specs=74 #number of specimens 
all_combined=array(dim=c(N,3,specs)) #3 is the columns of data we need (radii, azimuth, polar)

i=1
for (i in 1:specs)
{
  all_differences <- coordinates.difference(coordinates = MirroredAC[,,i],
                                            reference = manual_coords[,,i],
                                            type = "spherical",
                                            rounding = 9)
  
  all_combined[,,i]=all_differences[[1]]
  
  i=i+1
}


#55, 56, 57, 59, 60 are all missing data and should be zero 
#all_combined[1:18, 1:3, 1:74] <- c(0.000000, 0.000000, 0.000000)
#write.csv(all_combined, file = 'all_combined.csv')

radii=all_combined[,1,] #looking at the second column (usually x,y,z) but here it is the radii, aziumuth, and polar 
radii_mean=apply(radii, c(1), mean) #c(1) looking at the first column which is the radii 
#test=all_combined[[1]][,,1] #this is a test for specimen 1 to look at the differences between lms 
#test

radii=all_combined[,1,] #second column of whole dataset with just the radii [,1,]


#Looking at the average radii compared to specimen 21 (or an average specimen)
get.col.spectrum <- landvR::procrustes.var.plot(manual_coords[,,9], MirroredAC[,,9], col.val = radii_mean, col = colfunc)
#datcol2<-c(rep("black",66),get.col.spectrum)
#open3d()


################
#              #
#   T-test     #
#              #
################


#Test difference between sides in entire dataset - is the variance on the right higher than on the left+midline ?
side_test <- t.test(manual_coords, MirroredAC)
side_test



#############################
#                           #
#     MORPHOSPACE           #
#                           #
#############################

#Load species data .csv
#This should be a .csb with the specimen name, family, age, and any other data you want to plot by
Specimen_data_ALL=read.csv("D:\Checkpoint - ICVM\Species_data_asymmetry", header=TRUE, sep=",")
Specimen_data_ALL <- read.csv('Species_data_asymmetry_ALL.csv')

#Morphospace
PCA=plotTangentSpace(arranged_data, axis1=1, axis2=2, label = Specimen_data_ALL$species)

gp <- plotTangentSpace(arranged_data, 
                       groups = as.factor(paste(Specimen_data_ALL$suborder))) 

#PCA now has PC scores and variations 
library(ggplot2) #plot PCA
library(viridis)
library(ggfortify) # For prcomp 
#Manual plotting of data 
#To change the aethetics 
myplot = plotTangentSpace(arranged_data)
attributes(myplot) # shows extractable parts
PC.scores = myplot$pc.scores
df <- PC.scores
autoplot(prcomp(df)) #a rough and ready plot 

#Plotting suborder 
a <- autoplot(prcomp(df), data = Specimen_data_ALL, colour = 'suborder', size = 4)
a = a + scale_colour_manual(values=c("#D55E00", "#0072B2", "#73BFB8", "black")) 
a = a + theme_light()
a = a + labs(x="PC1 (38.5%)",y="PC2 (23.7%)", axes=FALSE, cex.lab = 10)
a = a + theme(axis.title.x = element_text(size = rel(1.15))) 
a = a + theme(axis.title.y = element_text(size = rel(1.15))) 
a = a +labs(color="suborder") ## change legend title
a


#Now read in the whole dataset 
#With lands 1-66 removed (we just want to look at the radii difference for the mirrored landmarks i.e. 67-123)
#This is for a PCA of the radii

library(factoextra)
#for prcomp (PCA)

#Compute PCA 
#Load landmark data 
landmarks <- read.csv('radii_X_PCA.csv')

# IF ERRORS #
#landmarks[,6] <- as.numeric(as.character(landmarks[,6]))
#omit.na(landmarks)
#any(is.na(landmarks)) #Want this to be FALSE 
#landmarks[is.na(landmarks)] <- 0.015850205 #Or whatever the value is 

myPr <- prcomp(landmarks[1:174,9:65], scale = TRUE) #pull out the columns of data which are nummerical only

#If the above isn't working, it's likely there is a numeric/factor problem
landmarks <- as.numeric(landmarks$X67)

plot(myPr, type = 'b') #boxplot to look at variation 
biplot(myPr, scale = 0) #explore the data 

#extract PC scores 
str(myPr)
myPr$x

summary(myPr) #look at PC variations 

Landradii <- cbind(landmarks, myPr$x[,1:2]) #pull out PC1 and PC2

library(ggplot2) #plot PCA
library(viridis)
library(ggfortify)
library(RColorBrewer)

#plot
a <- ggplot(Landradii, aes(PC1, PC2, col = suborder, fill = suborder, labels = TRUE)) + 
  geom_point(shape = 21) + 
  xlab('PC 1 (XX.X%)') + 
  ylab('XX.X%)') + 
  geom_point(aes(size = sum.radii)) +
  theme_classic() +
  #geom_text(aes(label= ID), hjust=1, vjust=2) - with numbers of specimens 
  #geom_text_repel(aes(label = test2$ID)) - with numbers of specimens repelled 
  
  a 
#a + scale_x_reverse() #if you want to reverse axis 
#a + scale_y_reverse() #if you want to reverse axis

a + scale_colour_viridis_d() + 
  scale_x_reverse()#viridis

#variations of the morphospace 
a + scale_colour_viridis_d(option = "inferno") #plasma
#geom_text repel to repel text 

a + scale_colour_viridis_d(option = "plasma")
a + scale_colour_brewer(palette="YlOrRd")
a

library(vegan)
library(ggplot2)
library(ggConvexHull)
install.packages('ggConvexHull')

#Adding a convex hull 
ggplot(Landradii,aes(x = PC1, y = PC2,col = frequency)) +
  geom_convexhull(alpha = 0.3,aes(fill = frequency)) + 
  xlab('PC 1 (XX.X%)') + 
  ylab('PC2 (XX.X%)') +
  geom_point() +
  theme_minimal() + 
  scale_x_reverse()
