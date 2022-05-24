
#========================================#
#      1. READ IN THE MANUAL LMS         #
#========================================#


library(Rvcg)
library(rgl)
library(Morpho)
library(rgl)
library(geomorph)
library(paleomorph)

###=========== LOADING DATA SET 1: WHOLE LANDMARKED SKULL 
#Read in landmarks manually placed on the whole skull 

ntaxa<-4 ## number of specimens (extant only) - NB can also put this in the code below (x,y,z) 
#data set .pts from Checkpoint

ptslist<-dir(pattern='.pts',recursive=T)
ptsarray<-array(dim=c(36,3,4)) #dim=c(number of landmarks and semilandmarks, number of dimensions, number of specimens)
for(i in 1:length(ptslist))
{
  ptsarray[,,i]<-as.matrix(read.table(file=ptslist[i],skip=2,header=F,sep="",row.names=1))
}

#Need the .plys for this 
#[3] stays the same 
dimnames(ptsarray)[3]<-list(
  substr(dir("./ply",pattern=".ply"),1,(nchar(dir("./ply",pattern=".ply"))-4)))
arraylm<-ptsarray #this is your array


##### MISSING LANDMARKS #########
#do this first and then rearrange landmarks 
#Landmarks need to be 'reshuffled' because we added 4 extra LMs on the nasal to be reflect that one midline landmark would 
#not work in the asymmetrical odontocetes 

arraylm[which(arraylm==9999)] <- NA
arraylm <- estimate.missing(arraylm,method="TPS")


#Check the LMs
text3d(arraylm[,,1], text = 1:36)

#Read in RHS and make a fake midline for mysticetes and any disarticulated archs 

LM1_bilat=arraylm[c(9,33),,] ## LM1 on left (5) and corresponding LM on the right side (69)
LM2_bilat=arraylm[c(10,34),,] ## LM2 on left (6) and corresponding LM on the right side (70)
LM3_bilat=arraylm[c(12,35),,] ## LM1 on left (5) and corresponding LM on the right side (69)
LM4_bilat=arraylm[c(13,36),,] 

#Check the dimensions  
#now want to find the midpoint of these pairs, 
#which is like finding the mean of each pair's position. 
#So to do this I summed each specimens x, then y then z coords, 
LM1_midline=colSums(LM1_bilat)/2
LM2_midline=colSums(LM2_bilat)/2
LM3_midline=colSums(LM3_bilat)/2
LM4_midline=colSums(LM4_bilat)/2


#Check

#Then just visually check that you're happy the code worked (IMPORTANT!):
## read in a specimen ply and plot these landmarks on it:
Pipa=ply2mesh(file="E:/Smithsonian Postdoc/Year 1/Mandible scans/LMs and curves/disarticulated/ply/Coronodon havensteini.ply")
shade3d(Pipa, col="white")
spheres3d(LM1_midline[,1]) #bracketed is specimen number
spheres3d(LM2_midline[,2])
spheres3d(LM3_midline[,3]) #bracketed is specimen number
spheres3d(LM4_midline[,4])
spheres3d(ptsarray[c(9,33),, 4], col = 'green') #check how these look with 5, 69, 6, 70
spheres3d(ptsarray[c(10,34),,4], col = 'red')
spheres3d(ptsarray[c(12,35),, 4], col = 'blue') #check how these look with 5, 69, 6, 70
spheres3d(ptsarray[c(13,36),,4], col = 'yellow')

#spheres3d(slidedlmsARCHS[c(58,112),,1], col = 'red')
#Then you add these two new midline landmarks to your landmark set:
Shape_data_with_bilats=abind::abind(arraylm, LM1_midline, LM2_midline, LM3_midline, LM4_midline, along=1)

#Check the midline LMs
spheres3d(Shape_data_with_bilats[c(1:36),,3], col = 'green', radius = 4)
spheres3d(Shape_data_with_bilats[c(37:40),,3], col = 'red', radius = 4)


#Delete superfluous LMs
arraylm_midline <- Shape_data_with_bilats[-c(9,10,12,13,33,34,35,36),,]

#Check 
text3d(arraylm_midline[c(1:32),,3], text = 1:32)




#Now rearrange so it's the same order as 

#You can use abind to bind datasets - this rearranges the data into the same order as the mysts (i.e. without the weird 1-120, 2 = 121, 19 - 122, 20 -123)
arranged_ptsarray=abind::abind(arraylm_midline[c(1:8),,], 
                             arraylm_midline[29:30,,],
                             arraylm_midline[11,,],
                             arraylm_midline[31:32,,],
                             arraylm_midline[c(14:32),,],
                             along= 1)

#Check 
text3d(arranged_ptsarray[c(1:32),,3], text = 1:32)

