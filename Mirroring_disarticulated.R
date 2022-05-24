
####################################################################
#                                                                  #
#     MIRROR YOUR SLID, RESAMPLED LMs                              #
#     NOTE: DIFFERENT CODE FOR SYMMETRIC AND ASYMMETRIC SPECIMENS  #
#                                                                  #
####################################################################

#========================================#
#      1. READ IN THE MANUAL LMS         #
#========================================#


###=========== LOADING DATA SET 1: WHOLE LANDMARKED SKULL 
#Read in landmarks manually placed on the whole skull 

ntaxa<-172 ## number of specimens (extant only) - NB can also put this in the code below (x,y,z) 
#data set .pts from Checkpoint

ptslist<-dir(pattern='.pts',recursive=T)
ptsarray<-array(dim=c(123,3,172)) #dim=c(number of landmarks and semilandmarks, number of dimensions, number of specimens)
for(i in 1:length(ptslist))
{
  ptsarray[,,i]<-as.matrix(read.table(file=ptslist[i],skip=2,header=F,sep="",row.names=1))
}

#Need the .plys for this 
#[3] stays the same 
dimnames(ptsarray)[3]<-list(
  substr(dir("./ply",pattern=".ply"),1,(nchar(dir("./ply",pattern=".ply"))-4)))
arraylm<-ptsarray #this is your array






#Read in RHS and make a fake midline for mysticetes and any disarticulated archs 

LM1_bilat=slidedlmsZARHI[c(5,69),,] ## LM1 on left (5) and corresponding LM on the right side (69)
LM2_bilat=slidedlmsZARHI[c(6,70),,] ## LM2 on left (6) and corresponding LM on the right side (70)

#Check the dimensions 
#now want to find the midpoint of these pairs, 
#which is like finding the mean of each pair's position. 
#So to do this I summed each specimens x, then y then z coords, 
LM1_midline=colSums(LM1_bilat)/2
LM2_midline=colSums(LM2_bilat)/2
#LM3_midline=colSums(LM3_bilat)/2


#Now you've created your 4 new mildine points, these should become 9, 10, 12, and 13 to correspond with odonts 












#called 'slidedlm_fake' because this data set still contains the additional midline LMs to help with alignment 
#removing the fake midline landmarks is carried out after mirroring (see below)

slidedlms_fake <- Shape_data_with_bilats
#slidedlms_fake[c(80:82,87:119),,]<-NA 

open3d();spheres3d(slidedlms_fake[,,1])
left.curves<-c(1:64)
left.lm <- c(1:37,39,41:47,50,52,53,57:60,62:66)
right.lm <- c(120, 121, 67:82, 122, 123, 83:119)

right.curves <- c(65:85)
left.curve.list<-unlist(my_curves$Curve.in[left.curves])
right.curve.list<-unlist(my_curves$Curve.in[right.curves])
leftside<-c(left.lm,left.curve.list) # RHS LMs+ RHS curves+all patch points
rightside<-c(right.lm, right.curve.list)
num.missing<-(length(leftside)-length(rightside)) # number of LMs to create= total RHS-current LHS LMs
blanks<-c((dim(slidedlms_fake)[1]+1):(dim(slidedlms_fake)[1]+num.missing))
# to fill in blanks from one row past the last current point, for the number of rows needed (num.missing)
rightside<-c(rightside,blanks)
add_col_or_row = function(x, n = 1, add_col = T, fill = 0)
{
  m1 = matrix(x, ncol = if(add_col) nrow(x) * ncol(x) else nrow(x), byrow = T)
  m2 = matrix(fill, nrow = if(add_col) dim(x)[3] else prod(dim(x)[-1]),
              ncol = if(add_col) nrow(x) * n else n)
  array(t(cbind(m1, m2)),
        c(nrow(x) + ((!add_col) * n), ncol(x) + (add_col * n), dim(x)[3]))
}
specimens2<-add_col_or_row(slidedlms_fake,n=num.missing,add_col=FALSE,fill=NA)
specimens2[c(80:82,87:119),,]<-NA
dimnames(specimens2)[3]<-dimnames(slidedlms_fake)[3]
#bilats <- read.csv('bilats.csv')
bilats<-cbind(leftside,rightside)
newarray<-mirrorfill(specimens2,l1=midline,l2=bilats)
dimnames(newarray)[3]<-dimnames(slidedlms)[3]
open3d();
spheres3d(newarray[,,3],radius=1.5)
spheres3d(newarray[bilats[,1],,1],col='red',radius=1.5)
spheres3d(newarray[bilats[,2],,1],col='blue',radius=1.5)
spheres3d(newarray[midline,,1], col = 'yellow', radius = 1.5)


midline<-as.integer(c(38,40,48,49,51,54,55,56,61,1449,1450)) 

spheres3d(final_mirrored_odonts[c(1:37,39,41:47,50,52,53,57:60, 62:66),,1],col='red',radius=3) #LHS landmarks - MIDLINE
spheres3d(final_mirrored_odonts[c(67:79,83:86,120:123),,1],col='red',radius=3) #RHS landmarks manually placed 
spheres3d(final_mirrored_odonts[c(80:82,87:119),,1],col='green',radius=3) #RHS landmarks mirrored 
spheres3d(final_mirrored_odonts[c(midline),,1],col='black',radius=3)
#spheres3d(newarray[c(116:119),,1],col='magenta',radius=1.5)
spheres3d(final_mirrored_odonts[c(124:1113),,1],col='green',radius=3) #LHS curves
spheres3d(final_mirrored_odonts[c(1114:1448),,1],col='brown',radius=3) #RHS curves manual 
spheres3d(final_mirrored_odonts[c(1449:1450),,1],col='yellow',radius=3) #fake midline LMs
spheres3d(final_mirrored_odonts[c(1451:2105),,1],col='green',radius=3) #RHS mirrored 


#Remove the extra fake landmarks first - numbers are different because of the extra curves on odonts 
#For odonts and archs 
final_mirrored_odonts=newarray[-c(1449:1450),,] 

#landmarks will now change 
spheres3d(final_mirrored_odonts[c(1:37,39,41:47,50,52,53,57:60, 62:66),,1],col='red',radius=1.5) #LHS 
spheres3d(final_mirrored_odonts[c(midline),,1],col='black',radius=1.5)
spheres3d(final_mirrored_odonts[c(67:123),,1],col='blue',radius=1.5) #RHS LMS
spheres3d(final_mirrored_odonts[c(124:1113),,1],col='green',radius=1.5) #LHS curves 
spheres3d(final_mirrored_odonts[c(1114:2103),,1],col='yellow',radius=1.5) #RHS curves

#Remove the double midline: curves 35, 45, 
final_mirrored_odonts=final_mirrored_odonts[-c(1629:1648, 1764:1783, 1824:1838, 1884:1903),,] #remove double midline 

#Final dataset 
spheres3d(final_mirrored_odonts[c(1:66),,2],col='red',radius=4)
spheres3d(final_mirrored_odonts[c(67:123),,2],col='red',radius=4)
spheres3d(final_mirrored_odonts[c(124:458),,2],col='blue',radius=4)
spheres3d(final_mirrored_odonts[c(459:1113),,2],col='yellow',radius=4)
spheres3d(final_mirrored_odonts[c(1114:1448),,2],col='blue',radius=4)
spheres3d(final_mirrored_odonts[c(1449:2028),,2],col='yellow',radius=4)

###################################################
#                                                 #
#  Make the same as mysticetes - correct order    #
#                                                 #
###################################################

#Slot in the landmarks 
#You can use abind to bind datasets - this rearranges the data into the same order as the mysts (i.e. without the weird 1-120, 2 = 121, 19 - 122, 20 -123)
arranged_odonts=abind::abind(final_mirrored_odonts[c(1:66),,], final_mirrored_odonts[120,,],
                             final_mirrored_odonts[c(67:68),,],
                             final_mirrored_odonts[121,,],
                             final_mirrored_odonts[c(69:82),,],
                             final_mirrored_odonts[122:123,,],final_mirrored_odonts[c(83:119),,],
                             final_mirrored_odonts[c(124:2028),,],
                             along= 1)
                             


#call it the final dataset again
final_mirrored_odonts <- arranged_odonts

#REMOVE THE EXTRA UNDESCRIBED SPECIMENS IF REQUIRED
#Save the undescribed specimens first 
final_mirrored_odonts_undescribed <- final_mirrored_odonts

save(final_mirrored_odonts_undescribed, file = 'final mirrored odonts undescribed2.R')


#Check
View(dimnames(final_mirrored_odonts)[[3]])
#Remove the undescribed 
final_mirrored_odonts <- final_mirrored_odonts[,,-c(84, 85, 86, 137)]     
#check 
View(dimnames(final_mirrored_odonts)[[3]])

save(final_mirrored_kogiids, file = 'final mirrored kogiids.R')

#Bind the 3 datasets 
final_dataset=abind::abind(final_mirrored_odonts, final_mirrored_mysts, final_mirrored_archs, along = 3)

#What order is this in? 
View(dimnames(final_dataset)[[3]])

#change the order to alphabetical for all 
final_dataset=final_dataset[,,sort(dimnames(final_dataset)[[3]])]
#check 
View(dimnames(final_dataset)[[3]])
save(final_dataset, file = 'final data set.R')
