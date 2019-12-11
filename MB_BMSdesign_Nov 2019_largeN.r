
library(MBHdesign)
library(spsurvey)
library(maptools)
library(rgeos)
library(raster)
library(car)
library(sf)
library(raster)
library(tidyverse)
library(sp)
library(rgdal)
library(lakemorpho)
library(dplyr)
library(car)
library(doBy)
library(gdata)
library(ggplot2)
library(nngeo)
library(reshape2)
library(viridis)
library(scales)


setwd("C:/Users/maderc/Documents/BMS/Manitoba_selection")
#load("C:\\Users\\vwilgs\\Documents\\BOREAL POSITION\\Projects\\Boreal Monitoring Strategy\\BMS design test\\Test_Design_largeN.RData")
#setwd("C:/Users/vwilgs/Documents/BOREAL POSITION/Projects/Saskatchewan Breeding Bird Atlas/Boreal Site Selection/")



###################################################################################################################################
##########################################          IMPORT THE SAMPLING FRAME            ##########################################
###################################################################################################################################


# Albers equal-area projection parameters based on above
#aea.proj <- "+proj=aea +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-96+x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m"
laea.proj <- "+proj=laea +lat_0=0 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs" # use this because it is lccc10 projection
NAD.proj <- "+proj=utm +zone=14 +ellps=GRS80 +datum=NAD83 +units=m +no_defs " #Manitoba hexagons are in this projection

hexagons<-readOGR("Mb_BOREAL_Hex_CM_UTM.shp", p4s= NAD.proj)
LCCC10 <- raster("C:/Users/maderc/Documents/BMS/Manitoba_selection/NA_LandCover_2010_25haMMU/NA_LandCover_2010_25haMMU.tif")
hexagons <- spTransform(hexagons, crs(LCCC10)) #transform to be in same projection as raster


#plot(hexagons, col=hexagons$Min_Cost) takes needless time
hexagons$cost_prob<-1/sqrt(hexagons$Min_Cost)


writeOGR(hexagons, "C:/Users/maderc/Documents/BMS/Manitoba_selection", "hexagons", driver = "ESRI Shapefile", overwrite = TRUE)

## now use the final_hexagons to clip lccc10 raster
# read in lccc10 file

###################------------    the following section is to create habitat selected MB.rast dataset.
##################------------     Can be removed after first successful creation of MB.rast after masking and removal of unwanted habitat classes
#################-------------     Temporarily replaced with the following load of dataset, because mask() takes around 45 minutes

MB.rast <- raster("C:/Users/maderc/Documents/BMS/Manitoba_selection/MB.rast.tif")
#plot(MB.rast)

#MB.rast <- crop(LCCC10,as(hexagons,"Spatial")) # LCCC raster for MB only
#MB.rast <- mask(MB.rast,hexagons,"Spatial", overwrite = TRUE)

# remove habitat classes not wanted in the model
#MB.rast <- calc(MB.rast, fun=function(x){x[x==19] <- NA; return(x)}) # remove snow and ice
#MB.rast <- calc(MB.rast, fun=function(x){x[x==18] <- NA; return(x)}) # water
#MB.rast <- calc(MB.rast, fun=function(x){x[x==17] <- NA; return(x)}) # urban
#MB.rast <- calc(MB.rast, fun=function(x){x[x==15] <- NA; return(x)}) # cropland
#MB.rast <- calc(MB.rast, fun=function(x){x[x==0] <- NA; return(x)}) # pixels not included in the cropping
#writeRaster(MB.rast,filename="MB.rast.tif",format="GTiff",overwrite=TRUE)

rm(LCCC10) # clear some memory

##############################################################################################################
############################ Now calculating Habitat Probabilities for ecoregion #############################


eco.list <- c("Aspen Parkland","Boreal Transition","Coastal Hudson Bay Lowland","Churchill River Upland",
              "Hudson Bay Lowland", "Hayes River Upland","Interlake Plain","Kazan River Upland", "Lake Manitoba Plain", 
              "Lac Seul Upland", "Lake of the Woods", "Mid-Boreal Lowland","Mid-Boreal Uplands","Maguse River Upland",
              "Selwyn Lake Upland" )

APhex<-subset(hexagons,REGION_NAM==eco.list[1])
BThex<-subset(hexagons,REGION_NAM==eco.list[2])
CHBLhex<-subset(hexagons,REGION_NAM==eco.list[3])
CRUhex<-subset(hexagons,REGION_NAM==eco.list[4])
HBLhex<-subset(hexagons,REGION_NAM==eco.list[5])
HRUhex<-subset(hexagons,REGION_NAM==eco.list[6])
ILPhex<-subset(hexagons,REGION_NAM==eco.list[7])
KRUhex<-subset(hexagons,REGION_NAM==eco.list[8])
LMPhex<-subset(hexagons,REGION_NAM==eco.list[9])
LSUhex<-subset(hexagons,REGION_NAM==eco.list[10])
LWhex<-subset(hexagons,REGION_NAM==eco.list[11])
MBLhex<-subset(hexagons,REGION_NAM==eco.list[12])
MBUhex<-subset(hexagons,REGION_NAM==eco.list[13])
MRUhex<-subset(hexagons,REGION_NAM==eco.list[14])
SLUhex<-subset(hexagons,REGION_NAM==eco.list[15])

AP<-crop(MB.rast,APhex)
AP <- mask(AP,APhex) # clip to the ecoregion
AP.freq <- as.data.frame(freq(AP)) # table of counts of pixels for each hab class in ecoregion
AP.freq <- subset(AP.freq,!value %in% NA)# remove NA from calculations
AP.freq$habprob <- 1/(nrow(AP.freq))/AP.freq$count # calculate probs for each class
AP.freq  

eco.old <- AP.freq[["value"]]
eco.new <- AP.freq[["habprob"]]
eco.mat <- matrix(c(rbind(eco.old,eco.new)),ncol=2,byrow=TRUE)
AP <- reclassify(AP,eco.mat) # new raster with pixel hab prob

BT<-crop(MB.rast,BThex)
BT <- mask(BT,BThex) # clip to the ecoregion
BT.freq <- as.data.frame(freq(BT)) # table of counts of pixels for each hab class in ecoregion
BT.freq <- subset(BT.freq,!value %in% NA)# remove NA from calculations
BT.freq$habprob <- 1/(nrow(BT.freq))/BT.freq$count # calculate probs for each class
BT.freq  
eco.old <- BT.freq[["value"]]
eco.new <- BT.freq[["habprob"]]
eco.mat <- matrix(c(rbind(eco.old,eco.new)),ncol=2,byrow=TRUE)
BT <- reclassify(BT,eco.mat) # new raster with pixel hab prob

CHBL<-crop(MB.rast,CHBLhex)
CHBL <- mask(CHBL,CHBLhex) # clip to the ecoregion
CHBL.freq <- as.data.frame(freq(CHBL)) # table of counts of pixels for each hab class in ecoregion
CHBL.freq <- subset(CHBL.freq,!value %in% NA)# remove NA from calculations
CHBL.freq$habprob <- 1/(nrow(CHBL.freq))/CHBL.freq$count # calculate probs for each class
CHBL.freq  

eco.old <- CHBL.freq[["value"]]
eco.new <- CHBL.freq[["habprob"]]
eco.mat <- matrix(c(rbind(eco.old,eco.new)),ncol=2,byrow=TRUE)
CHBL <- reclassify(CHBL,eco.mat) # new raster with pixel hab prob


CRU<-crop(MB.rast,CRUhex)
CRU <- mask(CRU,CRUhex) # clip to the ecoregion
CRU.freq <- as.data.frame(freq(CRU)) # table of counts of pixels for each hab class in ecoregion
CRU.freq <- subset(CRU.freq,!value %in% NA)# remove NA from calculations
CRU.freq$habprob <- 1/(nrow(CRU.freq))/CRU.freq$count # calculate probs for each class
CRU.freq
eco.old <- CRU.freq[["value"]]
eco.new <- CRU.freq[["habprob"]]
eco.mat <- matrix(c(rbind(eco.old,eco.new)),ncol=2,byrow=TRUE)
CRU <- reclassify(CRU,eco.mat) # new raster with pixel hab prob

HBL<-crop(MB.rast,HBLhex)
HBL <- mask(HBL,HBLhex) # clip to the ecoregion
HBL.freq <- as.data.frame(freq(HBL)) # table of counts of pixels for each hab class in ecoregion
HBL.freq <- subset(HBL.freq,!value %in% NA)# remove NA from calculations
HBL.freq$habprob <- 1/(nrow(HBL.freq))/HBL.freq$count # calculate probs for each class
HBL.freq  

eco.old <- HBL.freq[["value"]]
eco.new <- HBL.freq[["habprob"]]
eco.mat <- matrix(c(rbind(eco.old,eco.new)),ncol=2,byrow=TRUE)
HBL <- reclassify(HBL,eco.mat) # new raster with pixel hab prob

HRU<-crop(MB.rast,HRUhex)
HRU <- mask(HRU,HRUhex) # clip to the ecoregion
HRU.freq <- as.data.frame(freq(HRU)) # table of counts of pixels for each hab class in ecoregion
HRU.freq <- subset(HRU.freq,!value %in% NA)# remove NA from calculations
HRU.freq$habprob <- 1/(nrow(HRU.freq))/HRU.freq$count # calculate probs for each class
HRU.freq  

eco.old <- HRU.freq[["value"]]
eco.new <- HRU.freq[["habprob"]]
eco.mat <- matrix(c(rbind(eco.old,eco.new)),ncol=2,byrow=TRUE)
HRU <- reclassify(HRU,eco.mat) # new raster with pixel hab prob

ILP<-crop(MB.rast,ILPhex)
ILP <- mask(ILP,ILPhex) # clip to the ecoregion
ILP.freq <- as.data.frame(freq(ILP)) # table of counts of pixels for each hab class in ecoregion
ILP.freq <- subset(ILP.freq,!value %in% NA)# remove NA from calculations
ILP.freq$habprob <- 1/(nrow(ILP.freq))/ILP.freq$count # calculate probs for each class
ILP.freq  

eco.old <- ILP.freq[["value"]]
eco.new <- ILP.freq[["habprob"]]
eco.mat <- matrix(c(rbind(eco.old,eco.new)),ncol=2,byrow=TRUE)
ILP <- reclassify(ILP,eco.mat) # new raster with pixel hab prob

KRU<-crop(MB.rast,KRUhex)
KRU <- mask(KRU,KRUhex) # clip to the ecoregion
KRU.freq <- as.data.frame(freq(KRU)) # table of counts of pixels for each hab class in ecoregion
KRU.freq <- subset(KRU.freq,!value %in% NA)# remove NA from calculations
KRU.freq$habprob <- 1/(nrow(KRU.freq))/KRU.freq$count # calculate probs for each class
KRU.freq  

eco.old <- KRU.freq[["value"]]
eco.new <- KRU.freq[["habprob"]]
eco.mat <- matrix(c(rbind(eco.old,eco.new)),ncol=2,byrow=TRUE)
KRU <- reclassify(KRU,eco.mat) # new raster with pixel hab prob

LMP<-crop(MB.rast,LMPhex)
LMP <- mask(LMP,LMPhex) # clip to the ecoregion
LMP.freq <- as.data.frame(freq(LMP)) # table of counts of pixels for each hab class in ecoregion
LMP.freq <- subset(LMP.freq,!value %in% NA)# remove NA from calculations
LMP.freq$habprob <- 1/(nrow(LMP.freq))/LMP.freq$count # calculate probs for each class
LMP.freq  

eco.old <- LMP.freq[["value"]]
eco.new <- LMP.freq[["habprob"]]
eco.mat <- matrix(c(rbind(eco.old,eco.new)),ncol=2,byrow=TRUE)
LMP <- reclassify(LMP,eco.mat) # new raster with pixel hab prob

LSU<-crop(MB.rast,LSUhex)
LSU <- mask(LSU,LSUhex) # clip to the ecoregion
LSU.freq <- as.data.frame(freq(LSU)) # table of counts of pixels for each hab class in ecoregion
LSU.freq <- subset(LSU.freq,!value %in% NA)# remove NA from calculations
LSU.freq$habprob <- 1/(nrow(LSU.freq))/LSU.freq$count # calculate probs for each class
LSU.freq  

eco.old <- LSU.freq[["value"]]
eco.new <- LSU.freq[["habprob"]]
eco.mat <- matrix(c(rbind(eco.old,eco.new)),ncol=2,byrow=TRUE)
LSU <- reclassify(LSU,eco.mat) # new raster with pixel hab prob

LW<-crop(MB.rast,LWhex)
LW <- mask(LW,LWhex) # clip to the ecoregion
LW.freq <- as.data.frame(freq(LW)) # table of counts of pixels for each hab class in ecoregion
LW.freq <- subset(LW.freq,!value %in% NA)# remove NA from calculations
LW.freq$habprob <- 1/(nrow(LW.freq))/LW.freq$count # calculate probs for each class
LW.freq  

eco.old <- LW.freq[["value"]]
eco.new <- LW.freq[["habprob"]]
eco.mat <- matrix(c(rbind(eco.old,eco.new)),ncol=2,byrow=TRUE)
LW <- reclassify(LW,eco.mat) # new raster with pixel hab prob


MBL<-crop(MB.rast,MBLhex)
MBL<-mask(MBL,MBLhex)
MBL.freq <- as.data.frame(freq(MBL)) # table of counts of pixels for each hab class in ecoregion
MBL.freq <- subset(MBL.freq,!value %in% NA)# remove NA from calculations
MBL.freq$habprob <- 1/(nrow(MBL.freq))/MBL.freq$count # calculate probs for each class
MBL.freq
eco.old <- MBL.freq[["value"]]
eco.new <- MBL.freq[["habprob"]]
eco.mat <- matrix(c(rbind(eco.old,eco.new)),ncol=2,byrow=TRUE)
MBL <- reclassify(MBL,eco.mat) # new raster with pixel hab prob

MBU<-crop(MB.rast,MBUhex)
MBU<-mask(MBU,MBUhex)
MBU.freq <- as.data.frame(freq(MBU)) # table of counts of pixels for each hab class in ecoregion
MBU.freq <- subset(MBU.freq,!value %in% NA)# remove NA from calculations
MBU.freq$habprob <- 1/(nrow(MBU.freq))/MBU.freq$count # calculate probs for each class
MBU.freq
eco.old <- MBU.freq[["value"]]
eco.new <- MBU.freq[["habprob"]]
eco.mat <- matrix(c(rbind(eco.old,eco.new)),ncol=2,byrow=TRUE)
MBU <- reclassify(MBU,eco.mat) # new raster with pixel hab prob

MRU<-crop(MB.rast,MRUhex)
MRU<-mask(MRU,MRUhex)
MRU.freq <- as.data.frame(freq(MRU)) # table of counts of pixels for each hab class in ecoregion
MRU.freq <- subset(MRU.freq,!value %in% NA)# remove NA from calculations
MRU.freq$habprob <- 1/(nrow(MRU.freq))/MRU.freq$count # calculate probs for each class
MRU.freq
eco.old <- MRU.freq[["value"]]
eco.new <- MRU.freq[["habprob"]]
eco.mat <- matrix(c(rbind(eco.old,eco.new)),ncol=2,byrow=TRUE)
MRU <- reclassify(MRU,eco.mat) # new raster with pixel hab prob


SLU<-crop(MB.rast,SLUhex)
SLU<-mask(SLU,SLUhex)
SLU.freq <- as.data.frame(freq(SLU)) # table of counts of pixels for each hab class in ecoregion
SLU.freq <- subset(SLU.freq,!value %in% NA)# remove NA from calculations
SLU.freq$habprob <- 1/(nrow(SLU.freq))/SLU.freq$count # calculate probs for each class
SLU.freq
eco.old <- SLU.freq[["value"]]
eco.new <- SLU.freq[["habprob"]]
eco.mat <- matrix(c(rbind(eco.old,eco.new)),ncol=2,byrow=TRUE)
SLU <- reclassify(SLU,eco.mat) # new raster with pixel hab prob



HabProbraster <- raster::merge(AP,BT, CHBL, CRU, HBL, HRU,ILP,KRU,LMP,LSU,LW,MBL,MBU,MRU,SLU, tolerance = 0.1) #DOESN'T TAKE LONG
plot(HabProbraster)


########################  ---- up to this point updated for MB Dec 4  --CM

##############################################################################################################
############################# Add hab prob and cost prob to NL hexagon layer ##################################

#eco.hab.hex <- raster::extract(HabProbraster,hexagons,fun=sum,na.rm=TRUE,df=TRUE,sp=TRUE) # Takes a long time. new hexagon file with hab probs
#writeOGR(eco.hab.hex, "C:/Users/maderc/Documents/BMS/Manitoba_selection", "ecohabhex", driver = "ESRI Shapefile") #Wrote in new shapefile, temporarily replaced this step with load

eco.hab.hex<-readOGR("ecohabhex.shp", p4s= laea.proj)

hexagons$hab_prob<-eco.hab.hex@data$layer



###### CALCULATE THE INCLUSION PROBABILITY
###### NOTE 1: The variable Balanced is the habitat based inclusion probabilities for Balanced draw (weighted in Mahon parlance) 
###### NOTE 2: The variable invsqrCost is the inverse square root of the minimum cost (across differeint methods) to access the hexagon
###### NOTE 3: HOLDING HABITAT CONSTANT SO IT IS EASY TO TEST WHETHER COST IS HAVING ANY EFFECT

hexagons$p<-(hexagons$hab_prob*hexagons$cost_prob)/sum(hexagons$hab_prob*hexagons$cost_prob)

### dividing by sum puts on probability scale (i.e. all hexagons will sum to 1)

hexagons<-subset(hexagons, p!=0) #remove zeroes
plot(hexagons, col= hexagons$p)

#row.names(hexagons) <- seq(1,nrow(hexagons))  ### this will help for binding data later

rm(HabProbraster)
rm(eco.hab.hex)

##### Target sample sizes
#samplesize<-c(57,20,34,67,12,11,1,11)  # target sample sizes for priority squares
samplesize<-c(3,1,70,24,95,44,36,16,50,1,14,56,15,1,77)  # CODE TEST PLACEHOLDER #### target sample sizes for priority squares BASED on N=400  
Ecoregion<-c("Aspen Parkland","Boreal Transition","Coastal Hudson Bay Lowland","Churchill River Upland",
             "Hudson Bay Lowland", "Hayes River Upland","Interlake Plain","Kazan River Upland", "Lake Manitoba Plain", 
             "Lac Seul Upland", "Lake of the Woods", "Mid-Boreal Lowland","Mid-Boreal Uplands","Maguse River Upland",
             "Selwyn Lake Upland")
sample.size<-data.frame(samplesize, Ecoregion)
sample.size

### May need to add more oversamples for Boreal and Taiga Shield later



#sub<-subset(hexagons, LegacySite==1)
#table(sub$Ecoregion)
#sample.size$Legacy<- c(50,11,26,42,13,7,1,15)
#colnames(sample.size)<-c("Target","Ecoregion","Legacy")
#sample.size$samplesize<-sample.size$Target - sample.size$Legacy
#sample.size$samplesize<-ifelse(sample.size$samplesize< 1, 1,sample.size$samplesize)



#############################################################################################################################
###                                      
###                                      
###                                      		Select the ecoregion to draw from
###                                      
###                                      
#############################################################################################################################

### Total sample size target across ecoregions
N<- sum(sample.size$samplesize)


x<-hexagons$X
y<-hexagons$Y
X <- as.matrix(cbind(x, y))

# inclusion probabilties 
p <- as.matrix(hexagons$p)
p <- N * p / sum( p) #standardise to get n samples  




#############################################################################################################################
###                                      
###                                      
###                                      TEST THE COMPETING DESIGNS by drawing repeated random samples
###                                      
###                                      
#############################################################################################################################



writeSpatialShape(hexagons, "Sample_Frame")



frame <- hexagons


set.seed(3911548)  # Set seed 
# create mdcaty.  keeps variable it is based on in output in addition to mdcaty.  useful sometimes
# also want sum of inclusion probabilities to sum to expected sample size.
attframe <- read.dbf("Sample_Frame")
summary(attframe$p)


attframe$mdcaty <- N * attframe$p/sum(attframe$p)  
hexagons$mdcaty <- attframe$mdcaty
sum(attframe$mdcaty)



###  SPECIFY THE STRATIFIED SAMPLING DESIGN
oversample.size = 1 ### 2x oversample


################### RESTART HERE
# DESIGN 3: FULL BMS DESIGN -----------------------------------------------------

attframe$mdcaty <- N * attframe$p/sum(attframe$p)


#### STRATIFIED DESIGN WITH THE USE OF INCLUSION PROBABILITIES
# NO OVERSAMPLE because just testing the design with Panel one
Stratdsgn <- list("Aspen Parkland"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Aspen Parkland")$samplesize), over=0, seltype="Continuous"),
  "Boreal Transition"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Boreal Transition")$samplesize), over=0, seltype="Continuous"),
  "Coastal Hudson Bay Lowland"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Coastal Hudson Bay Lowland")$samplesize), over=0,seltype="Continuous"),
  "Churchill River Upland"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Churchill River Upland")$samplesize), over=0,seltype="Continuous"),
  "Hudson Bay Lowland"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Hudson Bay Lowland")$samplesize), over=0,seltype="Continuous"),
  "Hayes River Upland"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Hayes River Upland")$samplesize), over=0,seltype="Continuous"),
  "Interlake Plain"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Interlake Plain")$samplesize), over=0,seltype="Continuous"),
  "Kazan River Upland"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Kazan River Upland")$samplesize), over=0,seltype="Continuous"),
  "Lake Manitoba Plain"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Lake Manitoba Plain")$samplesize), over=0,seltype="Continuous"),
  "Lac Seul Upland"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Lac Seul Upland")$samplesize), over=0,seltype="Continuous"),
  "Lake of the Woods"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Lake of the Woods")$samplesize), over=0,seltype="Continuous"),
  "Mid-Boreal Lowland"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Mid-Boreal Lowland")$samplesize), over=0,seltype="Continuous"),
  "Mid-Boreal Uplands"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Mid-Boreal Uplands")$samplesize), over=0,seltype="Continuous"),
  "Maguse River Upland"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Maguse River Upland")$samplesize), over=0,seltype="Continuous"),
  "Selwyn Lake Upland"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Selwyn Lake Upland")$samplesize), over=0,seltype="Continuous")
  )
Stratdsgn

#  "Interlake Plain"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Interlake Plain")$samplesize), over=0,seltype="Continuous"),

#Stratdsgn <- list("Athabasca Plain"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Athabasca Plain")$samplesize), over=round(subset(sample.size, Ecoregion=="Athabasca Plain")$samplesize*oversample.size,0), seltype="Continuous"),
#  "Churchill River Upland"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Churchill River Upland")$samplesize), over=round(subset(sample.size, Ecoregion=="Churchill River Upland")$samplesize*oversample.size,0), seltype="Continuous"),
#  "Mid-Boreal Lowland"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Mid-Boreal Lowland")$samplesize), over=round(subset(sample.size, Ecoregion=="Mid-Boreal Lowland")$samplesize*oversample.size,0),seltype="Continuous"),
#  "Mid-Boreal Uplands"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Mid-Boreal Uplands")$samplesize), over=round(subset(sample.size, Ecoregion=="Mid-Boreal Uplands")$samplesize*oversample.size,0),seltype="Continuous"),
#  "Selwyn Lake Upland"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Selwyn Lake Upland")$samplesize), over=round(subset(sample.size, Ecoregion=="Selwyn Lake Upland")$samplesize*oversample.size,0),seltype="Continuous"),
#  "Tazin Lake Upland"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Tazin Lake Upland")$samplesize), over=round(subset(sample.size, Ecoregion=="Tazin Lake Upland")$samplesize*oversample.size,0),seltype="Continuous"),
#  "Interlake Plain"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Interlake Plain")$samplesize), over=round(subset(sample.size, Ecoregion=="Interlake Plain")$samplesize*oversample.size,0),seltype="Continuous"),
#  "Boreal Transition"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Boreal Transition")$samplesize), over=round(subset(sample.size, Ecoregion=="Boreal Transition")$samplesize*oversample.size,0),seltype="Continuous")
#)
#Stratdsgn


n.simulations=3 #### run 100 random draws
full.design.samp<-list() ###


#### this differs because it saves every iteration, other designs ONLY saving the last iteration
for (i in 1:n.simulations){
full.design.samp[i] <- grts(design=Stratdsgn,
DesignID="ET_ID",
type.frame="area",
src.frame="shapefile",
in.shape="Sample_Frame",
att.frame=attframe,
stratum="Ecoregion",
mdcaty="mdcaty",
shapefile=TRUE,
out.shape=paste0("Full_BMS_design", i)) 
}

#for (i in 1:n.simulations){
#full.design.samp[i] <- grts(design=Stratdsgn,
#DesignID="ET_ID",
#type.frame="area",
#src.frame="shapefile",
#in.shape="Sample_Frame",
#att.frame=attframe,
#stratum="Ecoregion",
#mdcaty="mdcaty",
#shapefile=TRUE,
#out.shape="Full_BMS_design") 
#}


#### get cost estimates from the draws (above)
cost.est<-rep (NA, length(full.design.samp))
for (i in 1:length(full.design.samp)){
  cost.est[[i]] <- sum((full.design.samp[[i]])$Cost)
}

sp.balance<-rep (NA, length(full.design.samp))
for (i in 1:length(full.design.samp)){
  sp.balance[[i]] <- as.matrix(unlist(spbalance(full.design.samp[[i]], frame, tess_ind = F, sbc_ind = T)$sbc$J_subp))[1,1] 
}



hab = list()
for (i in 1:length(full.design.samp)) {
    # ... make some data
    dat <- data.frame(VALUE_1= sum(full.design.samp[[i]]@data$VALUE_1),VALUE_2= sum(full.design.samp[[i]]@data$VALUE_2),VALUE_5= sum(full.design.samp[[i]]@data$VALUE_5),VALUE_6= sum(full.design.samp[[i]]@data$VALUE_6),VALUE_8= sum(full.design.samp[[i]]@data$VALUE_8),VALUE_10= sum(full.design.samp[[i]]@data$VALUE_10),VALUE_14= sum(full.design.samp[[i]]@data$VALUE_14),VALUE_15= sum(full.design.samp[[i]]@data$VALUE_15),VALUE_16= sum(full.design.samp[[i]]@data$VALUE_16),VALUE_17= sum(full.design.samp[[i]]@data$VALUE_17),VALUE_18= sum(full.design.samp[[i]]@data$VALUE_18))
    dat$i <- i  # maybe you want to keep track of which iteration produced it?
    hab[[i]] <- dat # add it to your list
}

habitat = do.call(rbind, hab)


Design.data.full.design.draw<-data.frame(Design=rep("Full BMS design", length(cost.est)), Cost=cost.est, Spatial.Balance= sp.balance)
Design.data.full.design.draw<-data.frame(Design.data.full.design.draw, habitat)
Design.data.full.design.draw
save.image("C:\\Users\\vwilgs\\Documents\\BOREAL POSITION\\Projects\\Boreal Monitoring Strategy\\BMS design test\\Test_Design_largeN.RData")


# DESIGN 4: SPATIAL BALANCE ONLY ------------------------------------------
attframe <- read.dbf("Sample_Frame")
summary(attframe$p)


attframe$mdcaty <- N * attframe$p/sum(attframe$p)  
sum(attframe$mdcaty)

# NO OVERSAMPLE because just testing the design with Panel one
Spatdsgn <- list("Athabasca Plain"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Athabasca Plain")$samplesize), over=0, seltype="Equal"),
  "Churchill River Upland"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Churchill River Upland")$samplesize), over=0, seltype="Equal"),
  "Mid-Boreal Lowland"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Mid-Boreal Lowland")$samplesize), over=0,seltype="Equal"),
  "Mid-Boreal Uplands"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Mid-Boreal Uplands")$samplesize), over=0,seltype="Equal"),
  "Selwyn Lake Upland"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Selwyn Lake Upland")$samplesize), over=0,seltype="Equal"),
  "Tazin Lake Upland"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Tazin Lake Upland")$samplesize), over=0,seltype="Equal"),
  "Boreal Transition"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Boreal Transition")$samplesize), over=0,seltype="Equal")
)
Spatdsgn

#  "Interlake Plain"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Interlake Plain")$samplesize), over=0,seltype="Equal"),

#Spatdsgn <- list("Athabasca Plain"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Athabasca Plain")$samplesize), over=round(subset(sample.size, Ecoregion=="Athabasca Plain")$samplesize*oversample.size,0), seltype="Equal"),
#  "Churchill River Upland"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Churchill River Upland")$samplesize), over=round(subset(sample.size, Ecoregion=="Churchill River Upland")$samplesize*oversample.size,0), seltype="Equal"),
#  "Mid-Boreal Lowland"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Mid-Boreal Lowland")$samplesize), over=round(subset(sample.size, Ecoregion=="Mid-Boreal Lowland")$samplesize*oversample.size,0),seltype="Equal"),
#  "Mid-Boreal Uplands"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Mid-Boreal Uplands")$samplesize), over=round(subset(sample.size, Ecoregion=="Mid-Boreal Uplands")$samplesize*oversample.size,0),seltype="Equal"),
#  "Selwyn Lake Upland"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Selwyn Lake Upland")$samplesize), over=round(subset(sample.size, Ecoregion=="Selwyn Lake Upland")$samplesize*oversample.size,0),seltype="Equal"),
#  "Tazin Lake Upland"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Tazin Lake Upland")$samplesize), over=round(subset(sample.size, Ecoregion=="Tazin Lake Upland")$samplesize*oversample.size,0),seltype="Equal"),
#  "Interlake Plain"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Interlake Plain")$samplesize), over=round(subset(sample.size, Ecoregion=="Interlake Plain")$samplesize*oversample.size,0),seltype="Equal"),
#  "Boreal Transition"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Boreal Transition")$samplesize), over=round(subset(sample.size, Ecoregion=="Boreal Transition")$samplesize*oversample.size,0),seltype="Equal")
#)
#Spatdsgn


n.simulations=100 #### run 100 random draws
spatial.design.samp<-list() ###

for (i in 1:n.simulations){
spatial.design.samp[i] <- grts(design=Spatdsgn,
DesignID="ET_ID",
type.frame="area",
src.frame="shapefile",
in.shape="Sample_Frame",
att.frame=attframe,
stratum="Ecoregion",
mdcaty="mdcaty",
shapefile=TRUE,
out.shape="Spatial_design") 
}



#### get cost estimates from the draws (above)
cost.est<-rep (NA, length(spatial.design.samp))
for (i in 1:length(spatial.design.samp)){
  cost.est[[i]] <- sum((spatial.design.samp[[i]])$Cost)
}

sp.balance<-rep (NA, length(spatial.design.samp))
for (i in 1:length(spatial.design.samp)){
  sp.balance[[i]] <- as.matrix(unlist(spbalance(spatial.design.samp[[i]], frame, tess_ind = F, sbc_ind = T)$sbc$J_subp))[1,1] 
}



hab = list()
for (i in 1:length(spatial.design.samp)) {
    # ... make some data
    dat <- data.frame(VALUE_1= sum(spatial.design.samp[[i]]@data$VALUE_1),VALUE_2= sum(spatial.design.samp[[i]]@data$VALUE_2),VALUE_5= sum(spatial.design.samp[[i]]@data$VALUE_5),VALUE_6= sum(spatial.design.samp[[i]]@data$VALUE_6),VALUE_8= sum(spatial.design.samp[[i]]@data$VALUE_8),VALUE_10= sum(spatial.design.samp[[i]]@data$VALUE_10),VALUE_14= sum(spatial.design.samp[[i]]@data$VALUE_14),VALUE_15= sum(spatial.design.samp[[i]]@data$VALUE_15),VALUE_16= sum(spatial.design.samp[[i]]@data$VALUE_16),VALUE_17= sum(spatial.design.samp[[i]]@data$VALUE_17),VALUE_18= sum(spatial.design.samp[[i]]@data$VALUE_18))
    dat$i <- i  # maybe you want to keep track of which iteration produced it?
    hab[[i]] <- dat # add it to your list
}

habitat = do.call(rbind, hab)



Design.data.spatial.draw<-data.frame(Design=rep("Spatial only", length(cost.est)), Cost=cost.est, Spatial.Balance= sp.balance)
Design.data.spatial.draw<-data.frame(Design.data.spatial.draw, habitat)
save.image("C:\\Users\\vwilgs\\Documents\\BOREAL POSITION\\Projects\\Boreal Monitoring Strategy\\BMS design test\\Test_Design_largeN.RData")
Design.data.spatial.draw


Design.data<-rbind(Design.data.cost.draw,Design.data.habitat.draw,Design.data.full.design.draw,Design.data.spatial.draw)
Design.data
save.image("C:\\Users\\vwilgs\\Documents\\BOREAL POSITION\\Projects\\Boreal Monitoring Strategy\\BMS design test\\Test_Design_largeN.RData")
write.csv(Design.data, file = "Design_comparison_largeN.csv") 

#library(car)
boxplot(Cost~Design, data=Design.data, ylab="Program Cost")
Design.data$factor<-factor(Design.data$Design, levels=c("Cost only", "Full BMS design", "Spatial only", "Habitat only"))

tiff(filename = "C:/Users/vwilgs/Documents/BOREAL POSITION/Projects/Boreal Monitoring Strategy/BMS Manuscript/SK_Cost_figure.tif",width = 174, height = 174, units = "mm",  res=600, compression = c("lzw"),restoreConsole = TRUE)
boxplot(Cost~factor, data=Design.data, ylab="Program Cost over Ten Years ($)", yaxt="n")
axis(2, at=axTicks(2), labels=sprintf("$%s", axTicks(2)))
dev.off()

tiff(filename = "C:/Users/vwilgs/Documents/BOREAL POSITION/Projects/Boreal Monitoring Strategy/BMS Manuscript/SK_spatial_balance_figure.tif",width = 174, height = 174, units = "mm",  res=600, compression = c("lzw"),restoreConsole = TRUE)
boxplot(Spatial.Balance~factor, data=Design.data, ylab="Spatial Balance")
dev.off()

tiff(filename = "C:/Users/vwilgs/Documents/BOREAL POSITION/Projects/Boreal Monitoring Strategy/BMS Manuscript/SK_Habitat_figure.tif",width = 174, height = 174, units = "mm",  res=600, compression = c("lzw"),restoreConsole = TRUE)
Boxplot(SSD~Design, data=Designs, ylab="Sum of Squared Differences", id.method="none")
dev.off()
