
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




setwd("C:/Users/vwilgs/Documents/BOREAL POSITION/Projects/Boreal Monitoring Strategy/BMS design test/")
load("C:\\Users\\vwilgs\\Documents\\BOREAL POSITION\\Projects\\Boreal Monitoring Strategy\\BMS design test\\Test_Design_largeN.RData")
#setwd("C:/Users/vwilgs/Documents/BOREAL POSITION/Projects/Saskatchewan Breeding Bird Atlas/Boreal Site Selection/")



###################################################################################################################################
##########################################          IMPORT THE SAMPLING FRAME            ##########################################
###################################################################################################################################


# Albers equal-area projection parameters based on above
#aea.proj <- "+proj=aea +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-96+x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m"
laea.proj <- "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs" # use this because it is lccc10 projection

hexagons<-readShapePoly("SK_Hexagons_w_Landcover_laea", proj4string= CRS(laea.proj))
hexagons

plot(hexagons, col=hexagons$Cost)
hexagons$cost_prob<-hexagons$invsqrCost


writePolyShape(hexagons, "hexagons")

## now use the final_hexagons to clip lccc10 raster
# read in lccc10 file
LCCC10 <- raster("C:/Users/vwilgs/Documents/BOREAL POSITION/Projects/Boreal Monitoring Strategy/BMS design test/NA_LandCover_2010_25haMMU.tif")
SK.rast <- crop(LCCC10,as(hexagons,"Spatial")) # LCCC raster for SK only
SK.rast <- mask(SK.rast,hexagons,"Spatial")

# remove habitat classes not wanted in the model
SK.rast <- calc(SK.rast, fun=function(x){x[x==19] <- NA; return(x)}) # remove snow and ice
SK.rast <- calc(SK.rast, fun=function(x){x[x==18] <- NA; return(x)}) # water
SK.rast <- calc(SK.rast, fun=function(x){x[x==17] <- NA; return(x)}) # urban
SK.rast <- calc(SK.rast, fun=function(x){x[x==15] <- NA; return(x)}) # cropland
SK.rast <- calc(SK.rast, fun=function(x){x[x==0] <- NA; return(x)}) # pixels not included in the cropping
writeRaster(SK.rast,filename="SK.rast.tif",format="GTiff",overwrite=TRUE)

rm(LCCC10) # clear some memory

##############################################################################################################
############################ Now calculating Habitat Probabilities for ecoregion #############################


eco.list <- c("Athabasca Plain","Boreal Transition","Churchill River Upland","Interlake Plain",
              "Mid-Boreal Lowland","Mid-Boreal Uplands","Selwyn Lake Upland",
              "Tazin Lake Upland")

APhex<-subset(hexagons,Ecoregion==eco.list[1])

BThex<-subset(hexagons,Ecoregion==eco.list[2])
CRUhex<-subset(hexagons,Ecoregion==eco.list[3])
IPhex<-subset(hexagons,Ecoregion==eco.list[4])
MBLhex<-subset(hexagons,Ecoregion==eco.list[5])
MBUhex<-subset(hexagons,Ecoregion==eco.list[6])
SLUhex<-subset(hexagons,Ecoregion==eco.list[7])
TLUhex<-subset(hexagons,Ecoregion==eco.list[8])

AP<-crop(SK.rast,APhex)
AP <- mask(AP,APhex) # clip to the ecoregion
AP.freq <- as.data.frame(freq(AP)) # table of counts of pixels for each hab class in ecoregion
AP.freq <- subset(AP.freq,!value %in% NA)# remove NA from calculations
AP.freq$habprob <- 1/(nrow(AP.freq))/AP.freq$count # calculate probs for each class
AP.freq  

eco.old <- AP.freq[["value"]]
eco.new <- AP.freq[["habprob"]]
eco.mat <- matrix(c(rbind(eco.old,eco.new)),ncol=2,byrow=TRUE)
AP <- reclassify(AP,eco.mat) # new raster with pixel hab prob

BT<-crop(SK.rast,BThex)
BT <- mask(BT,BThex) # clip to the ecoregion
BT.freq <- as.data.frame(freq(BT)) # table of counts of pixels for each hab class in ecoregion
BT.freq <- subset(BT.freq,!value %in% NA)# remove NA from calculations
BT.freq$habprob <- 1/(nrow(BT.freq))/BT.freq$count # calculate probs for each class
BT.freq  
eco.old <- BT.freq[["value"]]
eco.new <- BT.freq[["habprob"]]
eco.mat <- matrix(c(rbind(eco.old,eco.new)),ncol=2,byrow=TRUE)
BT <- reclassify(BT,eco.mat) # new raster with pixel hab prob


CRU<-crop(SK.rast,CRUhex)
CRU <- mask(CRU,CRUhex) # clip to the ecoregion
CRU.freq <- as.data.frame(freq(CRU)) # table of counts of pixels for each hab class in ecoregion
CRU.freq <- subset(CRU.freq,!value %in% NA)# remove NA from calculations
CRU.freq$habprob <- 1/(nrow(CRU.freq))/CRU.freq$count # calculate probs for each class
CRU.freq
eco.old <- CRU.freq[["value"]]
eco.new <- CRU.freq[["habprob"]]
eco.mat <- matrix(c(rbind(eco.old,eco.new)),ncol=2,byrow=TRUE)
CRU <- reclassify(CRU,eco.mat) # new raster with pixel hab prob

IP<-crop(SK.rast,IPhex)
IP<-mask(IP,IPhex)
IP.freq <- as.data.frame(freq(IP)) # table of counts of pixels for each hab class in ecoregion
IP.freq <- subset(IP.freq,!value %in% NA)# remove NA from calculations
IP.freq$habprob <- 1/(nrow(IP.freq))/IP.freq$count # calculate probs for each class
IP.freq
eco.old <- IP.freq[["value"]]
eco.new <- IP.freq[["habprob"]]
eco.mat <- matrix(c(rbind(eco.old,eco.new)),ncol=2,byrow=TRUE)
IP <- reclassify(IP,eco.mat) # new raster with pixel hab prob


MBL<-crop(SK.rast,MBLhex)
MBL<-mask(MBL,MBLhex)
MBL.freq <- as.data.frame(freq(MBL)) # table of counts of pixels for each hab class in ecoregion
MBL.freq <- subset(MBL.freq,!value %in% NA)# remove NA from calculations
MBL.freq$habprob <- 1/(nrow(MBL.freq))/MBL.freq$count # calculate probs for each class
MBL.freq
eco.old <- MBL.freq[["value"]]
eco.new <- MBL.freq[["habprob"]]
eco.mat <- matrix(c(rbind(eco.old,eco.new)),ncol=2,byrow=TRUE)
MBL <- reclassify(MBL,eco.mat) # new raster with pixel hab prob

MBU<-crop(SK.rast,MBUhex)
MBU<-mask(MBU,MBUhex)
MBU.freq <- as.data.frame(freq(MBU)) # table of counts of pixels for each hab class in ecoregion
MBU.freq <- subset(MBU.freq,!value %in% NA)# remove NA from calculations
MBU.freq$habprob <- 1/(nrow(MBU.freq))/MBU.freq$count # calculate probs for each class
MBU.freq
eco.old <- MBU.freq[["value"]]
eco.new <- MBU.freq[["habprob"]]
eco.mat <- matrix(c(rbind(eco.old,eco.new)),ncol=2,byrow=TRUE)
MBU <- reclassify(MBU,eco.mat) # new raster with pixel hab prob

SLU<-crop(SK.rast,SLUhex)
SLU<-mask(SLU,SLUhex)
SLU.freq <- as.data.frame(freq(SLU)) # table of counts of pixels for each hab class in ecoregion
SLU.freq <- subset(SLU.freq,!value %in% NA)# remove NA from calculations
SLU.freq$habprob <- 1/(nrow(SLU.freq))/SLU.freq$count # calculate probs for each class
SLU.freq
eco.old <- SLU.freq[["value"]]
eco.new <- SLU.freq[["habprob"]]
eco.mat <- matrix(c(rbind(eco.old,eco.new)),ncol=2,byrow=TRUE)
SLU <- reclassify(SLU,eco.mat) # new raster with pixel hab prob

TLU<-crop(SK.rast,TLUhex)
TLU<-mask(TLU,TLUhex)
TLU.freq <- as.data.frame(freq(TLU)) # table of counts of pixels for each hab class in ecoregion
TLU.freq <- subset(TLU.freq,!value %in% NA)# remove NA from calculations
TLU.freq$habprob <- 1/(nrow(TLU.freq))/TLU.freq$count # calculate probs for each class
TLU.freq
eco.old <- TLU.freq[["value"]]
eco.new <- TLU.freq[["habprob"]]
eco.mat <- matrix(c(rbind(eco.old,eco.new)),ncol=2,byrow=TRUE)
TLU <- reclassify(TLU,eco.mat) # new raster with pixel hab prob



HabProbraster <- raster::merge(AP,BT,CRU,IP,MBL,MBU,SLU,TLU, tolerance = 0.1)
plot(HabProbraster)


##############################################################################################################
############################# Add hab prob and cost prob to NL hexagon layer ##################################
eco.hab.hex <- raster::extract(HabProbraster,hexagons,fun=sum,na.rm=TRUE,df=TRUE,sp=TRUE) # new hexagon file with hab probs

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
samplesize<-c(113,40,68,133,23,22,1,22)  # target sample sizes for priority squares BASED on N=400
Ecoregion<-c("Mid-Boreal Uplands","Mid-Boreal Lowland","Athabasca Plain","Churchill River Upland","Tazin Lake Upland","Selwyn Lake Upland","Interlake Plain","Boreal Transition")
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


#set.seed(3911548)  # Set seed 
# create mdcaty.  keeps variable it is based on in output in addition to mdcaty.  useful sometimes
# also want sum of inclusion probabilities to sum to expected sample size.
attframe <- read.dbf("Sample_Frame")
summary(attframe$p)


attframe$mdcaty <- N * attframe$p/sum(attframe$p)  
sum(attframe$mdcaty)



###  SPECIFY THE STRATIFIED SAMPLING DESIGN
oversample.size = 1 ### 2x oversample


# DESIGN 1: COST ONLY -----------------------------------------------------
attframe$mdcaty <- N * attframe$cost_prob/sum(attframe$cost_prob) # Keep p standardized by N

# NO OVERSAMPLE because just testing the design with Panel one
Stratdsgn <- list("Athabasca Plain"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Athabasca Plain")$samplesize), over=0, seltype="Continuous"),
   "Churchill River Upland"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Churchill River Upland")$samplesize), over=0, seltype="Continuous"),
   "Mid-Boreal Lowland"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Mid-Boreal Lowland")$samplesize), over=0,seltype="Continuous"),
   "Mid-Boreal Uplands"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Mid-Boreal Uplands")$samplesize), over=0,seltype="Continuous"),
   "Selwyn Lake Upland"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Selwyn Lake Upland")$samplesize), over=round(subset(sample.size, Ecoregion=="Selwyn Lake Upland")$samplesize*oversample.size,0),seltype="Continuous"),
   "Tazin Lake Upland"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Tazin Lake Upland")$samplesize), over=0,seltype="Continuous"),
   "Boreal Transition"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Boreal Transition")$samplesize), over=0,seltype="Continuous")
)
Stratdsgn

n.simulations=100 #### run 100 random draws
cost.samp<-list() ###

for (i in 1:n.simulations){
cost.samp[i] <- grts(design=Stratdsgn,
DesignID="ET_ID",
type.frame="area",
src.frame="shapefile",
in.shape="Sample_Frame",
att.frame=attframe,
stratum="Ecoregion",
mdcaty="mdcaty",
shapefile=TRUE,
out.shape="Cost_design") 
}


cost.est<-rep (NA, length(cost.samp))
 for (i in 1:length(cost.samp)){
   cost.est[[i]] <- sum(cost.samp[[i]]$Cost)
}

sp.balance<-rep (NA, length(cost.samp))
 for (i in 1:length(cost.samp)){
   sp.balance[[i]] <- as.matrix(unlist(spbalance(cost.samp[[i]], frame, tess_ind = F, sbc_ind = T)$sbc$J_subp))[1,1] 
}


hab = list()

for (i in 1:length(cost.samp)) {
    # ... make some data
    dat <- data.frame(VALUE_1= sum(cost.samp[[i]]@data$VALUE_1),VALUE_2= sum(cost.samp[[i]]@data$VALUE_2),VALUE_5= sum(cost.samp[[i]]@data$VALUE_5),VALUE_6= sum(cost.samp[[i]]@data$VALUE_6),VALUE_8= sum(cost.samp[[i]]@data$VALUE_8),VALUE_10= sum(cost.samp[[i]]@data$VALUE_10),VALUE_14= sum(cost.samp[[i]]@data$VALUE_14),VALUE_15= sum(cost.samp[[i]]@data$VALUE_15),VALUE_16= sum(cost.samp[[i]]@data$VALUE_16),VALUE_17= sum(cost.samp[[i]]@data$VALUE_17),VALUE_18= sum(cost.samp[[i]]@data$VALUE_18))
    dat$i <- i  # maybe you want to keep track of which iteration produced it?
    hab[[i]] <- dat # add it to your list
}
habitat = do.call(rbind, hab)


Design.data.cost.draw<-data.frame(Design=rep("Cost only", length(cost.est)), Cost=cost.est, Spatial.Balance= sp.balance)
Design.data.cost.draw<-data.frame(Design.data.cost.draw, habitat)
Design.data.cost.draw
save.image("C:\\Users\\vwilgs\\Documents\\BOREAL POSITION\\Projects\\Boreal Monitoring Strategy\\BMS design test\\Test_Design_largeN.RData")


#### STRATIFIED DESIGN WITH THE USE OF INCLUSION PROBABILITIES


# DESIGN 2: HABITAT ONLY --------------------------------------------------

attframe$mdcaty <- N * attframe$hab_prob/sum(attframe$hab_prob)

# NO OVERSAMPLE because just testing the design with Panel one
Stratdsgn <- list("Athabasca Plain"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Athabasca Plain")$samplesize), over=0, seltype="Continuous"),
  "Churchill River Upland"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Churchill River Upland")$samplesize), over=0, seltype="Continuous"),
  "Mid-Boreal Lowland"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Mid-Boreal Lowland")$samplesize), over=0,seltype="Continuous"),
  "Mid-Boreal Uplands"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Mid-Boreal Uplands")$samplesize), over=0,seltype="Continuous"),
  "Selwyn Lake Upland"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Selwyn Lake Upland")$samplesize), over=round(subset(sample.size, Ecoregion=="Selwyn Lake Upland")$samplesize*oversample.size,0),seltype="Continuous"),
  "Tazin Lake Upland"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Tazin Lake Upland")$samplesize), over=0,seltype="Continuous"),
  "Boreal Transition"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Boreal Transition")$samplesize), over=0,seltype="Continuous")
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



n.simulations=100 #### run 100 random draws
hab.samp<-list() ###

for (i in 1:n.simulations){
hab.samp[i] <- grts(design=Stratdsgn,
DesignID="ET_ID",
type.frame="area",
src.frame="shapefile",
in.shape="Sample_Frame",
att.frame=attframe,
stratum="Ecoregion",
mdcaty="mdcaty",
shapefile=TRUE,
out.shape="Habitat_design") 
}
save.image("C:\\Users\\vwilgs\\Documents\\BOREAL POSITION\\Projects\\Boreal Monitoring Strategy\\BMS design test\\Test_Design_largeN.RData")



#### get cost estimates from the draws (above)
cost.est<-rep (NA, length(hab.samp))
for (i in 1:length(hab.samp)){
  cost.est[[i]] <- sum((hab.samp[[i]])$Cost)
}

sp.balance<-rep (NA, length(hab.samp))
for (i in 1:length(hab.samp)){
  sp.balance[[i]] <- as.matrix(unlist(spbalance(hab.samp[[i]], frame, tess_ind = F, sbc_ind = T)$sbc$J_subp))[1,1] 
}



hab = list()
for (i in 1:length(hab.samp)) {
    # ... make some data
    dat <- data.frame(VALUE_1= sum(hab.samp[[i]]@data$VALUE_1),VALUE_2= sum(hab.samp[[i]]@data$VALUE_2),VALUE_5= sum(hab.samp[[i]]@data$VALUE_5),VALUE_6= sum(hab.samp[[i]]@data$VALUE_6),VALUE_8= sum(hab.samp[[i]]@data$VALUE_8),VALUE_10= sum(hab.samp[[i]]@data$VALUE_10),VALUE_14= sum(hab.samp[[i]]@data$VALUE_14),VALUE_15= sum(hab.samp[[i]]@data$VALUE_15),VALUE_16= sum(hab.samp[[i]]@data$VALUE_16),VALUE_17= sum(hab.samp[[i]]@data$VALUE_17),VALUE_18= sum(hab.samp[[i]]@data$VALUE_18))
    dat$i <- i  # maybe you want to keep track of which iteration produced it?
    hab[[i]] <- dat # add it to your list
}

habitat = do.call(rbind, hab)



Design.data.habitat.draw<-data.frame(Design=rep("Habitat only", length(cost.est)), Cost=cost.est, Spatial.Balance= sp.balance)
Design.data.habitat.draw<-data.frame(Design.data.habitat.draw, habitat)
Design.data.habitat.draw
save.image("C:\\Users\\vwilgs\\Documents\\BOREAL POSITION\\Projects\\Boreal Monitoring Strategy\\BMS design test\\Test_Design_largeN.RData")

################### RESTART HERE
# DESIGN 3: FULL BMS DESIGN -----------------------------------------------------

attframe$mdcaty <- N * attframe$p/sum(attframe$p)


#### STRATIFIED DESIGN WITH THE USE OF INCLUSION PROBABILITIES
# NO OVERSAMPLE because just testing the design with Panel one
Stratdsgn <- list("Athabasca Plain"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Athabasca Plain")$samplesize), over=0, seltype="Continuous"),
  "Churchill River Upland"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Churchill River Upland")$samplesize), over=0, seltype="Continuous"),
  "Mid-Boreal Lowland"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Mid-Boreal Lowland")$samplesize), over=0,seltype="Continuous"),
  "Mid-Boreal Uplands"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Mid-Boreal Uplands")$samplesize), over=0,seltype="Continuous"),
  "Selwyn Lake Upland"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Selwyn Lake Upland")$samplesize), over=round(subset(sample.size, Ecoregion=="Selwyn Lake Upland")$samplesize*oversample.size,0),seltype="Continuous"),
  "Tazin Lake Upland"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Tazin Lake Upland")$samplesize), over=0,seltype="Continuous"),
  "Boreal Transition"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Boreal Transition")$samplesize), over=0,seltype="Continuous")
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


n.simulations=100 #### run 100 random draws
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
}


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
