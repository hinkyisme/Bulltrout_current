##For post-sim processing of Bull Trout simulations (CDMetaPop) and NMDS analysis
##post-sim processing code: J. Hinkle and J. Burkhart; NMDS code: M. Mims

##required packages:
install.packages("vegan")
install.packages("gstudio")
install.packages("pegas")
install.packages("plyr")
install.packages("scales")
install.packages("car")
library(vegan)
library(gstudio)
library(pegas)
library(plyr)
library(scales)
library(car)


##############################################
####Extracting and organizing genetic data####
##############################################


##set directory, scenario, and file information:

source('C:/Users/Owner/Desktop/Bull Trout Genetics/biostats.R', encoding = 'UTF-8')
setwd("C:/Users/Owner/Desktop/Bull Trout Genetics/BT_postSim processing")
data.dir <-"C:/Users/Owner/Desktop/Bull Trout Genetics/GenepopFiles_ExtTwoWay/"
scenario <- "ExtTwoWay"
species <- "_BT"
date <- "_20150701"
file.ext <- ".gen"


##set batch and generation information:

batch.no <-40
gen <- 150
batch <- paste("batchrun", batch.no, sep = "")
mc0 <- paste("mcrun0/", "genepop_ind", gen, sep = "")


##read in genetic data:

genmc0 <-read_population(paste(data.dir, batch, mc0, file.ext, sep = ""), type = "genepop")
genmc0$Pop <- genmc0$Population
    
##write strata data:

    for (id in 1:length(genmc0$Pop))
    {
      genmc0$Population[id] <- as.numeric(unlist(strsplit(genmc0$Pop[id], "Pop-")))[2]
    }
 
   
##read in watershed information (from M. Fuller):

wshed <-read.csv("C:/Users/Owner/Desktop/Bull Trout Genetics/WatershedIDsublevels.csv", header = T)
wat.df <-data.frame(Population = as.numeric(wshed$FID) + 1, watershed = wshed$watershed, extbar = wshed$extbarrsub, futbar = wshed$futbarrsub, ColID = wshed$ColorID, RGB = wshed$RGB)


##join watershed information with genetic data:

gen.mc0 <-join(x = genmc0, y = wat.df, by = "Population", type = "left")
gen.mc0$Population <- gen.mc0$trib


##OPTIONAL -- write to csv for external programs/analysis:

write.csv(x = gen.mc0, file = paste(data.dir, batch, mc0, ".csv", sep=""))


##covert to allele frequency matrix:

mv <- to_mv(gen.mc0, drop.allele = TRUE)


##OPTIONAL: write allele frequency matrix to csv:

write.csv(x = mv, file = paste(data.dir, batch, mc0, "_AlleleFreqMat.csv", sep=""))




##############################################
##Non-metric multidimensional scaling (NMDS)##
##############################################


##CHOOSE ONE -- drop rare alleles; 5% cutoff is most conservative & recommended:

mv.r <- drop.var(mv, min.po=1) # alleles dropped if in <1% of individuals
#mv.r <- drop.var(mv, min.fo=10) #alleles dropped if in <10 individuals


##scree plot to determine #dimensions (may take a few hours to run with >1000 individuals):

nmds.scree(mv.r, distance="bray", k=10, autotransform=FALSE, trymax=10)


##set K to # dimensions determined from scree plot (will take a few minutes):

gen.nmds <-metaMDS(mv.r, distance="bray", k=2, autotransform=FALSE, trymax=10)


##OPTIONAL -- quick plot:

plot(gen.nmds, "sites")

##create dataframe with NMDS scores, watershed data, and individual ID:

NMDSpops <- cbind(gen.mc0$Pop, gen.mc0$watershed, gen.mc0$extbar, gen.mc0$ColID, as.vector(gen.mc0$RGB), gen.nmds$points)
id <- rownames(NMDSpops)
NMDSpops2 <- cbind(id=id, NMDSpops)
NMDSpops <- data.frame(NMDSpops2)
names(NMDSpops)[names(NMDSpops)=="V6"] <- "RGB" #for following merge function
names(NMDSpops)[names(NMDSpops)=="V5"] <- "ColID" #for following merge function

##OPTIONAL -- write out NMDSpops:

write.csv(x = NMDSpops, file= paste(data.dir, batch, mc0, "_NMDSpops.csv", sep=""))


##ANOSIM (takes a while)

gen.d <- vegdist(mv.r, "bray")
gen.anosim <- anosim(gen.d, barrier)
summary(gen.anosim)
plot.anosim(gen.anosim)




##############################################
######Plotting options and figure export######
##############################################


##OPTIONAL -- read in NMDScol:

NMDSpops <- read.csv(paste(data.dir, batch, mc0, "_NMDSpops.csv", sep=""))
NMDSpops_sing <- read.csv(paste(data.dir, batch, mc0, "_NMDSpops_singles.csv", sep=""))

##plot prep:

NMDS1 <- as.vector(NMDSpops$MDS1)
NMDS2 <- as.vector(NMDSpops$MDS2)
wshedcol <- as.vector(NMDSpops$RGB)
barrier <- as.factor(NMDSpops$ColID)
NMDS1sing <- as.vector(NMDSpops_sing$MDS1)
NMDS2sing <- as.vector(NMDSpops_sing$MDS2)

#export at PNG:
png("NMDS_40.150.2.png", width=1500, height=1500, units='px', res=300)
par(mfrow=c(1,1), oma=c(0,0,0,0), mar=c(3,3,2,2))
plot(NMDS1, NMDS2, type="n", xlim=c(-1,1), ylim=c(-1,1), cex.axis=0.75,  xlab="", ylab="")
#points(NMDS1, NMDS2, pch=19, col="white")
points(NMDS1, NMDS2, pch=19, col=alpha(wshedcol, 0.5))
ordihull(gen.nmds, wshedcol, col="black")
points(NMDS1sing, NMDS2sing, pch=21, col="black", bg=c("#E1E1E1","#A2CF8C","#A2CF8C","#ABD1A7", "#ABD1A7"))
#text(0, 0.95, pos=1, "Existing one-way, early migrate", cex=1, col="grey")
#text(0.9, -0.9, pos=2, "150 Years", cex=1, col="grey")
dev.off()




##############################################
#########Additional plotting options##########
##############################################


##
##plot individuals with existing barrier group color, with outline:
plot(NMDS1, NMDS2, type="n", xlim=c(-1,1), ylim=c(-1,1))
points(NMDS1, NMDS2, pch=21, bg=wshedcol)


##plot individuals with existing barrier group color, without outline:
plot(NMDS1, NMDS2, type="n", xlim=c(-1,1), ylim=c(-1,1))
points(NMDS1, NMDS2, pch=19, col=wshedcol)


##plot individuals with existing barrier group, without outline, with convex hulls:
plot(NMDS1, NMDS2, type="n", xlim=c(-1,1), ylim=c(-1,1))
points(NMDS1, NMDS2, pch=19, col=wshedcol)
barrier <- as.vector(NMDSpops$extbar)
ordihull(gen.nmds, barrier)



##plot individuals with existing barrier group, without outline, transparent symbols:
plot(NMDS1, NMDS2, type="n", xlim=c(-1,1), ylim=c(-1,1))
dataEllipse(NMDS1, NMDS2, barrier, center.pch=FALSE, levels = 0.95, add=TRUE, lwd=0.75, col = c("grey", "grey"), plot.points=FALSE)
points(NMDS1, NMDS2, pch=19, col="white")
points(NMDS1, NMDS2, pch=19, col=alpha(wshedcol, 0.3))


##Export figure - UPDATE file and text label for each

png("NMDS_11.150.png", width=1500, height=3000, units='px', res=300)
par(mfrow=c(2,1), oma=c(0,0,0,0), mar=c(3,3,2,2))
plot(NMDS1, NMDS2, type="n", xlim=c(-1,1), ylim=c(-1,1), cex.axis=0.75,  xlab="", ylab="")
ordihull(gen.nmds, barrier, col="grey", lwd=0.75)
points(NMDS1, NMDS2, pch=19, col="white")
points(NMDS1, NMDS2, pch=19, col=alpha(wshedcol, 0.5))
barrier <- as.vector(NMDSpops$extbar)
text(-0.9, 0.9, pos=4, "Batch 11 Gen 150", cex=0.75)