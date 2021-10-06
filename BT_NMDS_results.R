library(vegan)
library(gstudio)
library(pegas)
library(plyr)
library(scales)
library(car)

#Importing Barrier runs:
data.dir <- "/Users/jamesonhinkle/Desktop/GenepopFiles_Riverine/"
#data.dir <- "C:/Users/Jacob/Desktop/RivGen/Meryl_Talk/Meryl_Runs_Genepop_Files/GenepopFiles_Barriers/"
scenario <- "Barriers"
species <- "_BT"
date <- Sys.Date()
file.ext <- ".gen"

#Generations to monitor
gen.seq <- c(150) #-1, 0, 1, 10, 150
batch.seq <- 3 #0:5

# read in watershed information
#wshed <- read.csv("C:/Users/Jacob/Desktop/RivGen/Meryl_Talk/PatchXYwatershed.csv", header=T)
wshed <- read.csv("PatchXYwatershed.csv", header=T)
wat.df <- data.frame(Population=as.numeric(wshed$FID)+1, POINT_X=wshed$POINT_X, POINT_Y=wshed$POINT_Y, trib=wshed$watershed)

for(i in 1:length(gen.seq)){
  for(j in 1:length(batch.seq)){
    batch.no <- batch.seq[j]
    gen <- gen.seq[i]
    
    batch <- paste("batchrun", batch.no, sep="")
    mc0 <- paste("mcrun0/", "genepop_ind", gen, sep="")
    mc1 <- paste("mcrun1/", "genepop_ind", gen, sep="")
    mc2 <- paste("mcrun2/", "genepop_ind", gen, sep="")
    
    # Read in Data 
    genmc0 <- read_population(paste(data.dir, batch, mc0, file.ext, sep=""), type="genepop")
    genmc1 <- read_population(paste(data.dir, batch, mc1, file.ext, sep=""), type="genepop")
    genmc2 <- read_population(paste(data.dir, batch, mc2, file.ext, sep=""), type="genepop")
    
    genmc0$Pop <- genmc0$Population
    genmc1$Pop <- genmc1$Population
    genmc2$Pop <- genmc2$Population
    
    #Add in strata metric
    for(id in 1:length(genmc0$Pop)){
      genmc0$Population[id] <- as.numeric(unlist(strsplit(genmc0$Pop[id], "Pop-")))[2]     #remove the pop- identifier in the genepop files
    }
    
    for(id in 1:length(genmc1$Population)){
      genmc1$Population[id] <- as.numeric(unlist(strsplit(genmc1$Pop[id], "Pop-")))[2]
    }
    
    for(id in 1:length(genmc2$Population)){
      genmc2$Population[id] <- as.numeric(unlist(strsplit(genmc2$Pop[id], "Pop-")))[2]
    }
    
    gen.mc0 <- join(x=genmc0, y=wat.df, by="Population", type="left")
    gen.mc1 <- join(x=genmc1, y=wat.df, by="Population", type="left")
    gen.mc2 <- join(x=genmc2, y=wat.df, by="Population", type="left")
    
    gen.mc0$Population <- gen.mc0$trib
    gen.mc1$Population <- gen.mc1$trib
    gen.mc2$Population <- gen.mc2$trib
  
    mv <- to_mv(gen.mc2, drop.allele = TRUE)
    ##############################################
    ##Non-metric multidimensional scaling (NMDS)##
    ##############################################
    
    
    ##CHOOSE ONE -- drop rare alleles; 5% cutoff is most conservative & recommended:
    
    mv.r <- drop.var(mv, min.po=1) # alleles dropped if in <1% of individuals
    #mv.r <- drop.var(mv, min.fo=10) #alleles dropped if in <10 individuals
    
    
    ##scree plot to determine #dimensions (may take a few hours to run with >1000 individuals):
    
    #nmds.scree(mv.r, distance="bray", k=10, autotransform=FALSE, trymax=10)
    
    
    ##set K to # dimensions determined from scree plot (will take a few minutes):
    
    gen.nmds <-metaMDS(mv.r, distance="bray", k=2, autotransform=FALSE, trymax=10) #try with K = 4
    
    ##create dataframe with NMDS scores, watershed data, and individual ID:
    
    NMDSpops <- cbind(gen.mc0$Pop, gen.mc0$watershed, gen.mc0$extbar, gen.mc0$ColID, as.vector(gen.mc0$RGB), gen.nmds$points)
    id <- rownames(NMDSpops)
    NMDSpops2 <- cbind(id=id, NMDSpops)
    NMDSpops <- data.frame(NMDSpops2)
    names(NMDSpops)[names(NMDSpops)=="V6"] <- "RGB" #for following merge function
    names(NMDSpops)[names(NMDSpops)=="V5"] <- "ColID" #for following merge function
    
    
  }
  
}

##OPTIONAL -- quick plot:

plot(gen.nmds, "sites")
