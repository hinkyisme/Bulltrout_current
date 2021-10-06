##For post-sim processing of Bull Trout simulations (CDMetaPop) and NMDS analysis
##post-sim processing code: J. Hinkle and J. Burkhart; NMDS code: M. Mims

#Load Packages
require(gstudio)
require(pegas)
require(plyr)
library(vegan)
library(scales)
library(car)


#Importing Data Files
  #Meryl
    #source('C:/Users/Owner/Desktop/Bull Trout Genetics/biostats.R', encoding = 'UTF-8')
    #setwd("C:/Users/Owner/Desktop/Bull Trout Genetics/BT_postSim processing")
    #data.dir <-"C:/Users/Owner/Desktop/Bull Trout Genetics/GenepopFiles_ExtTwoWay/"
    #wshed <- read.csv("UPDATE FILE PATH/WatershedIDsublevels.csv", header=T)              #Meryl, update file path

  #Jameson
    #source('UPDATE FILE PATH', encoding = 'UTF-8')                                   #Jameson, update file path
    #data.dir <- "/Users/jamesonhinkle/Desktop/Raw_Output07122015/"
    #wshed <- read.csv("WatershedIDsublevels_07122015.csv", sep = ",", header=T)         

  #Jake
    source("C:/Users/Jacob/Desktop/RivGen/Analyses_Fa2015/biostats.R", encoding = 'UTF-8')
    #data.dir <- "C:/Users/Jacob/Desktop/RivGen/Meryl_Talk/Meryl_Runs_Genepop_Files/GenepopFiles_Barriers/"
    data.dir <- "C:/Users/Jacob/Desktop/RivGen/Raw_Output07122015/"    #should only need to change batch number between scenario runs
    wshed <- read.csv("C:/Users/Jacob/Desktop/RivGen/Raw_Output07122015/WatershedIDsublevels_07122015.csv", header=T)   #update
  
  #Input/Output File Parameters
  species <- "_BT"
  date <- Sys.Date()
  file.ext <- ".gen"

#Generations to monitor
  #fol <- 0:11
  gen.seq <- c(-1,0,1,10,50,100,150)   #subset: c(-1,0,1,10,50,100,150); all gens c(-1,0,1,2,3,4,5,10,15,20,25,30,40,50,60,70,80,90,91,92,93,94,95,96,97,98,99,100,110,120,130,140,150,160,170,180,190,199)
  batch.seq <- 0:11                    #batch 0-3: riverine; batch 4-7: barriers; batch 8-11: panmictic

#Create a watershed dataframe for future plotting:
  wat.df <-data.frame(Population=as.numeric(wshed$FID)+1, watershed=wshed$watershed, 
                      extbar=wshed$ExistingBarriersublevelID, futbar=wshed$FutureBarriersublevelID,
                      ColID=wshed$ColorID, RGB=wshed$RGB)

#Create barrier object for ANOSIM 
  barrier <- as.factor(wat.df$ColID)

#Iterate across generations and batches to create NMDS data
  start.time = Sys.time() # Set timer for start of run
#for(i in 1:length(fol))                                          #ignore this line for now - used to loop through folders but JJB changed input file structure to remove 3rd for loop
  for(j in 1:length(gen.seq)){
    for(k in 1:length(batch.seq)){
      #Set up input parameters
        #fol.no <- fol[i]                                         #ignore this line for now - used to loop through folders but JJB changed input file structure to remove 3rd for loop
        batch.no <- batch.seq[k]
        gen <- gen.seq[j]
      
      #Create scenario object for output files
        if(batch.no >= 0 & batch.no <= 3){
          scenario <- "Riverine"
        } else if(batch.no >= 4 & batch.no <= 7) {
              scenario <- "Barriers"
            } else{
                scenario <- "Panmictic"
              }
         
      #Prepare input file objects  
        #batch <- paste("batch", fol.no, batch.no, sep="")        #ignore this line for now - used to loop through folders but JJB changed input file structure to remove 3rd for loop
        batch <- paste("batch", batch.no, sep="")   
        mc0 <- paste("mcrun0/", "genepop_ind", gen, sep="")
        mc1 <- paste("mcrun1/", "genepop_ind", gen, sep="")
        mc2 <- paste("mcrun2/", "genepop_ind", gen, sep="")
      
      #Read in Data 
        genmc0 <- read_population(paste(data.dir, batch, "/batchrun", batch.no, mc0, file.ext, sep=""), type="genepop")
        genmc1 <- read_population(paste(data.dir, batch, "/batchrun", batch.no, mc1, file.ext, sep=""), type="genepop")
        genmc2 <- read_population(paste(data.dir, batch, "/batchrun", batch.no, mc2, file.ext, sep=""), type="genepop")
        
        genmc0$Pop <- genmc0$Population
        genmc1$Pop <- genmc1$Population
        genmc2$Pop <- genmc2$Population
      
      #Write strata data
        for(id in 1:length(genmc0$Pop)){
          genmc0$Population[id] <- as.numeric(unlist(strsplit(genmc0$Pop[id], "Pop-")))[2]     #removes the pop- identifier in the genepop files
        }
        
        for(id in 1:length(genmc1$Population)){
          genmc1$Population[id] <- as.numeric(unlist(strsplit(genmc1$Pop[id], "Pop-")))[2]
        }
        
        for(id in 1:length(genmc2$Population)){
          genmc2$Population[id] <- as.numeric(unlist(strsplit(genmc2$Pop[id], "Pop-")))[2]
        }
      
      ##Join watershed information with genetic data:
        gen.mc0 <- join(x=genmc0, y=wat.df, by="Population", type="left")
        gen.mc1 <- join(x=genmc1, y=wat.df, by="Population", type="left")
        gen.mc2 <- join(x=genmc2, y=wat.df, by="Population", type="left")
        
        gen.mc0$Population <- gen.mc0$trib
        gen.mc1$Population <- gen.mc1$trib
        gen.mc2$Population <- gen.mc2$trib
      
      
      ##OPTIONAL -- write to csv for external programs/analysis:
        write.csv(x = gen.mc0, file = paste(data.dir, batch, mc0, ".csv", sep=""))
        write.csv(x = gen.mc1, file = paste(data.dir, batch, mc1, ".csv", sep=""))
        write.csv(x = gen.mc2, file = paste(data.dir, batch, mc2, ".csv", sep=""))
        
      ##Convert to allele frequency matrix:
        mv0 <- to_mv(gen.mc0, drop.allele = TRUE)
        mv1 <- to_mv(gen.mc1, drop.allele = TRUE)
        mv2 <- to_mv(gen.mc2, drop.allele = TRUE)
      
      ##OPTIONAL: write allele frequency matrix to csv:
        write.csv(x = mv0, file = paste(data.dir, batch, mc0, "_AlleleFreqMat.csv", sep=""))
        write.csv(x = mv1, file = paste(data.dir, batch, mc1, "_AlleleFreqMat.csv", sep=""))
        write.csv(x = mv2, file = paste(data.dir, batch, mc2, "_AlleleFreqMat.csv", sep=""))
      
      #-----------------------------------------------------------------------------------------------------------------
      #Non-metric multidimensional scaling (NMDS)
      #------------------------------------------------------------------------------------------------------------------
      
      ##CHOOSE ONE -- drop rare alleles; 5% cutoff is most conservative & recommended:
        mv.r0 <- drop.var(mv0, min.po=1) # alleles dropped if in <1% of individuals
        #mv.r0 <- drop.var(mv0, min.fo=10) #alleles dropped if in <10 individuals
        mv.r1 <- drop.var(mv1, min.po=1) # alleles dropped if in <1% of individuals
        #mv.r1 <- drop.var(mv1, min.fo=10) #alleles dropped if in <10 individuals
        mv.r2 <- drop.var(mv2, min.po=1) # alleles dropped if in <1% of individuals
        #mv.r2 <- drop.var(mv2, min.fo=10) #alleles dropped if in <10 individuals
      
      ##scree plot to determine #dimensions (may take a few hours to run with >1000 individuals):
        nmds.scree(mv.r0, distance="bray", k=10, autotransform=FALSE, trymax=10)
        nmds.scree(mv.r1, distance="bray", k=10, autotransform=FALSE, trymax=10)
        nmds.scree(mv.r2, distance="bray", k=10, autotransform=FALSE, trymax=10)
      
      ##set K to # dimensions determined from scree plot (will take a few minutes):
        gen.nmds0 <-metaMDS(mv.r0, distance="bray", k=2, autotransform=FALSE, trymax=10)
        gen.nmds1 <-metaMDS(mv.r1, distance="bray", k=2, autotransform=FALSE, trymax=10)
        gen.nmds2 <-metaMDS(mv.r2, distance="bray", k=2, autotransform=FALSE, trymax=10)
      
      
      ##OPTIONAL -- quick plot:
        plot(gen.nmds0, "sites")
        plot(gen.nmds1, "sites")
        plot(gen.nmds2, "sites")
      
      ##create dataframe with NMDS scores, watershed data, and individual ID:
        NMDSpops0 <- cbind(gen.mc0$Pop, gen.mc0$watershed, gen.mc0$extbar, 
                           gen.mc0$ColID, as.vector(gen.mc0$RGB), gen.nmds0$points)
        id <- rownames(NMDSpops0)
        NMDSpops0.2 <- cbind(id=id, NMDSpops0)
        NMDSpops0 <- data.frame(NMDSpops0.2)
        names(NMDSpops0)[names(NMDSpops0)=="V6"] <- "RGB"      #for following merge function
        names(NMDSpops0)[names(NMDSpops0)=="V5"] <- "ColID"    #for following merge function
        
        NMDSpops1 <- cbind(gen.mc1$Pop, gen.mc0$watershed, gen.mc1$extbar, 
                           gen.mc1$ColID, as.vector(gen.mc0$RGB), gen.nmds1$points)
        id <- rownames(NMDSpops1)
        NMDSpops1.2 <- cbind(id=id, NMDSpops1)
        NMDSpops1 <- data.frame(NMDSpops1.2)
        names(NMDSpops1)[names(NMDSpops1)=="V6"] <- "RGB"      #for following merge function
        names(NMDSpops1)[names(NMDSpops1)=="V5"] <- "ColID"    #for following merge function
        
        NMDSpops2 <- cbind(gen.mc2$Pop, gen.mc2$watershed, gen.mc2$extbar, 
                           gen.mc2$ColID, as.vector(gen.mc2$RGB), gen.nmds2$points)
        id <- rownames(NMDSpops2)
        NMDSpops2.2 <- cbind(id=id, NMDSpops2)
        NMDSpops2 <- data.frame(NMDSpops2.2)
        names(NMDSpops2)[names(NMDSpops2)=="V6"] <- "RGB"      #for following merge function
        names(NMDSpops2)[names(NMDSpops2)=="V5"] <- "ColID"    #for following merge function
        
      
      ##OPTIONAL -- write out NMDSpops:
        write.csv(x = NMDSpops0, file= paste(data.dir, batch, mc0, "_NMDSpops.csv", sep=""))
        write.csv(x = NMDSpops1, file= paste(data.dir, batch, mc1, "_NMDSpops.csv", sep=""))
        write.csv(x = NMDSpops2, file= paste(data.dir, batch, mc2, "_NMDSpops.csv", sep=""))
      
      
      ##ANOSIM (takes a while)
      #-------------------------------------------------------
        gen0.d <- vegdist(mv.r0, "bray")
        gen0.anosim <- anosim(gen0.d, barrier)    
        summary(gen0.anosim)
        plot.anosim(gen0.anosim)
        
        gen1.d <- vegdist(mv.r1, "bray")
        gen1.anosim <- anosim(gen1.d, barrier)
        summary(gen1.anosim)
        plot.anosim(gen1.anosim)
        
        gen2.d <- vegdist(mv.r2, "bray")
        gen2.anosim <- anosim(gen2.d, barrier)
        summary(gen2.anosim)
        plot.anosim(gen2.anosim)
    }
  }
#}                                                      #ignore this line for now - used to loop through folders but JJB changed input file structure to remove 3rd for loop
  end.time = Sys.time()                                                       #end timer
  round(difftime(end.time, start.time, units='mins'), dig = 2)                #prints elapsed time 


## JAMESON AND MERYL, I HAVE NOT PLAYED WITH ANYTHING BELOW THIS POINT YET *****************************************


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