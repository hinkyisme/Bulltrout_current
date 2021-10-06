##For post-sim processing of Bull Trout simulations (CDMetaPop) and NMDS analysis
##post-sim processing code: J. Hinkle and J. Burkhart; NMDS code: M. Mims

##Install most up to date gstudio
#install.packages("devtools")
require(devtools)
install_github("dyerlab/gstudio", ref = "develop", force=T)


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
    data.dir <- "/Users/jamesonhinkle/Desktop/6FocalScenarios/"
    wshed <- read.csv("WatershedIDsublevels_07122015.csv", sep = ",", header=T) 
    plot.dir <- "/Users/jamesonhinkle/Desktop/folders/vcu/landscape genetics/R code/landscape"

  #Jake
    #source("C:/Users/Jacob/Desktop/RivGen/Mani_Runs/biostats.R", encoding = 'UTF-8')
    #data.dir <- "C:/Users/Jacob/Desktop/RivGen/Mani_Runs/"    #should only need to change batch number between scenario runs   #used "Raw_Output07122015" before "Mani_Runs"
    #wshed <- read.csv("C:/Users/Jacob/Desktop/RivGen/Mani_Runs/WatershedIDsublevels_07122015.csv", header=T)   #update
    #plot.dir <- "C:/Users/Jacob/Desktop/RivGen/Mani_Runs/NMDS_plots/"
    
  # #Jake (Semlitsch Lab PC)
  #   source("C:/Users/semlitschlab/Jake/RivGen/Mani_Runs/biostats.R", encoding = 'UTF-8')
  #   data.dir <- "C:/Users/semlitschlab/Jake/RivGen/Mani_Runs/"    #should only need to change batch number between scenario runs
  #   wshed <- read.csv("C:/Users/semlitschlab/Jake/RivGen/Mani_Runs/WatershedIDsublevels_07122015.csv", header=T)   #update
  
  #Jake (Semlitsch Lab PC)
#     source("C:/Users/eggertlab/Jake/RivGen/Mani_Runs/biostats.R", encoding = 'UTF-8')
#     data.dir <- "C:/Users/eggertlab/Jake/RivGen/Mani_Runs/"    #should only need to change batch number between scenario runs
#     wshed <- read.csv("C:/Users/eggertlab/Jake/RivGen/Mani_Runs/WatershedIDsublevels_07122015.csv", header=T)   #update
#     
    
  #Input/Output File Parameters
  species <- "_BT"
  date <- Sys.Date()
  file.ext <- ".gen"

#Generations to monitor
  #fol <- 0:11
  gen.seq <- c(0)              #subset: c(-1,0,1,10,50,100,150); all gens c(-1,0,1,2,3,4,5,10,15,20,25,30,40,50,60,70,80,90,91,92,93,94,95,96,97,98,99,100,110,120,130,140,150,160,170,180,190,199)
  batch.seq <- 50#(1,2,25,26,49,50) #:11       #batch 0-3: riverine; batch 4-7: barriers; batch 8-11: panmictic
  mcrun.seq <- 0:9                           #starting and ending value for mcruns within each batch (e.g. mcrun0 through mcrun2)

#Create a watershed dataframe for future plotting:
  wat.df <-data.frame(Population=as.numeric(wshed$FID)+1, watershed=wshed$watershed, 
                      extbar=wshed$ExistingBarriersublevelID, futbar=wshed$FutureBarriersublevelID,
                      ColID=wshed$ColorID, RGB=wshed$RGB)

#Create barrier object for ANOSIM 
  barrier <- as.factor(wat.df$ColID)

#Iterate across generations and batches to create NMDS data
  start.time = Sys.time() # Set timer for start of run
  iter <- 1
  df.length <- rep(NA, length.out=(length(batch.seq)*length(gen.seq)))
  mc.df <- data.frame(batch=df.length, generation=df.length, mcrun=df.length)#, scenario=df.length)
  #for(i in 1:length(fol))                                          #ignore this line for now - used to loop through folders but JJB changed input file structure to remove 3rd for loop
    for(j in 1:length(batch.seq)){
      for(k in 1:length(gen.seq)){
        # Set up input parameters
          #fol.no <- fol[i]                                         #ignore this line for now - used to loop through folders but JJB changed input file structure to remove 3rd for loop
          batch.no <- batch.seq[j]
          gen <- gen.seq[k]
          mc.no <- sample(mcrun.seq, 1)
        
        # Create scenario object for output files
#           if(batch.no >= 0 & batch.no <= 3){
#             scenario <- "Riverine"
#           } else if(batch.no >= 4 & batch.no <= 7) {
#                 scenario <- "Barriers"
#               } else{
#                   scenario <- "Panmictic"
#                 }
        # Track the which mcrun we are analyzing for each scenario, batch, and generation.  
          #mc.df$scenario[iter] <- scenario
          mc.df$batch[iter] <- batch.no
          mc.df$generation[iter] <- gen
          mc.df$mcrun[iter] <- mc.no
          
        # Prepare input file objects  
          #batch <- paste("batch", fol.no, batch.no, sep="")        #ignore this line for now - used to loop through folders but JJB changed input file structure to remove 3rd for loop
          batch <- paste("batch", batch.no, sep="")   
          mc0 <- paste("mcrun", mc.no, "/batchrun0mcrun0/", "genepop_ind", gen, sep="")         ### HAD TO ADD AN EXTRA "batchrun0mcrun0" DUE TO NEW FORMAT IN THE ZIPPED STUFF. MAY REMOVE LATER!

        #Read in Data 
          gen.df <- read_population(paste(data.dir, batch, "/batchrun", batch.no, mc0, file.ext, sep=""), type="genepop")
          gen.df$Pop <- gen.df$Population
                  
        #Write strata data
          for(id in 1:length(gen.df$Pop)){
            gen.df$Population[id] <- as.numeric(unlist(strsplit(gen.df$Pop[id], "Pop-")))[2]     #removes the pop- identifier in the genepop files
          }
        
        ##Join watershed information with genetic data:
          gen.mc0 <- join(x=gen.df, y=wat.df, by="Population", type="left")
          gen.mc0$Population <- gen.mc0$trib
        
        
        ##OPTIONAL -- write to csv for external programs/analysis:
          write.csv(x = gen.mc0, file = paste(data.dir, batch, "/batchrun", batch.no, mc0, ".csv", sep=""))
          
        ##Convert to allele frequency matrix:
          mv0 <- to_mv(gen.mc0, drop.allele = TRUE)
        
        ##OPTIONAL: write allele frequency matrix to csv:
          write.csv(x = mv0, file = paste(data.dir, batch, "/batchrun", batch.no, mc0, "_AlleleFreqMat.csv", sep=""))
  
      #-----------------------------------------------------------------------------------------------------------------
      #Non-metric multidimensional scaling (NMDS)
      #------------------------------------------------------------------------------------------------------------------
        
        ##CHOOSE ONE -- drop rare alleles; 5% cutoff is most conservative & recommended:
          mv.r0 <- drop.var(mv0, min.po=1) # alleles dropped if in <1% of individuals
          #mv.r0 <- drop.var(mv0, min.fo=10) #alleles dropped if in <10 individuals

        ##scree plot to determine #dimensions (may take a few hours to run with >1000 individuals):
          nmds.scree(mv.r0, distance="bray", k=10, autotransform=FALSE, trymax=10)
        
        ##set K to # dimensions determined from scree plot (will take a few minutes):
          gen.nmds0 <-metaMDS(mv.r0, distance="bray", k=2, autotransform=FALSE, trymax=10)
        
        ##OPTIONAL -- quick plot:
          tiff(paste(data.dir, batch, "/batchrun", batch.no, mc0, "NMDS_quick_plot.tiff",sep=""),
               res=300, width=16.35, height=16.35, units="cm", compression=c("lzw"))   
            plot(gen.nmds0, "sites")
          dev.off()
        
        ##create dataframe with NMDS scores, watershed data, and individual ID:
          NMDSpops0 <- cbind(gen.mc0$Pop, gen.mc0$watershed, gen.mc0$extbar, 
                             gen.mc0$ColID, as.vector(gen.mc0$RGB), gen.nmds0$points)
          id <- rownames(NMDSpops0)
          NMDSpops0.2 <- cbind(id=id, NMDSpops0)
          NMDSpops0 <- data.frame(NMDSpops0.2)
          names(NMDSpops0)[names(NMDSpops0)=="V6"] <- "RGB"      #for following merge function
          names(NMDSpops0)[names(NMDSpops0)=="V5"] <- "ColID"    #for following merge function
          
        ##OPTIONAL -- write out NMDSpops:
          write.csv(x = NMDSpops0, file= paste(data.dir, batch, "/batchrun", batch.no, mc0, "_NMDSpops.csv", sep=""))
        
        ##OPTIONAL -- color plotting the NMDS 
          NMDS1 <- as.vector(NMDSpops0$MDS1)
          NMDS2 <- as.vector(NMDSpops0$MDS2)
          wshedcol <- as.vector(NMDSpops0$RGB)
          barrier <- as.factor(NMDSpops0$ColID)
          # NMDS1sing <- as.vector(NMDSpops_sing$MDS1)
          # NMDS2sing <- as.vector(NMDSpops_sing$MDS2)
          
          #export at PNG:
          png(paste(plot.dir, "NMDS_batch", batch.no, "_gen", gen, "_mc", mc.no,".png", sep=""), width=3000, height=3000, units='px', res=300)
          par(mfrow=c(1,1), oma=c(0,0,0,0), mar=c(3,3,2,2))
          plot(NMDS1, NMDS2, type="n", xlim=c(-1,1), ylim=c(-1,1), cex.axis=0.75,  xlab="", ylab="")
          ordihull(gen.nmds0, barrier, col="grey", lwd=0.75)
          points(NMDS1, NMDS2, pch=19, col="white")
          points(NMDS1, NMDS2, pch=19, col=alpha(wshedcol, 0.5))
          barrier <- as.vector(NMDSpops0$extbar)
          text(-0.9, 0.9, pos=4, paste("Batch ", batch.no, ", Generation ", gen, ", MC run", mc.no, sep=""), cex=0.75)
          dev.off()
          
        ##ANOSIM (takes a while)
        #-------------------------------------------------------
          #### THROWING AN ERROR AS OF 01 MARCH 2016 --- IGNORE FOR NOW
          
          gen0.d <- vegdist(mv.r0, "bray")
          gen0.anosim <- anosim(gen0.d, barrier)    
          # summary(gen0.anosim)
          # 
          # tiff(paste(data.dir, batch, "/batchrun", batch.no, mc0, "ANOSIM_quickplot.tiff",sep=""),
          #      res=300, width=16.35, height=16.35, units="cm", compression=c("lzw"))   
          #   plot.anosim(gen0.anosim)
          # dev.off()

        #update counter variable for mc.df output
        iter <- iter + 1
      }
    }
  #}                                                      #ignore this line for now - used to loop through folders but JJB changed input file structure to remove 3rd for loop
  end.time = Sys.time()                                                       #end timer
  round(difftime(end.time, start.time, units='mins'), dig = 2)                #prints elapsed time 

  write.csv(mc.df, paste(data.dir, "mcruns_used_for_NMDS", date, ".csv", sep=""))




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