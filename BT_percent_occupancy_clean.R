# -----------------------------------------------------------------------------
# Find Private Alleles in Bull Trout (BT) Riverscape Genetics Project
# 
# Script by Jake Burkhart (JJB) and Jameson Hinkle (JEH)
# -----------------------------------------------------------------------------

# Load packages & f(x)
# -----------------
require("gstudio")
require("plyr")
require("reshape2")
require("ggplot2")

rm(list = ls())

## Create multiplot function for plotting:
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {              ##MAKE SURE TO RUN THIS. OCCUPANCY PLOTS DON'T WORK OTHERWISE!
                library(grid)
                
                # Make a list from the ... arguments and plotlist
                plots <- c(list(...), plotlist)
                
                numPlots = length(plots)
                
                # If layout is NULL, then use 'cols' to determine layout
                if (is.null(layout)) {
                  # Make the panel
                  # ncol: Number of columns of plots
                  # nrow: Number of rows needed, calculated from # of cols
                  layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                                   ncol = cols, nrow = ceiling(numPlots/cols))
                }
                
                if (numPlots==1) {
                  print(plots[[1]])
                  
                } else {
                  # Set up the page
                  grid.newpage()
                  pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
                  
                  # Make each plot, in the correct location
                  for (i in 1:numPlots) {
                    # Get the i,j matrix positions of the regions that contain this subplot
                    matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
                    
                    print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                                    layout.pos.col = matchidx$col))
                  }
                }
              }
## -----------------


## Set-up Inputs Parameters: 
## -----------------
  # source("C:/Users/Jacob/Desktop/RivGen/Mani_Runs/private_alleles.r")  
  source("private_alleles.r")   ## put pathname where the "private_alleles.r" file is located within quotation marks

  data.dir <- "/Users/jamesonhinkle/Desktop/Existing_Rerun/"  ## change text within quotes to the pathname that houses all of the input data (batch folders)
  out.dir <- "/Users/jamesonhinkle/Desktop/folders/vcu/landscape genetics/R code/landscape/privateallele/"   ## change to the location in which you would like to output data for the run  ***CHANGE FOR EACH SURFACE***
  plot.dir <- out.dir                                ## location to output occupancy maps, occupancy rates plots, and private allele plots
  wshed <- read.csv("WatershedIDsublevels.csv", sep =",", header=T)  ## input watershed file to get barrier group information and RGB codes
  
  species <- "_BT"          
  date <- Sys.Date()
  file.ext <- ".gen"
  
  gen.seq <- c(0, 50, 150, 199) #50,150,199 #c(-1,0,1,2,3,4,5,10,15,20,25,30,40,50,60,70,80,90,91,92,93,94,95,96,97,98,99,100,110,120,130,140,150,160,170,180,190,199)
  batch.seq <- c(3,4)#0:11     #0-11 for full AD   ## 0-59 for riverine; 60-119 for existing barrs; 120-179 for future barriers  ***CHANGE FOR EACH SURFACE***
  mcrun.seq <- 0:9                          ## starting and ending value for mcruns within each batch (e.g. mcrun0 through mcrun2 --> 0:2)
  
  ## subset watershed data frame to only the columns of interest
  wat.df <-data.frame(Population=as.numeric(wshed$FID)+1, ## corrects FID's to corrspond with the patch IDs in popvars and genepop files
                      watershed=wshed$watershed,          ## holds watershed group information --> "Sullivan," "Le Clerc," "Tacoma," or "PendOreilleMS"
                      patch.x=wshed$POINT_X, patch.y=wshed$POINT_Y,
                      extbar=as.character(wshed$extbarrsub), futbar=as.character(wshed$futbarrsub), ##read these columns in as character vector so that they can be used as a strata field. won't allow strata fields using numeric vectors. 
                      ColID=wshed$ColorID, RGB=wshed$RGB)   ## houses hexadecimal values used for plotting
  
  barrier <- as.factor(wat.df$ColID)                                                     ##this line may or may not be used. could be holdover from the NMDS code -- JJB       
  popvars <- read.csv( "Patchvars_BT_FullAD_2K_EM_0001.csv", sep = ",", header=T)  ##read in popvars file, used to count the initial stocking locations
  popvars$Population <- popvars$Subpopulation      ##create join field 
  
  barr <- read.csv("Barriers20160615.csv", sep = ",", header = T)
  
  barr$bar.col <- ifelse(barr$Passibil_1 == 0, yes="#d95f0e", no=ifelse(barr$Passibil_1 == 33, yes= "#fe9929", no=ifelse(barr$Passibil_1 == 67, yes="#fec44f", no=NA)))
  
  wat.df <- join(wat.df, popvars, by="Population", type="left")  ##join the popvars and wat.df for easier reference later
# -----------------

  
## FOR RIVERINE AND EXISTING BARRIERS RUN THIS SECTION!  -->  REMEMBER TO CHANGE OUTPUT DIRECTORY!
## Map Occupancy, calculate occupancy rates, calculate and map private alleles across batches, mcruns, and generations 
# -----------------
  start.time = Sys.time() # Set timer for start of run
  
  #Set up output data frames
  df.length <- rep(NA, length.out=(length(batch.seq)))     ## Calculate output dataframe length
  occ.df <- data.frame(Population = wat.df$Population, patch.x=wat.df$patch.x, patch.y=wat.df$patch.y, watershed=wat.df$watershed)          ## Stores occupancy data

  ## Start calculations
  for(j in 1:length(batch.seq)){  ##loop over batches
    occ.count <- data.frame(Population=wat.df$Population)
    
    for(m in 1:length(mcrun.seq)){ ##loop over mc runs within batches
      count.df <- data.frame(Population=wat.df$Population)     ##create a dataframe to store 1/0 occupancy info per batch
      
      for(k in 1:length(gen.seq)){ ##loop over generations within mcruns
        # Set up input parameters
        batch.no <- batch.seq[j]
        gen <- gen.seq[k]
        mc.no <- mcrun.seq[m]                 
        
        # Prepare input file objects  
        batch <- paste("Batch", batch.no, sep="")   
        mc0 <- paste("mcrun", mc.no, "/batchrun0mcrun0/", "genepop_ind", gen, sep="")         ### HAD TO ADD AN EXTRA "batchrun0mcrun0" DUE TO NEW FORMAT IN THE ZIPPED STUFF. MAY REMOVE LATER!

       ## CHECK IF GENEPOP FILE EXISTS FOR GENERATION
        if(file.exists(paste(data.dir, batch, "/batchrun", batch.no, mc0, file.ext, sep = ""), type = "genepop")) { 
          #Read in Data
          gen.df <- read_population(paste0(data.dir, batch, "/batchrun", batch.no, mc0, file.ext), type="genepop")  ##read in genepop file
          gen.df$Pop <- gen.df$Population
          
          #Write strata data
          for(id in 1:length(gen.df$Pop)){
            gen.df$Population[id] <- as.numeric(unlist(strsplit(gen.df$Pop[id], "Pop-")))[2]     #removes the pop- identifier in the genepop files
          }
          
          ##Join watershed information with genetic data:
          gen.data <- join(x=gen.df, y=wat.df, by="Population", type="left")
          gen.data$Population <- as.factor(gen.data$Population)
          
          ##Count number inds per patch and 
          t.occ <- data.frame(Population = table(gen.data$Population))        ## count # inds per patch
          colnames(t.occ) <- c("Population", "N")                             ## rename columns for easy join

          t.occ$Population <- as.numeric(as.character(t.occ$Population))      ## 
          occ.data <- join(x=wat.df, y=t.occ, by="Population", type="left")
          # occ.data$
         
          ##Calculate Patch Occupancy
          count.df[,k+1] <- ifelse(is.na(occ.data$N), yes=0, no=ifelse(occ.data$N==0, yes=0, no=1))
          colnames(count.df)[k+1] <- paste0("Occ_gen", gen)
           
          # iter <- iter + 1
          
        } else {            ## Outputs place holder data and zero-occupancy plots for generations with missing genepop files
            print(paste0("file is missing, populations crashed at generation ", gen, " for batch ", batch.no, " MC run ", mc.no))
          
          ##Output Patch Occupancy
            count.df[, k+1] <- 0
            
            # iter <- iter + 1
          }  ## END IF/ELSE STATEMENT TO IGNORE CALCULATIONS ON GENERATIONS MISSING GENEPOP FILES
      
      } ##END GENERATION LOOP
      
      ## Count how many times a patch is occupied in the mcrun
      occ.count[, m+1] <- rowSums(count.df[,2:dim(count.df)[2]], na.rm=T) / length(gen.seq)
      colnames(occ.count)[m+1] <- paste0("MC", mc.no, "_perc_occ")
 
      ## OPTIONAL: output data frame for mcrun
      # write.csv(occ.count, paste(out.dir, "percent_occupancy_by_mcrun_batch_", batch.no, "_", date, ".csv", sep =""))
      
      ##Track progress
      print(paste("Batch", batch.no, ", mcrun", mc.no, "-", "time elapsed:", round(difftime(Sys.time(), start.time, units='mins'), dig = 2), "min"))  ##print progress
      
    } ## END MCRUN LOOP
    
    ## Summarize Occupancy Across MCRUNS
    occ.df[, j+4] <- rowSums(occ.count[,2:dim(occ.count)[2]], na.rm=T) / length(mcrun.seq)
    colnames(occ.df)[j+4] <- paste0("Batch", batch.no, "_perc_occ")
    occ.df$RGB <- ifelse(occ.df$watershed == "Sullivan", yes="#34ABFA", 
                         no=ifelse(occ.df$watershed == "Tacoma", yes="#8DD45D",
                                   no=ifelse(occ.df$watershed == "LeClerc", yes="#E3595D",
                                             no=ifelse(occ.df$watershed == "PendOreilleMS", yes="#9F44E9", NA))))
    occ.df$RGB <- ifelse(occ.df[,j+4]>0, yes=as.character(occ.df$RGB), no="#D3D3D3")
    
    ##Output MAPS <------------- NOTE THIS SECTION CAN BE TURNED OFF IF WE DON'T WANT PLOTS THAT SHOW ALL SAMPLING LOCATIONS!
    ##Graph the occupancy across the watershed
    
    occ.plot <- ggplot(wat.df, aes(x=patch.x, y=patch.y)) + 
      theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
      geom_point(data=wat.df, aes(x=patch.x, y=patch.y), size=3, color="light grey") +                    # plot all sites in grey
      geom_point(data=occ.df, aes(x=patch.x, y=patch.y), color=as.character(occ.df$RGB), size=(occ.df[,j+4]*4)) + geom_point(data=barr, size=3.5, pch=15, color=barr$bar.col, aes(x=POINT_X, y=POINT_Y)) +        
      labs(title=paste("Average Percent Occupied Patches - Batch ", batch.no, sep=""), x="Easting", y="Northing") + 
      theme(legend.position="none")
    
    # plot all occupied sites (>= 1 ind in genepop file)
    
    ##Save plot to file:  
     png(paste(plot.dir, "Percent_Occupancy_Batch", batch.no, "_", date ,".jpeg", sep=""), width=2000, height=3000, units='px', res=250)
      print(occ.plot)
     dev.off()
    
    ## Track Progress
    print(paste("Batch", batch.no, " Complete -", "time elapsed:", round(difftime(Sys.time(), start.time, units='mins'), dig = 2), "min"))  ##print progress
    print("____________________________________________________________")
    
  } ## END BATCH LOOP
  
  write.csv(occ.df, paste(out.dir, "occupancy_percentages_all_batches", date, ".csv", sep=""))
  end.time = Sys.time()                                                       #end timer
  round(difftime(end.time, start.time, units='mins'), dig = 2)                #prints elapsed time 
# -----------------
  
