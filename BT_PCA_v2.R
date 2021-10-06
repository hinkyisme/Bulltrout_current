# Load packages & f(x)
# -----------------
require("gstudio")
require("plyr")
require("reshape2")
require("ggplot2")

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

data.dir <- "/Users/jamesonhinkle/Desktop/6FocalScenarios/"
wshed <- read.csv("WatershedIDsublevels.csv", sep = ",", header=T) 
plot.dir <- "/Users/jamesonhinkle/Desktop/folders/vcu/landscape genetics/R code/landscape/abundance_plots"

#Input/Output File Parameters
species <- "_BT"
date <- Sys.Date()
file.ext <- ".gen"

#Generations to monitor
#fol <- 0:11
gen.seq <- c(0, 10, 50)              #subset: c(-1,0,1,10,50,100,150); all gens c(-1,0,1,2,3,4,5,10,15,20,25,30,40,50,60,70,80,90,91,92,93,94,95,96,97,98,99,100,110,120,130,140,150,160,170,180,190,199)
batch.seq <- 50#(1,2,25,26,49,50) #:11       #batch 0-3: riverine; batch 4-7: barriers; batch 8-11: panmictic
mcrun.seq <- 0:9                           #starting and ending value for mcruns within each batch (e.g. mcrun0 through mcrun2)

source("private_alleles.r")

                                                    ##this line may or may not be used. could be holdover from the NMDS code -- JJB       
popvars <- read.csv("Patchvars_BT_FullAD_2K_EM_0001.csv", sep = ",", header=T)  ##read in popvars file, used to count the initial stocking locations
popvars$Population <- popvars$Subpopulation ##create join field 

wat.df <-data.frame(Population=as.numeric(wshed$FID)+1, ## corrects FID's to corrspond with the patch IDs in popvars and genepop files
                    watershed=wshed$watershed,          ## holds watershed group information --> "Sullivan," "Le Clerc," "Tacoma," or "PendOreilleMS"
                    patch.x=wshed$POINT_X, patch.y=wshed$POINT_Y,
                    extbar=as.character(wshed$extbarrsub), futbar=as.character(wshed$futbarrsub), ##read these columns in as character vector so that they can be used as a strata field. won't allow strata fields using numeric vectors. 
                    ColID=wshed$ColorID, RGB=wshed$RGB)

barrier <- as.factor(wat.df$ColID)

wat.df <- join(wat.df, popvars, by="Population", type="left")  ##join the popvars and wat.df for easier reference later

start.time = Sys.time() # Set timer for start of run

#Set up output data frames
df.length <- rep(NA, length.out=(length(batch.seq)*length(mcrun.seq)*length(gen.seq)))                 ##calculate output dataframe length
occ.df <- data.frame(batch=df.length, mcrun=df.length, gen=df.length, pop.occ.total=df.length,
                     pop.occ.stock=df.length, bar.occ.total=df.length, bar.occ.stock=df.length)  ##stores occupancy and allelic diversity output for each generation
iter <- 1  ##sets up the writing index for the occ.df object..one day i shall sit down and figure out the indexing math rather than quick and dirty ways (JJB)

##Set up denominators for occupancy calculations
total.patches <- length(unique(wat.df$Population))                  ## Finds the total number of patches in your system
stock.patches <- sum(wat.df$N0 > 0)                                 ## Finds how many patches were stocked (all patches with K > 1)
total.bargroup <- length(unique(wat.df$extbar))                     ## total number EXISTING barrier groups
stock.bargroup <- sum(table(wat.df$N0>0, wat.df$extbar)[2,] > 0)    ## number stocked EXISTING barrier groups


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
    #mc.df$batch[iter] <- batch.no
    #mc.df$generation[iter] <- gen
    #mc.df$mcrun[iter] <- mc.no
    
    # Prepare input file objects  
    #batch <- paste("batch", fol.no, batch.no, sep="")        #ignore this line for now - used to loop through folders but JJB changed input file structure to remove 3rd for loop
    batch <- paste("batch", batch.no, sep="")   
    mc0 <- paste("mcrun", mc.no, "/batchrun0mcrun0/", "genepop_ind", gen, sep="")         ### HAD TO ADD AN EXTRA "batchrun0mcrun0" DUE TO NEW FORMAT IN THE ZIPPED STUFF. MAY REMOVE LATER!
    
    #Read in Data 
    if(file.exists(paste(data.dir, batch, "/batchrun", batch.no, mc0, file.ext, sep = ""), type = "genepop")) { 
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
    #write.csv(x = gen.mc0, file = paste(data.dir, batch, "/batchrun", batch.no, mc0, ".csv", sep=""))
    
    ##Convert to allele frequency matrix:
    mv0 <- to_mv(gen.mc0, drop.allele = TRUE)
    
    ##Run PCA on sampled mcruns and plot
    fit.pca <- princomp(mv, cor = TRUE)
    
    pred <- predict(fit.pca)
    
    names(gen.mc0)
    
    df <- data.frame(PC1 = pred[, 1], PC2 = pred[, 2], Pop = gen.mc0$ID, Trib = gen.mc0$trib, Patch = gen.mc0$Pop)
    
    p <- ggplot(df) + geom_point(aes(x = PC1, y = PC2, shape = Trib, color = Pop),size = 3, alpha = 0.75)
    
    } else{
      print("file does not exist")
    }
      
  }
  #write.csv(x = p, file = paste(data.dir, batch, "/batchrun", batch.no, mc0, "_PCA.tiff", sep=""))
  }

warning("file does not exist")
    
    ##OPTIONAL: write allele frequency matrix to csv:
    #write.csv(x = mv0, file = paste(data.dir, batch, "/batchrun", batch.no, mc0, "_AlleleFreqMat.csv", sep=""))