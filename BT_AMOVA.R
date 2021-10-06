# AMOVA script using gstudio
# runs AMOVA across some number of MCMC runs, and averages output value for a mean cluster coeffficitnt

require(gstudio)
require(pegas)
require(plyr)


# Manual, iterative runs because I gave up on life...
# -------------------------------------------------------------------------------------------------------------------------------
#Importing Barrier runs:
  data.dir <- "/Users/jamesonhinkle/Desktop/GenepopFiles_Barriers/"
  #data.dir <- "C:/Users/Jacob/Desktop/RivGen/Meryl_Talk/Meryl_Runs_Genepop_Files/GenepopFiles_Barriers/"
  scenario <- "Barriers"
  species <- "_BT"
  date <- Sys.Date()
  file.ext <- ".gen"

#Importing Riverine runs:
  #data.dir <- "C:/Users/semlitschlab/Jake/RivGen/Meryl_Talk/Meryl_Riverine_Runs/GenepopFiles_Riverine/"
  #data.dir <- "C:/Users/Jacob/Desktop/RivGen/Meryl_Talk/Meryl_Barrier_Runs/GenepopFiles_Riverine/"
  #scenario <- "Riverine"
  #species <- "_BT"
  #date <- Sys.Date()
  #file.ext <- ".gen"

#Generations to monitor
  gen.seq <- c(-1,0,1,10,150)
  batch.seq <- 0:5

# read in watershed information
  #wshed <- read.csv("C:/Users/Jacob/Desktop/RivGen/Meryl_Talk/PatchXYwatershed.csv", header=T)
  wshed <- read.csv("PatchXYwatershed.csv", header=T)
  wat.df <- data.frame(Population=as.numeric(wshed$FID)+1, POINT_X=wshed$POINT_X, POINT_Y=wshed$POINT_Y, trib=wshed$watershed)


  # n.mc <- mcrun[k]
  #batch.no <- 0  
  #batch.no <- 1
  #batch.no <- 2
  #batch.no <- 3                                
  #batch.no <- 4
  #batch.no <- 5
  #gen <- 150
  #gen <- 150                                

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
  
  
      # Genetic Distance Measures
        
      sink(paste(data.dir, scenario, "_batch", batch.no, "_gen", gen, "_AMOVA.txt", sep=""))
        print("Batch #")
        print(batch.no)
        print("Generation of Interest")
        print(gen)
  
        #AMOVA by population on 1st MC run
        D.amova.mc0 <- genetic_distance(gen.mc0, stratum="Population", mode="AMOVA")
        pops.mc0 <- as.factor(gen.mc0$Population)
        fit.mc0 <- amova(D.amova.mc0 ~ pops.mc0, nperm=0)
        print("fit.mc0 - AMOVA table for MC0")
        print(fit.mc0)
        #PhiST.mc0 <- fit.mc0$varcoef / (fit.mc0$varcoef + fit.mc0$varcomp[2,1])        #turn on if you DO run permutations
        PhiST.mc0 <- fit.mc0$varcoef / (fit.mc0$varcoef + fit.mc0$varcomp[2])          #turn on if you DON'T run permutations
        print("Ph1ST.mc0 - Proportion of Variance due to different populations (watersheds)")
        print(PhiST.mc0)
  
        
        #AMOVA by population on 2nd MC run
        D.amova.mc1 <- genetic_distance(gen.mc1, stratum="Population", mode="AMOVA")
        pops.mc1 <- as.factor(gen.mc1$Population)
        fit.mc1 <- amova(D.amova.mc1 ~ pops.mc1, nperm=0)
        print("fit.mc1 - AMOVA table for MC1")
        print(fit.mc1)
        #PhiST.mc1 <- fit.mc1$varcoef / (fit.mc1$varcoef + fit.mc1$varcomp[2,1])           #turn on if you DO run permutations
        PhiST.mc1 <- fit.mc1$varcoef / (fit.mc1$varcoef + fit.mc1$varcomp[2])             #turn on if you DON'T run permutations
        print("Ph1ST.mc1 - Proportion of Variance due to different populations (watersheds)")
        print(PhiST.mc1)                                                                         # a proportion of the difference (variance) due to individuals being in different populations
        
  
        # AMOVA by population on 3rd MC run
        D.amova.mc2 <- genetic_distance(gen.mc2, stratum="Population", mode="AMOVA")
        pops.mc2 <- as.factor(gen.mc2$Population)
        fit.mc2 <- amova(D.amova.mc2 ~ pops.mc2, nperm=0)      
        print("fit.mc2 - AMOVA table for MC2")
        print(fit.mc2)
        #PhiST.mc2 <- fit.mc2$varcoef / (fit.mc2$varcoef + fit.mc2$varcomp[2,1])         #turn on if you DO run permutations
        PhiST.mc2 <- fit.mc2$varcoef / (fit.mc2$varcoef + fit.mc2$varcomp[2])            #turn on if you DON'T run permutations
        print("Ph1ST.mc2 - Proportion of Variance due to different populations (watersheds)")
        print(PhiST.mc2)
        
  
        #Average variance
        AMOVA.mu <- mean(PhiST.mc0, PhiST.mc1, PhiST.mc2)
        print("AMOVA.mu - Average Proportion of Variance due to different populations (watersheds)")
        print(AMOVA.mu)        
  
      sink()
    }
  }

    #print(c(fit.mc0, PhiST.mc0, fit.mc1, PhiST.mc1, fit.mc2, PhiST.mc2, AMOVA.mu))
# ---------------------------------------------------------------






# IGNORE ALL THE STUFF BELOW 



# SLOW AS H-E-Double Hockey Sticks TRIPLE FOR LOOPS! Barf...
# ---------------------------------------------------------------------------------------------------------------------------
for(iter in 1:df.length){
  for(i in 1:n.gens){
    for(j in 1:n.batches){
      #for(k in 1:n.MCruns){
      #Prepare the Data
      batch.no <- batches[j]
      # n.mc <- mcrun[k]
      gen <- gens[i]
      
      batch <- paste("batchrun", batch.no, sep="")
      mc0 <- paste("mcrun0/", "genepop_ind", gen, sep="")
      mc1 <- paste("mcrun1/", "genepop_ind", gen, sep="")
      mc2 <- paste("mcrun2/", "genepop_ind", gen, sep="")
      
      # Read in Data 
      gen.mc0 <- read_population(paste(data.dir, batch, mc0, file.ext, sep=""), type="genepop")
      gen.mc1 <- read_population(paste(data.dir, batch, mc1, file.ext, sep=""), type="genepop")
      gen.mc2 <- read_population(paste(data.dir, batch, mc2, file.ext, sep=""), type="genepop")
      
      # Genetic Distance Measures
  
      #AMOVA by population on 1st MC run
      D.amova.mc0 <- genetic_distance(gen.mc0, stratum="Population", mode="AMOVA")
      pops.mc0 <- as.factor(gen.mc0$Population)
      fit.mc0 <- amova(D.amova.mc0 ~ pops.mc0, nperm=0)   
      #PhiST.mc0 <- fit.mc0$varcoef / (fit.mc0$varcoef + fit.mc0$varcomp[2,1])        #turn on if you DO run permutations
      PhiST.mc0 <- fit.mc0$varcoef / (fit.mc0$varcoef + fit.mc0$varcomp[2])          #turn on if you DON'T run permutations
      PhiST.mc0
      
      xMC0 <- D.amova.mc0[lower.tri(D.amova.mc0)]
      dfMC0 <- data.frame(vals=xMC0)
      ggplot(dfMC0) + geom_histogram(aes(x=vals), binwidth=1)
      
      
      #AMOVA by population on 2nd MC run
      D.amova.mc1 <- genetic_distance(gen.mc1, stratum="Population", mode="AMOVA")
      pops.mc1 <- as.factor(gen.mc1$Population)
      fit.mc1 <- amova(D.amova.mc1 ~ pops.mc1, nperm=0)
      #PhiST.mc1 <- fit.mc1$varcoef / (fit.mc1$varcoef + fit.mc1$varcomp[2,1])           #turn on if you DO run permutations
      PhiST.mc1 <- fit.mc1$varcoef / (fit.mc1$varcoef + fit.mc1$varcomp[2])             #turn on if you DON'T run permutations
      PhiST.mc1
      
      xMC1 <- D.amova.mc1[lower.tri(D.amova.mc1)]
      dfMC1 <- data.frame(vals=xMC1)
      ggplot(dfMC1) + geom_histogram(aes(x=vals), binwidth=1)
      
      
      # AMOVA by population on 3rd MC run
      D.amova.mc2 <- genetic_distance(gen.mc2, stratum="Population", mode="AMOVA")
      pops.mc2 <- as.factor(gen.mc2$Population)
      fit.mc2 <- amova(D.amova.mc2 ~ pops.mc2, nperm=0)      
      #PhiST.mc2 <- fit.mc2$varcoef / (fit.mc2$varcoef + fit.mc2$varcomp[2,1])         #turn on if you DO run permutations
      PhiST.mc2 <- fit.mc2$varcoef / (fit.mc2$varcoef + fit.mc2$varcomp[2])            #turn on if you DON'T run permutations
      PhiST.mc2
      
      xMC2 <- D.amova.mc1[lower.tri(D.amova.mc1)]
      dfMC2 <- data.frame(vals=xMC2)
      ggplot(dfMC2) + geom_histogram(aes(x=vals), binwidth=1)
      
      
      #Average variance
      AMOVA.mu <- mean(PhiST.mc0, PhiST.mc1, PhiST.mc2)
      AMOVA.mu
      
      # Write output to a dataframe
      AMOVA.results$batch[iter] <- batch.no
      AMOVA.results$generation[iter] <- gen
      AMOVA.results$PhiST.mc0[iter] <- PhiST.mc0
      AMOVA.results$PhiST.mc1[iter] <- PhiST.mc1
      AMOVA.results$PhiST.mc2[iter] <- PhiST.mc2
      AMOVA.results$AMOVA.mu[iter] <- AMOVA.mu
    }
  }
}
# ----------------------------------------------------------------------------------------------------------------

# IGNORE EVERYTHING BELOW HERE.
# -----------------------------------------------------------------------------------------------------------------------------


batch <- paste("batchrun", batch.no, sep="")
mc0 <- paste("mcrun0/", "genepop_ind", gen, sep="")
mc1 <- paste("mcrun1/", "genepop_ind", gen, sep="")
mc2 <- paste("mcrun2/", "genepop_ind", gen, sep="")

file.ext <- ".gen"

# Read in Data -- CHANGE THESE BACK TO MC1, MC2, MC3 LATER
gen.mc0 <- read_population(paste(data.dir, batch, mc0, file.ext, sep=""), type="genepop")
gen.mc1 <- read_population(paste(data.dir, batch, mc1, file.ext, sep=""), type="genepop")
gen.mc2 <- read_population(paste(data.dir, batch, mc2, file.ext, sep=""), type="genepop")


# Genetic Distance Measures
# -----
# Basic AMOVA distance matrix
D.amovaMC0 <- dist_amova(gen.mc0)
D.amovaMC1 <- dist_amova(gen.mc1)
D.amovaMC2 <- dist_amova(gen.mc2) 

low_d_muMC0 <- mean(D.amovaMC0[lower.tri(D.amovaMC0)])
low_d_muMC1 <- mean(D.amovaMC1[lower.tri(D.amovaMC1)])
low_d_muMC2 <- mean(D.amovaMC2[lower.tri(D.amovaMC2)])

#Plotting the histograms of AMOVA distances -- if you care about this?
xMC0 <- D.amovaMC0[lower.tri(D.amovaMC0)]
xMC1 <- D.amovaMC1[lower.tri(D.amovaMC1)]
xMC2 <- D.amovaMC2[lower.tri(D.amovaMC2)]

dfMC0 <- data.frame(valsMC0=xMC0)
dfMC1 <- data.frame(valsMC1=xMC1)
dfMC2 <- data.frame(valsMC2=xMC2)

ggplot(df) + geom_histogram(aes(x=valsMC0), binwidth=1)
ggplot(df) + geom_histogram(aes(x=valsMC1), binwidth=1)
ggplot(df) + geom_histogram(aes(x=valsMC2), binwidth=1)
# -----


#AMOVA by population on 1st MC run
D.amova.mc0 <- genetic_distance(gen.mc0, stratum="Population", mode="AMOVA")
pops.mc0 <- as.factor(gen.mc0$Population)
fit.mc0 <- amova(D.amova.mc0 ~ pops.mc0)

PhiST.mc0 <- fit.mc0$varcoef / (fit.mc0$varcoef + fit.mc0$varcomp[2,1])
PhiST.mc0


#AMOVA by population on 2nd MC run
D.amova.mc1 <- genetic_distance(gen.mc1, stratum="Population", mode="AMOVA")
pops.mc1 <- as.factor(gen.mc1$Population)
fit.mc1 <- amova(D.amova.mc1 ~ pops.mc1)

PhiST.mc1 <- fit.mc1$varcoef / (fit.mc1$varcoef + fit.mc1$varcomp[2,1])
PhiST.mc1


# AMOVA by population on 3rd MC run
D.amova.mc2 <- genetic_distance(gen.mc2, stratum="Population", mode="AMOVA")
pops.mc2 <- as.factor(gen.mc2$Population)
fit.mc2 <- amova(D.amova.mc2 ~ pops.mc2)

PhiST.mc2 <- fit.mc2$varcoef / (fit.mc2$varcoef + fit.mc2$varcomp[2,1])
PhiST.mc2


#Average variance
AMOVA.mu <- mean(PhisST.mc0, PhiST.mc1, PhiST.mc2)
AMOVA.mu