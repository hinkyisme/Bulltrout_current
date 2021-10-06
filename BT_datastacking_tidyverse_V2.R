### Stacking data frames from BT output with tidyverse for future
### downstream use (e.g. private alleles, patch occupancy, F, etc.)
### J. Hinkle
### 6/26/20

# Dependencies
require(gstudio)
require(tidyverse) # until further notice, use dplyr version 0.8.5, 1.0.0 does not bind_rows()
source("private_alleles.r")
#before beginning, and until Rodney updates gstudio, make sure to add the following in
#trace(".read_cdpop", edit=TRUE, where = read_population)
#idx <- grep(paste("L",i,"A",sep=""),locus_names) 
#locus_names <- grepl("L[[:digit:]]+A[[:digit:]]+$",col_names)

data.dir <- "/Users/jamesonhinkle/Desktop/BT_test/" 

species <- "_BT"          
date <- Sys.Date()
file.ext <- ".csv"
out.dir <-   "/Users/jamesonhinkle/Desktop/R_Output_results/"

# paste mig or res ahead of file path and riv or barrier at end of file path
gen.seq <-  c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 170, 175, 180, 185, 190, 195, 200) #  
batch.seq <- 0 #c(6,66) # c(1,2)#0:11     #0-11 for full AD   ## 0-59 for riverine; 60-119 for existing barrs; 120-179 for future barriers  ***CHANGE FOR EACH SURFACE***
mcrun.seq <- 0    #                       ## starting and ending value for mcruns within each batch (e.g. mcrun0 through mcrun2 --> 0:2)

# -----------------
start.time = Sys.time() # Set timer for start of run

#Set up output data frames
start.time = Sys.time()
#gen.tidy <- tibble() ## create empty gen.tidy data frame for storing genetic information
t.tidy <- tibble()



## Build tidyverse data frame with all data for resistance surface also, map_df()? look it up. broom:: for clean output table of model
for(j in 1:length(batch.seq)){ ## loop over batches (demographic scenarios; e.g., 2K, Full AD, 0.001 stray rate)
  
  for(m in 1:length(mcrun.seq)){ ##loop over mc runs within batches
    
    for(k in 1:length(gen.seq)){ ##loop over generations within mcruns
      # Set up input parameters
      batch.no <- batch.seq[j]
      gen <- gen.seq[k]
      mc.no <- mcrun.seq[m]    
      
      #batch <- paste("Batch", batch.no, sep="")   
      mc0 <- paste("batchrun", batch.no, "mcrun", mc.no, "/", "ind", gen, sep = "")
      
      if(file.exists(paste(data.dir, "mig_mp_canyon10", "/", mc0, file.ext, sep = ""))) {    #  batch.no,  batch,    ## check to see if genepop file exists --> if yes, read in file and process it
        
        # read in data and gather private alleles for these data (if you want it) 
        gen.df <- read_population(paste(data.dir, "mig_mp_canyon10", "/", mc0, file.ext, sep = ""), type="cdpop")         ## read in genepop file
        
        gen.df$PatchID <- as.factor(gen.df$PatchID)
        
        t <- priv.allele(gen.df, "PatchID")
        
        t_inc <- t %>% mutate(Batch = batch.seq[j], mcrun = mcrun.seq[m], gen = gen.seq[k])
        
        t.tidy <- t_inc %>% bind_rows(t.tidy)
        
       # }
        
        #stack all data into one data frame.  Be aware that running both private alleles and stacking
        #data is very computationally expensive.
        #gen.tidyish <- tibble(gen.df) 
        #gen.tidyishish <- inner_join(gen.tidyish, wat.df, by = "Population")
        ## convert genepop file to data frame
        
        #gen.tidy.inc <- gen.tidyish %>% mutate(Batch = batch.seq[j], mcrun = mcrun.seq[m], gen = gen.seq[k])              ## add batch, mcrun, generation information
        
        #gen.tidy <- gen.tidy.inc %>% bind_rows(gen.tidy, gen.tidy.inc)                                                                  ## make a tidyverse file to store ALL OF THE DATA!
      } ## end if statement 
      
      
      
      
      print(gen)
    } ##END GENERATION LOOP
    
    
    ##Track progress
    print(paste("Batch", batch.no, ", mcrun", mc.no, "-", "time elapsed:", round(difftime(Sys.time(), start.time, units='mins'), dig = 2), "min"))  
  } ## END MCRUN LOOP
  
  
  ## Track Progress
  print(paste("Batch", batch.no, " Complete -", "time elapsed:", round(difftime(Sys.time(), start.time, units='mins'), dig = 2), "min")); print("___________________________________")
} ## END BATCH LOOP

# create gen.df$PatchID <- as.factor(gen.df$PatchID) and run priv.allele in loop
# after read_population

## Looking at private allele tallys, read in files first
## get rid of useless weird, excel "ghost column" (X1)
## group_by() and tally()

Res_SO8_riv <- Res_SO8_riv %>% select(-X1)
Mig_SO6_bar <- Mig_SO6_bar %>% group_by(Stratum, Locus, Allele, mcrun, gen)
Mig_SO6_bar_tot <- Mig_SO6_bar %>% tally()
Mig_SO6_riv <- Mig_SO6_riv %>% group_by(Stratum, Locus, Allele, mcrun, gen)
Mig_SO6_riv_tot <- Mig_SO6_riv %>% tally()
Mig_SO8_bar <- Mig_SO8_bar %>% group_by(Stratum, Locus, Allele, mcrun, gen)
Mig_SO8_bar_tot <- Mig_SO8_bar %>% tally()
Mig_SO8_riv <- Mig_SO8_riv %>% group_by(Stratum, Locus, Allele, mcrun, gen)
Mig_SO8_riv_tot <- Mig_SO8_riv %>% tally()
Res_SO6_bar <- Res_SO6_bar %>% group_by(Stratum, Locus, Allele, mcrun, gen)
Res_SO6_bar_tot <- Res_SO6_bar %>% tally()
Res_SO6_riv <- Res_SO6_riv %>% group_by(Stratum, Locus, Allele, mcrun, gen)
Res_SO6_riv_tot <- Res_SO6_riv %>% tally()
Res_SO8_bar <- Res_SO8_bar %>% group_by(Stratum, Locus, Allele, mcrun, gen)
Res_SO8_bar_tot <- Res_SO8_bar %>% tally()
Res_SO8_riv <- Res_SO8_riv %>% group_by(Stratum, Locus, Allele, mcrun, gen)
Res_SO8_riv_tot <- Res_SO8_riv %>% tally()

# annotate column for scenario
Mig_SO8_bar_tot <- Mig_SO8_bar_tot %>% mutate(scenario = "Mig_SO8_bar")
Mig_SO8_riv_tot <- Mig_SO8_riv_tot %>% mutate(scenario = "Mig_SO8_riv")
Res_SO8_riv_tot <- Res_SO8_riv_tot %>% mutate(scenario = "Res_SO8_riv")
Res_SO6_bar_tot <- Res_S06_bar_tot %>% mutate(scenario = "Res_SO6_bar")
Res_SO6_bar_tot <- Res_SO6_bar_tot %>% mutate(scenario = "Res_SO6_bar")
Res_SO6_riv_tot <- Res_SO6_riv_tot %>% mutate(scenario = "Res_SO6_riv")
Mig_SO6_bar_tot <- Mig_SO6_bar_tot %>% mutate(scenario = "Mig_SO6_bar")
Mig_SO6_riv_tot <- Mig_SO6_riv_tot %>% mutate(scenario = "Mig_SO6_riv")

# bind_rows together of all scenarios and select for relevant columns
privatealleles_tally <- Mig_SO8_bar_tot %>% bind_rows(Mig_SO8_riv_tot, Res_SO8_riv_tot, Res_SO6_riv_tot, Res_SO6_bar_tot, Res_SO6_riv_tot, Mig_SO6_bar_tot, Mig_SO6_riv_tot)
privatealleles_tally <- Mig_SO8_bar_tot %>% select(Stratum, Locus, Allele, mcrun, gen, n, scenario)

# Need to bring in all raw files 
M6_riv <- read_csv("Mig_S06_riv_privatealleles.csv")
M6_bar <- read_csv("Mig_SO6_barrier_privatealleles.csv")
M8_bar <- read_csv("Mig_SO8_barrier_privatealleles.csv")
M8_riv <- read_csv("Mig_SO8_riv_privatealleles.csv")
R6_bar <- read_csv("Res_SO6_barrier_privatealleles.csv")
R6_riv <- read_csv("Res_SO6_riv_prviatealleles.csv")
R8_bar <- read_csv("Res_SO8_barrier_privatealleles.csv")
R8_riv <- read_csv("Res_SO8_riv_privatealleles.csv")
# mutate(SO = x, Migration = Res or Mig, distance = riv or bar) and group_by(mcrun, gen)
M6_riv <- M6_riv %>% mutate(S0 = 6, Migration = "Mig", distance = "riv")
M6_bar <- M6_bar %>% mutate(S0 = 6, Migration = "Mig", distance = "bar")
M8_bar <- M8_bar %>% mutate(S0 = 8, Migration = "Mig", distance = "bar")
M8_riv <- M8_riv %>% mutate(S0 = 8, Migration = "Mig", distance = "riv")
R6_bar <- R6_bar %>% mutate(S0 = 6, Migration = "Res", distance = "bar")
R6_riv <- R6_riv %>% mutate(S0 = 6, Migration = "Res", distance = "riv")
R8_bar <- R8_bar %>% mutate(S0 = 8, Migration = "Res", distance = "bar")
R8_riv <- R8_riv %>% mutate(S0 = 8, Migration = "Res", distance = "riv")
tib_all <- M6_riv %>% bind_rows(M6_bar, M8_bar, M8_riv, R6_bar, R6_riv, R8_bar, R8_riv)
tib_all$S0 <- as.factor(tib_all$S0)
tib_all$Migration <- as.factor(tib_all$Migration)
tib_all$distance <- as.factor(tib_all$distance)
tib_all <- tib_all %>% group_by(mcrun, gen)
# get separate mc's to plot mcruns by scenario of frequency of private alleles
mc0 <- tib_all %>% filter(mcrun == 0)
mc1 <- tib_all %>% filter(mcrun == 1)
mc2 <- tib_all %>% filter(mcrun == 2)
mc3 <- tib_all %>% filter(mcrun == 3)
mc4 <- tib_all %>% filter(mcrun == 4)
# plot frequency of private alleles through time.
p0 <- ggline(mc0, x = "gen", y = "Frequency", title = "Frequency of Private Alleles, Mc0", color = "S0", facet.by = "mcrun", add = "mean_se")
p0 + facet_grid(distance + Migration ~ .)
p1 <- ggline(mc1, x = "gen", y = "Frequency", title = "Frequency of Private Alleles, Mc1", color = "S0", facet.by = "mcrun", add = "mean_se")
p1 + facet_grid(distance + Migration ~ .)
p2 <- ggline(mc2, x = "gen", y = "Frequency", title = "Frequency of Private Alleles, Mc2", color = "S0", facet.by = "mcrun", add = "mean_se")
p2 + facet_grid(distance + Migration ~ .)
p3 <- ggline(mc3, x = "gen", y = "Frequency", title = "Frequency of Private Alleles, Mc3", color = "S0", facet.by = "mcrun", add = "mean_se")
p3 + facet_grid(distance + Migration ~ .)
p4 <- ggline(mc4, x = "gen", y = "Frequency", title = "Frequency of Private Alleles, Mc4", color = "S0", facet.by = "mcrun", add = "mean_se")
p4 + facet_grid(distance + Migration ~ .)

# maybe use this code (the face_grid() part) with juse one mce at a time.
# take this and add in all scenarios and get all tallys of all private alleles and give to
# Erin, Casey, and Jake.

