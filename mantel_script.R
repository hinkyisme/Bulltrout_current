## mantel_script
## Jameson Hinkle
## Takes a CDPOP output and does a mantel correlation by generation of the output in a 
## particular batch run and a corresponding Cost Distance matrix.  
## Variables needed:
## path (gdmatrix output)
## batch (which batch number)
## cd_path (where the CD matrices are stored)
## rho.hat (number of permutations for permuation test for p value for r)
## Output is the r correlation data frame, and p_all is p value data frame
## You can print k to track the permutation and make sure this is working correctly by 
## un-noting print(k).  Same with i and j.


path <- ("Max_Dispersal_Sims_GD/Ibd/")
batch <- ("batchrun1")
path_temp <- paste(path, batch, sep = "")
output <- data.frame()
cd_path <- ("Full_CDmat_andXY")
rho.hat <- numeric(9999)
p_all <- data.frame()


gd_files <- list.files(path_temp, pattern="Gdmatrix", full.names = TRUE)
cd_files <- list.files(cd_path, pattern="ibd_full.csv", full.names = TRUE)


for(i in 1:length(gd_files)){
  gd_file <- gd_files[i]
  gd.data <- read.csv(gd_file,header=F)
  gen <- gd.data[ lower.tri(gd.data) ]
  #unlist and correlate
  gen <- as.vector(gen)
  
  
  
  for( j in 1:length(cd_files)){
    cd_file <- cd_files[j]
    cd.data <- read.csv(cd_file,header=F)
    cd <- cd.data[lower.tri(cd.data)]
    #unlist and correlate
    cd <- as.vector(cd)
    rho <- cor(gen, cd)
    
  
   
    
    
    output[i,j] <- rho
    
    
  
  for (k in 1:9999) {
    
    rho.rand <- sample(cd, length(gen),replace=FALSE)
    rho.hat[k] <- cor( gen, rho.rand )
    #print(k)
    
    
  }
  
 
  N <- length(which(rho.hat >= rho))
  
  p <- N/(10000)
  p_all[i, j] <- p
  #print(j)
  
}

print(i)
}


write.csv(p_all, file = "pval.csv")
write.csv(output, file = "rho.csv")
#  mantel(gen ~ cd, permutations = 999)

