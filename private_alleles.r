# ------------------------------------------------------------------------------------
# J. Hinkle (JEH) and J. Burkhart (JJB)
# Identify and calculate frequency of private alleles in a dataset
# ------------------------------------------------------------------------------------

require(gstudio)
require(ggplot2)
require(plyr)

#private allele function
#before running, make sure stratum id is as.factor()
priv.allele <- function(data, stratum){
  freqs.strata <- frequencies(data, stratum = stratum)
  loci.names <- unique(freqs.strata$Locus)

  private.alleles <- data.frame()
  
  for(i in 1:length(loci.names)){
    df <- subset.data.frame(freqs.strata, subset=Locus==loci.names[i])
    
    test <- as.data.frame(table(df$Allele))
    colnames(test) <- c("Allele", "Allele_Count")
  
    join <- plyr::join(df, test, type="left", by="Allele")
  
    p <- subset.data.frame(join, subset=Allele_Count==1)
  
    private.alleles <- rbind(private.alleles, p)          
  }
  return(private.alleles)
}




