# need physical distance
phys <- read.csv(file = "xy_dams_full.csv", sep = ",", header = TRUE)
require(gstudio)
require(ggplot2)
require(vegan)
# get necessary packages
# read population in
gen <- read_population("genepopgrid0.gen", type = "genepop", sep = ",", header = TRUE)
# get genetic distance
gen.d <- genetic_structure(gen, stratum = "ID", mode = "Gst_prime", pairwise = TRUE)
# get physical distance
coords <- strata_coordinates(phys)
phys.d <- strata_distance(coords, stratum = "Stratum")
# see if points in matrices are out of order
cbind(rownames(gen.d), rownames(phys.d))
# if so, order physical distance to genetic distance
phys <- phys[order(phys$ID), ]
coords <- strata_coordinates(phys)
phys.d <- strata_distance(coords, stratum = "Stratum")
# plot out genetic distance vs. physical distance
df <- data.frame(P = phys.d[lower.tri(phys.d)], G = gen.d[lower.tri(gen.d)])
ggplot(df) + geom_point(aes(x = P, y = G)) + xlab("Physical Distance") + ylab("Genetic Distance")
# run mantel test
## still working out kinks on mantel test
