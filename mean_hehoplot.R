heho <- read.csv("BT_hetero.csv", sep = ",", header = T)
require(reshape2)
require(ggplot2)
heho <- melt(heho, id = "Year")
heho
p <- ggplot(heho, aes(x = Year, y = value, color = variable)) + geom_line() + geom_point()

