mv <- to_mv(gen.mc0, drop.allele = TRUE)

fit.pca <- princomp(mv, cor = TRUE)

pred <- predict(fit.pca)

names(gen.mc0)

df <- data.frame(PC1 = pred[, 1], PC2 = pred[, 2], Pop = gen.mc0$ID, Trib = gen.mc0$trib, Patch = gen.mc0$Pop)

p <- ggplot(df) + geom_point(aes(x = PC1, y = PC2, shape = Trib, color = Pop),size = 3, alpha = 0.75)
# must change color to match Matt's RGB codes.  Make a new vector in df that corresponds Pop to color and use that for color in p
# need vector of colors from Matt (I think)

summary(gen.mc0$trib)
