setwd("~/Dropbox/working/BD_pour_hist/codes/")
load("resultados_evolerror_new.Rdata")



df <- read.csv("datos2plot.csv", sep = "\t", dec = ",")
df$density <- factor(df$density, 
                     levels = unique(df$density),
                     labels = c("normal", "chi2", "mezcla1", 
                                "mezcla2", "bart", "triangular"))

library(ggplot2)
library(magrittr)
library(dplyr)

pdf("mise.pdf", width = 10)
ggplot(data = df, aes(x = size, y = MISE, col = method)) + 
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  scale_x_log10() + scale_y_log10() +
  facet_wrap(~ density, scales = "free_y") + theme_bw() +
  scale_color_discrete(name = "Method", 
                      breaks = c("H", "FP", "KDE", 
                                 "BagH", "BagFP", "BagKDE", "RASH"))
dev.off()

## Solo H, FP y KDE
#pdf("mise_simple.pdf", width = 10)
#dplyr::filter(df, method %in% c("H", "FP", "KDE")) %>%
#ggplot(aes(x = size, y = MISE, col = method)) + 
#  geom_point() +
#  geom_smooth(method = "lm", se = FALSE) +
#  scale_x_log10() + scale_y_log10() +
#  facet_wrap(~ density, scales = "free_y") + theme_bw() +
#  scale_color_discrete(name = "Method", 
#                       breaks = c("H", "FP", "KDE"))
#dev.off()



load(file = "resultados_mise-par.Rdata")
n <- rev(c(20, 50, 100, 200, 500, 1000, 2000, 5000, 10000))
df <- Reduce(rbind, res)
d <- rownames(df)
rownames(df) <- NULL
df <- data.frame(df, density = d, size = rep(n, each = 6))

df <- tidyr::gather(df, model = colnames(df)[1:7], "model", "MISE")

ggplot(data = df, aes(x = size, y = MISE, col = model)) + 
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  scale_x_log10() + scale_y_log10() +
  facet_wrap(~ density, scales = "free_y") + theme_bw() +
  scale_color_discrete(name = "model", 
                       breaks = c("H", "FP", "KDE", 
                                  "BagH", "BagFP", "BagKDE", "RASH"))

         
         
