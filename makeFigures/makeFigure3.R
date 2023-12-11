# make main figure

rm(list = ls())
# extrafont::loadfonts(device = "win")
library(tidyverse)
library(network)
library(ggthemes)
library(patchwork)
library(ggpubr)
library(viridis)
load("data/currentModels.rdata")
load("data/catAndVeg.RData")

nems <- c("actual", "resAdjacent", "resInner", "resMiddle", 
          "resOuter", "resExtreme")

distances <- c(3, 10, 25, 50, 150, 400)

getEffects <- function(modStack, diet = NULL){
  bestMod <- row.names((modStack$aic[1, ]))
  summ <- data.frame(summary(modStack[[bestMod]])$coefficients$cond)
  mood <- summ[row.names(summ) %in% nems, 1:2]
  colnames(mood) <- c("est", "stdE")
  mood$name <- row.names(mood)
  mood$dist <- distances[1:nrow(mood)]
  if(!is.null(diet)) mood$diet <- diet
  return(mood)
}

gg <- getEffects(generalsMods, "Generalists")
ss <- getEffects(specialsMods, "Specialists")

sg <- rbind(gg,ss)

ggplot(sg, aes(x = dist, y = exp(est))) +
  geom_errorbar(aes(ymin = exp(est - 1.96 * stdE), ymax = exp(est + 1.96 * stdE), group = diet), width = .1, position = position_dodge(.3)) +
  geom_point(aes(color = diet, shape = diet), size = 3, position = position_dodge(.3)) +
  scale_shape_manual(values = c(16,15), name = "Diet breadth") +
  stat_smooth(method = "lm", formula = y ~ poly(x, 2), se = F,
              aes(color = diet, linetype = diet), position = position_dodge(.3)) +
  scale_color_manual(values = c("#3A577E", "#6FBCEB"), name = "Diet breadth") +
  scale_linetype_manual(values = c("longdash", "solid"), name = "Diet breadth") +
  scale_y_continuous(breaks = c(.8, 1, 1.2, 1.4)) +
  scale_x_continuous(trans = "log", breaks = c(3, 10, 25, 50)) +
  theme_tufte() +
  theme(text = element_text(size=18, family = "Cambria"),
        # axis.title.x = element_blank(),
        # axis.ticks.x = element_blank(),
        # axis.text.x = element_blank(),
        axis.line = element_line(size = .2, color = "black"),
        axis.title.y = element_text(size=18,margin = unit(c(0, 5, 0, 0), "mm")),
        panel.grid.major.x = element_line(color = "grey", linetype = "dashed"),
        panel.grid.minor.x = element_line(color = "grey80", linetype = "dashed")) +
  ylab("Effect size (exp)") +
  xlab("\nMeasurement scale (m)") +
  geom_hline(yintercept = c(.8, 1.2, 1.4), linetype = "dashed", color = "grey80") +
  geom_hline(yintercept = 1, linetype = "dotted", color = "red")

ggsave("figures/mainFigure.png")
