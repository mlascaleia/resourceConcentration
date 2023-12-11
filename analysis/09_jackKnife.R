# the purpose of this script is to run the analysis for generalist and specialist models
# but each time I am going to leave out just one species

# initialize ####

rm(list = ls())
source("analysis/05_polyphaguosFunctions.R")
load("data/currentModels.rdata")

# make cat lists ####

specials <- rownames(electMatrix)[apply(electMatrix, 1, max) > .9]
specials <- c(specials, "MORREV", "CLADLI")
ggs <- c("MELACA", "NEMARE", "CROCNO", "IRIDEP", "LAMBFI", "MORRCO", "HIMEFI", "ACHADI", "ORTHXX", "EUTRCL", "ENNOSU", "ORGYLE")
all <- c(specials, ggs)

# run each model 
for(i in 1:length(specials)){
  assign(x = paste("not", specials[i], "mod", sep = "_"), 
         value = fullPolyphagousAnalysis(specials[!specials %in% specials[i]]))
  cat("completed ", specials[i], "\n", sep = "")
}

# run each model
for(i in 1:length(ggs)){
  assign(x = paste("not", ggs[i], "mod", sep = "_"), 
         value = fullPolyphagousAnalysis(ggs[!ggs %in% ggs[i]]))
  cat("completed ", ggs[i], "\n", sep = "")
}

# establish measurement distances
distances <- c(3, 10, 25, 50, 150, 400)

# extract specialist model effects
specialsJackKnife <- list(not_CLADLI_mod, not_EPIMHO_mod, not_MORREV_mod,
                          not_NOLACL_mod, not_NOLATR_mod, not_PSEUCO_mod,
                          not_PYREHE_mod, not_SPERPU_mod)

generalsJackKnife <- list(not_ACHADI_mod, not_CROCNO_mod, not_ENNOSU_mod,
                          not_EUTRCL_mod, not_HIMEFI_mod, not_IRIDEP_mod,
                          not_LAMBFI_mod, not_MELACA_mod, not_MORRCO_mod,
                          not_NEMARE_mod, not_ORGYLE_mod, not_ORTHXX_mod)

getEffects <- function(modStack, diet = NULL){
  bestMod <- row.names((modStack$aic[1, ]))
  summ <- data.frame(summary(modStack[["middle"]])$coefficients$cond)
  mood <- summ[row.names(summ) %in% nems, 1:2]
  colnames(mood) <- c("est", "stdE")
  mood$name <- row.names(mood)
  mood$dist <- distances[1:nrow(mood)]
  if(!is.null(diet)) mood$diet <- diet
  return(mood)
}

nems <- c("actual", "resAdjacent", "resInner", "resMiddle", 
          "resOuter", "resExtreme")

specialEffects <- lapply(specialsJackKnife, getEffects)
generalEffects <- lapply(generalsJackKnife, getEffects)

names(specialEffects) <- c("CLADLI", "EPIMHO", "MORREV", "NOLACL", "NOLATR", "PSEUCO", "PYREHE", "SPERPU")
names(generalEffects) <- c("ACHADI", "CROCNO", "ENNOSU",
                           "EUTRCL", "HIMEFI", "IRIDEP",
                           "LAMBFI", "MELACA", "MORRCO",
                           "NEMARE", "ORGYLE", "ORTHXX")

plotSpecials <- purrr::map_dfr(specialEffects, .f = mutate, .id = "Excluded")
plotGenerals <- purrr::map_dfr(generalEffects, .f = mutate, .id = "Excluded")

png("figures/jackKnife.png", height = 1200, width = 750) 

p1 <- ggplot(plotSpecials, aes(x = dist, y = exp(est))) +
  geom_point(aes(color = Excluded), size = 3) +
  stat_smooth(method = "lm", formula = y ~ poly(x, 2), se = F,
              aes(color = Excluded)) +
  scale_x_continuous(trans = "log", breaks = distances, labels = 2:7) +
  geom_errorbar(aes(ymin = exp(est - stdE), ymax = exp(est + stdE)), width = .1) +
  geom_hline(yintercept = 1, linetype = "dotted", color = "red") +
  theme_tufte() +
  theme(text = element_text(size=18),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.line.y = element_line(linewidth = .2, color = "black"),
        axis.line.x = element_blank(),
        axis.title.y = element_text(size=18,margin = unit(c(0, 5, 0, 0), "mm"))) +
  ylab(expression(effect~size)) +
  facet_wrap(~ Excluded)

p2 <- ggplot(plotGenerals, aes(x = dist, y = exp(est))) +
  geom_point(aes(color = Excluded), size = 3) +
  stat_smooth(method = "lm", formula = y ~ poly(x, 2), se = F,
              aes(color = Excluded)) +
  scale_x_continuous(trans = "log", breaks = distances, labels = 2:7) +
  geom_errorbar(aes(ymin = exp(est - stdE), ymax = exp(est + stdE)), width = .1) +
  geom_hline(yintercept = 1, linetype = "dotted", color = "red") +
  theme_tufte() +
  theme(text = element_text(size=18),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.line.y = element_line(linewidth = .2, color = "black"),
        axis.line.x = element_blank(),
        axis.title.y = element_text(size=18,margin = unit(c(0, 5, 0, 0), "mm"))) +
  ylab(expression(effect~size)) +
  facet_wrap(~ Excluded)

library(patchwork)
p1/p2

dev.off()








