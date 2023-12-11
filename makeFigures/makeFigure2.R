# make network figure

rm(list = ls())
library(tidyverse)
library(network)
library(igraph)
library(ggraph)
library(bipartite)
library(ggthemes)
library(patchwork)
library(ggpubr)
library(viridis)
load("data/electivity.RData")
load("data/catAndVeg.RData")

specials <- rownames(electMatrix)[apply(electMatrix, 1, max) > .9]
specials <- c(specials, "MORREV", "CLADLI")
ggs <- c("ORGYLE", "MELACA", "NEMARE", "CROCNO", "IRIDEP", "LAMBFI", "MORRCO", "HIMEFI", "ACHADI", "ORTHXX", "EUTRCL", "ENNOSU")
all <- c(specials, ggs)
electMatrix <- electMatrix[row.names(electMatrix) %in% all, ]


tc <- caterpillars %>%
  group_by(CatID, hostGroup) %>%
  summarise(total = n()) %>%
  filter(CatID %in% row.names(electMatrix) &
         hostGroup %in% colnames(electMatrix)) %>%
  pivot_wider(names_from = hostGroup, 
              values_from = total,
              values_fill = 0) %>%
  remove_rownames %>% 
  column_to_rownames(var="CatID") %>%
  as.matrix()

nem <- t(tc)

# nem <- sapply(1:nrow(tc), FUN = function(i) tc$total[i] * electMatrix[i, ])
# colnames(nem) <- row.names(electMatrix)
nem <- nem[, all ]
nem <- nem[c("LINBE", "CLEAL", "HAMVI", "ACERU", "KALLA", "TSUCA",
             "ACEBL", "AMEXX", "BETBL", "CARCA", "CARYA", "CASDE", 
             "CORYL", "FAGGR", "GAYLU", "ILEVE", "LYOLI", "OSTVI",
             "PINST", "QUEAL", "QUEBL", "QUEMO", "VACCO"), ]

plaCodes <- row.names(nem)

row.names(nem) <- c("L. benzoin", "C. alnifolia", "H. virginiana", "A. rubrum",
                    "K. latifolia", "T. canadensis", "Acer", "Amelanchier",
                    "Betula", "C. caroliniana", "Carya", "C. dentata",
                    "Corylus", "F. grandifolia", "Gaylussacia", "I. laevis",
                    "L. ligustrina", "O. virginiana", "P. strobus", "Q. alba", "Quercus",
                    "Q. montana", "Vaccinium")

plaNames <- row.names(nem)
tree_trans <- cbind.data.frame(hostGroup = plaCodes, name = plaNames)

nem <- t(nem)
codes <- row.names(nem)
row.names(nem) <- c("E. hortaria", "N. clethrae", "N. triquetrana", "P. costomaculana",
                    "P. hesperidago", "S. pustularia", "M. evicta", "C. limitaria",
                    "O. leucostigma", "M. canadaria", "N. resistaria", "C. normani",
                    "I. ephyraria", "L. fiscellaria", "M. confusa", "H. fidelis", "A. distincta",
                    "Orthosia", "E. clemataria", "E. subsignaria")
full_cat <- row.names(nem)

cat_trans <- cbind.data.frame(CatID = codes, name = full_cat)

png(height = 1200, width = 1500, file = "figures/network.png")

par(font = 3, family = "serif")

plotweb(nem, text.high.col = "grey5", text.low.col = "grey5", 
        bor.col.high = "grey5", bor.col.low = "grey5", 
        col.interaction = "grey80",
        method = "normal",
        col.low = ifelse(all %in% specials, "#6FBCEB","#3A577E"),
        col.high = "forestgreen",
        text.rot = 90,
        labsize = 4,
        low.y = .65, high.y = 1.25)

dev.off()

# plot b:
# caterpillar abundance rank curve w/mean electivity fill

# the most complicated way ever to do this... with math!!
# get number of hosts
nHost <- apply(electMatrix, 1, function(x){
  x <- x[x != 0]
  return(1/mean(x))
})

nHost

nHostFrame <- data.frame(CatID = names(nHost), nHosts = nHost)

catFig2 <- caterpillars %>%
  group_by(CatID) %>%
  filter(CatID %in% names(nHost)) %>%
  summarise(abundance = n()) %>%
  full_join(nHostFrame) %>%
  full_join(cat_trans) %>%
  arrange(-abundance)

p2 <- ggplot(data = catFig2, aes(x = 1:nrow(catFig2),
                           y = abundance)) +
  geom_line() +
  geom_point(aes(color = nHosts), size = 5) +
  theme_tufte() +
  scale_y_continuous(trans = "log",
                     breaks = c(30, 50, 100, 200, 300)) +
  scale_x_continuous(breaks = 1:nrow(catFig2),
                     labels = catFig2$name) +
  theme(text = element_text(size = 20),
        axis.line = element_line(color = "black", size = .1),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
        panel.grid.major = element_line(color = "grey90")
        ) +
  xlab("") +
  ylab("Total species abundance") +
  scale_color_viridis(option = "mako", name = "Number of\nhosts", direction = -1, end = .89) +
  theme(axis.text.x = element_text(face = "italic"))

e_total <- data.frame(sumE = rowSums(t(electMatrix))) %>%
  rownames_to_column("hostGroup") %>%
  filter(hostGroup %in% tree_trans$hostGroup)

b_total <- branches %>%
  group_by(hostGroup) %>%
  summarise(bTotal = sum(bSample)) %>%
  filter(hostGroup %in% tree_trans$hostGroup)

c_total <- caterpillars %>%
  filter(CatID %in% cat_trans$CatID) %>%
  group_by(hostGroup) %>%
  summarise(cTotal = n()) %>%
  filter(hostGroup %in% tree_trans$hostGroup) %>%
  full_join(b_total) %>%
  mutate(catDensity = cTotal/bTotal) %>%
  full_join(e_total) %>%
  full_join(tree_trans) %>%
  arrange(-cTotal)
  
p3 <- ggplot(data = c_total, aes(x = 1:nrow(c_total),
                           y = cTotal)) +
  geom_line() +
  geom_point(aes(col = log(sumE), size = catDensity)) +
  theme_tufte() +
  scale_y_continuous(trans = "log",
                     breaks = c(10, 25, 100, 300)) +
  scale_x_continuous(breaks = 1:nrow(c_total),
                     labels = c_total$name) +
  theme(text = element_text(size = 20),
        axis.line = element_line(color = "black", size = .1),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 16, face = "italic"),
        panel.grid.major = element_line(color = "grey90")
  ) +
  xlab("") +
  ylab("Total caterpillar abundance") +
  scale_size_continuous(name = "Caterpillars\nper branch") +
  scale_color_viridis(name = expression(paste("ln(",Sigma,hat(alpha)[i],")")))
p3

png("figures/fig2b2c.png", height = 1000, width = 800)
(p2/p3)
dev.off()


