# initialize ####

rm(list = ls())

# get everything needed
source("analysis/05_polyphaguosFunctions.R")

# establish specialists and generalist groups ####

specials <- c("EPIMHO", "NOLACL", "NOLATR", "PSEUCO", 
              "PYREHE", "SPERPU", "MORREV", "CLADLI")

generals <- c("MELACA", "NEMARE", "CROCNO", "IRIDEP", 
              "LAMBFI", "MORRCO", "HIMEFI", "ACHADI",
              "ORTHXX", "EUTRCL", "ENNOSU", "ORGYLE")
all <- c(specials, generals)

# run models ####

# specialsMods <- fullPolyphagousAnalysis(specials)
# summary(specialsMods$extreme)
# 
# generalsMods <- fullPolyphagousAnalysis(generals)
# summary(generalsMods$extreme)

# save models ####

# save(specialsMods, generalsMods, file = "data/currentModels.rdata")
load("data/currentModels.rdata")

# run other models ####

# each caterpillar individually ####
# this doesn't work very well, and probably shouldn't be done

# for(i in 1:length(all)){
#   assign(x = paste(all[i], "mod", sep = "_"), fullPolyphagousAnalysis(all[i]))
# }

# and another try with no generalist singletons ####
# gdat <- caterpillars[caterpillars$CatID %in% ggs, ]
# gtab <- as.matrix(table(gdat$CatID, gdat$hostGroup))
# sing <- apply(gtab, 1, function(x){
#   names(which(x == 1))
# })
# 
# em2 <- electMatrix
# for(i in names(sing)){
#   ts <- sing[[i]]
#   em2[i, ts] <- 0
# }

# gNoSingle <- fullPolyphagousAnalysis(ggs, catPreferences = em2, extremeOnly = T)
# save(gNoSingle, file = "noSingletonsModelGeneralists.rdata")

# summary(gNoSingle)

# lymantria model ####

# lymadiMods <- fullPolyphagousAnalysis("LYMADI")
