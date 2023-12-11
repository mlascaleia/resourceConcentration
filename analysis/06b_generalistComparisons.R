# initialize ####

rm(list = ls())
source("analysis/05_polyphaguosFunctions.R")
load("data/currentModels.rdata")

# load data ####

generals <- c("MELACA", "NEMARE", "CROCNO", "IRIDEP", 
              "LAMBFI", "MORRCO", "HIMEFI", "ACHADI",
              "ORTHXX", "EUTRCL", "ENNOSU", "ORGYLE")

# load("data/currentModels.rdata")

# run binary model ####

allSame <- ifelse(electMatrix > 0, 1, 0)

# genAllHostGood <- fullPolyphagousAnalysis(generals, catPreferences = allSame)

# make no preferences

monolith <- ifelse(is.na(electMatrix), 1, 1)

# make by hand (it breaks my catch-all functions)
aa <- apply(as.matrix(generals), 1, organizePolyphagous, vegetationDF = vegetation, 
            caterpillarDF = caterpillars,
            preferenceMatrix = monolith, preferenceMatrix.cat = allSame)
aa <- Map(cbind.data.frame, aa, cat = generals)
monolithData <- do.call("rbind", aa)

# run no preferences model

# genMonolith <- fullPolyphagousAnalysis(generals, catPreferences = monolith, supplyData = monolithData)

# # elimmed monolith pairs
# paste(colnames(allSame)[as.data.frame(which(allSame == 0, arr.ind = T))$col],
# row.names(allSame)[as.data.frame(which(allSame == 0, arr.ind = T))$row], sep = "_")

# save(genAllHostGood, genMonolith, file = "data/models/generalistComparisons.rdata")
load("data/models/generalistComparisons.rdata")

# jackknife to compare quality

withElectivityData <- generalsMods$data
hostsEqualData <- genAllHostGood$data
monolithData <- genMonolith$data

withElectivityData$dataType <- "electivity"
hostsEqualData$dataType <- "hostEquality"
monolithData$dataType <- "monolith"

allGenModDatas <- list(withElectivityData = withElectivityData,
                       hostsEqualData = hostsEqualData,
                       monolithData = monolithData)

allGenModDatas <- lapply(allGenModDatas, function(x){
  x$BlockID <- str_extract(x$PointID, "^..")
  return(x)
})

modelFormula <- catCount ~ offset(log(bSample)) + 
                           actual + 
                           resAdjacent + 
                           resInner + 
                           resMiddle + 
                           resOuter + 
                           resExtreme + 
                           Year + SurveyDate + 
                           (1|cat) + (1|SiteID/adjPoint/cat)

blocks <- unique(allGenModDatas$withElectivityData$BlockID)
# blocks <- 1:13

# for monolith only loop
monolithDataList <- list(allGenModDatas$monolithData)

# this takes a very, very long time to run

# eekMonolith <- lapply(blocks, function(b){
#   lapply(monolithDataList, function(d){
#     missingOneBlock <- d[!d$BlockID %in% b, ]
#     theOneBlockData <- d[d$BlockID %in% b, ]
#     oneJackKnife <- polyphagousAnalysis(organizedDF = missingOneBlock,
#                                         form = modelFormula)
#     bestModelName <- row.names(oneJackKnife$aic)[1]
#     valueOverNull <- oneJackKnife$aic["null", "dAIC"]
#     bestModelActual <- oneJackKnife[[bestModelName]]
#     predictedValues <- exp(predict(bestModelActual, newdata = theOneBlockData))
#     linearModel <- lm(theOneBlockData$catCount ~ predictedValues)
#     summaryDataFrame <- data.frame(dataType = d$dataType[1],
#                                    leftOutBlock = b,
#                                    bestModel = bestModelName,
#                                    linear_rSquared = summary(linearModel)$r.squared)
#     observedVsPredicted <- cbind.data.frame(theOneBlockData, prediction = predictedValues)
#     output <- list(summaryDataFrame = summaryDataFrame,
#                    observedVsPredicted = observedVsPredicted,
#                    jackKnifeModels = oneJackKnife,
#                    bestModel = bestModelActual,
#                    linearModel = linearModel)
#     return(output)
#   })
# })

# There's still one small issue, where monolith looks at more data, even in the predict
# fixed. New monolith input here

# for(i in 1:13){
#   eek[[i]][[3]] <- eekMonolith[[i]][[1]]
# }


# save(eek, file = "data/models/bigJackKnife.rdata")

# hey uhhhhh yikes at that file size
load("data/models/bigJackKnife.rdata")

pTypes <- c("Electivity", "Binary", "No specificity")

allDAIC <- lapply(eek, function(e){
  pp <- lapply(1:3, function(i){
    p <- e[[i]]$jackKnifeModels$aic["null", "dAIC"]
  })
  return(pp)
})

psuedoR_squared <- lapply(eek, function(e){
  # get obs vs pred from monolith's data, as monolith never beat the null
  # (this no longer works if monolith beats the null EVER)
  nullData <- e$monolithData$observedVsPredicted
  nullMod <- e$monolithData$bestModel
  nullDispersion <- summary(nullMod)$sigma
  nullPredicted <- nullData$prediction
  nullObserved <- nullData$catCount
  nullDeviance <- 2 * sum(dnbinom(x = nullObserved,
                                  size = nullDispersion,
                                  mu = nullPredicted,
                                  log = T))
  pseudoR_squared <- lapply(e, function(f, nullDeviance){
    data <- f$observedVsPredicted
    mod <- f$bestModel
    dispersion <- summary(mod)$sigma
    predicted <- data$prediction
    observed <- data$catCount
    deviance <- 2 * sum(dnbinom(x = observed,
                                    size = dispersion,
                                    mu = predicted,
                                    log = T))
    pseudoR_squared <- 1 - deviance/nullDeviance
    return(pseudoR_squared)
  }, nullDeviance = nullDeviance)
  return(pseudoR_squared)
})


dAICs <- data.frame(dAIC = unlist(allDAIC),
                    predictorType = rep(pTypes, 13),
                    predictedBlock = rep(blocks, each = 3),
                    linear_rSquared = unlist(psuedoR_squared))

dAICs$predictorType <- factor(dAICs$predictorType, levels = c("No specificity", 
                                                                      "Binary",
                                                                      "Electivity"))

summaryDF <- dAICs %>%
  group_by(predictorType) %>%
  summarise(mean_dAIC = mean(dAIC),
            stdErr_dAIC = sd(dAIC)/sqrt(n() - 1),
            maxAIC = max(dAIC),
            minAIC = min(dAIC),
            mean_rSquared = mean(linear_rSquared),
            stdErr_rSquared = sd(linear_rSquared)/sqrt(n() - 1)) %>%
  mutate(predictorType = factor(predictorType, levels = pTypes))

scaler <- summaryDF$mean_dAIC[1]/max(dAICs$linear_rSquared) + 100

set.seed(123)
dAICs$linear_rSquared[dAICs$predictorType %in% "No specificity"] <- rnorm(length(
  dAICs$linear_rSquared[dAICs$predictorType %in% "No specificity"]),
  0, 0.001)
r2.mod <- lmerTest::lmer(linear_rSquared ~ predictorType + (1|predictedBlock), data = dAICs)
summary(r2.mod)

# figure 4 ####

pp1 <- ggplot(data = summaryDF, aes(x = predictorType)) +
  geom_point(data = dAICs, aes(y = linear_rSquared), 
             position = position_jitter(width = .05),
             alpha = .35, color = "blue") +
  geom_pointrange(aes(y = mean_rSquared, 
                      ymax = mean_rSquared + 1.96 * stdErr_rSquared, 
                      ymin = mean_rSquared - 1.96 * stdErr_rSquared),
                  # position = position_nudge(x = -.2),
                  color = "darkblue", shape = 15) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  scale_y_continuous(name = "\nPseudo R-squared of prediction") +
  scale_x_discrete(name = "Host Specificity Data Used\n") +
  theme_tufte() +
  theme(axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 14)) +
  coord_flip()

labby <- ggplot() +
  ylab("Host Specificity Data Used") +
  theme_tufte() +
  theme(text = element_text(size = 20))

labby + pp1 + plot_layout(widths = c(1,40))

ggsave("figures/generalistTest.png", width = 8.18, height = 5.9)





