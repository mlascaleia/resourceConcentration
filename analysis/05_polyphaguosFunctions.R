# this script is where the bulk of the analysis functions are written
# it also loads the majority of the packages used
# make the functions for the polyphagous analysis
# initialize ####

# load electivity data, which includes cat data
load("data/electivity.RData")

# only valid plants allowed in the vegetation data
vegetation <- vegetation[vegetation$hostGroup %in% colnames(electMatrix),]

# load packages
library(lme4)
library(bbmle)
library(glmmTMB)
library(tidyr)
library(ggthemes)
library(ggallin)
library(future)
library(pscl)
library(MASS)
library(boot)
library(parallel)
library(patchwork)

# using the chesson matrices in modeling ####

hostPlantDensity <- function(csiq, vegl = vegetation, pref = electMatrix, scaled = T) {
  
  # get this caterpillar's preferences from the matrix
  thisCatPrefs <- pref[csiq,] 
  
  # only use vegetation that has a preference value
  # vegl <- vegl[vegl$hostGroup %in% colnames(pref), ] 
  
  # adjust the raw count values to reflect the caterpillars preferences
  # each count is the number of trees * the caterpillars electivity for that species
  nc <- apply(vegl, 1, function(v){
    host <- as.character(v["hostGroup"])
    count <- as.numeric(v["count"])
    multip <- as.numeric(thisCatPrefs[names(thisCatPrefs) %in% host])
    newCounts <- multip * count
    return(newCounts)
  })
  vegl$count <- unlist(nc) 
  
  # sum up the preferred count values at every SiteID (ignoring treeID, unlike the previous function)
  s <- lapply(unique(vegl$PointID), function(p) {
    some <- vegl$count[vegl$PointID %in% p]
    if (length(some) == 0)
      some <- 0
    r <- data.frame(PointID = p, adjacent = sum(some))
    return(r)
  })
  s <- do.call('rbind', s)

  # clean up the dataframe with Plot info
  v <- vegl %>%
    dplyr::select(BlockID, SiteID, PointID, ringID, inPlot) 
  
  v <- unique(v) # only one Plot info row per SiteID
  countData <- merge(s, v, by = "PointID") # merge the count data with the Plot info
  
  # pull out actual level because they're done differently
  countActuals <- countData[countData$ringID %in% "actual",]
  countData <- countData[!countData$ringID %in% "actual",]
  
  rings <- countData %>%
    group_by(SiteID, ringID) %>%
    dplyr::summarise(pc = sum(adjacent)/n(), .groups = 'drop') %>% # this is adjusted for how many Plots there are with the "/n()" operator
    spread(ringID, pc) %>%
    dplyr::select(-inner) # get the total sums for each ring, except the inner ring (which has exceptions)
  
  rings$extreme[is.na(rings$extreme)] <- 0 # saying that anything we didn't measure couldn't really possibly have had the psiq
  
  innervs23 <- countData %>%
    group_by(SiteID, inPlot) %>%
    summarise(inner = sum(adjacent), .groups = 'drop') %>%
    filter(inPlot %in% 2:3) %>%
    group_by(SiteID) %>%
    summarise(inner = sum(inner), .groups = 'drop')
  
  innervs34 <- countData %>%
    group_by(SiteID, inPlot) %>%
    summarise(inner = sum(adjacent), .groups = 'drop') %>%
    filter(inPlot %in% 3:4) %>%
    group_by(SiteID) %>%
    summarise(inner = sum(inner), .groups = 'drop') 
  
  innervs24 <- countData %>%
    group_by(SiteID, inPlot) %>%
    summarise(inner = sum(adjacent), .groups = 'drop') %>%
    filter(inPlot %in% c(2,4)) %>%
    group_by(SiteID) %>%
    summarise(inner = sum(inner), .groups = 'drop')
  
  p <- countData[countData$inPlot %in% 2,] %>% left_join(innervs34, by = "SiteID")
  q <- countData[countData$inPlot %in% 3,] %>% left_join(innervs24, by = "SiteID")
  s <- countData[countData$inPlot %in% 4,] %>% left_join(innervs23, by = "SiteID")
  
  psiqInners <- rbind(p, q, s)
  psiqInners$inner <- psiqInners$inner/2
  psiqAll <- left_join(psiqInners, rings, by = "SiteID") %>%
    dplyr::select(PointID, SiteID, adjacent, inner, middle, outer, extreme)

  # now add in actuals
  
  unScaled <- countActuals %>%
    mutate(cPlot = str_extract(PointID, "CP.$"),
           PointID = gsub("_CP.$","",PointID)) %>%
    rename(actual = adjacent) %>%
    left_join(psiqAll) %>%
    mutate(PointID = paste(PointID, cPlot, sep = "_")) %>%
    dplyr::select(PointID, SiteID, actual, adjacent, inner, middle, outer, extreme)
  
  # scale, make residuals, scale
  
  if(scaled){
    final <- unScaled %>%
      mutate(actual = scale(actual),
             adjacent = scale(adjacent), 
             inner = scale(inner),
             middle = scale(middle),
             outer = scale(outer),
             extreme = scale(extreme)) %>%
      mutate(resAdjacent = scale(resid(lm(adjacent ~ actual))), 
             resInner = scale(resid(lm(inner ~ adjacent))),
             resMiddle = scale(resid(lm(middle ~ inner))),
             resOuter = scale(resid(lm(outer ~ middle))),
             resExtreme = scale(resid(lm(extreme ~ outer))))
  } else {
    final <- unScaled %>%
      mutate(actual = actual,
             adjacent = adjacent, 
             inner = inner,
             middle = middle,
             outer = outer,
             extreme = extreme) %>%
      mutate(resAdjacent = scale(resid(lm(adjacent ~ actual))), 
             resInner = scale(resid(lm(inner ~ adjacent))),
             resMiddle = scale(resid(lm(middle ~ inner))),
             resOuter = scale(resid(lm(outer ~ middle))),
             resExtreme = scale(resid(lm(extreme ~ outer))))
  }
  
  return(final)
}


catDensity <- function(csiq, catl = caterpillars, pref = electMatrix, onlyPossibilities = T){
  
  # get the names of all trees this species eats at all
  pCode <- names(pref[csiq,][pref[csiq, ] > 0]) 
  
  # pull out all instances of any of these trees occurring
  catPlant <- catl[catl$hostGroup %in% pCode, ] 
  
  # find total number of branches sampled in each Plot for the host species in a given Year
  cc <- catPlant %>%                           
    group_by(PointID, Year, hostGroup, SiteID, SurveyDate) %>% 
    summarise(bSample = n_distinct(branchID), .groups = 'drop') 
    
  # get this caterpillar's preferences from the matrix
  thisCatPrefs <- pref[csiq,]
  
  # this apply loop multiplies the branches sampled by the caterpillar's preference
  nc <- apply(cc, 1, function(v){
    
    # get host group and branches sampled
    host <- as.character(v["hostGroup"])
    bcount <- as.numeric(v["bSample"])
    
    # get the preference of this cat species for this hostgroup
    multip <- as.numeric(thisCatPrefs[names(thisCatPrefs) %in% host])
    
    # multiply the branches sampled by the preference
    newBCounts <- bcount * multip
    return(newBCounts)
  })
  
  # replace branch count with our new pref values for this caterpillar
  # don't do this actually it's artefact city
  # cc$bSample <- nc 
  
  # aggregate all hostGroups together to get only one value for branches sampled per PointID
  if(onlyPossibilities){
    cct <- cc %>%
      group_by(PointID, Year) %>%
      mutate(bSample = sum(bSample)) %>%
      ungroup() %>%
      dplyr::select(-hostGroup) %>%
      distinct()
  } else {
    catl$pointYear <- paste(catl$PointID, catl$Year, sep = ".")
    cct <- cc %>%
      group_by(PointID, Year) %>%
      mutate(bSample = sum(bSample),
             pointYear = paste(PointID, Year, sep = ".")) %>%
      ungroup() %>%
      dplyr::select(-hostGroup) %>%
      distinct()
    zeroPointYears <- unique(catl$pointYear[!catl$pointYear %in% cct$pointYear])
    zeroPYDataFrame <- data.frame(PointID = gsub("\\..*$","", zeroPointYears),
                                  Year = as.integer(gsub("^.*\\.", "", zeroPointYears)),
                                  SiteID = str_extract(zeroPointYears, "(?<=_)[[:alpha:]]{3}(?=_)"),
                                  SurveyDate = 0,
                                  bSample = 0,
                                  pointYear = zeroPointYears)
    cct <- rbind(cct, zeroPYDataFrame) %>%
      arrange(Year, PointID)
  }
  
  # pull out any instance of this caterpillar being sampled
  csiqPlant <- catl[catl$CatID %in% csiq, ] 
  
  # aggregate by PointID and year
  cqp <- csiqPlant %>% 
    dplyr::select(PointID, Year, SiteID, SurveyDate) %>%
    group_by(PointID, Year) %>%
    mutate(catCount = n()) %>%
    ungroup() %>%
    distinct()

  # combine the amount of vegetation sampled, and the amount of caterpillars found
  cagg <- cqp %>%
    full_join(cct, by = c("PointID", "Year", "SiteID", "SurveyDate")) %>%
    arrange(Year, PointID) %>%
    replace_na(list(catCount = 0))

  return(cagg)
}

# pladf <- hostPlantDensity("NEMARE")
# cadf <- catDensity("NEMARE")

# merge plant and caterpillar data
polyphagousMerge <- function(pladf, cadf){
  suppressWarnings(oneCat <- pladf %>%
    left_join(cadf, by = c("PointID", "SiteID")) %>%
    mutate(adjPoint = gsub("_CP.$", "", PointID)) %>%
    filter(!is.na(bSample)))
  return(oneCat)
}

# finally, all the functions get put into one big organize function

organizePolyphagous <- function(caterpillarInQuestion, 
                                vegetationDF = vegetation, 
                                caterpillarDF = caterpillars,
                                preferenceMatrix = electMatrix,
                                preferenceMatrix.cat = electMatrix,
                                toScale = T,
                                possiblesOnly = T){
  # establish plant df
  plantDF <- hostPlantDensity(csiq = caterpillarInQuestion, 
                              vegl = vegetationDF, 
                              pref = preferenceMatrix,
                              scale = toScale)
  
  # establish caterpillar df
  catDF <- catDensity(csiq = caterpillarInQuestion, 
                      catl = caterpillarDF, 
                      pref = preferenceMatrix.cat,
                      onlyPossibilities = possiblesOnly)
  
  # combine
  pOrganized <- polyphagousMerge(pladf = plantDF, cadf = catDF)
  return(pOrganized)
}

# make the models which the data will be fit to
polyphagousAnalysis <- function(organizedDF, form = modForm){
  extreme  <- glmmTMB(formula = form, data = organizedDF, family = nbinom2,
                      control = glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS"), parallel = 40)) 
  print("completed 'extreme' model")
  outer <- update(extreme, ~ .-resExtreme -extreme)
  print("completed 'outer' model")
  middle <- update(outer, ~ .-resOuter -outer)
  print("completed 'middle' model")
  inner <- update(middle, ~ .-resMiddle -middle)
  print("completed 'inner' model")
  adjacent <- update(inner, ~ .-resInner -inner)
  print("completed 'adjacent' model")
  actual <- update(adjacent, ~ .-resAdjacent -adjacent)
  print("completed actual model")
  null <- update(actual, ~ .-actual)
  print("completed null model")
  ww <-
    list(
      actual = actual,
      adjacent = adjacent,
      inner = inner,
      middle = middle,
      outer = outer,
      extreme = extreme,
      null = null,
      aic = as.data.frame(AICtab(
        actual, adjacent, inner, middle, outer, extreme, null
      )),
      data = organizedDF
    )
  return(ww)
}

# bring it together ####
# This last function does all the work for me

fullPolyphagousAnalysis <- function(species, 
                                    vegSurvey = vegetation, 
                                    catSurvey = caterpillars,
                                    catPreferences = electMatrix, 
                                    useResiduals = T,
                                    randAnalysis = F,
                                    extremeOnly = F,
                                    supplyData = NULL){
  # make sure data has been entered correctly...
  if(length(species) == 1) 
    multipleSpecies <- F
  else
    multipleSpecies <- T
  if(!multipleSpecies){
    pOrganizedDF <- organizePolyphagous(caterpillarInQuestion = species, 
                                      vegetationDF = vegSurvey, 
                                      caterpillarDF = catSurvey,
                                      preferenceMatrix = catPreferences)
    pOrganizedDF$cat <- species
  }
  else{
    aa <- apply(as.matrix(species), 1, organizePolyphagous, vegetationDF = vegSurvey, 
                                              caterpillarDF = catSurvey,
                                              preferenceMatrix = catPreferences)
    aa <- Map(cbind.data.frame, aa, cat = species)
    pOrganizedDF <- do.call("rbind", aa)
  }
  if(randAnalysis){
    warning("data has been randomized")
    pOrganizedDF$catCount <- sample(pOrganizedDF$catCount)
    pOrganizedDF$bSample <- sample(pOrganizedDF$bSample)
  }
  if(!is.null(supplyData)) pOrganizedDF <- supplyData
  
  # make sure year is a factor
  pOrganizedDF$Year <- factor(pOrganizedDF$Year)
  contrasts(pOrganizedDF$Year) <- "contr.sum"
  print("data organized, starting models")
  
  # establish global model formula
  modForm <- catCount ~ offset(log(bSample)) + actual + resAdjacent + resInner + resMiddle + resOuter + resExtreme + Year + SurveyDate + (1|cat) + (1|SiteID/adjPoint/cat)

  # if the data has been entered correctly, do one of the following modeling methods:
  if(extremeOnly){
    extreme  <- glmmTMB(formula = modForm, data = pOrganizedDF, family = nbinom2,
                        control = glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS"),
                                                 parallel = 40)) 
    return(extreme)
  }
  else if(multipleSpecies & useResiduals)
    polyphagousAnalysis(pOrganizedDF, form = modForm)
  else if(multipleSpecies & !useResiduals)
    polyphagousAnalysis(organizedDF = pOrganizedDF, form = update(modForm, ~. -resAdjacent -resInner -resMiddle -resOuter -resExtreme + adjacent + inner + middle + outer + extreme))
  else if(!multipleSpecies & useResiduals)
    polyphagousAnalysis(organizedDF = pOrganizedDF, form = update(modForm, ~. -(1|cat) -(1|SiteID/adjPoint/cat) + (1|SiteID/adjPoint)))
  else 
    polyphagousAnalysis(organizedDF = pOrganizedDF, form = update(modForm, ~. -resAdjacent -resInner -resMiddle -resOuter -resExtreme + adjacent + inner + middle + outer + extreme - (1|cat) -(1|SiteID/PointID/cat) + (1|SiteID/PointID)))
    
}

# I honestly forget why this is needed but I'm afraid to delete it
asinh_trans <- function(){
  scales::trans_new(name = 'asinh', transform = function(x) asinh(x),
                    inverse = function(x) sinh(x))
}

# randomize data for articaft checking
sMods <- function(e){
  randHost <- t(apply(e, 1, sample))
  colnames(randHost) <- colnames(e)
  rHost <- fullPolyphagousAnalysis(catPreferences = randHost, randAnalysis = T)
  bb <- summary(rHost$extreme)$coefficients$cond[2:6,1:2]
  return(bb)
}


