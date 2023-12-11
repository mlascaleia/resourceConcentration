# This is the first file that should be run in the analysis folder #
# All scripts in this file should be run within the .RProj "CaterpillarHostDD" #
# some numbers will be missing - those have been deleted as they were deprecated
# and are ot used in the final analysis

# Read in caterpillar data and initialize packages ####

rm(list = ls())
library(tidyverse)
library(lubridate)

# Clean caterpillar data ####

# Read in caterpillar data from all years
cats <- read.csv("data/raw/CaterpillarSurveysAllYears-v1_2.csv")

# scrap data from 2015 because it was not collected with the same methodology of the other years
cats <- cats[!cats$Year %in% 2015,]

# There are many columns in the dataset that I don't need
# I find it easier to use when there are fewer unnecessary columns

cats$Page <- NULL
cats$LineNum <- NULL
cats$Enterer <- NULL
# cats$BlockID <- NULL
# cats$SiteID <- NULL
# cats$Plot <- NULL
cats$BranchDiam <- NULL
cats$NumLeaves <- NULL
cats$CatIDNum <- NULL
cats$RearID <- NULL
cats$EntererNotes <- NULL
cats$OpenRefineNotes <- NULL
cats$ID <- NULL

# then there are a few that I thought I might not need (still think so)
# but I have found a use for them

cats$Recorder <- NULL
cats$Team <- NULL
cats$PlantAbund <- NULL
cats$CatSpecies <- NULL
# cats$IndNum <- NULL
# cats$BranchNum <- NULL
# cats$BranchLength <- NULL
# cats$Year <- NULL
cats$BadHost <- NULL
cats$FieldNotes <- NULL
cats$X <- NULL
cats$CleaningNotes <- NULL
cats$HostAbund <- NULL
cats$BranchLength <- NULL

#saving this because it's useful

rawCats <- cats

# Make one record per caterpillar

# Create a record ID for each "bunch" of caterpillars found (same species found on the same branch) 
cats$bunch <- c(1001:(1000 + nrow(cats))) 
cats$Count <- as.numeric(cats$Count)

# get rid of some entirely blank rows, not sure why they exist
cats <- cats[!is.na(cats$Count), ]

# This line makes one row per caterpillar. It is a coincidence that it is very similar to the number of rows in "cats"
e <- cats[rep(row.names(cats), cats$Count), ]

# Add back in all branches where 0 caterpillars were found
z <- cats[cats$Count == 0, ]
longCats <- rbind(e,z)
rm(e,z)

# now sort them based on how they were originally listed
longCats <- longCats %>%
  arrange(bunch)

rm(cats)

# add in some new columns, and rename others

longCats$lineID <- 20001:(20000 + nrow(longCats))
colnames(longCats)[colnames(longCats) %in% "PlotID"] <- "PointID"
colnames(longCats)[colnames(longCats) %in% "HostID"] <- "PlantID"
longCats$PlantID <- as.character(longCats$PlantID)
longCats$branchID <- paste(longCats$PointID, longCats$PlantID, longCats$IndNum, longCats$BranchNum, longCats$Year, sep = "_")
longCats$eachCatID <- paste(longCats$branchID, longCats$CatID, longCats$lineID, sep = "_")
longCats$cplot <- str_extract(longCats$PointID, "(?<=_)CP[[:digit:]]$")

# remove some columns that I no longer need

longCats$Notes <- NULL
longCats$BranchNum <- NULL
# longCats$SurveyDate <- NULL
longCats$bunch <- NULL
longCats$lineID <- NULL
longCats$Count <- NULL
longCats$IndNum <- NULL


# make a data frame that shows the number of branches sampled per caterpillar plot, per year, per plant type
branchesPPY <- longCats %>%
  group_by(Year, PointID, PlantID) %>%
  summarise(bSample = n_distinct(branchID))

# now make the dates usable...
longCats$SurveyDate <- yday(as.Date(longCats$SurveyDate))  - min(yday(as.Date(longCats$SurveyDate)))

# and I can remove the final few columns I don't need
longCats$eachCatID <- NULL

# Vegetation data organization ####
vege <- read.csv("data/raw/VegetationSurveys.csv")

# I need the density of each tree species at each plot

# First, get rid of everything that's dead
vege <- vege[!vege$Dead %in% "yes",]

# select necessary columns
vege2 <- vege %>%
  dplyr::select(BlockID, SiteID, PointID, Tree) %>%
  group_by(PointID, Tree) %>%
  summarise(count = n()) 

vegetation <- full_join(vege, vege2) %>%
  dplyr::select(BlockID, SiteID, PointID, Tree, count)

rm(vege, vege2)

# Establish vegetation rings ####

vegetation$vegPlot <- str_extract(vegetation$PointID, "(?<=_)[[:digit:]]{1,2}(?=_(1|2)$)")
vegetation$vegPlot[is.na(vegetation$vegPlot)] <- str_extract(vegetation$PointID[is.na(vegetation$vegPlot)], "(2|3|4)$")

vegetation$ringID <- NA
vegetation$inPlot <- NA

vegetation$ringID[vegetation$vegPlot %in% c(2,3,4)] <- "inner"
vegetation$ringID[vegetation$vegPlot %in% c(5,6,7)] <- "middle"
vegetation$ringID[vegetation$vegPlot %in% c(8,9,10)] <- "outer"
vegetation$ringID[vegetation$vegPlot %in% c(11,12,13)] <- "extreme"

vegetation$inPlot[vegetation$vegPlot %in% c(2)] <- 2
vegetation$inPlot[vegetation$vegPlot %in% c(3)] <- 3
vegetation$inPlot[vegetation$vegPlot %in% c(4)] <- 4
vegetation$inPlot[vegetation$vegPlot %in% c(5,6,7)] <- NA
vegetation$inPlot[vegetation$vegPlot %in% c(8,9,10)] <- NA
vegetation$inPlot[vegetation$vegPlot %in% c(11,12,13)] <- NA

# error fixing, no longer used ####
# 
# c <- caterpillars %>%
#   
#   group_by(PointID, plantID, year) %>%
#   
#   filter(year == 2017) %>%
#   
#   summarise(plantCount = max(indNum)) %>%
#   
#   dplyr::select(PointID, plantID, plantCount, year) %>%
#   
#   full_join(caterpillars, by = c("plantID", "year", "PointID")) %>%
#   
#   rename(catCount = count) %>%
#   
#   mutate(plantCount = coalesce(plantCount.x, plantCount.y)) %>%
#   
#   dplyr::select(-plantCount.x, -plantCount.y) %>%
#   
#   ungroup() %>%
#   
#   dplyr::select(recorder, team, PointID, plantID, plantCount, indNum, catID, catCount, bSampled, year, block, site, plot, cplot, line)
# 
# # There's another issue in here where the plant abundance isn't consistent within a site in 2018 or 2019
# 
# c <- c %>% 
#   
#   arrange(PointID, year, plantID, plantCount) %>%
#   
#   mutate(line = 1:(nrow(c)))

# d <- c  %>%
#   
#   group_by(PointID, plantID, year) %>%
#   
#   summarise(tooMany = n_distinct(plantCount), minLine = min(line), maxLine = max(line), highestNum = max(plantCount)) %>%
#   
#   filter(tooMany > 1) %>%
#   
#   arrange(-tooMany) 
# 
# # sweet okay there's only 15 of them
# # I can do that in bulk here by always taking the highest value and applying it
# 
# for(i in 1:nrow(d)){
#   n <- d$minLine[i]
#   x <- d$maxLine[i]
#   c$plantCount[c$line %in% n:x] <- d$highestNum[i]
# } 

# check...
# d <- c  %>%
#   
#   group_by(PointID, plantID, year) %>%
#   
#   summarise(tooMany = n_distinct(plantCount), minLine = min(line), maxLine = max(line), highestNum = max(plantCount)) %>%
#   
#   filter(tooMany > 1)
# d is now empty, which means we're good!  

# make tree groups ####

vegetation$PlantID <- toupper(
  paste0(
    str_extract(pattern = "^...", string = vegetation$Tree),
    str_extract(pattern = "(?<= )..", string = vegetation$Tree)
  )
)

vegetation$Tree <- NULL

# some trees are hard to tell apart, and contain near-identical caterpillars
# this function lumps those together

# First, let's create a new column with just the same names 
treeGroups <- function(c){

  c$hostGroup <- as.character(c$PlantID)
  
  # All Hickories get group "CARYA"
  
  unique(c$PlantID[grepl("CAR..", c$PlantID)]) # this means it should inclue all these except CARCA
  c$hostGroup[c$PlantID %in% c('CARGL', 'CARTO', 'CARXX', 'CAROV')] <- "CARYA"
  
  # All Corylus get group "CORYL"
  
  unique(c$PlantID[grepl("COR..", c$PlantID)]) # there are dogwoods in here so only include CORCO and CORAM
  c$hostGroup[c$PlantID %in% c('CORCO', 'CORAM')] <- "CORYL"
  
  # Red + Black + Scarlet Oak # get code QUEBL (for black oak group)
  
  c$hostGroup[c$PlantID %in% c('QUERU', 'QUEVE', 'QUECO')] <- "QUEBL"
  
  # Sugar + Black Maple get code ACEBL for black maple group
  
  c$hostGroup[c$PlantID %in% c('ACESA', 'ACENI')] <- "ACEBL"
  
  # Black + Yellow Birch get BETBL for black birch group
  
  c$hostGroup[c$PlantID %in% c('BETAL', 'BETLE')] <- "BETBL"
  
  # Both Gaylussacia get GAYLU
  
  unique(c$PlantID[grepl("GAY..", c$PlantID)]) # GAYBA and GAYXX
  
  c$hostGroup[c$PlantID %in% c('GAYBA', 'GAYXX')] <- "GAYLU"
  
  return(c)

}


caterpillars <- treeGroups(c = longCats)
branches <- treeGroups(c = branchesPPY)
vegetation <- treeGroups(c = vegetation)
vegetation <- vegetation %>% distinct()

# eliminate contaminants ####
# see the end of 04 for details

# okay so I didn't catch everything at the end of four so I have to run it all MANUALLY

caterpillars <- caterpillars[!(caterpillars$CatID %in% "NOLATR" &
                              !caterpillars$hostGroup %in% "HAMVI"), ]

caterpillars <- caterpillars[!(caterpillars$CatID %in% "EPIMHO" &
                                 !caterpillars$hostGroup %in% "LINBE"), ]

caterpillars <- caterpillars[!(caterpillars$CatID %in% "PSEUCO" &
                                !caterpillars$hostGroup %in% "HAMVI"), ]

caterpillars <- caterpillars[!(caterpillars$CatID %in% "PYREHE" &
                                 !caterpillars$hostGroup %in% "HAMVI"), ]

caterpillars <- caterpillars[!(caterpillars$CatID %in% "SPERPU" &
                                 !caterpillars$hostGroup %in% "ACERU"), ]

caterpillars <- caterpillars[!(caterpillars$CatID %in% "MORREV" &
                                 !caterpillars$hostGroup %in% c("KALLA", "LYOLI")), ]

caterpillars <- caterpillars[!(caterpillars$CatID %in% "CLADLI" &
                                 !caterpillars$hostGroup %in% c("KALLA", "TSUCA")), ]

caterpillars$CatID[caterpillars$CatID %in% c("ORTHRU", "ORTHHI")] <- "ORTHXX" 

# I'm going to take one final stab at adding in the actual plots for vegetation density
# I will call this "actual"

actualPlots <- rawCats %>%
  
  # eliminate 2015
  filter(!is.na(Year)) %>%
  
  # rename PlantID column
  rename(PlantID = HostID) %>% 
  
  # make host groups
  treeGroups() %>% 
  
  # get rid of unneeded column
  dplyr::select(-PlantID) %>%
  
  # get total number of each tree type
  group_by(Year, PlotID, hostGroup) %>%
  summarise(individuals = max(IndNum)) %>%
  
  # reorder
  arrange(PlotID, hostGroup) %>%
  
  # get in wide format where each row is a caterpillar plot
  pivot_wider(id_cols = c(PlotID, hostGroup), names_from = Year, values_from = individuals, values_fill = 0) %>%
  
  # remove 0s
  filter(!is.na(hostGroup)) %>%
  
  # inconsistencies are averaged
  mutate(avg = ceiling((`2017` + `2018` + `2019`)/3)) %>%
  
  # put in the same format as other data for rbinding
  mutate(BlockID = str_extract(PlotID, "^.."),
         SiteID = str_extract(PlotID, "(?<=_).{3}(?=_)"),
         count = avg,
         vegPlot = str_extract(PlotID, "(?<=_)(2|3|4).*"),
         ringID = "actual",
         inPlot = NA,
         PlantID = "",
         PointID = PlotID,
         cPlot = str_extract(PlotID, "CP.")) %>%
  
  # reorganize
  ungroup() %>%
  dplyr::select(-`2017`, -`2018`, -`2019`, -avg, -PlotID)

vegetation$cPlot <- "None"

# combine
vegetation <- rbind(vegetation, actualPlots)

# adding in 0 caterpillar branches
caterpillars$CatID[caterpillars$CatID %in% ""] <- "NOCATS"

save(vegetation, caterpillars, branches,
     treeGroups, file = "data/catAndVeg.RData")

rm(list = ls())
