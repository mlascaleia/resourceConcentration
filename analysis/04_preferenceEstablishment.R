# The purpose of this script is to write the functions for analyzing the patterns of polyphagous caterpillars
# the idea is to look at how the densities of all their host plants affect their density
# Also, this will consider dietary preference by using a use-by-availability analysis

# initialize ####
rm(list = ls())
load("data/catAndVeg.RData")

# Much like I did host groups for trees, I need to make cat groups
# as well as eliminate some cat species that aren't real species

# eliminates arguments from data
elim <- function(cat, c = caterpillars){
  c <- c[!c$CatID %in% cat,]
  return(c)
}

# groups arguments in data
grup <- function(g, gg, c = caterpillars){
  c$CatID[c$CatID %in% g] <- gg
  return(c)
}

# unIDed Geometrids
caterpillars <- elim("GEOMET")

# unIDed microlepidoptera
caterpillars <- elim("MICROL")

# unIDed Gelechids
caterpillars <- elim("GELECH")

# unIDed Tortricids
caterpillars <- elim("TORTRI")

# unknowns 
caterpillars <- elim("UNKNXX")

# grouping Eupithecia spp.
caterpillars <- grup(c("EUPISW", "EUPIVA", "EUPIXX"), "EUPIXX")

# grouping Eupsilia spp.
caterpillars <- grup(c("EUPSSC","EUPSXX"), "EUPSXX")

# grouping Lithophane spp.
caterpillars <- grup(c("LITHAN", "LITHHB", "LITHIP", "LITHXX"), "LITHXX")

# grouping Orthosia spp.
caterpillars <- grup(c("ORTHRU","ORTHXX"), "ORTHXX")

# get plants w/ fewer than 20 caterpillars found
badGroups <- caterpillars %>%
  group_by(hostGroup, .drop = F) %>%
  summarise(cc = n()) %>%
  filter(cc < 20)

# get total numbers of branches sampled (of relevant tree types)
branches <- branches[!branches$hostGroup %in% badGroups$hostGroup, ]

b <- branches %>%
  distinct(PointID, hostGroup, Year, .keep_all = T) %>%
  group_by(hostGroup) %>%
  summarise(totalBranch = sum(bSample)) %>%
  filter(totalBranch > 20)

# update relevant tree types
rm(badGroups)
goodGroups <- b$hostGroup

# get total numbers of each cat (that is on a relevant tree-type)
caterpillars <- caterpillars[caterpillars$hostGroup %in% goodGroups,]

p <- caterpillars %>% 
  group_by(CatID) %>%
  summarise(totalCat = n()) %>%
  filter(totalCat > 19)

# here I am going to use Chesson's preference equation to determine food preference
# But I'm going to do it without food loss
# Instead: preference = (ri/ni)/summation(rj/nj)
# without depletion means that nj remains constant
# in this equation, ri is the number of food items of type in the diet
# ni is the total amount of food of type i in the environment
# If the consumer only eats one plant, its value for that i will be 1
# If it does not eat that type of plant, its value for that i will be 0

# In our situation, the amount of available food is the number of branches of the tree
# So, for example, in this experiment I know that 870 ACERU branches were sampled in our experiment
# I also know that 105 IRIDEP were found in our experiment, and not all on ACERU
# Let's try to find preference of IRIDEP on ACERU before trying to iterate for cats or plants
# Chesson methods w/ and w/o abundance ####
# 
# ni <- b$totalBranch[b$hostGroup == "ACERU"]
# ri <- sum(ppCat$catCount[ppCat$catID == "IRIDEP" & ppCat$plantID == "ACERU"])
# 
# prefs <- apply(b, 1 ,FUN = function(br, cid, ni,ri){
#   nj <- as.numeric(br["totalBranch"])
#   rj <- sum(ppCat$catCount[ppCat$catID == cid & ppCat$plantID == br["hostGroup"]])
#   smed <- rj/nj
#   return(smed)
# }, ni = ni, ri = ri, cid = "IRIDEP")
# (ri/ni)/(sum(prefs)) # This is the preference of IRIDEP on red maple

# now iterate for each cat/plant interaction
# each tree only has one ni, so there's really nothing to iterate there
# nis <- b$totalBranch # this is always the same across all cats

# This for loop performs a mathematical function from Chesson (1983)
# That function is a = (ri/ni)/Σ(rj/nj)
# in the context of this experiment, here are what all those values mean
# ri = the number of branches of tree species i on which caterpillar species c is found (note, c is not part of the formula)
# ni = the number of total branches of tree species i sampled
# rj = the number of branches of species j on which caterpillar species c is found
# nj = the number of total branches of species j
# note, the denominator is the sum of all the values of ni/ri for the caterpillar species
# a = the caterpillars "electivity" on a specific plant
# The goal of this loop is to get an a value for every caterpillar, for every tree species

# chessonNoAbund <- apply(b, 1, FUN = function(br){                               # this apply loop gives us every plant species to determine which plant we're looking at
#   piq <- as.character(br["hostGroup"])                                          # get the name of the plant we're looking at
#   oneTree <- apply(p, 1, FUN = function(pr, br, pic){                           # this apply loop esablishes which caterpillar species we're looking at. Now we have cat and plant
#     ciq <- as.character(pr["catID"])                                            # get the name of the caterpillar we're looking at
#     ni <- b$totalBranch[b$hostGroup == pic]                                     # get the total number of branches that were sampled for the plant species
#     ri <- length(ppCat$catCount[ppCat$catID == ciq & ppCat$hostGroup == pic])   # get the total number of branches of our plant that had our caterpillar on it
#     numer <- ri/ni                                                              # determine the numerator for this specific value of a
#     prefs <- apply(b, 1 ,FUN = function(br2, cid){                              # this apply loop allows us to do the summation (Σ) in the denominator, and requires us to iterate through the b data frame a within our iteration of the b data frame
#       nj <- as.numeric(br2["totalBranch"])                                      # get the total branches for species j (will do for all tree species)
#       bbb <- as.character(br2["hostGroup"])                                     # get the name of species j
#       rj <- length(ppCat$catCount[ppCat$catID == cid & ppCat$hostGroup == bbb]) # find all instances of our caterpillar on species j
#       smed <- rj/nj                                                             # get one value for the sigma in the denominator
#       return(smed)                                                              # return a vector that is the same length as the number of species we are looking at for the denominator
#     }, cid = ciq)                                                               # close apply loop (and make sure we're looking at the correct caterpillar)
#     denom <- sum(prefs)                                                         # establish the denominator by summing (Σ) the values in the prefs vector
#     pref <- (numer/denom)                                                       # get an output
#     return(pref)                                                                # return a for this cat/plant interaction
#   }, br = br, pic = piq)                                                        # close apply loop
# })                                                                              # close apply loop
# 
# colnames(chessonNoAbund) <- b$hostGroup # make the names in the matrix accurate
# rownames(chessonNoAbund) <- p$catID # make the rownames in the matrix accurate
# # chessonNoAbund <- round(chessonNoAbund, 3) # makes the matrix much easier to use


# Okay we have that now
# But it still doesn't take into account abundance on plants. I don't think that this works with the above equation
# I can try this out in two ways

# next, take into account when more than one caterpillar is found on any given branch, and scale accordingly
# There are other ways to do this, but let's do it within the framework of the previous apply for coding sake

cwa <- apply(b, 1, FUN = function(br){
  piq <- as.character(br["hostGroup"])
  oneTree <- apply(p, 1, FUN = function(pr, br, pic){
    ciq <- as.character(pr["CatID"])
    ni <- b$totalBranch[b$hostGroup == pic]
    ri <- length(unique(caterpillars$branchID[caterpillars$CatID == ciq & caterpillars$hostGroup == pic]))
    ab <- nrow(caterpillars[caterpillars$CatID == ciq & caterpillars$hostGroup == pic,])
    numer <- (ri/ni) * ab
      prefs <- apply(b, 1 ,FUN = function(br2, cid){
        nj <- as.numeric(br2["totalBranch"])
        bbb <- as.character(br2["hostGroup"])
        rj <- length(unique(caterpillars$branchID[caterpillars$CatID == cid & caterpillars$hostGroup == bbb]))
        smed <- rj/nj
        return(smed)
      }, cid = ciq)
    denom <- sum(prefs)
    pref <- (numer/denom)
    return(pref)
  }, br = br, pic = piq)
})

colnames(cwa) <- b$hostGroup
rownames(cwa) <- p$CatID

electMatrix <- (cwa/rowSums(cwa))

electivity <- as.data.frame(electMatrix) %>% 
  mutate(CatID = rownames(electMatrix)) %>%
  dplyr::select(CatID, everything()) %>%
  gather('species', 'electivity', -CatID)


#### finalize ####

save(electMatrix, vegetation, caterpillars, branches, file = "data/electivity.RData")
rm(list = ls())
