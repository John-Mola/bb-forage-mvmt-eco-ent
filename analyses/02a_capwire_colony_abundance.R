# Script to obtain Capwire estimates of Colony abundance


# PACKAGES ----------------------------------------------------------------

library(tidyverse)
library(beepr)
library(capwire)
library(foreach)


# CUSTOM FUNCTIONS --------------------------------------------------------

# Function to make lists with actual names (and not just index)
namedList <- function(...) {
  L <- list(...)
  snm <- sapply(substitute(list(...)),deparse)[-1]
  if (is.null(nm <- names(L))) nm <- snm
  if (any(nonames <- nm=="")) nm[nonames] <- snm[nonames]
  setNames(L,nm)
}

# Does everything I want except doesn't create a metadata column (Current workaround: Attach a comment to each inputted data frame)
list_colony_counter = function(listOfPops, maxPop, bootCount){
  
  require(tidyverse)
  require(foreach)
  require(capwire)
  
  foreach::foreach(i=listOfPops, .combine='rbind') %do% {
    
    sp_full = i %>% 
      group_by(ClusterIndex) %>% 
      summarise(n=n())
    
    captable=buildClassTable(sp_full$n)
    
    res.tirm <- fitTirm(data=captable, max.pop=maxPop)
    res.bootstrap <- bootstrapCapwire(res.tirm, bootstraps = bootCount)
    
    #df.name <- names(listOfPops[i])
    dFnm <- comment(i)
    
    data_frame(df.name = dFnm, ml.colony.num =  res.tirm[[3]], CI.lower= res.bootstrap[[2]][[1]], CI.upper = res.bootstrap[[2]][[2]])
    
  }}


# RAW DATA ----------------------------------------------------------------

sra_captures <- read_csv("./data_raw/sra_captures_post_COLONY.csv")

# WRANGLING ---------------------------------------------------------------

# filter to datasets for each year/species. vsw15 is "vosnesenskii sierra workers from 2015", bsw = bifarius yaddi yaddi yadda...
vsw15 <- filter(sra_captures, species == "vosnesenskii", year == 2015)
bsw15 <- filter(sra_captures, species == "bifarius", year == 2015)
bsw18 <- filter(sra_captures, species == "bifarius", year == 2018)

# this provides a commented "name" to each dataframe that goes into the list. 
comment(vsw15) <- "vsw15"
comment(bsw15) <- "bsw15"
comment(bsw18) <- "bsw18"

# saves a named list from the comments
list_all = namedList(vsw15, bsw15, bsw18)

# RUN CAPWIRE -------------------------------------------------------------

# run capwire with a maximum population size of 5000 and 95%CI from 1000 bootstraps (spoke with package developer and he said 10000 unlikely to give any better result and it takes fooorrreevveerrr)
capwire_out_all <- list_colony_counter(list_all,5000,1000) ; beepr::beep(sound="mario") 

write_csv(capwire_out_all, "./analyses_output/capwire_out_all_1kboot.csv")
