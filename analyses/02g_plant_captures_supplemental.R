# Rendering supplemental table for publication and submission data. 


# Supplemental Table of Captures and plant identity
# Desired columns: Plant Species, # B. vos, # B. bif 15, # B bif 18, Plant occurs in Forest or Meadow or Both, 


# PACKAGES ----------------------------------------------------------------

library(tidyverse)


# WRANGLING ---------------------------------------------------------------

sra_captures <- read_csv("./data_raw/sra_captures_post_COLONY.csv")
  
sra_captures %>% 
  select(plant.capture, species, site, year) %>% 
  group_by(species, year, plant.capture) %>% 
  tally() %>% 
  ungroup() %>% 
pivot_wider(names_from = c(species, year), values_from = n, values_fill = list(n = 0)) %>% 
  mutate(plant.capture = if_else(is.na(plant.capture), "Not recorded", plant.capture),
         total.captures = bifarius_2015+bifarius_2018+vosnesenskii_2015) %>% 
  select(plant.capture, vosnesenskii_2015, bifarius_2015, bifarius_2018, everything()) %>% 
  arrange(plant.capture)

#write_csv(df_plant_captures, "../supplemental/plant_captures_table.csv")

