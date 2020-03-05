# Recaptured bees from tagged vosnesenskii


# PACKAGES ----------------------------------------------------------------

library(tidyverse)


# RAW DATA ----------------------------------------------------------------

recaps <- read_csv("./data_raw/sra_recaptured_bees.csv")

vsw_col <- read_csv("./data_raw/sra_captures_post_COLONY.csv") %>% 
  filter(species == "vosnesenskii")


# WRANGLING ---------------------------------------------------------------

recap_bees = recaps %>% 
  select(unique.ID, lon, lat, site, date, plant.capture, recap.ID)

geno_bees = vsw_col %>% 
  select(unique.ID, lon, lat, bee.tag, site, date, plant.capture)

# match the two dataframes together and calculate the distance between the original observation and the recaptured location
recap_matched <- 
  inner_join(geno_bees, recap_bees, by = c("bee.tag"="recap.ID"), suffix=c(".1", ".2")) %>% 
  mutate(distance = sqrt((lon.1 - lon.2)^2+(lat.1 - lat.2)^2))



# DISPLAY RESULTS ---------------------------------------------------------

recap_table <- recap_matched %>% 
  select(site.1, site.2, distance) %>% 
  arrange(desc(distance))


# SAVE OUTPUT -------------------------------------------------------------

saveRDS(recap_table, "./analyses_output/recap_table.rds")
