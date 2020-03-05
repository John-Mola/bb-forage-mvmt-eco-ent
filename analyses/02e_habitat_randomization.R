##%######################################################%##
#                                                          #
####         HABITAT CONNECTIVITY RANDOMIZATION         ####
#                                                          #
##%######################################################%##


# This script recreates the randomization of family assignment across habitat types

# PACKAGES ----------------------------------------------------------------

library(tidyverse)
library(cowplot)

# CUSTOM STUFF ------------------------------------------------------------

# some colors

bif18_color <- "#31688EFF"
bif15_color <- "#31688EFF"
vos_color <- "#FDE725FF"

# function for reshuffling family membership (only for core sites)
type_shuffle_f_m = function(df) {
  
  # # Shuffle families
  df_shuff = df %>% mutate(FamilyID = sample(FamilyID))
  df_shuff$FamilyID = as.factor(df_shuff$FamilyID)

  # Calculate the types of sites a family is found in  
  n_type_col = df_shuff %>% 
    # this became redundant after rendering metadata, but kept because it works and I don't want to break it
    filter(site %in% forest_sites | site %in% meadow_sites) %>% 
    mutate(site.type = if_else(site %in% forest_sites, "forest", if_else(site %in% meadow_sites, "meadow", "distant"))) %>% 
    group_by(FamilyID, site.type) %>% 
    summarise(n=n()) %>% 
    spread(site.type, n) %>% 
    filter(FamilyID !="s") %>% 
    mutate(col.habs = case_when(
      forest >= 1 & is.na(meadow)  ~ "forest only",
      is.na(forest) & meadow >= 1  ~ "meadow only",
      forest >= 1 & meadow >= 1  ~ "forest and meadow",
      TRUE                                           ~ "other"
    ))
  
}


# RAW DATA ----------------------------------------------------------------

# capture data
sra_captures <- read_csv("./data_raw/sra_captures_post_COLONY.csv")
vsw15 <- filter(sra_captures, species == "vosnesenskii", year == 2015)
bsw15 <- filter(sra_captures, species == "bifarius", year == 2015)
bsw18 <- filter(sra_captures, species == "bifarius", year == 2018)

# site metadata
sites_metadata <- read_csv("./data_raw/sites_metadata.csv")



# WRANGLING ---------------------------------------------------------------


# pulling vectors of site names within categories
forest_sites <- sites_metadata %>% filter(habitat == "forest") %>% pull(site)
meadow_sites <- sites_metadata %>% filter(habitat == "meadow") %>% pull(site)
distant_sites <- sites_metadata %>% filter(habitat == "distant") %>% pull(site)

# set to do 1000 shuffles, but results don't really change when doing 10,000 (mostly it just plots cleaner this way)
n_shuffles = 1000



# VOSNESENSKII ------------------------------------------------------------

# calculating observed "crossings"
n_type_col_15 = vsw15 %>% 
  filter(FamilyID !="s", site %in% forest_sites | site %in% meadow_sites) %>% 
  mutate(site.type = if_else(site %in% forest_sites, "forest", if_else(site %in% meadow_sites, "meadow", "distant"))) %>% 
  group_by(FamilyID, site.type) %>% 
  summarise(n=n()) %>% 
  spread(site.type, n) %>% 
  mutate(col.habs = case_when(
    forest >= 1 & is.na(meadow) ~ "forest only",
    is.na(forest) & meadow >= 1 ~ "meadow only",
    forest >= 1 & meadow >= 1   ~ "forest and meadow",
    TRUE                                           ~ "other"
  )) %>%
  filter(col.habs == "forest and meadow")

# Find the randomized expectation
shuffled_types_15 <-lapply(1:n_shuffles, function(i) data.frame(iteration = i, type_shuffle_f_m(vsw15))) %>% 
  bind_rows() %>% 
  filter(col.habs == "forest and meadow")

# Summarise counts across iterations for easier plotting
sumR_sites_15 = shuffled_types_15 %>% 
  group_by(iteration, col.habs) %>%
  summarise(n=n()) %>% 
  filter(col.habs == "forest and meadow")

# Find mean for easier plotting
sumR_means_15 = sumR_sites_15 %>% 
  group_by(col.habs) %>% 
  summarise(mean_r = mean(n)) %>% 
  filter(col.habs == "forest and meadow")



# BIFARIUS 2015 ------------------------------------------------------------

# calculating observed "crossings"

n_type_col_bsw15 = bsw15 %>% 
  filter(FamilyID !="s", site %in% forest_sites | site %in% meadow_sites) %>% 
  mutate(site.type = if_else(site %in% forest_sites, "forest", if_else(site %in% meadow_sites, "meadow", "distant"))) %>% 
  group_by(FamilyID, site.type) %>% 
  summarise(n=n()) %>% 
  spread(site.type, n) %>% 
  mutate(col.habs = case_when(
    forest >= 1 & is.na(meadow) ~ "forest only",
    is.na(forest) & meadow >= 1 ~ "meadow only",
    forest >= 1 & meadow >= 1   ~ "forest and meadow",
    TRUE                                           ~ "other"
  )) %>%
  filter(col.habs == "forest and meadow")

# Find the randomized expectation
shuffled_types_bsw15 <-lapply(1:n_shuffles, function(i) data.frame(iteration = i, type_shuffle_f_m(bsw15))) %>% 
  bind_rows() %>% 
  filter(col.habs == "forest and meadow")

# Summarise counts across iterations for easier plotting
sumR_sites_bsw15 = shuffled_types_bsw15 %>% 
  group_by(iteration, col.habs) %>%
  summarise(n=n()) %>% 
  filter(col.habs == "forest and meadow")

# Find mean for easier plotting
sumR_means_bsw15 = sumR_sites_bsw15 %>% 
  group_by(col.habs) %>% 
  summarise(mean_r = mean(n)) %>% 
  filter(col.habs == "forest and meadow")


# BIFARIUS 2018 ------------------------------------------------------------

# calculating observed "crossings"

n_type_col_bsw18 = bsw18 %>% 
  filter(FamilyID !="s", site %in% forest_sites | site %in% meadow_sites) %>% 
  mutate(site.type = if_else(site %in% forest_sites, "forest", if_else(site %in% meadow_sites, "meadow", "distant"))) %>% 
  group_by(FamilyID, site.type) %>% 
  summarise(n=n()) %>% 
  spread(site.type, n) %>% 
  mutate(col.habs = case_when(
    forest >= 1 & is.na(meadow) ~ "forest only",
    is.na(forest) & meadow >= 1 ~ "meadow only",
    forest >= 1 & meadow >= 1   ~ "forest and meadow",
    TRUE                                           ~ "other"
  )) %>%
  filter(col.habs == "forest and meadow")

# Find the randomized expectation
shuffled_types_bsw18 <-lapply(1:n_shuffles, function(i) data.frame(iteration = i, type_shuffle_f_m(bsw18))) %>% 
  bind_rows() %>% 
  filter(col.habs == "forest and meadow")

# Summarise counts across iterations for easier plotting
sumR_sites_bsw18 = shuffled_types_bsw18 %>% 
  group_by(iteration, col.habs) %>%
  summarise(n=n()) %>% 
  filter(col.habs == "forest and meadow")

# Find mean for easier plotting
sumR_means_bsw18 = sumR_sites_bsw18 %>% 
  group_by(col.habs) %>% 
  summarise(mean_r = mean(n)) %>% 
  filter(col.habs == "forest and meadow")




# COMBINING AND PLOTTING THINGS -------------------------------------------


ntc_bif_18 <- n_type_col_bsw18 %>% 
  mutate(species = "B. bifarius", year = 2018)

ntc_bif_15 <- n_type_col_bsw15 %>% 
  mutate(species = "B. bifarius", year = 2015)

ntc_vos_15 <- n_type_col_15 %>% 
  mutate(species = "B. vosnesenskii", year = 2015)

ntc_combo <- bind_rows(ntc_bif_18, ntc_bif_15, ntc_vos_15)

sumR_bif_18 <- sumR_sites_bsw18 %>% 
  mutate(species = "B. bifarius", year = 2018)

sumR_bif_15 <- sumR_sites_bsw15 %>% 
  mutate(species = "B. bifarius", year = 2015)

sumR_vos_15 <- sumR_sites_15 %>% 
  mutate(species = "B. vosnesenskii", year = 2015)

sumR_combo <- bind_rows(sumR_bif_18, sumR_bif_15, sumR_vos_15)

habs_plot <- ggplot() +
  theme_classic() + 
  labs(x = "Species", y = "Count of colonies detected \n in both forest and meadow") +
  geom_jitter(data = sumR_combo, aes(x=species, y=n, fill=species), color = "black", stat="identity", alpha=0.05, position = position_jitterdodge(jitter.width = 0.4, dodge.width = 1)) +
  geom_point(data=ntc_combo, aes(x=species, fill = species), color="black", position = position_dodge(width = 1), shape =21, stat = "count", size = 5) +
  scale_fill_manual(values = c(bif18_color, vos_color)) +
  theme(legend.position = "none", axis.text.x = element_text(face = "italic"))

habs_plot

ggsave2("./figures/figure5.tiff", habs_plot, width = 3, height = 3, dpi = 600)


# SUMMARY STATS -----------------------------------------------------------


observed_crossings <- ntc_combo %>% 
  group_by(species) %>% 
  summarise(n=n())

randomized_crossings <- sumR_combo %>% 
  group_by(species) %>% 
  summarise(mean = mean(n), min = min(n), max = max(n), sd = sd(n))


# SAVE OUTPUT -------------------------------------------------------------

save(observed_crossings, randomized_crossings, file = "./analyses_output/habitat_crossings.Rdata")
