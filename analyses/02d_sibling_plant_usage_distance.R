##%######################################################%##
#                                                          #
####                SIBLING PLANT USAGE                 ####
#                                                          #
##%######################################################%##

# This script recreates the analysis of sibling plant usage across distance

# PACKAGES ----------------------------------------------------------------

library(tidyverse)
library(cowplot)
library(lme4)

# CUSTOM STUFF ------------------------------------------------------------

# some colors

bif18_color <- "#31688EFF"
bif15_color <- "#31688EFF"
vos_color <- "#FDE725FF"

# RAW DATA ----------------------------------------------------------------

sra_cap <- read_csv("./data_raw/sra_captures_post_COLONY.csv")


# WRANGLING ---------------------------------------------------------------

#VOSNESENSKii ----
vsw_15 <- filter(sra_cap, species == "vosnesenskii") %>% 
  rename(x=lon,y=lat)

# OBSERVED DISTANCEs

vsw_15$FamilyID = as.factor(vsw_15$FamilyID)

cm_ind_all_15 = vsw_15 %>% 
  filter(FamilyID !="s") %>% 
  mutate(occur = !FamilyID %in% FamilyID[duplicated(FamilyID)]
  ) %>% 
  filter(occur == F) %>% 
  dplyr::select(-occur)

# Joining the data frame to itself, to calculate all pairwise combinations of siblings. 
pairs_obs_15 = left_join(cm_ind_all_15, cm_ind_all_15, by = "FamilyID", suffix=c(".A", ".B")) %>% 
  filter(unique.ID.A != unique.ID.B) %>% 
  rowwise %>%
  mutate(name = toString(sort(c(unique.ID.A,unique.ID.B)))) %>% 
  distinct(name, .keep_all=T) %>% 
  mutate(distance = sqrt((x.A - x.B)^2+(y.A - y.B)^2))

#bifarius 15 ----
bsw_15 <- filter(sra_cap, species == "bifarius", year == 2015) %>% 
  rename(x=lon,y=lat)

# OBSERVED DISTANCEs

bsw_15$FamilyID = as.factor(bsw_15$FamilyID)

cm_ind_all_bsw15 = bsw_15 %>% 
  filter(FamilyID !="s") %>% 
  mutate(occur = !FamilyID %in% FamilyID[duplicated(FamilyID)]
  ) %>% 
  filter(occur == F) %>% 
  dplyr::select(-occur)

# Joining the data frame to itself, to calculate all pairwise combinations of siblings. 
pairs_obs_bsw15 = left_join(cm_ind_all_bsw15, cm_ind_all_bsw15, by = "FamilyID", suffix=c(".A", ".B")) %>% 
  filter(unique.ID.A != unique.ID.B) %>% 
  rowwise %>%
  mutate(name = toString(sort(c(unique.ID.A,unique.ID.B)))) %>% 
  distinct(name, .keep_all=T) %>% 
  mutate(distance = sqrt((x.A - x.B)^2+(y.A - y.B)^2))

#bifarius 18 ----
bsw_18 <- filter(sra_cap, species == "bifarius", year == 2018) %>% 
  rename(x=lon,y=lat)

# OBSERVED DISTANCEs

bsw_18$FamilyID = as.factor(bsw_18$FamilyID)

cm_ind_all_bsw18 = bsw_18 %>% 
  filter(FamilyID !="s") %>% 
  mutate(occur = !FamilyID %in% FamilyID[duplicated(FamilyID)]
  ) %>% 
  filter(occur == F) %>% 
  dplyr::select(-occur)

# Joining the data frame to itself, to calculate all pairwise combinations of siblings. 
pairs_obs_bsw18 = left_join(cm_ind_all_bsw18, cm_ind_all_bsw18, by = "FamilyID", suffix=c(".A", ".B")) %>% 
  filter(unique.ID.A != unique.ID.B) %>% 
  rowwise %>%
  mutate(name = toString(sort(c(unique.ID.A,unique.ID.B)))) %>% 
  distinct(name, .keep_all=T) %>% 
  mutate(distance = sqrt((x.A - x.B)^2+(y.A - y.B)^2))

# COMBINE THE YEARS AND SPECIES INTO A DATAFRAME ----
sib_sep_combo <- dplyr::select(pairs_obs_bsw15, unique.ID.A, unique.ID.B, plant.capture.A, plant.capture.B, date.A, date.B, distance, FamilyID) %>% 
  bind_rows(., dplyr::select(pairs_obs_bsw18, unique.ID.A, unique.ID.B, plant.capture.A, plant.capture.B, date.A, date.B, distance, FamilyID)) %>% 
  mutate(species = "bifarius") %>% 
  bind_rows(., dplyr::select(pairs_obs_15, unique.ID.A, unique.ID.B, plant.capture.A, plant.capture.B, date.A, date.B, distance, FamilyID)) %>% 
  mutate(diff_plant = if_else(plant.capture.A == plant.capture.B, 0, 1), species = if_else(!is.na(species), "bifarius", "vosnesenskii"),
         days_between = abs(date.A - date.B))


# RENDER OUT THE STATS ----------------------------------------------------

sib_sep_combo %>% 
  filter(distance < 1000) %>%
  glm(diff_plant ~ distance*species, data = ., family = "binomial") %>% 
  summary()

sib_usage_stats_summary <- summary(lmList(diff_plant ~ distance | species, data=filter(sib_sep_combo, distance < 1000), family = "binomial"))

sib_usage_vos_stats_summary <- sib_sep_combo %>% 
  filter(species == "vosnesenskii") %>% 
  glm(diff_plant ~ distance, data = ., family = "binomial") %>% 
  summary()

# PLOT --------------------------------------------------------------------


sep_plants_p1 <- sib_sep_combo %>% 
  filter(distance < 1000) %>% 
  ggplot(., aes(x = distance, y = diff_plant, color = species, fill = species)) +
  #geom_jitter(alpha = 0.7, shape = 21, color = "black", width = 0, height = 0.05) +
  geom_point(alpha = 0.4, shape = 21, color = "black", position=ggstance::position_dodgev(height=0.1)) +
  geom_smooth(method = "glm", 
              method.args = list(family = "binomial"), 
              se = FALSE) +
  scale_color_manual(values = c(bif18_color, vos_color), labels = c("B. bifarius", "B. vosnesenskii")) +
  scale_fill_manual(values = c(bif18_color, vos_color), labels = c("B. bifarius", "B. vosnesenskii")) +
  theme_classic() +
  labs(x = "", y = "Probability two siblings were \n captured on different floral species") +
  theme(legend.position = c(0.8, 0.3), legend.title = element_blank(), legend.text = element_text(face = "italic"))

#sep_plants_p1

sep_plants_p2 <- sib_sep_combo %>% 
  filter(species == "vosnesenskii") %>% 
  ggplot(., aes(x = distance, y = diff_plant, color = species, fill = species)) +
  geom_point(alpha = 0.4, shape = 21, color = "black") +
  geom_smooth(method = "glm", 
              method.args = list(family = "binomial"), 
              se = FALSE) +
  scale_color_manual(values = c(vos_color)) +
  scale_fill_manual(values = c(vos_color)) +
  theme_classic() +
  labs(x = "", y = "") +
  theme(legend.position = "none", axis.text.y = element_blank())



sep_plot_left <- ggdraw(add_sub(plot_grid(sep_plants_p1, sep_plants_p2, rel_widths = c(1.2,1), labels = c("A", "B")), "Sibling separation distance (m)", vpadding=grid::unit(0,"lines"), vjust = -1))

sep_plot_left

ggsave2("./figures/figure4.tiff", plot = sep_plot_left, width = 9, height = 4, dpi = 600)


# REPEATING THE ANALYSIS WITHIN DAYS --------------------------------------

# within 3 days

sib_usage_3day_summary <- summary(lmList(diff_plant ~ distance | species, data=filter(sib_sep_combo, distance < 1000, days_between <= 3), family = "binomial"))

# within 1 to 20 days
for (i in 1:20) {

  summary(lmList(diff_plant ~ distance | species, data=filter(sib_sep_combo, distance < 1000, days_between <= i), family = "binomial")) %>% 
  print()
}



# SAVE OUTPUTS ------------------------------------------------------------

save(sib_usage_stats_summary, sib_usage_vos_stats_summary, sib_usage_3day_summary, file = "./analyses_output/sib_usage_objects.Rdata")
