##%######################################################%##
#                                                          #
####   HABITAT CONNECTIVITY VISITATION DISSIMILARITY    ####
#                                                          #
##%######################################################%##

# This script recreates the habitat connectivity by visitation dissimilarity analysis

# PACKAGES ----------------------------------------------------------------

library(tidyverse)
library(broom)
library(cowplot)
library(vegan)

# CUSTOM STUFF ------------------------------------------------------------

# some colors

bif18_color <- "#31688EFF"
bif15_color <- "#31688EFF"
vos_color <- "#FDE725FF"


# RAW DATA ----------------------------------------------------------------

# capture data
sra_captures <- read_csv("./data_raw/sra_captures_post_COLONY.csv")
vsw15 <- filter(sra_captures, species == "vosnesenskii", year == 2015)
bsw15 <- filter(sra_captures, species == "bifarius", year == 2015)
bsw18 <- filter(sra_captures, species == "bifarius", year == 2018)

# site metadata
sites_metadata <- read_csv("./data_raw/sites_metadata.csv")



# WRANGLING ----------------------------------------------------------------

# rather inefficient way of making matrix, but works
vsw_matrix <- vsw15 %>% 
  select(unique.ID, site, date, plant.capture, FamilyID) %>% 
  filter(FamilyID != "s") %>% 
  select(unique.ID, FamilyID, site) %>% 
  mutate(JM01 = if_else(site=="JM01", T, F),
         JM02 = if_else(site=="JM02", T, F),
         JM03 = if_else(site=="JM03", T, F),
         JM04 = if_else(site=="JM04", T, F),
         JM05 = if_else(site=="JM05", T, F),
         JM06 = if_else(site=="JM06", T, F),
         JM07 = if_else(site=="JM07", T, F),
         JM10 = if_else(site=="JM10", T, F),
         JMX = if_else(site=="JMX", T, F),
         JMY = if_else(site=="JMY", T, F),
         SH = if_else(site=="SH", T, F),
         IN = if_else(site=="IN", T, F))  %>% 
  group_by(FamilyID) %>% 
  summarise_if(is.logical, sum) %>% 
  gather(site, count, -FamilyID) %>% 
  filter(site != "RM", site != "IN", site != "SH") %>% 
  spread(FamilyID, count) %>% 
  column_to_rownames(var = "site")

nmds.distance <- sites_metadata %>% 
  filter(year == 2015) %>% 
  mutate(k=1) %>% 
  full_join(.,., by = "k") %>% 
  filter(site.x != site.y) %>%
  mutate(dist = sqrt((x.center.x - x.center.y)^2 + (y.center.x - y.center.y)^2), #nmds.distance = sqrt((NMDS1.x - NMDS1.y)^2 + (NMDS2.x - NMDS2.y)^2), 
         site.separation = case_when(
           habitat.x == "meadow" & habitat.y == "meadow" ~ "meadow only",
           habitat.x == "meadow" & habitat.y == "forest" ~ "mixed",
           habitat.x == "forest" & habitat.y == "meadow" ~ "mixed",
           habitat.x == "forest" & habitat.y == "forest" ~ "forest"
         ))


vsw_d_matrix <- vegdist(vsw_matrix, method = "jaccard")

vsw_dissim <- tidy(vsw_d_matrix) %>% 
  rename(site.y = item1, site.x = item2, dis.distance = distance)

nmds.dis.distance <- inner_join(nmds.distance, vsw_dissim, by = c("site.x" = "site.x", "site.y" = "site.y"))  %>% 
  mutate(
  hab.cross = case_when(
    site_direction.y == site_direction.x  ~ "same",
    site_direction.y !=  site_direction.x ~ "different"
  ))


nmds_plot <- nmds.dis.distance %>% 
  ggplot(., aes(x = dist, y = dis.distance)) +
  geom_point(shape = 21, color = "black", aes(fill = hab.cross)) +
  geom_smooth(method = "lm", color = "black", size = 0.5, se = F) +
  theme_classic() +
  scale_fill_viridis_d(option = "A", labels = c("Different", "Same")) +
  labs(x ="", y = "Jaccard dissimilarity index", fill = "Vegetation type \n between sites", title = "B. vosnesenskii") +
  theme(legend.position = "none", plot.title = element_text(face = "italic"), axis.text.x = element_blank())



bsw15_matrix <- bsw15 %>% 
  select(unique.ID, site, date, plant.capture, FamilyID) %>% 
  filter(FamilyID != "s") %>% 
  select(unique.ID, FamilyID, site) %>% 
  mutate(#JM01 = if_else(site=="JM01", T, F),
         JM02 = if_else(site=="JM02", T, F),
         JM03 = if_else(site=="JM03", T, F),
         JM04 = if_else(site=="JM04", T, F),
         #JM05 = if_else(site=="JM05", T, F),
         JM06 = if_else(site=="JM06", T, F),
         JM07 = if_else(site=="JM07", T, F),
         JM10 = if_else(site=="JM10", T, F),
         #JMX = if_else(site=="JMX", T, F),
         #JMY = if_else(site=="JMY", T, F),
         SH = if_else(site=="SH", T, F),
         IN = if_else(site=="IN", T, F))  %>% 
  group_by(FamilyID) %>% 
  summarise_if(is.logical, sum) %>% 
  gather(site, count, -FamilyID) %>% 
  filter(site != "RM", site != "IN", site != "SH") %>% 
  spread(FamilyID, count) %>% 
  column_to_rownames(var = "site")


bsw15_d_matrix <- vegdist(bsw15_matrix, method = "jaccard")

bsw15_dissim <- tidy(bsw15_d_matrix) %>% 
  rename(site.y = item1, site.x = item2, dis.distance = distance)

bsw15_nmds.dis.distance <- inner_join(nmds.distance, bsw15_dissim, by = c("site.x" = "site.x", "site.y" = "site.y"))  %>% 
  mutate(
    hab.cross = case_when(
      site_direction.y == site_direction.x  ~ "same",
      site_direction.y !=  site_direction.x ~ "different"
    ))


nmds_plot_bsw15 <- bsw15_nmds.dis.distance %>% 
  ggplot(., aes(x = dist, y = dis.distance)) +
  geom_point(shape = 21, color = "black", aes(fill = hab.cross)) +
  geom_smooth(method = "lm", color = "black", size = 0.5, se = F) +
  theme_classic() +
  scale_fill_viridis_d(option = "A", labels = c("Different", "Same")) +
  labs(x ="", y = "Jaccard dissimilarity index", fill = "Vegetation type \n between sites", title = "B. vosnesenskii") +
  theme(legend.position = "none", plot.title = element_text(face = "italic"), axis.text.x = element_blank())










bsw_matrix <- bsw18 %>% 
  select(unique.ID, site, date, plant.capture, FamilyID) %>% 
  filter(FamilyID != "s") %>% 
  select(unique.ID, FamilyID, site) %>% 
  mutate(mdw_north = if_else(site=="mdw_north", T, F),
         mdw_south = if_else(site=="mdw_south", T, F),
         forest_north = if_else(site=="forest_north", T, F),
         forest_south = if_else(site=="forest_south", T, F),
         crk.x = if_else(site=="crk.x", T, F),
         north_east = if_else(site=="north_east", T, F),
         sh = if_else(site=="sh", T, F))  %>% 
  group_by(FamilyID) %>% 
  summarise_if(is.logical, sum) %>% 
  gather(site, count, -FamilyID) %>% 
  #filter(site != "RM", site != "IN", site != "SH") %>% 
  #filter(site != "RM") %>% 
  spread(FamilyID, count) %>% 
  column_to_rownames(var = "site")



bsw.distance <- sites_metadata %>% 
  filter(year == 2018) %>% 
  mutate(k=1) %>% 
  full_join(.,., by = "k") %>%
  filter(site.x != site.y) %>%
  mutate(dist = sqrt((x.center.x - x.center.y)^2 + (y.center.x - y.center.y)^2),
         site.separation = case_when(
           habitat.x == "meadow" & habitat.y == "meadow" ~ "meadow only",
           habitat.x == "meadow" & habitat.y == "forest" ~ "mixed",
           habitat.x == "forest" & habitat.y == "meadow" ~ "mixed",
           habitat.x == "forest" & habitat.y == "forest" ~ "forest"
         ), 
         hab.cross = if_else(site.separation == "mixed", "different", "same"))

# dissimilarity scores ----------------------------------------------------

bsw_d_matrix <- vegdist(bsw_matrix, method = "jaccard")

bsw_dissim <- tidy(bsw_d_matrix) %>% 
  rename(site.y = item1, site.x = item2, dis.distance = distance)

bsw.dis.distance <- inner_join(bsw.distance, bsw_dissim, by = c("site.x" = "site.x", "site.y" = "site.y"))


bif_dissim <- bind_rows(bsw15_nmds.dis.distance, bsw.dis.distance) %>% 
  filter(dist < 1000)

nmds_plot_bsw <- bif_dissim %>% 
  ggplot(., aes(x = dist, y = dis.distance)) +
  geom_point(shape = 21, color = "black", aes(fill = hab.cross)) +
  #geom_jitter(shape = 21, color = "black", width= 0.5) +
  #geom_smooth(method = "lm", color = "black", size = 0.5, se = F) +
  theme_classic() +
  scale_fill_viridis_d(option = "A", labels = c("Heterogenous", "Homogenous")) +
  labs(x ="Linear distance between sites (m)", y = "Jaccard dissimilarity index", fill = "Vegetation type \n between sites", title = "B. bifarius") +
  theme(legend.position = c(0.8, 0.2), legend.background = element_rect(color = "black"), plot.title = element_text(face = "italic"))

nmds_plot_bsw


nmds_plot_combined <- plot_grid(nmds_plot +
  coord_cartesian(xlim = c(100, 800), ylim = c(0.5, 1)), nmds_plot_bsw +
  coord_cartesian(xlim = c(100, 800), ylim = c(0.5, 1)), nrow = 2, labels = c("A", "B"), rel_heights = c(0.95, 1))

nmds_plot_combined
ggsave("./figures/figure6.tiff", nmds_plot_combined, width = 5, height = 9, dpi = 600)

nmds_dist_vos_summary <- summary(lm(data = nmds.dis.distance, dis.distance ~ dist))
nmds_dist_vos_int_summary <- summary(lm(data = nmds.dis.distance, dis.distance ~ dist * hab.cross))
anova(lm(data = nmds.dis.distance, dis.distance ~ dist * hab.cross))

nmds_dist_bif_summary <- summary(lm(data = bif_dissim, dis.distance ~ dist))
nmds_dist_bif_int_summary <- summary(lm(data = bif_dissim, dis.distance ~ dist * hab.cross))
anova(lm(data = bif_dissim, dis.distance ~ dist * hab.cross))


# SAVE OUTPUT -------------------------------------------------------------

save(nmds_dist_vos_int_summary, nmds_dist_bif_int_summary, nmds_dist_bif_summary, nmds_dist_vos_summary, file = "./analyses_output/nmds_distance_objects.Rdata")
