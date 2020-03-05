##%######################################################%##
#                                                          #
####        PAIRWISE SIBLING SEPARATION DISTANCE        ####
####       AND COLONY SPECIFIC FORAGING DISTANCE        ####
#                                                          #
##%######################################################%##

# This script is used to calculate sibling separation distances using relative frequency, area-adjusted measurement, and find the colony-specific foraging distance. 

# I apologize for the state of the code, but it runs and is reproducible and will allow one to see how we generate the results of the manuscript. If you have any questions please contact me at jmola@usgs.gov

# PACKAGES ----------------------------------------------------------------

library(tidyverse)
library(cowplot)
library(lme4)


# CUSTOM FUNCTIONS --------------------------------------------------------

#despite the fantastic naming, this calculates colony-specific foraging distance from a dataframe and a minimum family size (mfs), I usually use 2 but it was nice to explore and make sure it wasn't sensitive to that. 
sib_separation = function(df, mfs)
{
  
  df$date = as.Date(df$date, "%m/%d/%y")
  
  df$FamilyID = as.factor(df$FamilyID)
  
  cm_ind = df %>% 
    filter(FamilyID !="s")
  
  cm_fam = cm_ind %>% 
    group_by(FamilyID) %>% 
    summarise(x.fam.center = mean(x), y.fam.center = mean(y), n = n()) %>% 
    filter(n >= mfs)
  
  cm_combo = inner_join(cm_ind,cm_fam, by=c("FamilyID"))
  
  cm_filt = cm_combo %>% 
    mutate(ctr.dist = sqrt((x - x.fam.center)^2+(y - y.fam.center)^2
    )) 
  
  cm_filt
  
  
}

# some colors

bif18_color <- "#31688EFF"
bif15_color <- "#31688EFF"
vos_color <- "#FDE725FF"

# RAW DATA ----------------------------------------------------------------

sra_cap <- read_csv("./data_raw/sra_captures_post_COLONY.csv")


# VOSNESENSKII ---------------------------------------------------------------

# set annulus breaks to 100m
break_dist <- 100

vsw_15 <- filter(sra_cap, species == "vosnesenskii") %>% 
  rename(x=lon,y=lat)

# CALCULATING ALL POSSIBLE PAIRWISE DISTANCEs
# this chunk creates a column and then joins it to itself and calculates all of the distances between every observation
pairs_poss_15_incomplete = vsw_15 %>% 
  mutate(k=1) %>% 
  full_join(.,., by = "k") %>% 
  filter(unique.ID.x != unique.ID.y) %>%
  mutate(dist = sqrt((x.x - x.y)^2 + (y.x - y.y)^2)) %>%
  dplyr::select(-k)

IDx = pairs_poss_15_incomplete$unique.ID.x # be careful with this chunk because it could be redundant with other names
IDy = pairs_poss_15_incomplete$unique.ID.y 
length = pairs_poss_15_incomplete$dist 
df <- data.frame(IDx,IDy,length) 
df <- data.frame(t(apply(df, 1, sort))) 
df <- unique(df) 
IDx = df$X2 
IDy = df$X3 
length = as.numeric(paste(df$X1)) 
pairs_poss_15 <- data.frame(IDx,IDy,length)

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


# Save the observed counts to an object ----

observed_counts15 = hist(pairs_obs_15$distance, breaks=seq(0,6500,break_dist), plot = F)$counts


# Calculate frequency annulus 215 ----
df_annulus_vsw = data_frame(m = hist(pairs_poss_15$length, breaks=seq(0,6500,break_dist), plot = F)$counts, n = observed_counts15, r1 = seq(0,6500-break_dist,break_dist), r2 = seq(break_dist,6500,break_dist)) %>% 
  mutate(area = (r2^2-r1^2),
         #adj_area = (seq(10,0.02, by = -0.2)/10)*area,
         fi = (n/m)*area) %>% 
  na.omit() %>% 
  mutate(sum_fi = fi/sum(fi), 
         nm = n/m,
         nm = if_else(is.na(nm), 0, nm)) %>% 
  filter(m > 1)

# Model and predicted
vsw_distance_model <- df_annulus_vsw %>% 
  glm(nm~r1, data =., family = "binomial", weights = m)

predicted_range <- data.frame(r1 = seq(0,6500,break_dist))

predicted_data_vsw <- data.frame(y = predict(vsw_distance_model, newdata = predicted_range, type = "response"), r1 = seq(0,6500,break_dist))


# Calculate colony-specific foraging distances
vsw_dist <- sib_separation(vsw_15,2)


# BIFARIUS 2015 -----------------------------------------------------------

bsw_15 <- filter(sra_cap, species == "bifarius", year == 2015) %>% 
  rename(x=lon,y=lat)

# CALCULATING ALL POSSIBLE PAIRWIE DISTANCEs
pairs_poss_bsw15_incomplete = bsw_15 %>% 
  mutate(k=1) %>% 
  full_join(.,., by = "k") %>% 
  filter(unique.ID.x != unique.ID.y) %>%
  mutate(dist = sqrt((x.x - x.y)^2 + (y.x - y.y)^2)) %>%
  dplyr::select(-k)

IDx = pairs_poss_bsw15_incomplete$unique.ID.x # be careful with this chunk because it could be redundant with other names
IDy = pairs_poss_bsw15_incomplete$unique.ID.y 
length = pairs_poss_bsw15_incomplete$dist 
df <- data.frame(IDx,IDy,length) 
df <- data.frame(t(apply(df, 1, sort))) 
df <- unique(df) 
IDx = df$X2 
IDy = df$X3 
length = as.numeric(paste(df$X1)) 
pairs_poss_bsw15 <- data.frame(IDx,IDy,length)


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


# Save the observed counts ----

observed_countsbsw15 = hist(pairs_obs_bsw15$distance, breaks=seq(0,6500,break_dist), plot = F)$counts


# Calculate frequency annulus 215 ----
df_annulus_bsw15 = data_frame(m = hist(pairs_poss_bsw15$length, breaks=seq(0,6500,break_dist), plot = F)$counts, n = observed_countsbsw15, r1 = seq(0,6500-break_dist,break_dist), r2 = seq(break_dist,6500,break_dist)) %>% 
  mutate(area = (r2^2-r1^2),
         # no longer used
         #adj_area = (seq(10,0.02, by = -0.2)/10)*area,
         fi = (n/m)*area) %>% 
  na.omit() %>% 
  mutate(sum_fi = fi/sum(fi), 
         nm = n/m,
         nm = if_else(is.na(nm), 0, nm)) %>% 
  filter(m > 1)


# Model and predicted
bsw15_distance_model <- df_annulus_bsw15 %>% 
  glm(nm~r1, data =., family = "binomial")

predicted_range <- data.frame(r1 = seq(0,6500,break_dist))

predicted_data_bsw15 <- data.frame(y = predict(bsw15_distance_model, newdata = predicted_range, type = "response"), r1 = seq(0,6500,break_dist))


# colony-specific foraging distance (despite the _great_ function name...)
bsw15_dist <- sib_separation(bsw_15,2)



# BIFARIUS 2018 -----------------------------------------------------------

bsw_18 <- filter(sra_cap, species == "bifarius", year == 2018) %>% 
  rename(x=lon,y=lat)

# CALCULATING ALL POSSIBLE PAIRWIE DISTANCEs
pairs_poss_bsw18_incomplete = bsw_18 %>% 
  mutate(k=1) %>% 
  full_join(.,., by = "k") %>% 
  filter(unique.ID.x != unique.ID.y) %>%
  mutate(dist = sqrt((x.x - x.y)^2 + (y.x - y.y)^2)) %>%
  dplyr::select(-k)

IDx = pairs_poss_bsw18_incomplete$unique.ID.x # be careful with this chunk because it could be redundant with other names
IDy = pairs_poss_bsw18_incomplete$unique.ID.y 
length = pairs_poss_bsw18_incomplete$dist 
df <- data.frame(IDx,IDy,length) 
df <- data.frame(t(apply(df, 1, sort))) 
df <- unique(df) 
IDx = df$X2 
IDy = df$X3 
length = as.numeric(paste(df$X1)) 
pairs_poss_bsw18 <- data.frame(IDx,IDy,length)

# EMPIRICAL DISTANCEs

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


# Save the observed counts ----

observed_countsbsw18 = hist(pairs_obs_bsw18$distance, breaks=seq(0,6500,break_dist), plot = F)$counts


# Calculate annulus 2018 ----
df_annulus_bsw18 = data_frame(m = hist(pairs_poss_bsw18$length, breaks=seq(0,6500,break_dist), plot = F)$counts, n = observed_countsbsw18, r1 = seq(0,6500-break_dist,break_dist), r2 = seq(break_dist,6500,break_dist)) %>% 
  mutate(area = (r2^2-r1^2),
         #adj_area = (seq(10,0.02, by = -0.2)/10)*area,
         fi = (n/m)*area) %>% 
  na.omit() %>% 
  mutate(sum_fi = fi/sum(fi), 
         nm = n/m,
         nm = if_else(is.na(nm), 0, nm)) %>% 
  filter(m > 1)

# Model and predicted
bsw18_distance_model <- df_annulus_bsw18 %>% 
  glm(nm~r1, data =., family = "binomial", weights = m)

predicted_range <- data.frame(r1 = seq(0,6500,break_dist))

predicted_data_bsw18 <- data.frame(y = predict(bsw18_distance_model, newdata = predicted_range, type = "response"), r1 = seq(0,6500,break_dist))

# Model and predicted
bsw_combined_dist <- inner_join(select(df_annulus_bsw18, m, n, r1, r2), select(df_annulus_bsw15, m, n, r1, r2), by = c("r1", "r2")) %>% 
  mutate(area = (r2^2-r1^2),
         m = m.x + m.y,
         n = n.x + n.y) %>% 
  select(r1, r2, m, n, area) %>% 
  mutate(fi = (n/m)*area) %>% 
  na.omit() %>% 
  mutate(sum_fi = fi/sum(fi), 
         nm = n/m,
         nm = if_else(is.na(nm), 0, nm)) %>% 
  filter(m > 1)

bsw_distance_model <- bsw_combined_dist %>% 
  glm(nm~r1, data =., family = "binomial", weights = m)

predicted_range <- data.frame(r1 = seq(0,6500,break_dist))

predicted_data_bsw <- data.frame(y = predict(bsw_distance_model, newdata = predicted_range, type = "response"), r1 = seq(0,6500,break_dist))

# sivakoff Model and predicted

bsw_sivakoff_model <- bsw_combined_dist %>% 
  glm(sum_fi~r1, data =., family = "binomial", weights = m)

predicted_bsw_sivakoff <- data.frame(y = predict(bsw_sivakoff_model, newdata = predicted_range, type = "response"), r1 = seq(0,6500,break_dist))


# Colony specific foraging distance
bsw18_dist <- sib_separation(bsw_18,2)



# COMBINING FURTHER WRANGLING  ------------------------------------

bsw18_col_dist <- bsw18_dist %>% 
  group_by(FamilyID) %>% 
  summarise(col.ctr.dist = mean(ctr.dist)) %>% 
  mutate(species = "bifarius")

bsw15_col_dist <- bsw15_dist %>% 
  group_by(FamilyID) %>% 
  summarise(col.ctr.dist = mean(ctr.dist)) %>% 
  mutate(species = "bifarius")

vos_col_dist <- vsw_dist %>% 
  group_by(FamilyID) %>% 
  summarise(col.ctr.dist = mean(ctr.dist)) %>% 
  mutate(species = "vosnesenskii")

combined_col_dist <- bind_rows(bsw18_col_dist, bsw15_col_dist, vos_col_dist)

# I think this may be redundant with some improvements made but it works
pairs_obs_15$species = "vosnesenskii"
pairs_obs_bsw18$species = "bifarius"
pairs_obs_bsw15$species = "bifarius"

# PLOTS -------------------------------------------------------------------


# for adding bars to plot aesthetics
p2df1 <- data.frame(a = c(1.79, 1.79,2.21,2.21), b = c(0.28,0.29,0.29,0.28))

combined_sib_dist <- bind_rows(select(pairs_obs_15, distance, species), select(pairs_obs_bsw15, distance, species), select(pairs_obs_bsw18, distance, species))

col_distance_plot <- ggplot(combined_col_dist, aes(x = species, y = col.ctr.dist, fill = species)) +
  #geom_violin() +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(alpha = 0.7, shape = 21, color = "black") +
  scale_fill_manual(labels = c("B. bifarius", "B. vosnesenskii"), values = c(bif18_color, vos_color)) +
  #geom_hline(yintercept = 362) +
  labs(x = "Species", y = "Colony-specific foraging distance (m)", fill = "") +
  theme_classic() +
  theme(legend.position = "none", axis.text.x = element_text(face = "italic")) +
  annotate("segment", x = 1, xend = 2, y = 2000, yend = 2000, lwd = 0.3, color = "black") +
  annotate("text", x=1.5, y=2010, label="***", size=4, color = "black")


bsw_combined_annulus <- bsw_combined_dist %>% 
  mutate(species = "bifarius")

vos_annulus <- df_annulus_vsw %>% 
  mutate(species = "vosnesenskii")

combined_annulus <- bind_rows(bsw_combined_annulus, vos_annulus)





# Relative distances
combined_rel_sib_plot <- ggplot(data=combined_annulus, aes(x=r1, y=nm)) + 
  geom_point(aes(fill = species), shape = 21, color = "black") + 
  #geom_smooth(color="black", se = F, method="lm", formula = y ~ exp(x), size = 0.5) +
  #geom_smooth(method="glm", method.args=list(family=binomial(link="log"))) +
  geom_path(data = predicted_data_bsw, aes(x = r1, y = y), color = bif18_color) +
  geom_path(data = predicted_data_vsw, aes(x = r1, y = y), color = vos_color) +
  scale_fill_manual(labels = c("B. bifarius", "B. vosnesenskii"), values = c(bif18_color, vos_color)) +
  theme_classic() +
  labs(x = "", y = "\n Relative frequency") +
  theme(legend.position = c(0.8, 0.7), legend.background = element_rect(color = "black"), legend.title = element_blank(), legend.text = element_text(face = "italic"), axis.text.x = element_blank(), axis.text.y = element_text(size = 10)) +
  coord_cartesian(xlim = c(0, 5000))

# Sivakoff distance
combined_area_adjust_plot <- ggplot(data=combined_annulus, aes(x=r1, y=sum_fi)) + 
  geom_point(shape = 21, color = "black", aes(fill = species)) + 
  #geom_path(data = predicted_bsw_sivakoff, aes(x = r1, y = y), color = bif18_color) +
  scale_fill_manual(labels = c("B. bifarius", "B. vosnesenskii"), values = c(bif18_color, vos_color)) +
  theme_classic() +
  labs(x = "Pairwise sibling separation distance (m)", y = "Area-adjusted \n relative frequency") +
  theme(legend.position = "none") +
  coord_cartesian(xlim = c(0, 5000))


left_plots <- plot_grid(combined_rel_sib_plot, combined_area_adjust_plot, ncol = 1, labels = c("A", "B"))

distance_grid <- plot_grid(left_plots, col_distance_plot, rel_widths = c(2,1), scale = c(1,1), labels = c("", "C"))

# the used figure 
distance_grid

ggsave2(filename = "./figures/figure3.tiff", plot = distance_grid, width = 9, height = 5, dpi = 600)


# OUTPUT ASSOCIATED STATS -------------------------------------------------------

#sibling separation distance test
sib_sep_summary_table <- combined_sib_dist %>% 
  group_by(species) %>% 
  summarise(median = median(distance), q1 = quantile(distance, probs = 0.25), q3 = quantile(distance, probs = 0.75), n = n())

shapiro.test(combined_sib_dist$distance)
sib_sep_wilcox <- wilcox.test(data = combined_sib_dist, distance ~ species)

# colony specific foraging distance test
col_dist_summary_table <- combined_col_dist %>% 
  group_by(species) %>% 
  summarise(median = median(col.ctr.dist), q1 = quantile(col.ctr.dist, probs = 0.25), q3 = quantile(col.ctr.dist, probs = 0.75), n = n())

shapiro.test(combined_col_dist$col.ctr.dist)
col_dist_wilcox <- wilcox.test(data = combined_col_dist, col.ctr.dist ~ species)

# interaction slope
combined_annulus %>% 
  glm(data = ., nm ~ r1*species, family = "binomial", weights = m) %>% 
  summary()

#by species relative frequency
rel_freq_stats_summary <- summary(lmList(nm ~ r1 | species, data=combined_annulus, family = "binomial", weights = m))
#by species area-adjusted
area_adj_stats_summary <- summary(lmList(sum_fi ~ r1 | species, data=combined_annulus, family = "binomial"))





# DID WE UNDERSAMPLE BIFARIUS? --------------------------------------------



subsample_prop_individuals = function(to_sample, set_size) {

  df_samp = to_sample %>%
    add_row(unique.ID = 1:600, FamilyID = "s") %>%
    sample_n(size = nrow(set_size))
  #df_shuff$FamilyID = as.factor(df_shuff$FamilyID)

  df_samp$FamilyID = as.factor(df_samp$FamilyID)

  cm_ind = df_samp %>%
    filter(FamilyID !="s")

  cm_fam = cm_ind %>%
    group_by(FamilyID) %>%
    summarise(x.fam.center = mean(x), y.fam.center = mean(y), n = n()) %>%
    filter(n >= 2)

  cm_combo = inner_join(cm_ind,cm_fam, by=c("FamilyID"))

  cm_filt = cm_combo %>%
    mutate(ctr.dist = sqrt((x - x.fam.center)^2+(y - y.fam.center)^2
    )) %>%
    group_by(FamilyID) %>%
    summarise(col.ctr.dist = mean(ctr.dist))



}

#set_size = nrow(bsw_col18)
n_shuffles = 10000


# Find the randomized expectation
subsampled_vsw_distance18 <-lapply(1:n_shuffles, function(i) data.frame(iteration = i, subsample_prop_individuals(vsw_15, bsw_18))) %>%
  bind_rows()



med_subsampled_Sum_summary18 <- subsampled_vsw_distance18 %>%
  mutate(species = "vosnesenskii subsampled") %>%
  group_by(iteration, species) %>%
  summarise(med.ctr.it = median(col.ctr.dist)) %>%
  summarise(med.ctr = median(med.ctr.it),
            count = n(),
            count25 = sum(med.ctr > 25)) %>%
  summarise(prop25 = sum(count25)/n_shuffles, total25 = sum(count25), mean = mean(med.ctr), sd = sd(med.ctr), median = median(med.ctr))

subsample_summary <- med_subsampled_Sum_summary18

bsw_col_dist_combined <- bind_rows(bsw15_col_dist, bsw18_col_dist)

# sub_sample_plot18 <- subsampled_vsw_distance18 %>%
#   mutate(species = "vosnesenskii subsampled") %>%
#   group_by(iteration, species) %>%
#   summarise(med.ctr.it = median(col.ctr.dist)) %>%
#   ggplot(., aes(x = species, y = med.ctr.it)) +
#   geom_jitter(alpha = 0.1, width = 0.3) +
#   geom_hline(yintercept = c(median(bsw_col_dist_combined$col.ctr.dist),median(vos_col_dist$col.ctr.dist), pull(med_subsampled_Sum_summary18, var = mean)
#   ), linetype = c(1,1,2), color = c(bif18_color, vos_color, vos_color), size = 1.2) +
#   annotate("label", y = c(median(bsw_col_dist_combined$col.ctr.dist),median(vos_col_dist$col.ctr.dist), pull(med_subsampled_Sum_summary18, var = median)), x=1, label = c("bifarius observed", "vosnesenskii observed","vosnesenskii mean subsample"), color = "black", alpha = 0.8) +
#   labs(x = "", y = "Median Colony Foraging Distance (m)") +
#   theme_classic()
# 
# sub_sample_plot18

subsampled_histogram <- subsampled_vsw_distance18 %>%
  mutate(species = "vosnesenskii subsampled") %>%
  group_by(iteration, species) %>%
  summarise(med.ctr.it = median(col.ctr.dist)) %>%
  ggplot(., aes(x = med.ctr.it)) +
  #geom_density(fill = "grey80") +
  geom_histogram(binwidth = 10, fill = "grey50") +
  geom_vline(xintercept = c(median(bsw_col_dist_combined$col.ctr.dist),median(vos_col_dist$col.ctr.dist), pull(med_subsampled_Sum_summary18, var = median)), linetype = c(1,1,2), color = c(bif18_color, vos_color,vos_color), size = 1) +
  #annotate("label", x = c(median(bsw_col_dist_combined$col.ctr.dist),median(vos_col_dist$col.ctr.dist), pull(med_subsampled_Sum_summary18, var = median)), y=2000, label = c("bifarius observed", "vosnesenskii observed","vosnesenskii mean subsample"), color = "black", alpha = 0.8) +
  #annotate("text", x = c(median(bsw_col_dist_combined$col.ctr.dist) + 20,median(vos_col_dist$col.ctr.dist))+20, y=2000, label = c("bifarius observed", "vosnesenskii observed"), color = "black", alpha = 0.8, angle = 90) +
  labs(x = "Median Colony Foraging Distance (m)", y = "Count") +
  theme_classic()

subsampled_histogram

ggsave2("./figures/supplemental_figure.tiff", plot = subsampled_histogram, width = 8, height = 4, dpi = 600)


# SAVE OUTPUTS -------------------------------------------------------------

save(sib_sep_summary_table, sib_sep_wilcox, col_dist_summary_table, col_dist_wilcox, rel_freq_stats_summary, area_adj_stats_summary, subsample_summary, file = "./analyses_output/sibsep_foragingdist_objects.Rdata")

