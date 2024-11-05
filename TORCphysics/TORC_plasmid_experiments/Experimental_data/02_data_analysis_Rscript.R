# Libraries ---------------------------------------------------------------
library(tidyverse)
library(rstatix)
library(ggpubr)

# Import data -------------------------------------------------------------

data <- read_tsv("M_MT_MIT_S_ST_SIT_pTORC1_49_50_63_64_65_overnight_culture_24h_measurements.tsv")

# Blank subtraction
data <- data %>% 
  group_by(date) %>% 
  mutate(A600 = A600-mean(A600[code=="blank"])) %>% 
  mutate(A700 = A700-mean(A700[code=="blank"])) %>% 
  mutate(mhYFP = mhYFP-mean(mhYFP[code=="blank"])) %>% 
  mutate(mRaspberry = mRaspberry-mean(mRaspberry[code=="blank"]))

# Fluorescence/ A600
data$mhYFP_by_A600 <- data$mhYFP/data$A600
data$mRas_by_A600 <- data$mRaspberry/data$A600

data <- filter(data, code !="blank")
write_tsv(data, "blank_subtracted_data_export.tsv")

# Leave out pTORC66 since it's not fluorescent
data <- filter(data, code != "pTORC66")

# Check number of replicates
temp <- data %>% 
  group_by(bacterium, strain, plasmid) %>% 
  count()



# Themes for plotting -----------------------------------------------------

theme_01 <- theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(),
                  axis.line = element_line(),
                  axis.text = element_text(size = 12),
                  axis.title = element_text(size = 15),
                  strip.background = element_blank(),
                  strip.text = element_text(size = 12),
                  legend.key = element_blank(),
                  legend.title = element_blank(),
                  legend.text = element_text(size = 15),
                  legend.position = "bottom")

theme_02 <- theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(),
                  axis.line = element_line(),
                  axis.text = element_text(size = 12),
                  axis.text.x = element_text(size = 9),
                  axis.title = element_text(size = 15),
                  strip.background = element_blank(),
                  strip.text = element_text(size = 10),
                  legend.key = element_blank(),
                  legend.title = element_blank(),
                  legend.text = element_text(size = 15),
                  legend.position = "bottom")


# For more colors info: https://r-graph-gallery.com/ggplot2-color.html
# Scale colors picked from https://htmlcolorcodes.com/


####### E. coli data analysis ---------------------------------------------------
data_mg1655 <- data[data$bacterium == "Escherichia coli K12 MG1655",]

data_mg1655$strain <- gsub("FRT ", "FRT\n", data_mg1655$strain)

code_order <- c("pTORC1", "pTORC49", "pTORC50", "pTORC65", "pTORC63", "pTORC64")
data_mg1655$code <- factor(data_mg1655$code, levels = code_order)

promoter_order <- c("None", 
                    "PleuWT.1min mhYFP", "PleuWT.1minRBS mhYFP", "PleuWT.1 mhYFP",
                    "PleuWT.1min mRaspberry", "PleuWT.1minRBS mRaspberry")
data_mg1655$promoter <- factor(data_mg1655$promoter, levels = promoter_order)

strain_order <- c("WT", "ΔtopA::cat", "ΔlacIZYA::FRT\nΔtopA::cat")
data_mg1655$strain <- factor(data_mg1655$strain, levels = strain_order)


# Pairwise tests ----------------------------------------------------------

data_mg1655_mhYFP_t_tests <- data_mg1655 %>% 
  group_by(bacterium, promoter) %>% 
  #tukey_hsd(mhYFP_by_A600~strain) + # maybe not very conservative
  t_test(mhYFP_by_A600~strain, p.adjust.method = "BH") %>% 
  add_y_position(fun = "max", step.increase = 0.02)
data_mg1655_mhYFP_t_tests$strain <- NA

data_mg1655_mRaspberry_t_tests <- data_mg1655 %>% 
  group_by(bacterium, promoter) %>% 
  #tukey_hsd(mRas_by_A600~strain) # maybe not very conservative
  t_test(mRas_by_A600~strain, p.adjust.method = "BH") %>% 
  add_y_position(fun = "max", step.increase = 0.02)
data_mg1655_mRaspberry_t_tests$strain <- NA


# Leaving out the promoter with RBS ---------------------------------------

data_mg1655_noRBS <- data_mg1655[-grep("RBS", data_mg1655$promoter),]

data_mg1655_noRBS_mhYFP_t_tests <- data_mg1655_noRBS %>% 
  group_by(bacterium, promoter) %>% 
  #tukey_hsd(mhYFP_by_A600~strain) # maybe not very conservative
  t_test(mhYFP_by_A600~strain, p.adjust.method = "BH") %>% 
  add_y_position(fun = "max", step.increase = 0.02)
data_mg1655_noRBS_mhYFP_t_tests$strain <- NA

data_mg1655_noRBS_mRaspberry_t_tests <- data_mg1655_noRBS %>% 
  group_by(bacterium, promoter) %>% 
  #tukey_hsd(mRas_by_A600~strain) # maybe not very conservative
  t_test(mRas_by_A600~strain, p.adjust.method = "BH") %>% 
  add_y_position(fun = "max", step.increase = 0.02)
data_mg1655_noRBS_mRaspberry_t_tests$strain <- NA

# Plotting ----------------------------------------------------------------


data_mg1655_mhYFP_plot_all_promoters <- ggplot(data_mg1655, aes(x = strain, y = mhYFP_by_A600, fill = strain)) +
  facet_grid(~bacterium+promoter) +
  stat_summary(fun = "mean", geom = "bar", show.legend = F) +
  stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1), 
               geom = "errorbar", width = 0.4) +
  geom_jitter(width = 0.2, show.legend = F) +
  xlab(NULL) +
  ylab(~F[mhYFP]/A[600]) +
  scale_fill_manual(values = c("#b7950b", "#f1c40f", "#f7dc6f")) +
  theme_02 +
  #stat_compare_means() # non-parameteric
  stat_compare_means(method = "anova", label.x = 1.5, label.y = 1.1e5) + # parametric 
  stat_pvalue_manual(data_mg1655_mhYFP_t_tests,
                     tip.length = 0.01)
ggsave("data_mg1655_mhYFP_plot_all_promoters.png", data_mg1655_mhYFP_plot_all_promoters, width = 22, height = 5, units = "in")

data_mg1655_mRas_plot_all_promoters <- ggplot(data_mg1655, aes(x = strain, y = mRas_by_A600, fill = strain)) +
  facet_grid(~bacterium+promoter) +
  stat_summary(fun = "mean", geom = "bar", show.legend = F) +
  stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1), 
               geom = "errorbar", width = 0.4) +
  geom_jitter(width = 0.2, show.legend = F) +
  xlab(NULL) +
  ylab(~F[mRaspberry]/A[600]) +
  scale_fill_manual(values = c("#76448a", "#9b59b6", "#c39bd3")) +
  theme_02 + 
  #stat_compare_means() # non-parameteric
  stat_compare_means(method = "anova", label.x = 1.5, label.y = 6e4) + # parametric
  stat_pvalue_manual(data_mg1655_mRaspberry_t_tests,
                     tip.length = 0.01)
ggsave("data_mg1655_mRas_plot_all_promoters.png", data_mg1655_mRas_plot_all_promoters, width = 20, height = 5, units = "in")


# Leaving out the promoter with RBS

data_mg1655_mhYFP_plot_no_RBS_promoters <- ggplot(data_mg1655_noRBS, aes(x = strain, y = mhYFP_by_A600, fill = strain)) +
  facet_grid(~bacterium+promoter) +
  stat_summary(fun = "mean", geom = "bar", show.legend = F) +
  stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1), 
               geom = "errorbar", width = 0.4) +
  geom_jitter(width = 0.2, show.legend = F) +
  xlab(NULL) +
  ylab(~F[mhYFP]/A[600]) +
  scale_fill_manual(values = c("#b7950b", "#f1c40f", "#f7dc6f")) +
  theme_02 +
  #stat_compare_means() # non-parameteric
  stat_compare_means(method = "anova", label.x = 1.5, label.y = 2.8e4) + # parametric 
  stat_pvalue_manual(data_mg1655_noRBS_mhYFP_t_tests,
                     tip.length = 0.01)
ggsave("data_mg1655_mhYFP_plot_no_RBS_promoters.png", data_mg1655_mhYFP_plot_no_RBS_promoters, width = 15, height = 5, units = "in")


data_mg1655_mRas_plot_no_RBS_promoters <-ggplot(data_mg1655_noRBS, aes(x = strain, y = mRas_by_A600, fill = strain)) +
  facet_grid(~bacterium+promoter) +
  stat_summary(fun = "mean", geom = "bar", show.legend = F) +
  stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1), 
               geom = "errorbar", width = 0.4) +
  geom_jitter(width = 0.2, show.legend = F) +
  xlab(NULL) +
  ylab(~F[mRaspberry]/A[600]) +
  scale_fill_manual(values = c("#76448a", "#9b59b6", "#c39bd3")) +
  theme_02 + 
  #stat_compare_means() # non-parameteric
  stat_compare_means(method = "anova", label.x = 1.5, label.y = 1e4) + # parametric
  stat_pvalue_manual(data_mg1655_noRBS_mRaspberry_t_tests,
                     tip.length = 0.01)
ggsave("data_mg1655_mRas_plot_no_RBS_promoters.png", data_mg1655_mRas_plot_no_RBS_promoters, width = 11, height = 5, units = "in")



####### Salmonella data analysis ------------------------------------------------

data_sl1344 <- data[data$bacterium == "Salmonella enterica Typhimurium SL1344",]

data_sl1344$strain <- gsub("FRT ", "FRT\n", data_sl1344$strain)
data_sl1344$strain <- gsub("lacI", "\nlacI", data_sl1344$strain)

code_order <- c("pTORC1", "pTORC49", "pTORC50", "pTORC65", "pTORC63", "pTORC64", "pTORC66")
data_sl1344$code <- factor(data_sl1344$code, levels = code_order)

promoter_order <- c("None", 
                    "PleuWT.1min mhYFP", "PleuWT.1minRBS mhYFP", "PleuWT.1 mhYFP",
                    "PleuWT.1min mRaspberry", "PleuWT.1minRBS mRaspberry",
                    "PleuWT.1 mRaspberry")
data_sl1344$promoter <- factor(data_sl1344$promoter, levels = promoter_order)

strain_order <- c("WT", "ΔSL1483::\nlacIMG1655-FRT\nΔtopA::cat", "ΔtopA::cat")
data_sl1344$strain <- factor(data_sl1344$strain, levels = strain_order)


# Pairwise tests ----------------------------------------------------------

data_sl1344_mhYFP_t_tests <- data_sl1344 %>% 
  group_by(bacterium, promoter) %>% 
  #tukey_hsd(mhYFP_by_A600~strain) + # maybe not very conservative
  t_test(mhYFP_by_A600~strain, p.adjust.method = "BH") %>% 
  add_y_position(fun = "max", step.increase = 0.02)
data_sl1344_mhYFP_t_tests$strain <- NA

data_sl1344_mRaspberry_t_tests <- data_sl1344 %>% 
  group_by(bacterium, promoter) %>% 
  #tukey_hsd(mRas_by_A600~strain) # maybe not very conservative
  t_test(mRas_by_A600~strain, p.adjust.method = "BH") %>% 
  add_y_position(fun = "max", step.increase = 0.02)
data_sl1344_mRaspberry_t_tests$strain <- NA


# Leaving out the promoter with RBS ---------------------------------------

data_sl1344_noRBS <- data_sl1344[-grep("RBS", data_sl1344$promoter),]

data_sl1344_noRBS_mhYFP_t_tests <- data_sl1344_noRBS %>% 
  group_by(bacterium, promoter) %>% 
  #tukey_hsd(mhYFP_by_A600~strain) # maybe not very conservative
  t_test(mhYFP_by_A600~strain, p.adjust.method = "BH") %>% 
  add_y_position(fun = "max", step.increase = 0.02)
data_sl1344_noRBS_mhYFP_t_tests$strain <- NA

data_sl1344_noRBS_mRaspberry_t_tests <- data_sl1344_noRBS %>% 
  group_by(bacterium, promoter) %>% 
  #tukey_hsd(mRas_by_A600~strain) # maybe not very conservative
  t_test(mRas_by_A600~strain, p.adjust.method = "BH") %>% 
  add_y_position(fun = "max", step.increase = 0.02)
data_sl1344_noRBS_mRaspberry_t_tests$strain <- NA

# Plotting ----------------------------------------------------------------


data_sl1344_mhYFP_plot_all_promoters <- ggplot(data_sl1344, aes(x = strain, y = mhYFP_by_A600, fill = strain)) +
  facet_grid(~bacterium+promoter) +
  stat_summary(fun = "mean", geom = "bar", show.legend = F) +
  stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1), 
               geom = "errorbar", width = 0.4) +
  geom_jitter(width = 0.2, show.legend = F) +
  xlab(NULL) +
  ylab(~F[mhYFP]/A[600]) +
  scale_fill_manual(values = c("#b7950b", "#f1c40f", "#f7dc6f")) +
  theme_02 +
  #stat_compare_means() # non-parameteric
  stat_compare_means(method = "anova", label.x = 1.5, label.y = 1.6e5) + # parametric 
  stat_pvalue_manual(data_sl1344_mhYFP_t_tests,
                     tip.length = 0.01)
ggsave("data_sl1344_mhYFP_plot_all_promoters.png", data_sl1344_mhYFP_plot_all_promoters, width = 22, height = 5, units = "in")

data_sl1344_mRas_plot_all_promoters <- ggplot(data_sl1344, aes(x = strain, y = mRas_by_A600, fill = strain)) +
  facet_grid(~bacterium+promoter) +
  stat_summary(fun = "mean", geom = "bar", show.legend = F) +
  stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1), 
               geom = "errorbar", width = 0.4) +
  geom_jitter(width = 0.2, show.legend = F) +
  xlab(NULL) +
  ylab(~F[mRaspberry]/A[600]) +
  scale_fill_manual(values = c("#76448a", "#9b59b6", "#c39bd3")) +
  theme_02 + 
  #stat_compare_means() # non-parameteric
  stat_compare_means(method = "anova", label.x = 1.5, label.y = 4.5e4) + # parametric
  stat_pvalue_manual(data_sl1344_mRaspberry_t_tests,
                     tip.length = 0.01)
ggsave("data_sl1344_mRas_plot_all_promoters.png", data_sl1344_mRas_plot_all_promoters, width = 20, height = 5, units = "in")


# Leaving out the promoter with RBS

data_sl1344_mhYFP_plot_no_RBS_promoters <- ggplot(data_sl1344_noRBS, aes(x = strain, y = mhYFP_by_A600, fill = strain)) +
  facet_grid(~bacterium+promoter) +
  stat_summary(fun = "mean", geom = "bar", show.legend = F) +
  stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1), 
               geom = "errorbar", width = 0.4) +
  geom_jitter(width = 0.2, show.legend = F) +
  xlab(NULL) +
  ylab(~F[mhYFP]/A[600]) +
  scale_fill_manual(values = c("#b7950b", "#f1c40f", "#f7dc6f")) +
  theme_02 +
  #stat_compare_means() # non-parameteric
  stat_compare_means(method = "anova", label.x = 1.5, label.y = 5.3e4) + # parametric 
  stat_pvalue_manual(data_sl1344_noRBS_mhYFP_t_tests,
                     tip.length = 0.01)
ggsave("data_sl1344_mhYFP_plot_no_RBS_promoters.png", data_sl1344_mhYFP_plot_no_RBS_promoters, width = 15, height = 5, units = "in")


data_sl1344_mRas_plot_no_RBS_promoters <-ggplot(data_sl1344_noRBS, aes(x = strain, y = mRas_by_A600, fill = strain)) +
  facet_grid(~bacterium+promoter) +
  stat_summary(fun = "mean", geom = "bar", show.legend = F) +
  stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1), 
               geom = "errorbar", width = 0.4) +
  geom_jitter(width = 0.2, show.legend = F) +
  xlab(NULL) +
  ylab(~F[mRaspberry]/A[600]) +
  scale_fill_manual(values = c("#76448a", "#9b59b6", "#c39bd3")) +
  theme_02 + 
  #stat_compare_means() # non-parameteric
  stat_compare_means(method = "anova", label.x = 1.5, label.y = 1.8e4) + # parametric
  stat_pvalue_manual(data_sl1344_noRBS_mRaspberry_t_tests,
                     tip.length = 0.01)
ggsave("data_sl1344_mRas_plot_no_RBS_promoters.png", data_sl1344_mRas_plot_no_RBS_promoters, width = 11, height = 5, units = "in")

