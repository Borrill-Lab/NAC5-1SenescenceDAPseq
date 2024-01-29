# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#'31.05.21
#'Draw graphs for qPCR data
#'Here, NAC-5A qPCR experiment with individual samples
#'
#'01.06.21
#'Add theme and graphs for Plant Theme seminar
#'
#'03/06/21
#'Use for NAP-A pooled samples run 25/05/21. Simplified plot section.
#'
#'21/09/21
#'Use for NAP-A individual samples run 15/09/21 and 16/09/21.
#'
#'28/10/21
#'Use data from Pfaffl 2001 method
#'
#'26/09/22
#'Specifically plot individual plants that were used for the next experiment
#'
#'29/09/22
#'Omit CTA8.9
#'
#'28/07/23
#'Make pretty graphs for thesis
#'
#'05/12/23
#'Corrections: use new data from QPCR_run_me_NAPA_2023-12-05.R and QPCR_run_me_NAC5A_2023-12-05.R
#'
#'06/12/23
#'Paper: use alternative data w/o 8.26 and paper themes for NAC5
#'
#'23/01/24
#'Omit 8.3. Highlight key points.
#'New plan, new script: run only 5.4 and 8.2 without faceting. Show expression data categorised by copy number
#'
#'26/01/24
#'Minimum font size 8pt
#'
#'Built under R 4.0.5
#'
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#'SOURCE section
#'Packages and source file QPCR_functions.R
library(readr); library(readxl)
library(tidyverse)
library(reshape2)
library(lubridate)
setwd("C:/Users/evansc/OneDrive - Norwich BioScience Institutes/NAC transgenics/NAC expression R")
source("QPCR_functions.R")

#'EDIT section
#'Choose file names and parameters
setwd("C:/Users/evansc/OneDrive - Norwich BioScience Institutes/NAC transgenics/NAC expression R")

gene = "NAC-5"
source('C:/Users/evansc/OneDrive - Norwich BioScience Institutes/NAC transgenics/NAC expression R/Overexpression_paper_graphs_source_2023_12_06.R', echo=TRUE)


#Date is the date raw data were run through QPCR_run_me.R
date = "2023-12-06"

genotype_key <- read_excel("NAPA_genotype_key.xlsx")
genotype_key2 <- read_excel("C:/Users/evansc/OneDrive - Norwich BioScience Institutes/NAC transgenics/NAC expression R/NAC-5A_genotype_key.xlsx")
genotype_key <- rbind(genotype_key, genotype_key2)

NAP_analysis<- read_csv(paste("qPCR_NAPA_individual_analysis_pfaffl_", date, ".csv", sep = "")) #Pfaffl method fold_change

NAC5_analysis<- read_csv(paste("qPCR_NAC5A_individual_wo8.26_analysis_pfaffl_", date, ".csv", sep = "")) #Pfaffl method fold_change


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#TIDY

#Filter and add Line_name and Parental Line_name data 
analysis_combined <- full_join(NAP_analysis, NAC5_analysis)

analysis_filtered_v1 <- analysis_combined %>%
  filter(Count.x >=2 & Count.y >=2) %>%
  mutate(Sample_2 = str_remove(Sample, "CTA"), Sample_3 = str_remove(Sample, "CTA")) %>% #Don't want to lose sample column
  tidyr::separate(Sample_2, into = c("Line_name", "PlantNum"), sep = "-") %>%
  left_join(genotype_key, by = "Line_name")

#Load copy number data
copy_number <- read_excel("C:/Users/evansc/OneDrive - Norwich Bioscience Institutes/NAC transgenics/2021-11-30 Taqman copy number assay/2021-12-17 Taqman copy number analysis edited.xlsx", 
                          skip = 7)

copy_number <- copy_number %>%
  dplyr::select(Sample, `2^-∆∆Ct`, `Copy number`)

analysis_filtered <- analysis_filtered_v1 %>%
  left_join(copy_number, by = c("Sample_3" = "Sample")) %>%
  mutate(`Copy number` = ifelse(Parental_genotype == "control", 0, `Copy number`))

analysis_kept <- analysis_filtered %>%
  filter(Sample %in% c("CTA5.4-9", "CTA5.4-2", "CTA8.9-7", "CTA8.9-1", "CTA8.2-3")) #"CTA8.3-3"


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#Some more plot ideas for NAC5

new_plot_labels <- c(`0.8.con1` = "8.con1 (WT)",
                     `0.5.4`  = "5.4 (WT)",
                     `1.5.4` = "5.4 (Het)",
                     `2.5.4` = "5.4 (Hom)",
                     `4.8.2` = "8.2 (Multi)")

analysis_filtered <- analysis_filtered %>%
  filter(Sample != "CTA8.con1-1") %>%
  filter(Line_name %in% c("5.4", "8.2", "8.con1")) %>%
  mutate(Kept = Sample %in% c("CTA5.4-9", "CTA5.4-2", "CTA8.2-3"),
         `Copy number` = replace_na(`Copy number`, "4"))

new_plot_CENAC57 <- analysis_filtered %>%
  filter(Primer.y == "CENAC5-7") %>%
  ggplot(aes(x = interaction(`Copy number`, Line_name), y = fold_change, fill = Kept)) +
  geom_dotplot(binaxis = "y", width = 0.2, stackdir = "center") +
  my_theme +
  scale_fill_manual(values = c("#ffffff", "blue") )+
  scale_x_discrete(limits = c("0.8.con1", "0.5.4", "1.5.4", "2.5.4", "4.8.2"), labels = new_plot_labels) +
  scale_y_continuous(breaks = seq(0, 8, 2)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), legend.position = "none") +
  labs(x = "Line name \n(Genotype)", y = str_wrap("Relative transcript level of NAC5-A1", width = 25))
new_plot_CENAC57

new_plot_CENAC53 <- analysis_filtered %>%
  filter(Primer.y == "CENAC5-3" & !is.na(`Copy number`)) %>%
  ggplot(aes(x = interaction(`Copy number`, Line_name), y = fold_change, fill = Kept)) +
  geom_dotplot(binaxis = "y", width = 0.2, stackdir = "center") +
  my_theme +
  scale_fill_manual(values = c("#ffffff", "blue") )+
  scale_x_discrete(limits = c("0.8.con1", "0.5.4", "1.5.4", "2.5.4", "4.8.2"), labels = new_plot_labels) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), legend.position = "none") +
  labs(x = "Line name \n(Genotype)", y = str_wrap("Relative transcript level of NAC5-1, all homoeologs", width = 25))
new_plot_CENAC53

new_plot_CENAC55 <- analysis_filtered %>%
  filter(Primer.y == "CENAC5-5" & !is.na(`Copy number`)) %>%
  ggplot(aes(x = interaction(`Copy number`, Line_name), y = fold_change, fill = Kept)) +
  geom_dotplot(binaxis = "y", width = 0.2, stackdir = "center") +
  my_theme +
  scale_fill_manual(values = c("#ffffff", "blue") )+
  scale_x_discrete(limits = c("0.8.con1", "0.5.4", "1.5.4", "2.5.4", "4.8.2"), labels = new_plot_labels) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), legend.position = "none") +
  labs(x = "Line name \n(Genotype)", y = str_wrap("Relative transcript level of 3*FLAG-NAC5-1 construct", width = 26))
new_plot_CENAC55

setwd("C:/Users/evansc/OneDrive - Norwich BioScience Institutes/NAC transgenics/NAC expression R")
save_small_svg(new_plot_CENAC57, "new_plot_CENAC5-7", gene, ratio = 3/4)
save_small_svg(new_plot_CENAC53, "new_plot_CENAC5-3", gene, ratio = 3/4)
save_small_svg(new_plot_CENAC55, "new_plot_CENAC5-5", gene, ratio = 3/4)

save_panel_svg(ggarrange(new_plot_CENAC55, new_plot_CENAC57, new_plot_CENAC53, ncol = 1), plot_name = "Expression_pfaffl_copynumber", gene = "NAC5",
               n_panel_cols = 3, n_panel_rows = 3, ratio = 3/4)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #