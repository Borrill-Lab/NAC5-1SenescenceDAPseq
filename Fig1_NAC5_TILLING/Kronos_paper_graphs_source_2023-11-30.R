#'NAC5 Kronos paper graphs SOURCE
#'
#'26/09/2022
#'Build graphs for Kronos data for Crop Gen Seminar on 03/10/2022
#'Based on Kronos_exploratory_2022-09-07
#'
#'06/10/2022
#'Build graphs for JIC ASM poster
#'New themes.
#'
#'09/03/2023
#'Added local_file_path
#'
#'30/11/2023
#'Build graphs for NAC5 paper. Themes only, doesn't load data.
#'
#'26/01/2024
#'Minimum of 8pt font requested. Because graphs have been built at double size, this actually equates to 16pt in the script.
#'
#'
#'This SOURCE script loads packages and data.
#'
#'Syntax rules: variables snake_case, Column Names in Caps, use pipe %>% where beneficial.
#'Built under R 4.0.5 (wait no, this computer only has 4.0.4)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#LOAD
require(car) #Anova for non-orthogonal data
require(emmeans) #marginal means estimates for when p-values are not enough
require(arm)
require(lubridate)
require(tidyverse) #Load tidyverse last to avoid masking
require(ggpubr) #for ggarrange
theme_set(theme_bw())

local_file_path <- stringr::str_extract(getwd(), ".*OneDrive\\s?.?\\s?Norwich\\s?Bio[sS]cience\\s?Institutes/")
setwd(paste(local_file_path, "NAC transgenics/2021-03 NAC phenotyping", sep = ""))
source("Functions_for_phenotype_data.R")
setwd(paste(local_file_path, "NAC transgenics/2021-11-10 NAC5 NAP Kronos greenhouse/Figures", sep = ""))
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#THEME
#Smaller font sizes for a paper
#minimum of 16pt required
axis_title_size = 16
axis_text_size = 16 

my_theme <- theme_linedraw() +
  theme(panel.grid = element_blank(),
        line = element_line(colour = "#000000"),
        axis.text = element_text(colour = "#000000"),
        strip.background = element_rect(fill = "#cccccc"),
        strip.text = element_text(colour = "#000000"),
        legend.position = "right") +
  theme(axis.title = element_text(size = axis_title_size), plot.title = element_text(size = axis_title_size), legend.title = element_text(size=axis_title_size),
        axis.text = element_text(size = axis_text_size), legend.text = element_text(size = axis_text_size), strip.text = element_text(size = axis_text_size),
        panel.grid = element_blank(), legend.position = "right")

theme_set(my_theme)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#LABELS
genotype_lab <- paste(gene, "Genotype")
leaf_senescence_lab <- "Days to 25% leaf senescence"
peduncle_senescence_lab <- "Days to 100% peduncle senescence"
SPAD_lab <- "Flag leaf chlorophyll content (SPAD)"
genotype_limits = c("Y:Y.Y:Y",
                    "X:X.Y:Y",
                    "Y:Y.X:X",
                    "X:X.X:X")
genotype_names <- c(
  `X:X.X:X` = "aabb",
  `Y:Y.X:X` = "AAbb",
  `X:X.Y:Y` = "aaBB",
  `Y:Y.Y:Y` = "AABB"
)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#Scales for ggplot2
genotype_scale_fill <- scale_fill_manual(limits = c("Y:Y.Y:Y",
                                                             "X:X.Y:Y",
                                                             "Y:Y.X:X",
                                                             "X:X.X:X"),
                                         values = c("Y:Y.Y:Y" = "#ffffffff",
                                         "X:X.Y:Y" = "#deebf7ff",
                                         "Y:Y.X:X" = "#9ecae1ff",
                                         "X:X.X:X" = "#3182bdff"), 
                                         labels = genotype_names)
genotype_scale_col <- scale_colour_manual(limits = c("Y:Y.Y:Y",
                                                     "X:X.Y:Y",
                                                     "Y:Y.X:X",
                                                     "X:X.X:X"),
                                          values = c("Y:Y.Y:Y" = "#000000ff",
                                                   "X:X.Y:Y" = "#deebf7ff",
                                                   "Y:Y.X:X" = "#9ecae1ff",
                                                   "X:X.X:X" = "#3182bdff"), 
                                        labels = genotype_names)
genotype_scale_col_2 <- scale_colour_manual(values = c("Y:Y.Y:Y" = "#000000ff",
                                                     "X:X.X:X" = "#3182bdff"), 
                                          labels = genotype_names)
genotype_scale_x <- scale_x_discrete(limits = c("Y:Y.Y:Y",
                                                "X:X.Y:Y",
                                                "Y:Y.X:X",
                                                "X:X.X:X"), labels = genotype_names)
genotype_scale_shape <- scale_shape(limits = c("X:X.X:X",
                                               "Y:Y.X:X",
                                               "X:X.Y:Y",
                                               "Y:Y.Y:Y"),
                                    labels = genotype_names)
date_scale <- scale_x_date(date_breaks = "2 days")
percent_scale <- scale_y_continuous(limits = c(0,100), breaks = seq(0,100, 20))
week_scale <- scale_x_continuous(breaks = seq(0,9,1))
subgenome_scale <- scale_colour_manual(values=c("A"="#3182bdff","B"="#62c2acff","D"="#e5f58dff"))
spad_scale <- scale_y_continuous(breaks = c(0, 20, 40, 60), limits = c(0, 60))

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

variable_list <- c("Leaf_senescence_DAH", "TT_SPAD30","Dur_Leaf_Sen", "AUC_Leaf_Sen",
                   "Peduncle_senescence_DAH", "Heading_date", "Height", "Tiller_1",
                   "Seed_num", "Seeds_per_tiller", "Weight", "Weight_dry", "TGW", "TGW_dry",
                   "O_length", "O_width", "O_area",
                   "`Predicted Moisture %`",          "`Predicted Protein Dry basis %`",
                   "`Predicted Starch Dry basis %`",  "`Predicted NDF Dry basis %`",
                   "`Predicted Dry gluten As is %`", "`Predicted Hardness -`",
                   "Protein_yield", "Protein_TGW"
)

variable_list_unquote <- c("Leaf_senescence_DAH", "TT_SPAD30","Dur_Leaf_Sen", "AUC_Leaf_Sen",
                           "Peduncle_senescence_DAH", "Heading_date", "Height", "Tiller_1",
                           "Seed_num", "Seeds_per_tiller", "Weight", "Weight_dry", "TGW", "TGW_dry",
                           "O_length", "O_width", "O_area",
                           "Predicted Moisture %",          "Predicted Protein Dry basis %",
                           "Predicted Starch Dry basis %",  "Predicted NDF Dry basis %",
                           "Predicted Dry gluten As is %", "Predicted Hardness -",
                           "Protein_yield", "Protein_TGW"
)

variable_list_truncated <- c(str_trunc(variable_list_unquote, width = 12, side = "right", ellipsis = "-"))

variable_labels <- c(
  Leaf_senescence_DAH = "Days to 25% leaf senescence",
  TT_SPAD30 = "Days to SPAD=30",
  Dur_Leaf_Sen = "Duration of leaf senescence",
  AUC_Leaf_Sen = "AUC of leaf senescence",
  Peduncle_senescence_DAH = "Days to 100% peduncle senescence",
  Heading_date = "Heading date",
  Height = "Height (mm)",
  Tiller_1 = "Main tiller number",
  Seed_num = "Grain number",
  Seeds_per_tiller = "Grain number per tiller",
  Weight = "Grain mass (g)",
  Weight_dry = "Grain mass dry basis (g)",
  TGW = "Thousand grain weight (g)",
  TGW_dry = "Thousand grain weight dry basis (g)",
  O_length = "Average grain length (mm)",
  O_width = "Average grain width (mm)",
  O_area = "Average grain area (mm2)",
  `Predicted Moisture %` = "Predicted moisture (%)",
  `Predicted Protein Dry basis %` = "Predicted protein dry basis (%)",
  `Predicted Starch Dry basis %` = "Predicted starch dry basis (%)",
  `Predicted NDF Dry basis %` = "Predicted NDF dry basis (%)",
  `Predicted Dry gluten As is %` = "Predicted dry gluten (%)",
  `Predicted Hardness -` = "Predicted hardness",
  Light = "Distance from light",
  SPAD = "Flag leaf chlorophyll content (SPAD)",
  Protein_yield = "Predicted protein mass per plant (g)",
  Protein_TGW = "Protein mass per thousand grain (g)"
)

#No more than 20 characters, to ensure it fits on one line
variable_labels_short <- c(
  Leaf_senescence_DAH = "Days to 25% leaf sen",
  TT_SPAD30 = "Days to SPAD=30",
  Dur_Leaf_Sen = "Duration leaf sen",
  AUC_Leaf_Sen = "AUC leaf sen",
  Peduncle_senescence_DAH = "Days to peduncle sen",
  Heading_date = "Heading date",
  Height = "Height (mm)",
  Tiller_1 = "Main tiller num",
  Seed_num = "Grain num",
  Seeds_per_tiller = "Grain num / tiller",
  Weight = "Grain mass (g)",
  Weight_dry = "Grain mass DB (g)",
  TGW = "TKW (g)",
  TGW_dry = "TKW DB (g)",
  O_length = "Grain length (mm)",
  O_width = "Grain width (mm)",
  O_area = "Grain area (mm2)",
  `Predicted Moisture %` = "Moisture (%)",
  `Predicted Protein Dry basis %` = "Protein DB (%)",
  `Predicted Starch Dry basis %` = "Starch DB (%)",
  `Predicted NDF Dry basis %` = "NDF DB (%)",
  `Predicted Dry gluten As is %` = "Gluten DB (%)",
  `Predicted Hardness -` = "Hardness",
  Light = "Distance from Light",
  SPAD = "Chlorophyll content",
  Protein_yield = "Protein (g/plant)",
  Protein_TGW = "Protein (mg/grain)"
)

#Prepare manual breaks for the "panel sizes the same" fudge
#These are prepared based on the range of the Kronos NAC5 data, and to ensure break labels are
#'nice' numbers where possible
#'start at 0 where possible
#'consistent between similar variables
#If manual breaks are used, the SPACE between the last data point and the edge of the plot will be larger and not be uniform
manual_breaks <- c(
  Leaf_senescence_DAH = list(seq(0, 80, 20)),
  TT_SPAD30 = list(seq(0, 80, 20)),
  Dur_Leaf_Sen = list(seq(0, 40, 10)),
  AUC_Leaf_Sen = list(seq(0, 4000, 1000)),
  Peduncle_senescence_DAH = list(seq(20, 80, 10)),
  Heading_date = list(seq(50, 70, 5)),
  Height = list(seq(450, 850, 100)),
  Tiller_1 = list(seq(0, 15, 5)),
  Seed_num = list(seq(0, 600, 100)),
  Seeds_per_tiller = list(seq(0, 60, 10)),
  Weight = list(seq(0, 40, 2)),
  Weight_dry = list(seq(0, 40, 2)),
  TGW = list(seq(0, 80, 20)),
  TGW_dry = list(seq(0, 80, 20)),
  O_length = list(seq(6.5, 8.5, 0.5)),
  O_width = list(seq(2.5, 4.5, 0.5)),
  O_area = list(seq(16, 26, 2)),
  `Predicted Moisture %` = list(seq(14, 16, 1)),
  `Predicted Protein Dry basis %` = list(seq(10, 18, 2)),
  `Predicted Starch Dry basis %` = list(seq(65, 75, 5)),
  `Predicted NDF Dry basis %` = list(seq(0, 25, 5)),
  `Predicted Dry gluten As is %` = list(seq(8, 16, 2)),
  `Predicted Hardness -` = list(seq(40, 120, 20)),
  Protein_yield = list(seq(0, 4, 0.5)),
  Protein_TGW = list(seq(0, 15, 5))
)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#Automate standard plot options

#ratio 4:3
save_small_svg <- function(plot, plot_name, gene, ratio = (4/3)){
  ggsave(
    filename = paste(gene, "_", plot_name, "_paper_", lubridate::today(), ".svg", sep = ""),
    plot,
    width = 80,
    height = 80 / ratio,
    units = "mm",
    device = "svg"
  )
}

#ratio 16:9
save_big_svg <- function(plot, plot_name, gene, ratio = (16/9)){
  ggsave(
    filename = paste(gene, "_", plot_name, "_paper_", lubridate::today(), ".svg", sep = ""),
    plot,
    width = 160,
    height = 160 / ratio,
    units = "mm",
    device = "svg"
  )
}

save_pdf_1 <- function(plot, plot_name, gene, ratio = (16/9)){
  ggsave(
    filename = paste(gene, "_", plot_name, "_paper_", lubridate::today(), ".pdf", sep = ""),
    plot,
    width = 160,
    height = 160 / ratio,
    units = "mm",
    device = "pdf"
  )
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #