#'Overexpression paper graphs SOURCE
#'This SOURCE script loads packages and defines themes
#'
#'04/12/2023
#'Build possible graphs for paper on NAC5
#'Based on Overexpression_thesis_graphs_source and Kronos_paper_graphs_source
#'
#'23/01/2023
#'Add seeds per tiller
#'
#'26/01/2024
#'Minimum 8pt font required, effectively 16pt font in this script
#'
#'Syntax rules: variables snake_case, Column Names in Caps, use pipe %>% where beneficial.
#'Built under R 4.0.4
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
setwd(paste(local_file_path, "NAC transgenics/2021-11-10 NAC5 NAP overexpression CER/Figures", sep = ""))

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#THEME
#Smaller font sizes for a paper

#16pt font required
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

#AXIS SCALES
genotype_lab <- c(null = "Null control", hom = paste("Homozygous ", gene, "A", sep  =""), multi = paste("Multi-copy ", gene, "A", sep=""))


line_name_limits <- c("5.4_null", "5.4_hom", "8.9_null", "8.9_hom", "8.23_null", "8.23_hom", "8.2_multi", "8.3_multi")
pair_limits <- c("5.4", "8.9", "8.23", "8.2", "8.3")

genotype_scale_fill <- scale_fill_manual(values=c("hom"="#fecc5cff", "null"="#ffffffff"), limits = c("null", "hom"), labels = genotype_lab)
genotype_scale_col <- scale_colour_manual(values=c("hom"="#fecc5cff", "null"="#000000ff"), limits = c("null", "hom"), labels = genotype_lab)
genotype_scale_shape <- scale_shape(limits = c("null", "hom"), labels = genotype_lab)

genotype_scale_fill_multi <- scale_fill_manual(values=c("hom"="#fecc5cff", "null"="#ffffffff", "multi" = "#e31a1cff"))
genotype_scale_col_multi <- scale_colour_manual(values=c("hom"="#fecc5cff", "null"="#000000ff", "multi" = "#e31a1cff"), limits = c("null", "hom", "multi"))
genotype_scale_shape_multi <- scale_shape(limits = c("null", "hom", "multi"))

date_scale <- scale_x_date(date_breaks = "2 days")
percent_scale <- scale_y_continuous(limits = c(0,100), breaks = seq(0,100, 20))
week_scale <- scale_x_continuous(breaks = seq(0,9,1))
spad_scale <- scale_y_continuous(breaks = c(0, 20, 40, 60), limits = c(0, 60), expand = expansion(0, 0))

#Additional scales for use with re-organised Combined_data
genotype_scale_fill_category <- scale_fill_manual(values=c("hom"="#fecc5cff", "null"="#ffffffff", "multi" = "#e31a1cff", "combo1" = "#c0c0c0ff"))
genotype_scale_col_category <- scale_colour_manual(values=c("hom"="#fecc5cff", "null"="#000000ff", "multi" = "#e31a1cff", "combo1" = "#c0c0c0ff"),
                                                limits = c("null", "hom", "multi"))
genotype_scale_shape_category <- scale_shape(limits = c("null", "hom", "multi", "combo1", "combo2"))

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#VARIABLES

#Transgenics have Marvin data but no NIR data
#Adjusted length width and area due to 2 Marvin machines being used
variable_list <- c("Leaf_senescence_DAH", "TT_SPAD30","Dur_Leaf_Sen", "AUC_Leaf_Sen",
                   "Peduncle_senescence_DAH", "Heading_date", "Tiller_1",
                   "Seed_num", "Seeds_per_tiller", "Weight", "TGW",
                   "O_length_adj", "O_width_adj", "O_area_adj",
                   "Flag_leaf_length"
)

#only needs defining if any variables have special characters escaped by backticks
variable_list_unquote <- variable_list

variable_list_truncated <- c(str_trunc(variable_list_unquote, width = 12, side = "right", ellipsis = "-"))

#It doesn't matter if the variable labels contain unused variables
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
  O_length_adj = "Average grain length (mm)",
  O_width_adj = "Average grain width (mm)",
  O_area_adj = "Average grain area (mm2)",
  `Predicted Moisture %` = "Predicted moisture (%)",
  `Predicted Protein Dry basis %` = "Predicted protein dry basis (%)",
  `Predicted Starch Dry basis %` = "Predicted starch dry basis (%)",
  `Predicted NDF Dry basis %` = "Predicted NDF dry basis (%)",
  `Predicted Dry gluten As is %` = "Predicted dry gluten (%)",
  `Predicted Hardness -` = "Predicted hardness",
  Light = "Distance from light",
  SPAD = "Flag leaf chlorophyll content (SPAD)",
  Protein_yield = "Predicted protein mass per plant (g)",
  Protein_TGW = "Protein mass per thousand grain (g)",
  Flag_leaf_length = "Flag leaf length (mm)"
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
  O_length_adj = "Grain length (mm)",
  O_width_adj = "Grain width (mm)",
  O_area_adj = "Grain area (mm2)",
  `Predicted Moisture %` = "Moisture (%)",
  `Predicted Protein Dry basis %` = "Protein DB (%)",
  `Predicted Starch Dry basis %` = "Starch DB (%)",
  `Predicted NDF Dry basis %` = "NDF DB (%)",
  `Predicted Dry gluten As is %` = "Gluten DB (%)",
  `Predicted Hardness -` = "Hardness",
  Light = "Distance from Light",
  SPAD = "Chlorophyll content",
  Protein_yield = "Protein (g/plant)",
  Protein_TGW = "Protein (mg/grain)",
  Flag_leaf_length = "Flagleaf length (mm)"
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
  Weight = list(seq(0, 40, 10)),
  Weight_dry = list(seq(0, 40, 10)),
  TGW = list(seq(0, 80, 20)),
  TGW_dry = list(seq(0, 80, 20)),
  O_length = list(seq(6.5, 8.5, 0.5)),
  O_width = list(seq(2.5, 4.5, 0.5)),
  O_area = list(seq(16, 26, 2)),
  O_length_adj = list(seq(6.5, 8.5, 0.5)),
  O_width_adj = list(seq(2.5, 4.5, 0.5)),
  O_area_adj = list(seq(16, 26, 2)),
  `Predicted Moisture %` = list(seq(14, 16, 1)),
  `Predicted Protein Dry basis %` = list(seq(10, 18, 2)),
  `Predicted Starch Dry basis %` = list(seq(65, 75, 5)),
  `Predicted NDF Dry basis %` = list(seq(0, 25, 5)),
  `Predicted Dry gluten As is %` = list(seq(8, 16, 2)),
  `Predicted Hardness -` = list(seq(40, 120, 20)),
  Protein_yield = list(seq(0, 4, 0.5)),
  Protein_TGW = list(seq(0, 15, 5)),
  Flag_leaf_length = list(seq(0, 400, 100))
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

save_panel_svg <- function(plot, plot_name, gene, n_panel_cols, n_panel_rows, ratio = (16/9)){
  ggsave(
    filename = paste(gene, "_", plot_name, "_", lubridate::today(), ".svg", sep = ""),
    plot,
    width = 180,
    height = 20 + (160/n_panel_cols) / ratio * n_panel_rows,
    units = "mm",
    device = "svg"
  )
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #