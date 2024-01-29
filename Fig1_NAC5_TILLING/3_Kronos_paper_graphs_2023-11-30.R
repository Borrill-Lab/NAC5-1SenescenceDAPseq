#'Kronos paper graphs
#'WARNING this script creates a very large number of boxplots
#'
#'26/09/2022
#'Build graphs for Kronos data for Crop Gen Seminar on 03/10/2022
#'Based on Kronos_exploratory_2022-09-07
#'
#'06/10/2022
#'Adapt graphs for ASM poster
#'new sizes
#'
#'26/10/2022
#'Adapt graphs for AAB seminar
#'Best of both Crop Gen seminar and ASM poster
#'
#'30/11/2023
#'Adapt graphs for NAC5 paper.
#'Load analysis scripts from 30/11/2023 (after tidying source and data_preparation scripts)
#'Remove unnecessary script
#'updated boxplotting to avoid copy-paste based on NAM2_figures_for_thesis_2023-07-31
#'updated aes() and vars() following vignette("ggplot2-in-packages")
#'
#'19/01/2024
#'SPAD plot with single and double mutant lines together
#'
#'26/01/2024
#'Minimum of 8pt font requested. Set x-axis at 90 degrees and increase ratio so it fits
#'
#'Syntax rules: variables snake_case, Column Names in Caps, use pipe %>% where beneficial.
#'Built under R 4.0.5 (wait no, this computer only has 4.0.4)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#EDIT GENE AND DATE TO SUIT OCCASION
#date is the date that Kronos_data_preparation.R was run to generate analysis files
gene = "NAC5"
date = "2023-11-30"
print(paste("Gene:", gene, "Date:", date))
#

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#LOAD SOURCE SCRIPTS AND DATA FILES
local_file_path <- stringr::str_extract(getwd(), ".*OneDrive\\s?.?\\s?Norwich\\s?Bio[sS]cience\\s?Institutes/")

#THEMES, LABELS, SCALES, COLOUR SCHEMES, VARIABLE LISTS
source(paste(local_file_path, "NAC transgenics/2021-11-10 NAC5 NAP Kronos greenhouse/Kronos_paper_graphs_source_2023-11-30.R", sep = ""), echo = TRUE)

#DATA FILES
setwd(paste(local_file_path, "NAC transgenics/2021-11-10 NAC5 NAP Kronos greenhouse/2022-03-07_clean_data", sep = ""))

#Read in analysis files; exclude heterozygous genotypes
Kronos_SPAD_analysis <- read_csv(paste("Kronos_SPAD_yday_", gene, "_", date, ".csv", sep= "")) %>%
  filter(Genotype %in% c("X:X.X:X","Y:Y.X:X","X:X.Y:Y","Y:Y.Y:Y")) %>%
  mutate(Block = as.factor(Block), Genotype = as.factor(Genotype))
Kronos_analysis <- read_csv(paste("Kronos_phenotyping_yday_", gene, "_", date, ".csv", sep= "")) %>%
  filter(Genotype %in% c("X:X.X:X","Y:Y.X:X","X:X.Y:Y","Y:Y.Y:Y")) %>%
  mutate(Block = as.factor(Block), Genotype = as.factor(Genotype))

setwd(paste(local_file_path, "NAC transgenics/2021-11-10 NAC5 NAP Kronos greenhouse/Figures", sep = ""))

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#SPAD PLOTS

setwd(paste(local_file_path, "NAC transgenics/2021-11-10 NAC5 NAP Kronos greenhouse/Figures", sep = ""))

# #Example of testing individual timepoints with a Wilcoxon non-parametric test
# compare <- ggpubr::compare_means(SPAD ~ Genotype,
#                       data = filter(Kronos_SPAD_analysis, Genotype %in% c("X:X.Y:Y","Y:Y.Y:Y")),
#                       group.by = "Week",
#                       method = "wilcox.test",
#                       p.adjust.method = "fdr",
#                       paired = FALSE
# )
# View(compare)

#Plot SPAD plots annotated with Wilcoxon test (NB significance stars use non-adjusted p-values)
#This section should be turned into a function really
Averages_plot_wilcox <- ggline(data = filter(Kronos_SPAD_analysis,  Genotype %in% c("X:X.X:X","Y:Y.Y:Y")), 
       x = "Day_approx", y = "SPAD", add = "mean_se",
       color = "Genotype", shape = "Genotype") +
  stat_compare_means(aes(group = Genotype), label = "p.signif",
                     label.y = 58, method = "wilcox.test", size = axis_text_size / .pt) +
  genotype_scale_col +
  genotype_scale_shape +
  spad_scale +
  labs(y = str_wrap(variable_labels["SPAD"], width = 25), x = "Days after Heading") +
  my_theme +
  theme(legend.position = "none")
Averages_plot_wilcox

#SAVE PDFs as svg files are problematic for this plot
save_pdf_1(Averages_plot_wilcox, "Averages_plot_wilcox_aabb", gene = gene)
save_big_svg(Averages_plot_wilcox, "Averages_plot_wilcox_aabb", gene = gene)

#Also make this plot for the single mutants
Averages_plot_wilcox <- ggline(data = filter(Kronos_SPAD_analysis,  Genotype %in% c("Y:Y.X:X","Y:Y.Y:Y")), 
                               x = "Day_approx", y = "SPAD", add = "mean_se",
                               color = "Genotype", shape = "Genotype", size = 0.3, point.size = 0.15) +
  stat_compare_means(aes(group = Genotype), label = "p.signif",
                     label.y = 58, method = "wilcox.test", size = axis_text_size / .pt) +
  genotype_scale_col +
  genotype_scale_shape +
  spad_scale +
  labs(y = str_wrap(variable_labels["SPAD"], width = 25), x = "Days after Heading") +
  my_theme +
  theme(legend.position = "none")
Averages_plot_wilcox

save_pdf_1(Averages_plot_wilcox, "Averages_plot_wilcox_bb", gene = gene)
save_big_svg(Averages_plot_wilcox, "Averages_plot_wilcox_bb", gene = gene)

Averages_plot_wilcox <- ggline(data = filter(Kronos_SPAD_analysis,  Genotype %in% c("X:X.Y:Y","Y:Y.Y:Y")), 
                               x = "Day_approx", y = "SPAD", add = "mean_se",
                               color = "Genotype", shape = "Genotype", size = 0.3, point.size = 0.15) +
  stat_compare_means(aes(group = Genotype), label = "p.signif",
                     label.y = 58, method = "wilcox.test", size = axis_text_size / .pt) +
  genotype_scale_col +
  genotype_scale_shape +
  spad_scale +
  labs(y = str_wrap(variable_labels["SPAD"], width = 25), x = "Days after Heading") +
  my_theme +
  theme(legend.position = "none")
Averages_plot_wilcox

save_pdf_1(Averages_plot_wilcox, "Averages_plot_wilcox_aa", gene = gene)
save_big_svg(Averages_plot_wilcox, "Averages_plot_wilcox_aa", gene = gene)


#Solely to generate a COLOUR LEGEND
Averages_plot_wilcox <- ggline(data = Kronos_SPAD_analysis,
                               x = "Day_approx", y = "SPAD", add = "mean_se",
                               color = "Genotype", shape = "Genotype", group = "Genotype") +
  stat_compare_means(aes(group = Genotype), label = "p.signif",
                     label.y = 58, method = "wilcox.test", size = axis_text_size / .pt) +
  genotype_scale_col +
  genotype_scale_shape +
  spad_scale +
  labs(y = str_wrap(variable_labels["SPAD"], width = 25), x = "Days after Heading") +
  my_theme +
  theme(legend.position = "right")
Averages_plot_wilcox

save_pdf_1(Averages_plot_wilcox, "Averages_plot_wilcox_legend", gene = gene)
save_big_svg(Averages_plot_wilcox, "Averages_plot_wilcox_legend", gene = gene)

#Alternative version with all lines but no stats
#re format
#Use standard ggplot approach (can't do Wilcox tests but this is fine)
Averages_plot <- ggplot(Kronos_SPAD_analysis, aes(x=Week, y=SPAD, group = Genotype, fill=Genotype, color = Genotype, shape = Genotype))
Averages_plot <- Averages_plot +
  stat_summary(fun.data = "mean_se", geom=("errorbar"), width=0.1, lwd = 0.8) +
  stat_summary(fun.data = "mean_se", geom="line", lwd = 0.8) +
  stat_summary(fun.data = "mean_se", geom="point", size=2, color = "black") +
  scale_x_continuous(breaks = seq(0, 9, 1), labels = ~(.x*7)) +
  genotype_scale_fill +
  genotype_scale_col +
  scale_shape_manual(limits = c("Y:Y.Y:Y",
                                "X:X.Y:Y",
                                "Y:Y.X:X",
                                "X:X.X:X"),
                     labels = genotype_names,
                     values=c(21, 22, 23, 24)) +
  spad_scale +
  labs(y = str_wrap(variable_labels["SPAD"], width = 25), x = "Days after Heading") +
  my_theme + theme(legend.position = "none")

Averages_plot
opt = "opt3"
save_pdf_1(Averages_plot, paste("Averages_plot", opt, sep = "_"), gene = gene)
save_big_svg(Averages_plot, paste("Averages_plot", opt, sep = "_"), gene = gene)



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#BOXPLOTS

setwd(paste(local_file_path, "NAC transgenics/2021-11-10 NAC5 NAP Kronos greenhouse/Figures", sep = ""))

#MODELS

#Collect plots
plot_list <- list()
list_index = 1

#Set model parameters
model_formula <- c(" ~ Light + Block + Genotype", " ~ Light + Block + Cross*Genotype")
n_variables <- as.integer(c(1,2))
explanatory <- c("Genotype", "Cross*Genotype")
facet_by <- c(NA, "Cross")

#Iterate through variables
#This was not made into a function because it relies on many variables/parameters 
for(mdl in 1:length(model_formula)){
  
  #Plot everything at once
  for(i in 1:length(variable_list)){
    #quick sanity check
    if(is.integer(n_variables[mdl]) &
       is.character(model_formula[mdl]) &
       is.character(explanatory[mdl])){
      print("#### Input parameters exist, generating boxplots ####")
    }else{
      next()
    }
    print(variable_list[i])
    
    #Prepare linear model
    linear_model <- lm(formula(paste(variable_list[i], model_formula[mdl], sep = "")), data = Kronos_analysis)
    
    print(paste(variable_list[i], model_formula[mdl], sep = ""))
    # #Check linear model assumptions
    # par(mfrow=c(1,2))
    # plot(linear_model, which = c(1,2))
    # par(mfrow=c(1,1))
    
    #ANOVA and "compact letter display" labels
    labels <- plot_anova_and_return_labels_emm(linear_model, explanatory = explanatory[mdl], alpha = 0.05, variables = n_variables[mdl])
    #axes
    ymin <- min(Kronos_analysis[variable_list_unquote[i]], na.rm = TRUE)
    ymax <- max(Kronos_analysis[variable_list_unquote[i]], na.rm = TRUE)
    label_y <- ymax + 0.12*(ymax - ymin)
    axis_ymin <- ymin
    axis_ymax <- ymax + 0.18*(ymax - ymin)
    
    #final boxplot
    #str_wrap reduced to 20 because plots are smaller relative to text
    plot <- ggplot(Kronos_analysis,
                   aes(x = Genotype, y = .data[[variable_list_unquote[i]]], fill = Genotype)
    ) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(position = position_jitter(height = 0)) +
      genotype_scale_x +
      genotype_scale_fill +
      labs(y = str_wrap(variable_labels[variable_list_unquote[i]], width = 20), x = NULL, fill = genotype_lab) +
      geom_label(
        data = labels,
        x = labels[, (n_variables[mdl])],
        label = labels[, (1 + n_variables[mdl])],
        y = label_y,
        fill = "white",
        colour ="black", #specify colours or it will inherit from aes() and throw an error
        size = axis_text_size / .pt, #axis_text_size defined in paper_graphs_source
        label.size = 0
      ) +
      theme(legend.position = "none")
    
    #optional faceting
    if(!is.na(facet_by[mdl])){
      plot1 <- plot + facet_wrap(vars(.data[[facet_by[mdl]]]), ncol = 2) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
        coord_cartesian(ylim = c(axis_ymin, axis_ymax))
      print(plot1)
      save_big_svg(plot1, str_c("Boxplot", facet_by[mdl], variable_list_truncated[i], sep = "_"), gene, ratio = (4/3))
    }else{
      plot1 <- plot + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
        coord_cartesian(ylim = c(axis_ymin, axis_ymax))
      print(plot1)
      save_small_svg(plot1, str_c("Boxplot", variable_list_truncated[i], sep = "_"), gene, ratio = (1))
    }
    
    #plot_short demonstrates alternative formatting options
    #this should result in uniform panel sizes (due to uniform y-axis depth)
    #and works with a narrower aspect ratio
    breaks <- as.vector(unlist(manual_breaks[variable_list_unquote[i]]))
    plot_short <- plot +
      labs(y = str_wrap(variable_labels_short[variable_list_unquote[i]], width = 20), x = NULL, fill = genotype_lab) +
      scale_y_continuous(breaks = breaks, limits = c(min(breaks), max(breaks))) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
    
    if(!is.na(facet_by[mdl])){
      plot_short <- plot_short + facet_wrap(vars(.data[[facet_by[mdl]]]), ncol = 2)
      print(plot_short)
      save_big_svg(plot_short, str_c("Boxplot", "short", facet_by[mdl], variable_list_truncated[i], sep = "_"), gene, ratio = (4/3))
    }else{
      print(plot_short)
      save_small_svg(plot_short, str_c("Boxplot", "short", variable_list_truncated[i], sep = "_"), gene, ratio = 1)
    }
    
    #collect plots into a list (required for combining plots in R with ggarrange)
    plot_list[[list_index]] <- plot
    list_index <- list_index + 1
  }
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
