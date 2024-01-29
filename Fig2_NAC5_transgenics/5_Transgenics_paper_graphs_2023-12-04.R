#'Overexpression paper graphs
#'
#'04/12/2023
#'Build possible graphs for paper on NAC5
#'
#'05/12/2023
#'Remove 8.23 because a couple plants are missing and it doesn't come with expression data
#'
#'22/01/2023
#'Remove 8.3, run Averages_plot_wilcox with Combined_data
#'
#'23/01/2023
#'Just show 5.4 and 8.2, without facets, similar to Kronos_paper_graphs
#'
#'26/01/2024
#'Minimum 8pt font required
#'
#'Syntax rules: variables snake_case, Column Names in Caps, use pipe %>% where beneficial.
#'Built under R 4.0.4
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#EDIT GENE AND DATE TO SUIT OCCASION
#date is the date the scripts were run through Overexpression_data_preparation.R
gene = "NAC5"
#gene = "NAP"
date = "2024-01-23"

source('C:/Users/evansc/OneDrive - Norwich BioScience Institutes/NAC transgenics/2021-11-10 NAC5 NAP overexpression CER/2022-08-12_scripts/Transgenics_paper_graphs_source_2023-12-04.R', echo=TRUE)

setwd(paste(local_file_path, "NAC transgenics/2021-11-10 NAC5 NAP Overexpression CER/2022-08-12_clean_data", sep = ""))

Overexpression_SPAD_analysis <- read_csv(paste("Overexpression_SPAD_yday_", gene, "_", date, ".csv", sep= "")) %>%
  mutate(Block = as.factor(Block))
Overexpression_analysis <- read_csv(paste("Overexpression_phenotyping_yday_", gene, "_", date, ".csv", sep= "")) %>%
  mutate(Block = as.factor(Block))

#Combined_data has some custom columns and filtering for the paper
#Warning: col_types specification must be updated if the number of data columns is changed
#Key purpose of this is to assign some variables as factors
Combined_data <- read_csv(paste("Combined_control_data_", gene, "_", date, ".csv", sep= ""),
                          col_types = "cfcffdddddddddddddddddddddddddddcldclcldddddddddccddddddddddTcdddfff") %>%
  filter(Line_prefix %in% c("5.4", "8.2"))

Combined_SPAD_data <- read_csv(paste("Combined_control_SPAD_data_", gene, "_", date, ".csv", sep= "")) %>%
  filter(Line_prefix %in% c("5.4", "8.2")) %>%
  mutate(Genotype = as.factor(Genotype))

setwd(paste(local_file_path, "NAC transgenics/2021-11-10 NAC5 NAP Overexpression CER/Figures", sep = ""))

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#SPAD PLOT
setwd(paste(local_file_path, "NAC transgenics/2021-11-10 NAC5 NAP Overexpression CER/Figures", sep = ""))

#This plot isn't compatible with facet so I'll plot individually
make_wilcox_plot <- function(data){
  Averages_plot_wilcox <- ggline(data = data, 
                                 x = "Day_approx", y = "SPAD", add = "mean_se",
                                 color = "Genotype", shape = "Genotype"
                                 ) +
    stat_compare_means(aes(group = Genotype), label = "p.signif",
                       label.y = 57, method = "wilcox.test", size = 16 / .pt, hide.ns = TRUE) +
    spad_scale +
    genotype_scale_fill_category +
    genotype_scale_col_category +
    genotype_scale_shape_multi +
    labs(y = str_wrap(variable_labels["SPAD"], width = 25), x = "Days after Heading") +
    my_theme +
    theme(legend.position = "none") +
    coord_cartesian(ylim = c(0, 65))
  return(Averages_plot_wilcox)
}

#SAVE SAVE SAVE

plots <- list()

pair <- "5.4"
Averages_plot_wilcox <- make_wilcox_plot(filter(Combined_SPAD_data, Pair == pair)) +
  ggtitle(pair)
save_small_svg(Averages_plot_wilcox, plot_name = paste("SPAD_wilcox", pair, sep = "_"), gene = gene, ratio = (4/3))
plots[[1]] <- Averages_plot_wilcox

pair <- "8.9"
Averages_plot_wilcox <- make_wilcox_plot(filter(Combined_SPAD_data, Pair == pair)) +
  ggtitle(pair)
save_small_svg(Averages_plot_wilcox, plot_name = paste("SPAD_wilcox", pair, sep = "_"), gene = gene, ratio = (4/3))
plots[[2]] <- Averages_plot_wilcox

pair <- "8.2"
Averages_plot_wilcox <- make_wilcox_plot(filter(Combined_SPAD_data, Pair == pair)) +
  ggtitle(pair)
save_small_svg(Averages_plot_wilcox, plot_name = paste("SPAD_wilcox", pair, sep = "_"), gene = gene, ratio = (4/3))
plots[[3]] <- Averages_plot_wilcox

pair <- "8.2"
Averages_plot_wilcox <- make_wilcox_plot(filter(Combined_SPAD_data, Line_name %in% c("5.4_null", "8.2_multi"))) +
  ggtitle(pair)
# save_small_svg(Averages_plot_wilcox, plot_name = paste("SPAD_wilcox", pair, sep = "_"), gene = gene, ratio = (4/3))
plots[[3]] <- Averages_plot_wilcox

ggarrange(plots[[1]], plots[[3]], plots[[2]])

#SAVE PANELS
save_panel_svg(ggarrange(plots[[1]], plots[[3]], plots[[2]]), plot_name = "SPAD_wilcox_all", gene = gene,
               n_panel_cols = 2, n_panel_rows = 2, ratio = (4/3))

#What happens if I plot them all at once
#This one has custom scaling
#sortof works except legend doesn't show line type, would have to set manually

#this version relies on only one line name per genotype (careful)
SPAD_plot <- ggline(data = filter(Overexpression_SPAD_analysis,  Line_name %in% c("5.4_null", "5.4_hom", "8.2_multi") & Shelf == "Shelf 2"), 
       x = "Day_approx", y = "SPAD", add = "mean_se",
       color = "Genotype", group = "Genotype", shape = "Genotype") +
  spad_scale +
  genotype_scale_col_multi +
  genotype_scale_shape_multi +
  labs(y = str_wrap(variable_labels["SPAD"], width = 25), x = "Days after Heading") +
  my_theme +
  theme(legend.position = "none") +
  coord_cartesian(ylim = c(0, 65))
SPAD_plot

# save_small_svg(SPAD_plot, plot_name = paste("SPAD_plot", "small", sep = "_"), gene = gene, ratio = (16/9))
save_big_svg(SPAD_plot, plot_name = paste("SPAD_plot", "big", sep = "_"), gene = gene, ratio = (16/9))

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# check straight wilcox tests

#this code makes a table comparing all line names to a reference group for all variables using wilcoxon test
Combined_data <- Combined_data %>%
  filter(Category %in% c("hom", "null", "multi")) %>%
  mutate(Seeds_per_tiller = Seed_num/Tiller_1) 

argument <- paste("c(",  paste(variable_list, collapse = ", "), ")", " ~ Line_name", sep = "")

pvalue_table <- ggpubr::compare_means(formula(argument), data = Combined_data,
                      method = "wilcox.test", ref.group = "5.4_null")

write_csv(pvalue_table, paste("Transgenics_wilcox_pvalues_", today(), ".csv", sep = ""))

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#BOXPLOTS
#CLASSIC boxplots using Overexpression_analysis

setwd(paste(local_file_path, "NAC transgenics/2021-11-10 NAC5 NAP Overexpression CER/Figures", sep = ""))

#MODELS

#Collect plots
plot_list <- list()
list_index = 1

#Set model parameters. Can be a vector; ensure all vector are the same length
model_formula <- c(" ~ Shelf * Line_name")
n_variables <- as.integer(c(2))
explanatory <- c("Shelf*Line_name")
facet_by <- c("Shelf")

#Iterate through variables
#This was not made into a function because it relies on many variables/parameters 
for(mdl in 1:length(model_formula)){
  
  #Plot everything at once
  for(i in 1:length(variable_list)){
    print(variable_list[i])
    #quick sanity check
    if(is.integer(n_variables[mdl]) &
       is.character(model_formula[mdl]) &
       is.character(explanatory[mdl])){
      print("#### Input parameters exist, generating boxplots ####")
    }else{
      print(paste("Error: Input parameter(s) are incorrect class for model", mdl, "- skipping"))
      next()
    }
    
    #Prepare linear model
    linear_model <- lm(formula(paste(variable_list[i], model_formula[mdl], sep = "")), data = Overexpression_analysis)
    
    print(paste(variable_list[i], model_formula[mdl], sep = ""))
    #Check linear model assumptions
    par(mfrow=c(1,2))
    plot(linear_model, which = c(1,2))
    par(mfrow=c(1,1))
    
    #ANOVA and "compact letter display" labels
    labels <- plot_anova_and_return_labels_emm(linear_model, explanatory = explanatory[mdl], alpha = 0.05, variables = n_variables[mdl])
    #axes
    ymin <- min(Overexpression_analysis[variable_list_unquote[i]], na.rm = TRUE)
    ymax <- max(Overexpression_analysis[variable_list_unquote[i]], na.rm = TRUE)
    label_y <- ymax + 0.12*(ymax - ymin)
    axis_ymin <- ymin
    axis_ymax <- ymax + 0.18*(ymax - ymin)
    
    #final boxplot
    #str_wrap reduced to 20 because plots are smaller relative to text
    plot <- ggplot(Overexpression_analysis,
                   aes(x = Line_name, y = .data[[variable_list_unquote[i]]], fill = Genotype)
    ) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(position = position_jitter(height = 0)) +
      scale_x_discrete(limits = line_name_limits) +
      genotype_scale_fill_multi +
      labs(y = str_wrap(variable_labels[variable_list_unquote[i]], width = 20), x = NULL, fill = Line_name_lab) +
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
        coord_cartesian(ylim = c(axis_ymin, axis_ymax))
      print(plot1)
      save_big_svg(plot1, str_c("Boxplot", "ox1", facet_by[mdl], variable_list_truncated[i], sep = "_"), gene, ratio = (16/9))
    }else{
      plot1 <- plot + 
        coord_cartesian(ylim = c(axis_ymin, axis_ymax))
      print(plot1)
      save_small_svg(plot1, str_c("Boxplot", "ox1", variable_list_truncated[i], sep = "_"), gene, ratio = (4/3))
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
      save_big_svg(plot_short, str_c("Boxplot", "ox1", "short", facet_by[mdl], variable_list_truncated[i], sep = "_"), gene, ratio = (4/3))
    }else{
      print(plot_short)
      save_small_svg(plot_short, str_c("Boxplot", "ox1", "short", variable_list_truncated[i], sep = "_"), gene, ratio = 1)
    }
    
    #collect plots into a list (required for combining plots in R with ggarrange)
    plot_list[[list_index]] <- plot
    list_index <- list_index + 1
  }
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#BOXPLOTS
#DIFFERENT boxplots using Combined_data
#And paired Wilcoxon tests

setwd(paste(local_file_path, "NAC transgenics/2021-11-10 NAC5 NAP Overexpression CER/Figures", sep = ""))

#Collect plots
plot_list <- list()
list_index = 1

#Set model parameters. Can be a vector; ensure all vector are the same length
model_formula <- c("")
n_variables <- as.integer(c(2))
explanatory <- c("")
facet_by <- c(NA)

custom_limits <- c("5.4_null", "5.4_hom", "8.2_multi")
custom_labels <- c(`5.4_null` = "5.4\n(WT)",
                   `5.4_hom` = "5.4\n(Hom)",
                   `8.2_multi` = "8.2\n(Multi)")

#remove unnecessary stuff again
Combined_data <- Combined_data %>%
  filter(Category %in% c("hom", "null", "multi")) %>%
  mutate(Seeds_per_tiller = Seed_num/Tiller_1) #that one extra variable

#Iterate through variables
#This was not made into a function because it relies on many variables/parameters

for(mdl in 1:length(model_formula)){
  
  #Plot everything at once
  for(i in 1:length(variable_list)){
    print(variable_list[i])
    #quick sanity check
    if(is.integer(n_variables[mdl]) &
       is.character(model_formula[mdl]) &
       is.character(explanatory[mdl])){
      print("#### Input parameters exist, generating boxplots ####")
    }else{
      print(paste("Error: Input parameter(s) are incorrect class for model", mdl, "- skipping"))
      next()
    }
    
    #axes
    ymin <- min(Combined_data[variable_list_unquote[i]], na.rm = TRUE)
    ymax <- max(Combined_data[variable_list_unquote[i]], na.rm = TRUE)
    label_y <- ymax + 0.12*(ymax - ymin)
    axis_ymin <- ymin
    axis_ymax <- ymax + 0.18*(ymax - ymin)
    
    #final boxplot
    #str_wrap reduced to 20 because plots are smaller relative to text
    plot <- ggplot(Combined_data,
                   aes(x = Line_name, y = .data[[variable_list_unquote[i]]], fill = Genotype)
    ) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(position = position_jitter(height = 0)) +
      genotype_scale_fill_multi +
      labs(y = str_wrap(variable_labels[variable_list_unquote[i]], width = 20), x = NULL) + 
      scale_x_discrete(limits = custom_limits, labels = custom_labels) +
      # stat_compare_means(aes(group = Binary), label = "p.signif",
      #                    label.y = label_y, method = "wilcox.test", size = axis_text_size / .pt, label.x.npc = "centre") +
      theme(legend.position = "none")
    
    #optional faceting
    if(!is.na(facet_by[mdl])){
      plot1 <- plot + facet_wrap(vars(factor((.data[[facet_by[mdl]]]), levels = pair_limits)), ncol = 5) +
        coord_cartesian(ylim = c(axis_ymin, axis_ymax))
      print(plot1)
      save_big_svg(plot1, str_c("Boxplot", "ox2", facet_by[mdl], variable_list_truncated[i], sep = "_"), gene, ratio = (16/9))
    }else{
      plot1 <- plot + 
        coord_cartesian(ylim = c(axis_ymin, axis_ymax))
      print(plot1)
      save_small_svg(plot1, str_c("Boxplot", "ox2", variable_list_truncated[i], sep = "_"), gene, ratio = (4/3))
    }
    
    #plot_short tries alternative formatting / saving options
    #this should result in uniform panel sizes (due to uniform y-axis depth)
    #and works with a narrower aspect ratio
    breaks <- as.vector(unlist(manual_breaks[variable_list_unquote[i]]))
    plot_short <- plot +
      #labs(y = str_wrap(variable_labels_short[variable_list_unquote[i]], width = 20), x = NULL, fill = genotype_lab) +
      coord_cartesian(ylim = c(axis_ymin, axis_ymax)) #+
      #scale_y_continuous(breaks = breaks, limits = c(min(breaks), max(breaks))) +
      #theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
    
    if(!is.na(facet_by[mdl])){
      plot_short <- plot_short + facet_wrap(vars(factor((.data[[facet_by[mdl]]]), levels = pair_limits)), ncol = 5)
      print(plot_short)
      #save_big_svg(plot_short, str_c("Boxplot", "ox2", "short", facet_by[mdl], variable_list_truncated[i], sep = "_"), gene, ratio = (4/3))
      ggsave(
        filename = paste(gene, "_", str_c("Boxplot", "ox2", "short", facet_by[mdl], variable_list_truncated[i], sep = "_"), "_paper_", lubridate::today(), ".svg", sep = ""),
        plot_short,
        width = 120,
        height = 120 / (4/3),
        units = "mm",
        device = "svg"
      )
    }else{
      print(plot_short)
      save_small_svg(plot_short, str_c("Boxplot", "ox2", "short", variable_list_truncated[i], sep = "_"), gene, ratio = 1)
    }
    
    #collect plots into a list (required for combining plots in R with ggarrange)
    plot_list[[list_index]] <- plot
    list_index <- list_index + 1
  }
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
