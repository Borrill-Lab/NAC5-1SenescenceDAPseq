#'Overexpression phenotyping data preparation
#'21/10/2021
#'Open the Excel file, fill in blank columns
#'calculate dates as "days after sowing" or "days after heading"
#'and save as a .csv file.
#'
#'07/03/2021
#'Adapt to use with Kronos NAC5 NAP data (incomplete). Remove unnecessary lines.
#'Using sowing date as column rather than fixed.
#'
#'29/03/2021
#'Use with Kronos_phenotyping_2022-03-28.xlsx
#'I have manually added the genotypes to the spreadsheet (so will need to sort this again when I get the final datafile)
#'
#'08/04/2022
#'Calculate metrics
#'Linear interpolation method Using the r function stats::approx
#'Based on NAM_calculate_metrics_2022-03-17
#'
#'27/06/2022
#'Use full dataset from 2022-04-28
#'
#'12/08/2022
#'Copy Kronos_data_preparation_2022-06-27
#'Use NAC5 NAP overexpression CER data
#'
#'05/09/2022
#'Join Marvin data
#'
#'29/09/2022
#'Use properly unique and adjusted Marvin data
#'
#'01/12/2023
#'Prepare data for NAC5 paper
#'Make a version of the data with rows for 'combined' controls associated with multi-copy lines
#'
#'23/01/2024
#'Make a version of SPAD table with 'combined' controls and variable Day_approx. Omit 8.3.
#'
#'This is the main script.
#'
#'Syntax rules: variables snake_case, Column Names in Caps, use pipe %>% where beneficial.
#'Built under R 4.0.4
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#'SOURCE
#'Packages and source files
require(readxl); require(readr); require(lubridate); require(MESS) #for auc()
require(tidyverse)

#Run twice, once for each gene
gene = "NAC5"
# gene= "NAP"

setwd("C:/Users/evansc/OneDrive - Norwich BioScience Institutes/NAC transgenics/2021-11-10 NAC5 NAP overexpression CER")
# setwd("/Users/u1984449/OneDrive - Norwich BioScience Institutes/NAC transgenics/2021-11-10 NAC5 NAP overexpression CER")

#34 columns vs 33 columns
if(gene == "NAC5"){
Overexpression_phenotyping <- read_excel("2022-02-14_phenotyping_data/NAC5_NAP_CER_phenotyping_2022-04-28.xlsx",
                                 sheet = paste(gene),
                                 col_types = c("numeric", "numeric", "text", 
                                 "numeric", "numeric", "numeric", 
                                 "date", "numeric", "date", "date", 
                                 "numeric", "numeric", "date", "numeric", 
                                 "date", "numeric", "date", "numeric", 
                                 "date", "numeric", "date", "numeric", 
                                 "date", "numeric", "date", "numeric", 
                                  "numeric", "numeric", "skip", "numeric", 
                                 "numeric", "text", "date", "numeric"))
} else {
  Overexpression_phenotyping <- read_excel("2022-02-14_phenotyping_data/NAC5_NAP_CER_phenotyping_2022-04-28.xlsx",
                                           sheet = paste(gene),
                                           col_types = c("numeric", "numeric", "text", 
                                                         "numeric", "numeric", "numeric", 
                                                         "date", "numeric", "date", "date", 
                                                         "numeric", "numeric", "date", "numeric", 
                                                         "date", "numeric", "date", "numeric", 
                                                         "date", "numeric", "date", "numeric", 
                                                         "date", "numeric", "date", "numeric", 
                                                         "numeric", "skip", "numeric", 
                                                         "numeric", "text", "date", "numeric"))
}

Marvin_data <- read_csv("2022-08-12_clean_data/Overexpression_Marvin_unique_2022-09-29.csv")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#'EDIT
#'Key parameters

sowing_date <- lubridate::ymd('2021-12-22')

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#'TIDY
#'Explore the file and reformat data

summary(Overexpression_phenotyping)

#Neaten various columns, add new columns, and remove photo columns
if(gene == "NAC5"){
  Overexpression_phenotyping_clean <- Overexpression_phenotyping %>%
    mutate(
      SPAD_8 = replace_na(Date_8, 0)
    )
} else if(gene == "NAP"){
  Overexpression_phenotyping_clean <- Overexpression_phenotyping %>%
    mutate(SPAD_8 = replace_na(SPAD_8, 0)) %>%
    rename(Date_6 = SPAD_52)
}

Overexpression_phenotyping_clean <- Overexpression_phenotyping_clean %>%
  arrange(Index) %>%
  mutate(Heading_date = lubridate::ymd(Heading_date),
         Tiller_2 = replace_na(Tiller_2, 0),
         Date_8 = Date_7 + days(7),
         Gene_name = gene,
         snapped_pre_ped_sen = snapped < Peduncle_senescence | !is.na(snapped) & is.na(Peduncle_senescence),
         Shelf = as.factor(paste("Shelf", ((Block -1) %/% 5) + 1)),
         Mildew = grepl("mildew", Comment)) %>%
  tidyr::extract(Line_name, "Genotype", regex = "_([a-z]+)$", remove = FALSE)

not_needed <- c("photo", "Photo", "Column") #remove photo comments and unlabelled columns
Overexpression_phenotyping_clean <- Overexpression_phenotyping_clean %>%
  select(!contains(not_needed))

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#CALCULATE
#Calculate all dates as days since SOWING
#Using lubridate
#Also calculate Leaf and Peduncle Senescence as Days after Heading (DAH)
#Done slightly differently to NAC-5A data where heading was set as 'day of the year'
days.after.sowing <- function(test_date, sowing_date){
  return(as.numeric(interval(start = sowing_date, end = test_date), "days"))
} 

Overexpression_phenotyping_yday <- Overexpression_phenotyping_clean %>%
  mutate(Leaf_senescence_DAH = days.after.sowing(Leaf_senescence, sowing_date = Heading_date),
         Peduncle_senescence_DAH = days.after.sowing(Peduncle_senescence, sowing_date = Heading_date),
         Date_0 = Heading_date,
         across(c(starts_with("Date")),
                ~ days.after.sowing(.x, sowing_date = Heading_date)),
         across(c(Heading_date, Leaf_senescence, Peduncle_senescence),
                ~ days.after.sowing(.x, sowing_date = sowing_date)))



Long_data <- Overexpression_phenotyping_yday %>%
  rename(SPAD_0 = SPAD_heading) %>%
  pivot_longer(cols = starts_with("Date")|starts_with("SPAD"), names_pattern = "(.*)_([0-9])$", names_to = c("Category","Week")) %>%
  pivot_wider(names_from = Category, values_from = value)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#AUC NORMALISATION TEST

#Mean senescence values are dataset-specific as they depend on how many data are collected.
#Therefore if I want to compare between sowing dates I'll need to normalise to a range in days after sowing or heading

#normalise to a range in days after sowing or heading
Long_data %>% group_by(Week) %>%
  summarise(max(Date),
            min(Date),
            mean(Date)
  )
#End AUC at minimum scoring endpoint - ensures interpolation only, no extrapolation.
#Start point is less important as additional score points are 0 so do not add to total auc
auc_limit <- Long_data %>% filter(Week == max(Week)) %>%
  summarise(
    min_DAH = floor(min(Date))
  )


#Calculate metrics from SPAD table
NA_count <- Long_data %>%
  filter(!is.na(Plant_num)) %>%
  group_by(Line_name, Plant_num) %>%
  arrange(-(Date)) %>%
  summarise(
    NA_count = sum(is.na(SPAD))
  )

Summary_Days_after_heading <- inner_join(Long_data, NA_count, by = c("Line_name", "Plant_num")) %>%
  filter(!is.na(Plant_num) & NA_count < 4) %>%
  group_by(Line_name, Plant_num) %>%
  arrange(-(Date)) %>%
  summarise(
    TT_SPAD10 = approx(SPAD, Date, c(10), ties = "ordered")$y,
    TT_SPAD30 = approx(SPAD, Date, c(30), ties = "ordered")$y,
    TT_SPAD40 = approx(SPAD, Date, c(40), ties = "ordered")$y,
    Dur_Leaf_Sen = TT_SPAD10 - TT_SPAD40,
    AUC_Leaf_Sen = auc(Date, SPAD, type = "linear")
  )

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#JOIN
#
#Join calculated values AND Marvin data

Overexpression_phenotyping_yday_v2 <- left_join(Overexpression_phenotyping_yday, NA_count, by = c("Line_name", "Plant_num"))

Overexpression_phenotyping_yday_v3 <- left_join(Overexpression_phenotyping_yday_v2, Summary_Days_after_heading, by = c("Line_name", "Plant_num"))

Overexpression_phenotyping_yday_v4 <- Overexpression_phenotyping_yday_v3 %>%
  mutate(ID = paste(Line_name, Plant_num, sep = "_")) %>%
  left_join(Marvin_data, by = c("ID"))

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#OUTLIERS

#EDIT: keep snapped_pre_ped_sen
Overexpression_outliers <- Overexpression_phenotyping_yday_v4 %>%
  filter(is.na(Plant_num)|grepl("damage", Comment)|grepl("tagged", Comment))

View(Overexpression_outliers)

Overexpression_analysis <- Overexpression_phenotyping_yday_v4 %>%
  anti_join(Overexpression_outliers, by = "Index")

Overexpression_SPAD_analysis <- Long_data %>%
  anti_join(Overexpression_outliers, by = "Index") %>%
  mutate(Day_approx = as.numeric(Week) * 7)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#GENERATE COMBINED CONTROLS DATA for the paper

#To add additional binning columns, Category and Binary
Category_table <- data.frame(Genotype = c("hom", "multi", "null", "null"),
                             Category = c("hom", "multi", "null", "combo1"),
                             Binary = c("transgenic", "transgenic", "WT", "WT")
)

#Join category table with data; this will create 3 rows for each original null row (80 obs -> 140 obs)
#Set one copy of each null for its original partner and one copy for each of 8.2_multi and 8.3_multi
#Remove Shelf 1 as these data were lower quality
Combined_control_data <- Overexpression_analysis %>%
  filter(Shelf == "Shelf 2") %>%
  full_join(Category_table, relationship = "many-to-many", by = "Genotype") %>%
  tidyr::extract(Line_name, into = "Line_prefix", regex = "([0-9.]{3,5})_", remove = FALSE) %>%
  mutate(Pair = ifelse(Category == "combo1", "8.2",
                          Line_prefix))

Combined_control_SPAD_data <- Overexpression_SPAD_analysis %>%
  filter(Shelf == "Shelf 2") %>%
  full_join(Category_table, relationship = "many-to-many", by = "Genotype") %>%
  tidyr::extract(Line_name, into = "Line_prefix", regex = "([0-9.]{3,5})_", remove = FALSE) %>%
  mutate(Pair = ifelse(Category == "combo1", "8.2",
                       Line_prefix))

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#WRITE

setwd("./2022-08-12_clean_data")
# setwd("/Users/u1984449/OneDrive - Norwich BioScience Institutes/NAC transgenics/2021-11-10 NAC5 NAP Overexpression greenhouse/2022-03-07_clean_data")
write_csv(Overexpression_phenotyping_clean, file = paste('Overexpression_phenotyping_clean_', gene, "_", today(), '.csv', sep = ""))
write_csv(Overexpression_analysis, file = paste('Overexpression_phenotyping_yday_', gene, "_", today(), '.csv', sep = ""))
write_csv(Overexpression_SPAD_analysis, file = paste('Overexpression_SPAD_yday_', gene, "_", today(), '.csv', sep = ""))

write_csv(Combined_control_data, file = paste('Combined_control_data_', gene, "_", today(), '.csv', sep = ""))
write_csv(Combined_control_SPAD_data, file = paste('Combined_control_SPAD_data_', gene, "_", today(), '.csv', sep = ""))

