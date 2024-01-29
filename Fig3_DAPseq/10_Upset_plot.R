# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#'Upset plot
#'Build a table and an Upset plot to visualise key overlaps
#'
#'04/01/2024
#'Exploratory UpSetR plots
#'Jaccard index tests to assess similarity between sets
#'
#'09/01/2024
#'Final UpsetR plot
#'Update jaccard test to use all genes
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

#LOAD PACKAGES
require(tidyverse)
require(readxl)
require(qvalue)

#UpSetR used for Upset plots
#Conway et al. (2017) doi:10.1093/bioinformatics/btx364
if(!require(UpSetR)){
  install.packages("UpSetR")
  library(UpsetR)
}

#Jaccard used for Jaccard tests
#Chung, N. C., Miasojedow, B., Startek, M., & Gambin, A. (2019). Jaccard/Tanimoto similarity test and estimation methods for biological presence-absence data. BMC Bioinformatics, 20(Suppl 15), 644. https://doi.org/10.1186/s12859-019-3118-5 
if(!require(jaccard)){
  install.packages("jaccard")
  library(jaccard)
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#LOAD DATA

#Catherine's DAP-seq and other datasets
NAC5_DAP_seq_targets_vs_published <- read_excel("U:/DAP-seq analysis/Target gene analysis/NAC5_DAP-seq_targets_vs_published_2024-01-04.xlsx", 
                                                           sheet = "Targets_vertical", skip = 2)
#GENIE3 targets
GENIE3_targets <- read_csv("U:/DAP-seq analysis/GENIE3 targets/NAC_downstream_targets.csv")

#Metadata
metadata <- read_excel("U:/DAP-seq analysis/Target gene analysis/NAC5_DAP-seq_targets_vs_published_2024-01-04.xlsx", 
                       sheet = "Samples", skip = 9)
#Gene list
transcript_to_gene_refseqv1_1 <- read_csv("U:/DAP-seq analysis/Metadata/transcript_to_gene_refseqv1.1.csv")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#TIDY DATA
#Requires matrix data frame with 1 row per predicted target gene and one column per dataset

#Use key to correspond datasets with an alphabetic index
key <- metadata[c("Index", "GeneID", "Experiment", "Gene", "Homoeolog", "Alternative label")] %>%
  filter(!(Index %in% c("K", "L", "X"))) %>%
  mutate(Label = paste(Experiment, Homoeolog, sep = "_"))

#One column per dataset, one column per gene
#Remove sets K and L (drought treatment from Mao 2022)
DAPseq_v2 <- NAC5_DAP_seq_targets_vs_published %>%
  pivot_longer(names_to = "Index", values_to = "to.gene", cols = everything()) %>%
  filter(!is.na(to.gene) & !(Index %in% c("K", "L")))

#GENIE3 data already has 1 column per dataset, one column per gene
#Convert Gene IDs from v1.0 to v1.1 to match other data; these Gene IDs are directly equivalent (only difference being some LC genes removed in v1.1)

GENIE3_v2 <- GENIE3_targets %>%
  mutate("GeneID" = str_replace(TF, "01G", "02G"),
         "to.gene" = str_replace(to.gene, "01G", "02G")
  ) %>%
  left_join(filter(key, Experiment == "Harrington_2020"), by = c("GeneID")) %>%
  select(Index, to.gene)
  
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#BUILD MATRIX
#Set value to one, pivot to one column per dataset, set values for missing combinations to 0
#Unique removes a few instances where 2 peaks called for same gene or two rice genes had same wheat ortholog
targets <- rbind(DAPseq_v2, GENIE3_v2)

target_matrix <- targets %>%
  arrange(Index) %>%
  unique() %>%
  mutate(value = 1) %>%
  pivot_wider(id_cols = to.gene, names_from = Index) %>%
  mutate(across(everything(), ~replace_na(.x, 0))) %>%
  as.data.frame() #change tibble to data frame to make upset() function happy

with_overlap <- targets %>%
  group_by(to.gene) %>%
  summarise(n_overlaps = n()) %>%
  filter(n_overlaps > 1)

target_matrix_with_overlap <- target_matrix %>%
  inner_join(with_overlap, by = join_by(to.gene))

setwd("U:/DAP-seq analysis/Target gene analysis")
write_csv(target_matrix, paste("upset_matrix_", today(), ".csv", sep = ""))

write_csv(target_matrix_with_overlap, paste("upset_matrix_with_overlap", today(), ".csv", sep = ""))

genes <- transcript_to_gene_refseqv1_1 %>%
  select(gene) %>%
  unique() %>%
  filter(!grepl("LC", gene))
# produces 107892 lines, matching the expected 107891 high-confidence genes

#create target matrix containing ALL HC genes
target_matrix <- targets %>%
  arrange(Index) %>%
  unique() %>%
  mutate(value = 1) %>%
  pivot_wider(id_cols = to.gene, names_from = Index) %>%
  full_join(genes, by = join_by(to.gene == gene)) %>%
  mutate(across(everything(), ~replace_na(.x, 0))) %>%
  as.data.frame() #change tibble to data frame to make upset() function happy

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#Make UpSetR plot

#Automatically picks the largest sets in standard settings
upset(target_matrix, nsets = 8, number.angles = 30, point.size = 3.5, line.size = 2, 
      mainbar.y.label = "Intersections", sets.x.label = "Target genes per dataset", 
      text.scale = c(1.3, 1.3, 1, 1, 2, 0.75)
      )

#Consider sets from within a specific experiment
focus <- key %>% filter(Experiment %in% c("Harrington_2020")) %>% select(Index) %>% unlist() %>% as.vector()

#OR Consider sets for a specific gene
focus <- key %>% filter(Gene == "NAC5-1") %>% select(Index) %>% unlist() %>% as.vector()

#Plot
upset(target_matrix, sets = focus, number.angles = 30, point.size = 3.5, line.size = 2, 
      mainbar.y.label = "Intersections", sets.x.label = "Target genes per dataset", 
      text.scale = c(1.3, 1.3, 1, 1, 2, 0.75), keep.order = TRUE,
      set.metadata = list(data = key, plots = list(list(type = "text", 
                                                        column = "Label", assign = 10))
                          ),
      
)

upset(target_matrix, sets = focus, number.angles = 30, point.size = 3.5, line.size = 2, 
      mainbar.y.label = "Intersections", sets.x.label = "Target genes per dataset", 
      text.scale = c(1.3, 1.3, 1, 1, 2, 0.75), keep.order = TRUE,
      set.metadata = list(data = key, plots = list(list(type = "text", 
                                                        column = "Label", assign = 10))
      ),
      order.by = "freq"
)


#Highest degree interactions across entire set
upset(target_matrix, nsets = 21, number.angles = 30, point.size = 3.5, line.size = 2, 
      mainbar.y.label = "Intersections", sets.x.label = "Target genes per dataset", 
      text.scale = c(1.3, 1.3, 1, 1, 2, 0.75), keep.order = TRUE,
      set.metadata = list(data = key, plots = list(list(type = "text", 
                                                        column = "Label", assign = 10))
      ),
      mb.ratio = c(0.5, 0.5),
      order.by = "degree"
)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#Figure for paper

focus <- key %>% filter(Gene == "NAC5-1") %>% select(Index) %>% unlist() %>% as.vector() %>% rev()
focus <- rev(c("A", "B", "N", "O", "P", "I", "J", "M"))


upset_plot <- upset(target_matrix, sets = focus, number.angles = 30, point.size = 3.5, line.size = 2, 
                    mainbar.y.label = "Genes per intersection", sets.x.label = "Genes per dataset", 
                    text.scale = c(2, 2, 2, 1, 2, 1.5),
                    keep.order = TRUE,
                    mb.ratio = c(0.6, 0.4),
                    set.metadata = list(data = key, plots = list(
                      list(type = "text", column = "Label", assign = 10),
                      list(type="matrix_rows", column = "Experiment",
                           colors = c("Evans" = "#66c2a5", "Mao_2022" = "white", "Chung_2018" = "#fc8d62", "Harrington_2020" = "#8da0cb")
                      )
                    ))
)

upset_plot

# save manually because it isn't a ggplot

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# OVERLAPS
n_samples <- dim(target_matrix)[2]-1

output_matrix <- matrix(nrow = n_samples, ncol = n_samples)
output_jaccard <- matrix(nrow = n_samples, ncol = n_samples)

for(i in 1:n_samples){
  for(j in 1:n_samples){
    x <- target_matrix[,i+1]
    y <- target_matrix[,j+1]
    intersection <- sum(x==1 & y==1)
    union <- sum(x==1 | y==1)
    jaccard_manual <- ifelse(union>0, intersection / union, 0)
    jaccard_auto <- jaccard(x, y) #these two lines should be the same
    # print(paste(intersection, union, jaccard_manual, jaccard_auto))
    output_matrix[i,j] <- intersection
    output_jaccard[i,j] <- jaccard_auto
  }
}


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# JACCARD TESTS

#run jaccard test
#use target matrix containing ALL HC genes
output <- jaccard.test.pairwise(t(target_matrix[,-1]), method = "mca", accuracy = 1e-05, compute.qvalue = FALSE)

#fix issue with qval function
rm_na <- !is.na(output$pvalues)
p <- output$pvalues / max(output$pvalues[rm_na])

#find qvalues
#lambda selected based on section of histogram with "baseline" frequency of p-values
#algorithm based on assuming that p-values from the true null hypothesis will follow a certain distribution
hist(x = output$pvalues, breaks = 20)

output_q <- qvalue(p, lambda = seq(0.1, 0.4, 0.05))
print(output_q$pi0)

plot(x = output$statistics, y = output_q$qvalues) + lines(x = c(0, 1), y = c(0, 1))

# cbind(key, output$statistics) #double check that target_matrix is sorted correctly before doing this
setwd("U:/DAP-seq analysis/Target gene analysis")
write.csv(cbind(key, output_matrix), paste("overlap_", "matrix_", today(), ".csv", sep = ""))
write.csv(cbind(key, output_jaccard), paste("jaccard_", "rawstatistics_", today(), ".csv", sep = ""))
write.csv(cbind(key, output$statistics), paste("jaccard_", "statistics_", today(), ".csv", sep = ""))
write.csv(cbind(key, output$pvalues), paste("jaccard_", "pvalues_", today(), ".csv", sep = ""))
write.csv(cbind(key, output$expectation), paste("jaccard_", "expectation_", today(), ".csv", sep = ""))
write.csv(cbind(key, output_q$qvalues), paste("jaccard_", "qvalues_", today(), ".csv", sep = ""))

# #Examples and comments
# set.seed(1234)
# x = c(rbinom(100,1,.5), rep(0, 1000))
# y = sample(x, 1100, replace = FALSE) #not co-occurring with x
# z = c(x[0:50], rbinom(50,1,.5), rep(0, 1000)) #partially co-occurring with x
# jaccard(x,y)
# jaccard.ev(x,y)
# jaccard.test(x, y, method = "mca", accuracy = 1e-05)
# jaccard.test(x, y, method = "bootstrap")
# jaccard(x,z)
# jaccard.ev(x,z)
# jaccard.test(x, z, method = "mca", accuracy = 1e-05)
# jaccard.test(x, z, method = "bootstrap")
# 
# #pairwise function is offered, effectively this is a for() loop around jaccard.test which I could recode myself if needed
# #this runs with matrix with sets as rows and elements as columns, i.e. datasets as rows and target genes as columns
# test <- matrix(rbinom(500,1,.5), ncol = 100, nrow = 5) #not co-occurring
# 
# output <- jaccard.test.pairwise(test, method = "mca", accuracy = 1e-05, compute.qvalue = TRUE)
# #$statistics centred jaccard statistic (centred against expected values)
# #$pvalues raw p values
# #$expectation expected raw jaccard statistic
# #which scales depending on size of input sets and overall proportion of 1 vs 0 in input sets
# #$qvalues output of function qvalue::qvalue()
# #$qvalues $qvalues q values corrected by false discovery rate fdr
# #$qvalues $pi0 An estimate of the proportion of null p-values.
# 
# output_q <- qvalue::qvalue(output$pvalues)
# #using qvalue() requires q values to be spread across the full distribution from <0.05 to >0.95
# #pi0est throws an error if this is not the case
# #using qvalue_trunc() fixes this by normalizing p values before calling pi0est
# #however this function is only available on github and not on the one I have installed
# #therefore I have added this step manually

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#VISUALISE JACCARD OUTPUT

#Note the jaccard test is two-tailed
#and reports both sets with more cc-occurrence and less co-occurrence as significant
high <- which(output$statistics > 0)

plot(x = output$statistics, y = output$pvalues) + lines(x = c(0, 0), y = c(0, 1))
plot(x = output$expectation, y = output$pvalues)
plot(x = output_matrix, y = output$pvalues, xlim = c(0,40))

plot(x = output$statistics[high], y = output$pvalues[high])  + lines(x = c(0, 0), y = c(0, 1))
plot(x = output$expectation[high], y = output$pvalues[high])
plot(x = output_matrix[high], y = output$pvalues[high])

heatmap(output_matrix)
heatmap(replace_na(output$statistics, 0), symm = TRUE, Rowv = NA)
heatmap(replace_na(output_q$qvalues, 0), symm = TRUE, Rowv = NA)
