library(readr)
library(ggpicrust2)
library(tibble)
library(tidyverse)
library(ggprism)
library(patchwork)

# If you want to analyze KEGG pathway abundance instead of KO within the pathway, turn ko_to_kegg to TRUE.
# KEGG pathways typically have more explainable descriptions.

# Load metadata as a tibble
# data(metadata)
metadata <- read_delim(
  "ps1a.txt",
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) 
str(metadata) #Treatment.Type is character

# Load KEGG pathway abundance
# data(kegg_abundance)
kegg_abundance <- ko2kegg_abundance("picrust2_out_pipeline8125/KO_metagenome_out/pred_metagenome_unstrat.tsv") 
dim(kegg_abundance) #262 x 37
#write.csv(kegg_abundance,"kegg_abundance8625.csv")

# Perform pathway differential abundance analysis (DAA) using ALDEx2 method
# Please change group to "your_group_column" if you are not using example dataset
daa_results_df <- pathway_daa(abundance = kegg_abundance, metadata = metadata, group = "Treatment.Type", daa_method = "edgeR", select = NULL, reference = NULL) 
colnames(daa_results_df) #"feature"    "method"     "group1"     "group2"     "p_values"   "adj_method" "p_adjust"
head(daa_results_df)
tail(daa_results_df)
dim(daa_results_df) #262 x 7
#daa_results_df_sig <- daa_results_df %>% filter(p_adjust < 0.05)
#dim(daa_results_df_sig) #102 x 7


# Filter results for ALDEx2_Welch's t test method
# Please check the unique(daa_results_df$method) and choose one
daa_sub_method_results_df <- daa_results_df[daa_results_df$method == "edgeR", ]
dim(daa_sub_method_results_df) #262 x 7 
head(daa_sub_method_results_df)

daa_sub_method_results_df1<-daa_sub_method_results_df[1:50,]
dim(daa_sub_method_results_df1) #50x7
daa_sub_method_results_df2<-daa_sub_method_results_df[51:102,]
dim(daa_sub_method_results_df2) #52x7
daa_sub_method_results_df3<-daa_sub_method_results_df[103:193,]
dim(daa_sub_method_results_df3) #91x7
daa_sub_method_results_df4<-daa_sub_method_results_df[194:262,]
dim(daa_sub_method_results_df4) #69x7

# Annotate pathway results using KO to KEGG conversion
library("KEGGREST"); packageVersion("KEGGREST") #1.44.1 #make sure that ggpicrust2 is fully updated.
daa_annotated_sub_method_results_df1 <- ggpicrust2::pathway_annotation(pathway = "KO", daa_results_df = daa_sub_method_results_df1, ko_to_kegg = TRUE)
dim(daa_annotated_sub_method_results_df1) #17 x 11
#head(daa_annotated_sub_method_results_df1)

daa_annotated_sub_method_results_df2 <- pathway_annotation(pathway = "KO", daa_results_df = daa_sub_method_results_df2, ko_to_kegg = TRUE)
dim(daa_annotated_sub_method_results_df2) #24 x 11
#head(daa_annotated_sub_method_results_df2)

daa_annotated_sub_method_results_df3 <- pathway_annotation(pathway = "KO", daa_results_df = daa_sub_method_results_df3, ko_to_kegg = TRUE)
dim(daa_annotated_sub_method_results_df3) #31x 11

daa_annotated_sub_method_results_df4 <- pathway_annotation(pathway = "KO", daa_results_df = daa_sub_method_results_df4, ko_to_kegg = TRUE)
dim(daa_annotated_sub_method_results_df4) #29 x 11

daa_annotated_sub_method_results_df<-rbind(daa_annotated_sub_method_results_df1,daa_annotated_sub_method_results_df2,daa_annotated_sub_method_results_df3,daa_annotated_sub_method_results_df4)
dim(daa_annotated_sub_method_results_df) #101 x 11
str(daa_annotated_sub_method_results_df)
AnnotatedSig05KOtoKeGGS <- apply(daa_annotated_sub_method_results_df,2,as.character)
#write.csv(AnnotatedSig05KOtoKeGGS,"AnnotatedSig05KOtoKeGGS_8725.csv")

#Generate pathway error bar plot
#Please change Group to metadata$your_group_column if you are not using example dataset
#this function is ugly. I feel like I can make a figure outside of "pathway_errorbar" on my own computer.

# max(daa_annotated_sub_method_results_df$p_adjust) #max= 0.04851148
# daa_annotated_sub_method_results_df[56,]
# anno.fix1<-daa_annotated_sub_method_results_df[-56,]
# anno.fix1[90,]
# aano.fix2<-anno.fix1[-90,]
# min(aano.fix2$p_values)
# sel<-c("ko04062", "ko04111", "ko05412", "ko04650", "ko00311", "ko04142", "ko00604", "ko05142", "ko00909", "ko00513", "ko04973", "ko00510", "ko05110", "ko04974", "ko00905", "ko05140", "ko00592", "ko00590", "ko00062", "ko04662", "ko03070", "ko04370", "ko04730", "ko04740", "ko05120", "ko05322", "ko04660", "ko04380", "ko00940", "ko00941", "ko00943", "ko00100", "ko00945", "ko01056", "ko04145", "ko04360", "ko05215", "ko05214", "ko01055", "ko00902", "ko00531", "ko04916", "ko00633", "ko00364", "ko05130", "ko04114", "ko04914", "ko00430", "ko00130", "ko04330", "ko03050", "ko00361", "ko00603", "ko04270", "ko03015", "ko00523", "ko05144", "ko04020", "ko00120", "ko00965", "ko05414", "ko04614", "ko05131", "ko00196", "ko00785", "ko05220", "ko02030", "ko04744", "ko04745", "ko04530", "ko00522", "ko04612", "ko04621", "ko04971", "ko04970", "ko00511", "ko04972", "ko00232", "ko04350", "ko03022", "ko04512", "ko04630", "ko05100", "ko05410", "ko04514", "ko04320", "ko00906", "ko00901", "ko04520", "ko04510", "ko04012", "ko04010", "ko00944", "ko04640", "ko04310", "ko04912", "ko05160", "ko04144", "ko04720", "ko04722")
# p <- pathway_errorbar(abundance = kegg_abundance, daa_results_df = aano.fix2, Group = metadata$Treatment.Type, p_values_threshold = 0.05, order = "pathway_class", select = sel, ko_to_kegg = TRUE, p_value_bar = TRUE, colors = NULL, x_lab = "pathway_name")



# If you want to analyze EC, MetaCyc, and KO without conversions, turn ko_to_kegg to FALSE.
# Load metadata as a tibble

# Load KO abundance as a data.frame
# data(ko_abundance)
ko_abundance <- read.delim("picrust2_out_pipeline8125/KO_metagenome_out/pred_metagenome_unstrat.tsv")
head(ko_abundance)

# Perform pathway DAA using edgeR method
# Please change column_to_rownames() to the feature column if you are not using example dataset
# Please change group to "your_group_column" if you are not using example dataset
# EXPLORATORY but likely going to use the picrust2_edgeR analyses in another script.
ko_daa_results_df <- pathway_daa(abundance = ko_abundance %>% column_to_rownames("function."), metadata = metadata, group = "Treatment.Type", daa_method = "edgeR", select = NULL, reference = NULL)
head(ko_daa_results_df)
dim(ko_daa_results_df) #7194 x 7

# Filter results for ALDEx2_Kruskal-Wallace test method
ko_daa_sub_method_results_df <- ko_daa_results_df[ko_daa_results_df$method == "edgeR", ]
dim(ko_daa_sub_method_results_df)
min(ko_daa_sub_method_results_df$p_adjust)
ko_daa_sub_method_results_df <- ko_daa_sub_method_results_df %>% filter(p_adjust < 0.05)
dim(ko_daa_sub_method_results_df) #4336 x 7

# Annotate pathway results without KO to KEGG conversion
ko_daa_annotated_sub_method_results_df <- pathway_annotation(pathway = "KO", daa_results_df = ko_daa_sub_method_results_df, ko_to_kegg = FALSE)
dim(ko_daa_annotated_sub_method_results_df) #4336 x 8 
colnames(ko_daa_annotated_sub_method_results_df) #feature"     "method"      "group1"      "group2"      "p_values"    "adj_method"  "p_adjust"    "description"
head(ko_daa_annotated_sub_method_results_df)

# Generate pathway error bar plot
# Please change column_to_rownames() to the feature column
# Please change Group to metadata$your_group_column if you are not using example dataset
ko_p <- pathway_errorbar(abundance = ko_abundance %>% column_to_rownames("function."), daa_results_df = ko_daa_annotated_sub_method_results_df, Group = metadata$Treatment.Type, p_values_threshold = 0.05, order = "group",
                      select = ko_daa_annotated_sub_method_results_df %>% arrange(p_adjust) %>% slice(1:30) %>% dplyr::select(feature) %>% pull(), 
                      ko_to_kegg = FALSE, 
                      p_value_bar = TRUE, 
                      colors = NULL, 
                      x_lab = "description")

# Workflow for MetaCyc Pathway and EC

# Load MetaCyc pathway abundance and metadata
#EXPLORATORY ONLY NoOT OFFICIAL ANALYSES
metacyc_abundance<-read.delim("picrust2_out_pipeline8125/pathways_out/path_abun_unstrat.tsv")
dim(metacyc_abundance) #413 x 38
colnames(metacyc_abundance)

# Perform pathway DAA using LinDA method
# Please change column_to_rownames() to the feature column if you are not using example dataset
# Please change group to "your_group_column" if you are not using example dataset
metacyc_daa_results_df <- pathway_daa(abundance = metacyc_abundance %>% column_to_rownames("pathway"), metadata = metadata, group = "Treatment.Type", daa_method = "edgeR")
dim(metacyc_daa_results_df) #413 x 7
colnames(metacyc_daa_results_df)
head(metacyc_daa_results_df)

# Annotate MetaCyc pathway results without KO to KEGG conversion
metacyc_daa_annotated_results_df <- pathway_annotation(pathway = "MetaCyc", daa_results_df = metacyc_daa_results_df, ko_to_kegg = FALSE)
head(metacyc_daa_annotated_results_df)
dim(metacyc_daa_annotated_results_df) #413 x 8
metacyc_daa_annotated_results_df_sig <- metacyc_daa_annotated_results_df %>% filter(p_adjust < 0.05)
dim(metacyc_daa_annotated_results_df_sig) #178 x 8

# Generate pathway error bar plot
# Please change column_to_rownames() to the feature column
# Please change Group to metadata$your_group_column if you are not using example dataset
myc<-pathway_errorbar(abundance = metacyc_abundance %>% column_to_rownames("pathway"), daa_results_df = metacyc_daa_annotated_results_df_sig, Group = metadata$Treatment.Type, ko_to_kegg = FALSE, p_values_threshold = 0.05, order = "group", select = metacyc_daa_annotated_results_df_sig %>%
                   arrange(p_adjust) %>%
                   slice(1:30) %>%
                   select("feature") %>% pull(), p_value_bar = TRUE, colors = NULL, x_lab = "description")

# Generate pathway heatmap
# Please change column_to_rownames() to the feature column if you are not using example dataset
# Please change group to "your_group_column" if you are not using example dataset
feature_with_p_0.05 <- metacyc_daa_results_df %>% filter(p_adjust < 0.05)
dim(feature_with_p_0.05) #178 x 7
library(ggh4x) #requires ggplot 3.4.2 or greater
pathway_heatmap(abundance = metacyc_abundance %>% filter(pathway %in% feature_with_p_0.05$feature) %>% column_to_rownames("pathway"), metadata = metadata, group = "Treatment.Type")

# Generate pathway PCA plot
# Please change column_to_rownames() to the feature column if you are not using example dataset
# Please change group to "your_group_column" if you are not using example dataset
metacyc_abundance1<- metacyc_abundance %>% column_to_rownames("pathway")
head(metacyc_abundance1)
pathway_pca(abundance = metacyc_abundance1, metadata = metadata, group = "Treatment.Type")

# Run pathway DAA for multiple methods
# Please change column_to_rownames() to the feature column if you are not using example dataset
# Please change group to "your_group_column" if you are not using example dataset
methods <- c("ALDEx2", "DESeq2", "edgeR")
dddddaa_results_list <- lapply(methods, function(method) {
  pathway_daa(abundance = metacyc_abundance %>% column_to_rownames("pathway"), metadata = metadata, group = "Treatment.Type", daa_method = method)
})

# Compare results across different methods
comparison_results <- compare_daa_results(daa_results_list = daa_results_list, method_names = c("ALDEx2_Welch's t test", "ALDEx2_Wilcoxon rank test", "DESeq2", "edgeR"))

#SAVE THE WHOLE ENVIRONMENT FOR NEW SCRIPTS
#save(list=ls(),file="ggpicrust2.rda")
