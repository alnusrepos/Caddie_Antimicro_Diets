############################################################
############################################################
############################################################
# Metabolomics Analyses (Single R Script)
############################################################
############################################################
############################################################

library(tidyverse)
library(vegan)
library(mixOmics)
library(dplyr)
library(stringr)
library(tibble)
library(caret)
library(ggplot2)
library(ggthemes)
library(pairwiseAdonis)

set.seed(33)

###############################################################################
# 1. Load data (modify to your file paths)
###############################################################################
# feature_table: rows = samples, columns = Feature_XXXX
# annotations:   rows = features, column "Feature"
# metadata:      rows = samples, sample metadata

library(dplyr)

# setwd
setwd("/Users/jrdickey/Documents/UCSD_EBE_PostDoc/Mentee_Materials/Dahlia_Loomis/Revisions/MetaBolomics/Files/RSave")

#load 

load(file="dat_ALL_C.rda")
load("an_final.rda")
load("md.rda")

dat_C <- dat_ALL_C %>% 
  rownames_to_column("filename") %>% 
  left_join(., md[, c("filename", "Type")]) %>% 
  dplyr::filter(Type %in% c("RTE_Decomp", "SAMPLE_POWDER_METAB","Fresh_Leaves")) %>% #remove pools, keep just samples
  dplyr::select(-Type) %>% 
  column_to_rownames("filename") %>% 
  dplyr::select_if(~sum(.) > 0) %>% 
  dplyr::select_if(~sum(0 != .) > 1) 

dim(dat_C) #365 x 7128

#rm(dat_ALL_C)

#update metadata and annotation objects too
an_final <- an_final %>% 
  dplyr::filter(Feature %in% colnames(dat_C))

md <- md %>% 
  dplyr::filter(filename %in% rownames(dat_C))

feature_table <- read.csv("2026-04-21ft_tab.csv", row.names = 1)

annotations <- read.csv("2026-04-21anno_table.csv")
dim(annotations)
annotations <- annotations[,2:132] #remove X column
annotations <- annotations %>% dplyr::select(sort(names(.))) #alphabetically sort
annotations <- annotations %>% dplyr::select(-sus_deltamz) #remove this column

metadata <- read.csv("2026-02-11_metadata.csv")
dim(metadata) #365 by 12
metadata <- metadata[,2:11] #remove X column

#find Internal Standard, sulfachlorpyrazidine 
length(colnames(feature_table)) #7128
head(colnames(feature_table))
which(colnames(feature_table) == "Feature_14924") #interger 0 (oop!)
"Feature_14924" %in% colnames(feature_table) #FALSE
"Feature_15989" %in% colnames(feature_table) #FALSE
"Feature_14924" %in% colnames(dat_C) #FALSE
"Feature_15989" %in% colnames(dat_C)

#combine data_samp_metabext, data_samp_fresh, data_samp_decomp
#briefdat<-rbind(data_samp_metabext,data_samp_fresh,data_samp_decomp)
load("briefdat.rda")

briefdat_new <- briefdat %>%
  arrange(SampleID)
briefdat_new$SampleID==rownames(feature_table) #nice

is_ft_save<-cbind(briefdat_new$SampleID,briefdat_new$Feature_14924)
#write.csv(is_ft_save,file="is_ft.csv")

load("2026-04-16_CGMDiet_feature_annotations1.rda")
#an_final<-read.csv(file="2026-04-16_CGMDiet_feature_annotations1.csv") #load in file created on line 745/746 to get full annotation to prune IS row.

which(an_final$Feature == "Feature_14924") #4586
is_ann_save<-an_final[4586,]
#write.csv(is_ann_save,file="is_ann_save.csv")

is_ft <- read.csv("is_ft.csv")
is_ft <- is_ft[,2:3]
colnames(is_ft)<-c("SampleID","Feature_14924")
rownames(is_ft) <- is_ft$SampleID
is_ft <- is_ft[rownames(is_ft) %in% rownames(feature_table), ]

is_ann <- read.csv("is_ann_save.csv")
is_ann <- is_ann[,2:131]
is_ann <- is_ann %>% dplyr::select(sort(names(.))) #alphabetically sort
colnames(is_ann) == colnames(annotations) #all true
rownames(is_ann) <- "Feature_14924" #"Feature_14413" #upon re-running this doesnt match anymore. Need to make sure the IS isn't some other feature number. 

# ft_ans<-read.csv(file="2026-04-16_CGMDiet_feature_annotations1.csv")
# 
# fbmn_IS <- ft_ans %>% 
#   dplyr::filter(str_detect(Combined_Name, regex("sulf", ignore_case = TRUE))) %>% 
#   distinct(Combined_Name, .keep_all = TRUE) 
# dim(fbmn_IS)
# 
# fbmn_IS[1:5, c("Scan", "Combined_Name", "row.m.z", "row.retention.time")]
# fbmn_IS[1:5, c("Scan", "Combined_Name")]
# fbmm_ISsimp <- fbmn_IS %>% dplyr::select("Scan", "Combined_Name", "row.m.z", "row.retention.time") %>% as_tibble()
# fbmm_ISsimp. #doesnt appear so.   

#join internal standard data with clean annotation file and with feature table.
rownames(annotations) <- annotations$Feature
annotations <- rbind(annotations,is_ann[1:130])

fbmn_IS <- annotations %>% 
  dplyr::filter(str_detect(Combined_Name, regex("sulf", ignore_case = TRUE))) %>% 
  distinct(Combined_Name, .keep_all = TRUE)
fbmn_IS[1:12, c("Scan", "Combined_Name", "row.m.z", "row.retention.time")] #sulfadimethoxine is annotated, but i believe these arent TRUE IS based on mass and retention time. 

rownames(feature_table) == rownames(is_ft)
feature_table <- cbind(feature_table,is_ft$Feature_14924)
colnames(feature_table)[7000:7129]
colnames(feature_table)[7129] <- "Feature_14924"

###############################################################################
# 2. Check alignment of sample names
###############################################################################
metadata <- metadata %>% dplyr::arrange(filename)

if(!identical(rownames(feature_table), metadata$filename)){
  stop("Sample names in feature_table rows do NOT match metadata$SampleID.")
} #sanity check, if anything pops up you'll need to fix it

###############################################################################
# 3. Internal Standard Normalization (replace Feature_14413 if needed)
###############################################################################
IS_feature <- "Feature_14924"

if(!(IS_feature %in% colnames(feature_table))){
  stop(paste("Internal standard", IS_feature, "NOT found in feature_table."))
} #sanity check

# Create info_feature from annotations
info_feature <- annotations %>%
  dplyr::select(
    Feature,               # corresponds to row ID / Feature_ID
    row.m.z,          # corresponds to row m/z
    row.retention.time,         # corresponds to retention time
    correlation.group.ID   # correlation group ID
  ) %>%
  dplyr::rename(
    mz = row.m.z,
    RT = row.retention.time,
    Corr_ID = correlation.group.ID
  )

info_feature$Feature <- as.character(info_feature$Feature)

info_feature_complete <- info_feature %>%
  left_join(annotations, by = "Feature") #OK just renaming

#Create dataframe for analysis
data <- feature_table %>%
  rownames_to_column("SampleID") %>% 
  arrange(SampleID) %>% 
  distinct(SampleID, .keep_all = TRUE)

setdiff(data$SampleID, metadata$filename) #sanity check, should be zero
setdiff(metadata$filename, data$SampleID) #sanity check, should be zero

metadata_metabolomics <- data %>%
  dplyr::select(SampleID) %>%
  left_join(metadata, by = c("SampleID" = "filename"))

metadata_metabolomics <- metadata_metabolomics %>%
  mutate(
    hfa = case_when(
      Deploy == "None" ~ "None",
      Origin == Deploy ~ "Home",
      Origin != Deploy ~ "Away"
    )
  ) #creating HFA column

metadata_metabolomics$Day_Hfa <- paste(metadata_metabolomics$Day,
                                       metadata_metabolomics$hfa,
                                       sep = "_") #creating combined column "Day15_Away", "Day15_Home"

# ---- Internal standard normalization adapted for your data ----
# internal standard feature id
IS_feature <- "Feature_14924"
IS_name <- "sulfachlorpyridazine"  #"human" name for clarity

# sanity check: IS exists in the data
if(!IS_feature %in% colnames(data)){
  stop(paste0("Internal standard '", IS_feature, "' not found in 'data' columns."))
}

# create a small IS table (SampleID + IS value + optional sample order from metadata)
table_IS <- data %>%
  dplyr::select(SampleID, all_of(IS_feature)) %>%
  dplyr::rename(!!IS_name := all_of(IS_feature))

library(scales)
###############################################################################
# 1. New step: filter out features seen as instument noise (%% of internal standard height)
###############################################################################
str(data) 
rm(list = setdiff(ls(), c("data", "metadata","info_feature_complete","table_IS")))

#get the global intensity distribution
# remove SampleID
mat <- as.matrix(data[ , -1]) #365 7129
rownames(mat)<-data$SampleID
colnames(mat)

# vector mat# vector of non-zero intensities
nonzero <- mat[mat > 0]

summary(nonzero)
#Min.  1st Qu.   Median      Mean     3rd Qu.     Max. 
# 47      5303     18736    157302     64953 112350600

quantile(nonzero, probs = c(0.001, 0.01, 0.05, 0.1))
#   0.1%        1%        5%       10% 
# 115.4855  282.4131  1231.3724 2032.4350 

#examine log distribution
hist(log10(nonzero), breaks = 100)

#check out internal standard intensity
summary(table_IS$sulfachlorpyridazine)
#Min.    1st Qu.  Median  Mean   3rd Qu.   Max. 
#0    1947104     2192674   2188132   4416563

#what is one percent of the internal standard mean?
median(table_IS$sulfachlorpyridazine) * 0.01 # 21,926.74
mean(table_IS$sulfachlorpyridazine) * 0.01 # 21,881.32 (just curious)

#this is not an insignificant amount ^^^
quantile(nonzero, 0.5) #18,735.66 #that's over 50% removal
quantile(nonzero, 0.1) #2032.435

###############################################################################
# 2. Check alignment of sample names
###############################################################################
if(!identical(rownames(mat), metadata$filename)){
  stop("Sample names in feature_table rows do NOT match metadata$SampleID.")
} #sanity check

###############################################################################
# 3. Internal Standard Normalization (replace Feature_14413 if needed)
###############################################################################
#i want to combine rownames(filtered_mat) with filtered_mat as SampleID
dim(mat) #365 x 7129

data1 <- data.frame(
  SampleID = rownames(mat),
  mat,
  check.names = FALSE
)

metadata_metabolomics <- data1 %>%
  dplyr::select(SampleID) %>%
  left_join(metadata, by = c("SampleID" = "filename"))

metadata_metabolomics <- metadata_metabolomics %>%
  mutate(
    hfa = case_when(
      Deploy == "None" ~ "None",
      Origin == Deploy ~ "Home",
      Origin != Deploy ~ "Away"
    )
  ) #creating HFA column

metadata_metabolomics$Day_Hfa <- paste(metadata_metabolomics$Day,
                                       metadata_metabolomics$hfa,
                                       sep = "_") #creating combined column "Day15_Away", "Day15_Home"

IS_feature <- "Feature_14924"
IS_name <- "sulfachlorpyridazine"  #"human" name for clarity

# sanity check: IS exists in the data
if(!(IS_feature %in% colnames(data1))){
  stop(paste("Internal standard", IS_feature, "not found in feature_table."))
} #sanity check

#recreate a small IS table (SampleID + IS value + optional sample order from metadata)
table_IS <- data1 %>%
  dplyr::select(SampleID, all_of(IS_feature)) %>%
  dplyr::rename(!!IS_name := all_of(IS_feature))

# compute CV (coefficient of variation) for QC
sd_IS <- sd(table_IS[[IS_name]], na.rm = TRUE) #390276.044325366
mean_IS <- mean(table_IS[[IS_name]], na.rm = TRUE) #2188132.23534247
cv_IS <- sd_IS / mean_IS #0.178360355933554

# detect outliers using 1.5 * IQR rule
Q1 <- quantile(table_IS[[IS_name]], 0.25, na.rm = TRUE)
Q3 <- quantile(table_IS[[IS_name]], 0.75, na.rm = TRUE)
IQRv <- Q3 - Q1
IS_lower_bound <- Q1 - 1.5 * IQRv
IS_upper_bound <- Q3 + 1.5 * IQRv

IS_outliers <- table_IS %>%
  filter(.data[[IS_name]] < IS_lower_bound | .data[[IS_name]] > IS_upper_bound)

if(nrow(IS_outliers) > 0){
  message("Found ", nrow(IS_outliers), " IS outlier(s). They will be excluded from IS mean calculation.")
} #Found 9 IS outlier(s). They will be excluded from IS global mean calculation.

# create cleaned IS table (without outliers)
table_IS_cleaned <- table_IS %>%
  filter(!(SampleID %in% IS_outliers$SampleID))

# compute global mean of IS from cleaned samples
global_IS_mean <- mean(table_IS_cleaned[[IS_name]], na.rm = TRUE) #2156635.80234604

data_with_IS <- data1 %>%
  left_join(table_IS_cleaned, by = "SampleID") 

# If some samples had IS outlier (and thus were removed from table_IS_cleaned),
# they will have NA for the IS. Decide how to handle: here we keep them and normalize using global_IS_mean
# (alternatively you could remove them entirely or impute)

missing_IS_count <- sum(is.na(data_with_IS[[IS_name]]))

if(missing_IS_count > 0){
  message("Hi, hello: ", missing_IS_count, " sample(s) have missing IS after outlier removal. Those samples will be normalized using global IS mean.")
  # fill with global mean so they still get normalized
  data_with_IS[[IS_name]] <- ifelse(is.na(data_with_IS[[IS_name]]), global_IS_mean, data_with_IS[[IS_name]])
}

# Perform normalization across the entire dataset
# Identify feature columns to normalize. Assumes feature columns start with "Feature_"
feature_cols <- grep("^Feature_", colnames(data_with_IS), value = TRUE)

# Normalize: for each feature column, divide by the sample IS value and multiply by global_IS_mean
data_IS_norm <- data_with_IS %>%
  mutate(across(all_of(feature_cols),
                ~ .x / .data[[IS_name]] * global_IS_mean)) #347 x 6350

# compute IS_scaled (for QC) = IS / global mean
data_IS_norm <- data_IS_norm %>%
  mutate(IS_scaled = .data[[IS_name]] / global_IS_mean) #347 6351 (note new column added)

# drop the IS column (keep IS_scaled if you want QC)
data_IS_norm <- data_IS_norm %>%
  dplyr::select(-all_of(IS_name)) 

IS_scaled_info<-data_IS_norm$IS_scaled

#I'm going to take the opportunity here to remove Feature_14413 and IS_scaled as these don't need to be in my data frames anymore, normalization is complete.
data_IS_norm <- data_IS_norm %>%
  dplyr::select(-all_of(c("Feature_14924", "IS_scaled")))

# final message
message("Internal-standard normalization complete. Normalized dataset: ", nrow(data_IS_norm), " samples × ", length(feature_cols), " features (plus SampleID).") #Sanity check; Internal-standard normalization complete. Normalized dataset: 365 samples × 7129 features (plus SampleID).

###############################################################################
# Step next: Filter features based on median feature intensity 
###############################################################################

#subset to just dry samples
dim(data_IS_norm) #365 x 7129
dim(info_feature_complete) #7129 x 133
dim(metadata) #365 x 10
unique(metadata$Type) # "SAMPLE_POWDER_METAB" "Fresh_Leaves"        "RTE_Decomp"  
unique(metadata$Origin) # "LH1"  "ELK"  "HOK1" "HOK2" "SEK1" "SEK2" 
# maybe useful line filter(Origin %in% c("LH1", "ELK")) %>%

# Get filenames you want to keep
keep_ids <- metadata %>%
  dplyr::filter(Origin %in% c("LH1", "ELK")) %>%
  pull(filename)

# Filter main data
data_IS_norm_powder <- data_IS_norm %>%
  filter(SampleID %in% keep_ids)

# Check dimensions
dim(data_IS_norm_powder) #16 x 7129

# remove column sums = 0 from data_IS_norm_powder and remove it from info_feature_complete
data_IS_norm_powder1<- data_IS_norm_powder %>%
  dplyr::select(
    where(~ !is.numeric(.x) || sum(.x, na.rm = TRUE) > 0)
  )

dim(data_IS_norm_powder1) #16 x 5433, beautiful and seems close to the 5106 that was  pulled for 19 samples of dried A. rubra leaf litter. 

# Compute median intensity per feature
str(data_IS_norm_powder1[1:5,1:5])
rownames(data_IS_norm_powder1)<-data_IS_norm_powder1$SampleID
dim(data_IS_norm_powder1)
data_IS_norm_powder2 <- data_IS_norm_powder1[,2:5433]

#implement prevalence threshold. 

#the design for this experiment is 10 samples per site, but since some did not have enough leaf powder for metabolomics and we had one sample fail during injection - we are at 16 samples. Half of which is 8, where 4 is 25%. Thus, I think we can still use 4 samples as a threshold for noise. 

#for the HFA experiment we had 2 rivers, 4 sites with 5 trees each. n = 20, and used 4 as well, which is basically saying if it isnt present in practically one site, toss it. In this study, since we only had 2 sites: it would be if it isnt present in at least half of the trees on one river (or approx. 25% of the data). 

prevalence_counts <- colSums(data_IS_norm_powder2 > 0)
keep_prev <- prevalence_counts >= 4
data_prev_filtered <- data_IS_norm_powder2[, keep_prev, drop = FALSE]
cat("Features retained (prevalence ≥4):", sum(keep_prev), "\n") #Features retained (prevalence ≥4): 4586

#are many features still near or at zero for median?
feature_medians_prev4 <- apply(data_IS_norm_powder2[, prevalence_counts >= 4], 2, median)
summary(feature_medians_prev4)
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0     3034    18835   269597    81216 60983560

#distribution of max intensity features with median = 0 after pervalence >= 4 filtering
zero_median <- feature_medians_prev4 == 0
zero_max <- apply(
  data_IS_norm_powder2[, prevalence_counts >= 4][, zero_median],
  2,
  max
)

summary(zero_max)
#    Min.   1st Qu.   Median    Mean   3rd Qu.     Max. 
#.  264.7   5768.6   12217.3  29936.9  28902.0  801636.1

quantile(zero_max, c(0.9, 0.95, 0.99))
#. 90%       95%       99% 
#66739.11 109056.25 256436.54

#filter based on feature intensity now. 
intensity_threshold <- 20000  # adjust as desired, set minimum intensity threshold (IS-normalized, not log-transformed)

feature_medians_filtered <- apply(data_prev_filtered, 2, max)  # compute per-feature median (or max) across samples

keep_features_intensity <- feature_medians_filtered > intensity_threshold # Identify features to keep (median above threshold)

data_intensity_filtered <- data_prev_filtered[, keep_features_intensity] # Subset data

cat("Features retained after intensity filter:", sum(keep_features_intensity), "\n") # Summary, Retained 3736 for this data set. 
cat("Features removed:", sum(!keep_features_intensity), "\n") #Removed 850 for this data set. 

summary(apply(data_intensity_filtered, 2, median))
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0     9309    32746   330604     114658     60983560

#Filter annotation matrix to match filtered feature table
str(data_intensity_filtered)
str(info_feature_complete)
kept_features <- colnames(data_intensity_filtered)
length(kept_features)   # should be 3736.

info_feature_filtered_s2 <- info_feature_complete[
  info_feature_complete$Feature %in% kept_features,
]

nrow(info_feature_filtered_s2) #3736

#match order in case it changed somehow
info_feature_filtered_s2 <- info_feature_filtered_s2[
  match(kept_features, info_feature_filtered_s2$Feature),
]

#check
all(info_feature_filtered_s2$Feature == colnames(data_intensity_filtered))
setdiff(kept_features, info_feature_filtered_s2$Feature)

###check superclass retention
table_original <- table(info_feature_complete$Sirius_NPC.superclass)
table_filtered <- table(info_feature_filtered_s2$Sirius_NPC.superclass)

retention_df <- data.frame(
  superclass = names(table_original),
  original_n = as.numeric(table_original),
  retained_n = as.numeric(table_filtered[names(table_original)]),
  retention_prop = as.numeric(table_filtered[names(table_original)]) /
    as.numeric(table_original)
)

retention_df[is.na(retention_df$retained_n), "retained_n"] <- 0
retention_df$retention_prop[is.na(retention_df$retention_prop)] <- 0 #If plant secondary metabolite classes are disproportionately removed, you may need to relax max intensity threshold
str(retention_df)
print(retention_df)
#write.csv(retention_df,paste0(Sys.Date(),"4prev_20KmaxIntensityFiltering.csv"))

#pull out diarylheptanoids and ellagitannins to manually verify
Diaryl_target_features <- info_feature_filtered_s2 %>%
  filter(Sirius_NPC.superclass == "Diarylheptanoids") %>%
  pull(Feature)

head(Diaryl_target_features)
length(Diaryl_target_features) #90

#write.csv(Diaryl_target_features,file="CGM_Diaryl_target_features.csv")

#now do ellagetannnins
parse_mq <- function(x) {
  if (is.na(x) || x == "") return(NA_real_)
  
  # Split on "AND" and convert to numeric
  vals <- str_split(x, " AND ", simplify = TRUE)
  nums <- suppressWarnings(as.numeric(vals))
  
  # Return max of available values
  max(nums, na.rm = TRUE)
}

mass_check <- function(mass_str, min_mass = 600) {
  if (is.na(mass_str)) return(FALSE)
  
  masses <- str_split(mass_str, " AND ", simplify = TRUE)
  nums <- suppressWarnings(as.numeric(masses))
  
  any(nums >= min_mass, na.rm = TRUE)
}

#pull broadly
ellagitannins_df <- info_feature_filtered_s2 %>%
  filter(if_any(everything(),
                ~ grepl("ellag", .x, ignore.case = TRUE)))

#prioritize based on Sirius name
ellag_prioritized <- ellagitannins_df %>%
  mutate(
    sirius_flag = grepl("ellag|ellagic|ellagitannin", Sirius_name, ignore.case = TRUE) |
      grepl("tannin", Sirius_NPC.superclass, ignore.case = TRUE),
    fbmn_flag   = grepl("ellag|ellagic|tannin", FBMN_Compound_Name, ignore.case = TRUE)
  ) %>%
  arrange(desc(sirius_flag), desc(fbmn_flag))

ellag_filtered <- ellag_prioritized %>%
  
  # Clean MQ scores
  mutate(
    CMMC_MQScore_clean = map_dbl(CMMC_MQScore, parse_mq),
    FBMN_MQScore_clean = map_dbl(FBMN_MQScore, parse_mq),
    
    Combined_MQScore = pmax(CMMC_MQScore_clean,
                            FBMN_MQScore_clean,
                            na.rm = TRUE)
  ) %>%
  
  mutate(
    Combined_MQScore = ifelse(!is.finite(Combined_MQScore),
                              NA,
                              Combined_MQScore)
  ) %>%
  
  # Rowwise needed for per-row mass_check
  rowwise() %>%
  filter(mass_check(Sirius_ionMass, min_mass = 300)) %>%
  ungroup() %>%
  
  mutate(
    confidence = case_when(
      !is.na(Combined_MQScore) & Combined_MQScore >= 0.90 ~ "high",
      !is.na(Combined_MQScore) & Combined_MQScore >= 0.70 ~ "medium",
      TRUE ~ "low"
    )
  ) %>%
  
  arrange(desc(Combined_MQScore))

dim(ellag_filtered) #not sure how useful this one is. 
dim(ellag_prioritized) #95, number should drop greatly after verification

#write.csv(ellag_prioritized, file= "CGM_ellag_prioritized.csv")

#other useful table
library(dplyr)

#feature verification and filtering happens here.
diaryl_features<-read.csv(file="CGM_Diaryl_target_features.csv",header=F)
diaryl_keeps<- diaryl_features %>% 
  filter(V2 == "1")
diaryl_keeps1<-diaryl_keeps$V1
length(diaryl_keeps1) #59

powder_df<-data_intensity_filtered #contains SampleID and Origin
length(Diaryl_target_features) #90
length(diaryl_keeps1) #59
length(intersect(diaryl_keeps1, Diaryl_target_features)) #59
newft_tormv<- setdiff(Diaryl_target_features, diaryl_keeps1) #the object order matters here
length(newft_tormv) #31

# Remove newft_tormv (the 38 rejected features) from both objects
data_intensity_filtered_clean <- data_intensity_filtered[
  , !colnames(data_intensity_filtered) %in% newft_tormv]

info_feature_filtered_s2_v1 <- info_feature_filtered_s2[
  !info_feature_filtered_s2$Feature %in% newft_tormv, ]

# not_diaryl <- dry_PLSDA_VIP.df_diaryl$ID[
#   is.na(dry_PLSDA_VIP.df_diaryl$Sirius_NPC.superclass) | 
#     dry_PLSDA_VIP.df_diaryl$Sirius_NPC.superclass != "Diarylheptanoids"] #for this step the diaryl plsda vip needs to be set to 0.0 to retain all feature numbers for median value and all features PLSDA. This is so clunky sorry future Jonathan. 

#rename ft from not_diaryl
#info_feature_filtered_s2_v1$Sirius_NPC.superclass[info_feature_filtered_s2_v1$Feature %in% not_diaryl] <- "Diarylheptanoids"

#now do the ellagitannins and ellagitannin derivatives (i.e., ellagic acid and HHDP frags)
ellagic_data<-read.csv(file = "candidate_ellagitannins.csv") 
ellagic_keeps<- ellagic_data %>% 
  filter(Keep_hypothesis == "1") #42
ellagic_keeps1<-ellagic_keeps$Feature
length(ellagic_keeps1) #42

length(intersect(ellagic_keeps1, ellag_prioritized$Feature)) #42
misID_tormv<- setdiff(ellag_prioritized$Feature, ellagic_keeps1) #the object order matters here
length(misID_tormv) #53

data_intensity_filtered_clean_v2 <- data_intensity_filtered_clean[
  , !colnames(data_intensity_filtered_clean) %in% misID_tormv]

info_feature_filtered_s2_v2 <- info_feature_filtered_s2_v1[
  !info_feature_filtered_s2_v1$Feature %in% misID_tormv, ]

############################################################
############################################################
############################################################

# 1. Gather feature intensities into long format
data_long <- data_intensity_filtered_clean_v2 %>%
  # Keep SampleID if it exists; else assume just numeric columns
  tibble::rownames_to_column(var = "SampleID") %>%
  tidyr::pivot_longer(
    cols = -SampleID,
    names_to = "Feature",
    values_to = "Intensity"
  )

# 2. Join with annotation table
data_annotated <- data_long %>%
  left_join(info_feature_filtered_s2_v2 %>% 
              dplyr::select(Feature, Sirius_NPC.superclass),
            by = "Feature")

# 3. Summarize number of features and median intensity per superclass
summary_by_superclass_simple <- data_annotated %>%
  filter(!is.na(Sirius_NPC.superclass) & Sirius_NPC.superclass != "") %>%
  
  # First: compute per-feature mean & median intensity
  group_by(Sirius_NPC.superclass, Feature) %>%
  summarize(
    mean_feature_intensity   = mean(Intensity, na.rm = TRUE),
    median_feature_intensity = median(Intensity, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  
  # Then: summarize across features within superclass
  group_by(Sirius_NPC.superclass) %>%
  summarize(
    n_features = n(),
    mean_intensity   = mean(mean_feature_intensity, na.rm = TRUE),
    median_intensity = median(median_feature_intensity, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  
  arrange(desc(n_features))
dim(summary_by_superclass_simple) #59 x 4
print(summary_by_superclass_simple,n=59)
sum(summary_by_superclass_simple$n_features) #3323
#write.csv(summary_by_superclass_simple,file="10June26_4prev_20Kmax_summary_by_superclass.csv")

###############################################################################
# Step 4: Prepare sample matrix
###############################################################################

# Identify features with missing or empty superclass
features_missing_superclass <- info_feature_filtered_s2_v2 %>%
  filter(is.na(Sirius_NPC.superclass) | Sirius_NPC.superclass == "") %>%
  pull(Feature) #329

# Filter info_feature_filtered_s2 to keep only annotated features
info_feature_filtered_s2_clean <- info_feature_filtered_s2_v2 %>%
  filter(!Feature %in% features_missing_superclass) #3323 x 133 

# Filter data_intensity_filtered to remove those features (columns)
data_intensity_filtered_clean <- data_intensity_filtered_clean_v2 %>%
  dplyr::select(-all_of(features_missing_superclass)) #16 x 3323

#filter out superclasses that are not PSMs given analyses from Dickey et al. Nat Commun
unique(info_feature_filtered_s2_clean$Sirius_NPC.superclass)

#"Nucleosides", "Small peptides", "Fatty Acids and Conjugates", "Saccharides", "Oligopeptides",  "β-lactams", "Aminosugars and aminoglycosides" , "Fatty acyls", "Fatty esters" , "Macrolides", "Glycerophospholipids", "Polyols", "Spingolipids", "Terphenyls" ,"Eicosanoids" , "Glycerolipids", "Steroids", "Polyethers", "Tetramate alkaloids"

#"Tropolones" (not quite convinced, seems pretty rare tbh. Reported in Cyprus)                      

psm_classes <- c("Pseudoalkaloids (transamidation)", "Phenolic acids (C6-C1)", "Histidine alkaloids" , "Coumarins" , "Nicotinic acid alkaloids", "Amino acid glycosides", "Anthranilic acid alkaloids", "Ornithine alkaloids" ,   "Monoterpenoids" ,"Chromanes" , "Fatty amides" , "Cyclic polyketides", "Peptide alkaloids","Phenylpropanoids (C6-C3)", "Flavonoids" ,"Naphthalenes", "Lignans", "Lysine alkaloids" , "Linear polyketides", "Tyrosine alkaloids", "Guanidine alkaloids", "Fatty acyl glycosides",  "Tryptophan alkaloids" ,"Isoflavonoids", "Styrylpyrones" ,"Phloroglucinols" ,"Polycyclic aromatic polyketides", "Diterpenoids" ,"Sesquiterpenoids" ,"Diarylheptanoids" ,"Apocarotenoids" ,"Phenylethanoids (C6-C2)" ,"Octadecanoids" ,  "Aromatic polyketides" ,"Triterpenoids", "Stilbenoids", "Carotenoids (C40)", "Sesterterpenoids" ,"Meroterpenoids")

info_feature_filtered_s2_clean_2 <- info_feature_filtered_s2_clean %>%
  mutate(metabolite_type = ifelse(Sirius_NPC.superclass %in% psm_classes, "PSM", "Other"))

info_feature_filtered_s2_psms <- info_feature_filtered_s2_clean_2 %>%
  filter(metabolite_type %in% "PSM") #1871 x 134

fts_psm<-info_feature_filtered_s2_psms$Feature
length(fts_psm) #1871

data_intensity_filtered_clean_psm <- data_intensity_filtered_clean %>%
  dplyr::select(all_of(fts_psm)) #16 x 1871

#save point
#save(list=ls(),file="workingdir.rda")
#setwd("~/Documents/UCSD_EBE_PostDoc/Mentee_Materials/Dahlia_Loomis/Revisions/MetaBolomics/Files/RSave")
#load(file="workingdir.rda")

#renaming
data_numeric <- data_intensity_filtered_clean_psm

# Log-transform to stabilize variance
data_log <- log(data_numeric + 1)

# Re-add SampleID for metadata joining
data_log$SampleID <- rownames(data_numeric)

# Join with metadata
data_log.df <- data_log %>%
  left_join(metadata_metabolomics, by = "SampleID") #16 x  (joined up with metadata)

dim(data_log.df) #16 x 1883
unique(data_log.df$Origin)#Type: SAMPLE_POWDER_METAB" and Origin: LH1 anfd ELK
unique(data_log.df$Type)

#redundant steps but good for renaming. 
dry_pca.df <- data_log %>%
  left_join(metadata_metabolomics, by = "SampleID") %>%
  filter(Type == "SAMPLE_POWDER_METAB") %>%
  dplyr::select(SampleID, where(is.numeric)) #16 x 1873 (numeric plus 2 columns: sampleid and sample order)

colnames(dry_pca.df)[1865:1872] #Sample_Order
colnames(dry_pca.df)[1:5] #SampleID

dry_is.log.df <- data_log.df %>%
  filter(Type == "SAMPLE_POWDER_METAB") %>%
  dplyr::select(SampleID, Origin, where(is.numeric)) #16 x 1874

colnames(dry_is.log.df)[1865:1873] #Sample_Order
colnames(dry_is.log.df)[1:5] #SampleID, Origin

dim(dry_pca.df) #16 1873

#remove columns that sum to zero. 
dry_pca.df1 <- dry_pca.df %>%
  dplyr::select(
    where(~ !is.numeric(.x) || sum(.x, na.rm = TRUE) > 0)
  ) #16 x 1873

colnames(dry_pca.df1)[1865:1872] #Sample_order
colnames(dry_pca.df1)[1:5] #SampleID

dry_pca.df1[1:5,1:5]
dry_pca.df1[1:5,1865:1872]

dry_is.log.df1 <- dry_is.log.df %>%
  dplyr::select(
    where(~ !is.numeric(.x) || sum(.x, na.rm = TRUE) > 0)
  )

dry_is.log.df1[1:5,1:5] #sampleID and origin here 
dry_is.log.df1[1:5,1865:1873] #sample_order at 3326

dim(dry_pca.df1) #16 x 1873
rowSumsDryAllft<-rowSums(dry_pca.df1[,2:1871]) #remove sample-order and SampleID

dim(dry_is.log.df1) #16 x 1873; contains Sample ID and Origin and sample_order

dry_PLSDA_sample_ids <- dry_is.log.df$SampleID

# Near-zero variance removal (keep feature names)
dry_PLSDA_nzv <- caret::nearZeroVar(dry_is.log.df1 %>% dplyr::select(-SampleID, -Origin,-Sample_Order), names = TRUE)
message("Removing ", length(dry_PLSDA_nzv), " near-zero variance features.") #removing 0 features (theyre removed above, this is just done for renaming and verification).

#Run PLS-DA with Origin as the Y factor
pls_X <- dry_is.log.df1 %>%
  dplyr::select(-SampleID, -Origin, -Sample_Order) %>%
  dplyr::select(-all_of(dry_PLSDA_nzv)) #16 x 1871

pls_X[1:5,1:5]
pls_X[1:5,1865:1870]

pls_Y <- factor(dry_is.log.df1$Origin)

# final check: rows match length of Y
if(nrow(pls_X) != length(pls_Y)) stop("Mismatch between X rows and Y length.")

dry_PLSDA <- mixOmics::plsda(
  X = pls_X,
  Y = pls_Y,
  ncomp = 2,
  scale = TRUE #running PLSDA right here 
)

#Extract PLS-DA scores and reattach metadata
dry_PLSDA_scores <- as.data.frame(dry_PLSDA$variates$X) %>%
  rownames_to_column("tmp_rowname") %>%    # mixOmics rownames map to original row order
  mutate(SampleID = dry_PLSDA_sample_ids,
         Origin = as.character(pls_Y)) %>%
  left_join(metadata_metabolomics, by = "SampleID")

#Outlier detection (Mahalanobis) on comp1 & comp2 — identify, but do NOT drop
plsda_mahal <- mahalanobis(dry_PLSDA_scores[, c("comp1", "comp2")],
                           colMeans(dry_PLSDA_scores[, c("comp1", "comp2")]),
                           cov(dry_PLSDA_scores[, c("comp1", "comp2")]))

threshold <- qchisq(0.90, df = 4) #re-examine the degrees of freedom here. Should it be 2?
outlier_idx <- which(plsda_mahal > threshold)
message("Identified ", length(outlier_idx), " PLS-DA outlier(s) by Mahalanobis (90% cutoff).") #0 outlier

# custom site colors
site_colors <- c("black","#C4CBCA")
dry_PLSDA_scores$Origin.x[dry_PLSDA_scores$Origin.x=="LH1"]<-"Little Hoko River"
dry_PLSDA_scores$Origin.x[dry_PLSDA_scores$Origin.x=="ELK"]<-"Pysht River"
dry_PLSDA_scores$Origin.x<-factor(dry_PLSDA_scores$Origin.x,levels=c("Little Hoko River","Pysht River"))

dry_PLSDA_scores$Origin.x<-as.character(dry_PLSDA_scores$Origin.x)
dry_PLSDA_scores$Origin.x[dry_PLSDA_scores$Origin.x=="Little Hoko River"]<-"Non-local"
dry_PLSDA_scores$Origin.x[dry_PLSDA_scores$Origin.x=="Pysht River"]<-"Local"
dry_PLSDA_scores$Origin.x<-factor(dry_PLSDA_scores$Origin.x,levels=c("Local","Non-local"))

# ---- Build Plot ----
leaf_PLSDA.plot <- ggplot(dry_PLSDA_scores,
                          aes(x = comp1, y = comp2)) +
  scale_x_continuous(limits=c(-30,30),breaks=c(-30, -15, 0, 15, 30))+
  scale_y_continuous(limits=c(-30,30),breaks=c(-30, -15, 0, 15, 30))+
  
  # ELLIPSES (68% CI)
  stat_ellipse(
    aes(group = Origin.x, color = Origin.x),
    type = "t",
    level = 0.68,
    linewidth = 0.8,
    linetype = 1,
    show.legend = FALSE
  ) +
  
  # POINTS
  geom_point(
    aes(fill = Origin.x, shape = Origin.x),
    size = 3.5,
    alpha = 0.9,
    color = "#119822",
    stroke = 2
  ) +
  
  # MANUAL SCALES
  scale_fill_manual(name = "Origin River", values = site_colors) +
  scale_color_manual(name = "Origin River", values = site_colors) +
  scale_shape_manual(name = "Origin River", values = c("Non-local" = 21, "Local" = 21)
  ) +
  
  # THEMING
  theme_bw(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
  ) +
  
  # LABELS
  labs(
    x = paste0("PLS-DA Comp 1 (", round(dry_PLSDA$prop_expl_var$X[2] * 100, 1), "%)"),
    y = paste0("PLS-DA Comp 2 (", round(dry_PLSDA$prop_expl_var$X[1] * 100, 1), "%)")
  ) #all features dried leaves

leaf_PLSDA.plot

#metrics
dry_PLSDA$prop_expl_var$X 
dry_PLSDA$prop_expl_var$Y 
cumsum(dry_PLSDA$prop_expl_var$X) 
cumsum(dry_PLSDA$prop_expl_var$Y) 
table(dry_PLSDA$Y)

perf_PLSDA <- perf(dry_PLSDA, 
                   validation = "Mfold",  
                   folds = 3,             # adjust to what you used
                   nrepeat = 100,   
                   progressBar = TRUE)

cumsum(dry_PLSDA$prop_expl_var$X) #cumsum comp2
perf_PLSDA$error.rate$BER #max.dist comp2
perf_PLSDA$error.rate$overall #max.dist comp2

#get the loadings
dry_PLSDA_loadings <- mixOmics::plotLoadings(
  dry_PLSDA, 
  plot = FALSE, 
  contrib = "max")

#convert into a df
dry_PLSDA_loadings.df <- as.data.frame(dry_PLSDA_loadings$X) %>%
  rownames_to_column("Feature") %>%
  dplyr::select(Feature, GroupContrib)

#generate VIP scores
dry_PLSDA_VIP <- as.data.frame(mixOmics::vip(dry_PLSDA))

#filter for scores > 1.2 (anything > 1 is contributing to group separation)
dry_PLSDA_filtered <- dplyr::filter(dry_PLSDA_VIP, dry_PLSDA_VIP$comp1 > 1.0) #set threshold to zero for main figure

dry_PLSDA_filtered$ID <- rownames(dry_PLSDA_filtered)

#select comp 1 (similar to PC axis for PCA) and comp 2
dry_PLSDA_select <- dry_PLSDA_filtered %>% dplyr::select(ID, comp1, comp2)

#create new VIP df
dry_PLSDA_VIP.df <- dry_PLSDA_select %>% 
  left_join(dry_PLSDA_loadings.df, by = c("ID" = "Feature")) %>%
  left_join(info_feature_filtered_s2_clean, by = c("ID" = "Feature")) %>% 
  arrange(desc(comp1))

dim(dry_PLSDA_VIP.df) #VIP > 0: 1871 x 136; VIP > 1: 684 x 136

#write.csv(dry_PLSDA_VIP.df,file="leaf_PLSDA_VIP.df.csv")

#filter for the top 20 VIP features
dry_PLSDA_20_VIP <- dry_PLSDA_VIP.df %>%
  head(20)

#filter NA / aggregate with annotations / arrange by median VIP
dry_PLSDA_aggregated_data <- dry_PLSDA_VIP.df %>%
  filter(!is.na(Sirius_NPC.superclass)) %>%
  dplyr::group_by(Sirius_NPC.superclass) %>%
  dplyr::summarise(median_VIP = median(comp1, na.rm = TRUE)) %>%
  arrange(desc(median_VIP))

#plot median VIP score by superclass
dry_PLSDA_VIP.plot <- dry_PLSDA_aggregated_data %>%
  ggplot(aes(x = reorder(Sirius_NPC.superclass, median_VIP), y = median_VIP)) +
  geom_point(size = 3, color = "black") + 
  coord_flip() +  
  labs() +
  xlab("") +
  ylab("Median VIP Score") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 30, face = "bold"),  
    axis.title.x = element_text(size = 20),  
    axis.title.y = element_text(size = 20),  
    axis.text.x = element_text(size = 12),  
    axis.text.y = element_text(size = 12)
  )

dry_PLSDA_VIP.plot

#VIP figure for total vips in pls-da. 
vip_df <- dry_PLSDA_VIP %>%
  as.data.frame() %>%
  rownames_to_column("FeatureID")

vip_df <- vip_df %>%
  left_join(dry_PLSDA_VIP.df, by = c("FeatureID" = "ID")) #OK for this figure I went up to the top and changed threshold to 0.0, but for other analyses it will need to be re-set to 1.2. DONT FORGET. 

dim(vip_df)
colnames(vip_df)

vip_superclass <- vip_df %>%
  filter(!is.na(Sirius_NPC.superclass)) %>%
  dplyr::group_by(Sirius_NPC.superclass) %>%
  dplyr::summarise(median_VIP = median(comp1.x, na.rm = TRUE),n_features = n()) %>%
  mutate(metabolite_type = ifelse(Sirius_NPC.superclass %in% psm_classes, "PSM", "Other"))

dim(vip_superclass)
colnames(vip_superclass)
vip_superclass$n_features

vip_superclass <- vip_superclass %>%
  mutate(
    superclass_label = paste0(Sirius_NPC.superclass, " (", n_features, ")"))

# Assign green to PSMs, black to others
label_colors <- setNames(
  ifelse(vip_superclass$metabolite_type == "PSM", "black", "black"),
  vip_superclass$superclass_label
)

vip_plot <- vip_superclass %>%
  ggplot(aes(
    x = reorder(superclass_label, median_VIP),
    y = median_VIP,
    fill = metabolite_type
  )) +
  
  geom_point(
    size = 4.2,
    shape = 21,
    color = "black",
    stroke = 0.8
  ) +
  
  coord_flip() +
  
  geom_hline(
    yintercept = 1.0,
    linetype = "dashed",
    linewidth = 0.7
  ) +
  
  scale_fill_manual(
    values = c(
      "PSM" = "#2bc016",
      "Other" = "grey65"
    )
  ) +
  
  scale_y_continuous(
    limits = c(0, 1.75),
    breaks = seq(0, 1.75, 0.25)
  ) +
  
  labs(
    y = "Median VIP score",
    x = "",
    fill = "Metabolite class"
  ) +
  
  theme_bw(base_size = 14) +
  theme(
    # Axis text color for y-axis based on PSM status
    axis.text.y = element_text(
      size = 12,
      color = label_colors[levels(reorder(vip_superclass$superclass_label, vip_superclass$median_VIP))]
    ),
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    
    # Light horizontal gridlines
    panel.grid.major.y = element_line(color = "grey85", linewidth = 0.5),
    panel.grid.minor.y = element_line(color = "grey85", linewidth = 0.25)
  )

vip_plot

#test 
vip_superclass <- dry_PLSDA_VIP.df %>%
  filter(!is.na(Sirius_NPC.superclass)) %>%
  group_by(Sirius_NPC.superclass) %>%
  summarise(
    median_VIP = median(comp1, na.rm = TRUE),
    min_VIP = min(comp1, na.rm = TRUE),
    max_VIP = max(comp1, na.rm = TRUE),
    n_features = n(),
    n_vip1 = sum(comp1 > 1, na.rm = TRUE)
  )

vip_superclass <- vip_superclass %>%
  mutate(superclass_label = paste0(Sirius_NPC.superclass, " (", n_features, ")"))

vip_superclass <- vip_superclass %>%
  mutate(metabolite_type = ifelse(Sirius_NPC.superclass %in% psm_classes, "PSM", "Other"))

vip_plot <- vip_superclass %>%
  ggplot(aes(
    x = reorder(superclass_label, median_VIP),
    y = median_VIP,
    fill = metabolite_type
  )) +
  
  # VIP range within superclass
  geom_linerange(
    aes(ymin = min_VIP, ymax = max_VIP),
    linewidth = 0.8,
    color = "grey40"
  ) +
  
  # Median VIP point
  geom_point(
    size = 4,
    shape = 21,
    color = "black",
    stroke = 0.6
  ) +
  
  coord_flip() +
  
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    linewidth = 0.7
  ) +
  
  scale_fill_manual(
    values = c(
      "PSM" = "forestgreen",
      "Other" = "grey70"
    )
  ) +
  
  scale_y_continuous(
    limits = c(0, 3.5),
    breaks = seq(0, 3.5, 0.25)
  ) +
  
  labs(
    y = "Median VIP score",
    x = "",
    fill = "Metabolite class"
  ) +
  theme_bw(base_size = 14) +
  theme(
    # Axis text color for y-axis based on PSM status
    axis.text.y = element_text(
      size = 12,
      color = label_colors[levels(reorder(vip_superclass$superclass_label, vip_superclass$median_VIP))]
    ),
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    
    # Light horizontal gridlines
    panel.grid.major.y = element_line(color = "grey85", linewidth = 0.5),
    panel.grid.minor = element_blank()
  )

vip_plot

vip_superclass <- vip_superclass %>%
  mutate(prop_vip_gt1 = n_vip1 / n_features)

#vip_plot <- vip_superclass %>%
vip_plot1 <- vip_superclass %>%
  ggplot(aes(
    x = reorder(superclass_label, median_VIP),
    y = median_VIP,
    fill = metabolite_type
  )) +
  
  # VIP >1 proportion bar (scaled for visual alignment)
  geom_col(
    aes(y = prop_vip_gt1),
    width = 0.65,
    alpha = 0.40,
    color = NA
  ) +
  
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    linewidth = 0.7
  ) +
  
  # Median VIP point — size differentiated by metabolite type
  geom_point(
    aes(size = metabolite_type),
    shape = 21,
    color = "black",
    stroke = 0.8
  ) +
  
  coord_flip(clip = "off") +
  
  scale_fill_manual(
    values = c(
      "PSM" = "#2bc016",
      "Other" = "grey70"
    )
  ) +
  
  scale_size_manual(
    values = c(
      "PSM" = 4.5,
      "Other" = 4.5
    )
  ) +
  
  scale_y_continuous(
    limits = c(0, 1.75),
    breaks = seq(0, 1.75, 0.25),
    expand = expansion(mult = c(0.02, 0.02)),
    sec.axis = sec_axis(
      ~ . ,
      name = "Proportion VIP > 1",
      breaks = c(0, 0.25, 0.5, 0.75, 1.0),
      labels = c("0", "0.25", "0.50", "0.75", "1.0")
    )
  ) +
  
  labs(
    y = "Median VIP score",
    x = "Superclass (SIRIUS)",
    fill = "Metabolite class"
  ) +
  
  theme_bw(base_size = 14) +
  
  theme(
    legend.position = "none",
    axis.text.y = element_text(
      size = 12,
      color = label_colors[
        levels(reorder(vip_superclass$superclass_label, vip_superclass$median_VIP))
      ]
    ),
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.x.top = element_text(size = 12, color = "grey50",hjust = 0.30),
    
    panel.grid.major.y = element_line(color = "grey85", linewidth = 0.5),
    panel.grid.minor.y = element_line(color = "grey85", linewidth = 0.25)
  )

vip_plot1 #ok now go back up and threshold to >1

#save(list=ls(),file="Apr30_2026_CaddieMetabolomics_analyses.rda")

#filter
all_vips<-dry_PLSDA_VIP.df$ID
all_vip_table <- pls_X[,colnames(pls_X) %in% all_vips]
dim(all_vip_table) #16 x 684

#rowsum intensities
all_VIPintensities <- rowSums(all_vip_table) 

#getting a weighted score instead
vip_features <- intersect(colnames(all_vip_table), dry_PLSDA_VIP.df$ID)

length(vip_features) #45
vip_matrix <- all_vip_table[, vip_features]

vip_df <- dry_PLSDA_VIP.df[match(vip_features, dry_PLSDA_VIP.df$ID), ]

vip_weights <- vip_df$comp1 / sum(vip_df$comp1) #weighted vip scores. 
vip_weights2 <- vip_df$comp2 / sum(vip_df$comp1) #weighted vip scores. 
length(vip_weights)

str(vip_matrix)
vip_matrix <- as.matrix(vip_matrix)
str(vip_weights)

#Per sample VIP-weighted Score
vip_score <- as.vector(vip_matrix %*% vip_weights) #A weighted sum of all VIP features, Essentially: “how much does this sample express the features that most drive group separation?”
vip_score_comp2 <- as.vector(vip_matrix %*% vip_weights2)
vip_total_weighted_score<-vip_score+vip_score_comp2

meta_y<-dry_is.log.df[,1:2]

meta_xy_all<-cbind(meta_y,all_VIPintensities,vip_score,vip_score_comp2,vip_total_weighted_score)
str(meta_xy_all)

hist(meta_xy_all$all_VIPintensities)
hist(meta_xy_all$vip_score)
hist(vip_score_comp2)
hist(vip_total_weighted_score)

meta_xy_all$Origin[meta_xy_all$Origin=="LH1"]<-"Little Hoko River"
meta_xy_all$Origin[meta_xy_all$Origin=="ELK"]<-"Pysht River"
meta_xy_all$Origin<-factor(meta_xy_all$Origin,levels=c("Little Hoko River","Pysht River"))

Anova(lm(meta_xy_all$all_VIPintensities~meta_xy_all$Origin)) #summed intensities of VIPs
Anova(lm(meta_xy_all$vip_score~meta_xy_all$Origin)) #weighted summed VIP scores.
Anova(lm(meta_xy_all$vip_score_comp2~meta_xy_all$Origin))
Anova(lm(meta_xy_all$vip_total_weighted_score~meta_xy_all$Origin))

vp1<-aov(meta_xy_all$all_VIPintensities~meta_xy_all$Origin)
vp1b<-TukeyHSD(vp1)
vp1b

multcompLetters4(vp1, vp1b)
#$`meta_xy$Origin`
#Pysht River Little Hoko River 
#"a"               "a" 

all_relative_abundVIP<-ggplot(data=meta_xy_all, aes(x=Origin,y=all_VIPintensities,fill=Origin,shape=Origin))+
  geom_boxplot(aes(fill=Origin),alpha=0.5, outlier.size=0) +
  geom_point(aes(group=Origin, fill=Origin), color = "#119822", stroke = 1.5, size=3.0, position=position_jitterdodge(jitter.width=1.25)) +
  #facet_grid(~ Origin_River)+
  scale_shape_manual(values = c(21,21,21,21,21)) +
  scale_fill_manual(values = c(site_colors)) + # Boxplot fill color
  # scale_y_continuous(label=scales::comma,limits=c(400,500),breaks=c(400,420,440,460,480,500))+
  theme_classic() +
  labs(y="Cumulative Feature Abundance\n of all PSMs",x="Origin Site") +
  theme(legend.position = "none")

all_relative_abundVIP_final<-all_relative_abundVIP + theme(axis.text.x = element_text(size = 12)) + theme(text = element_text(size = 15),legend.text.align = 0)

all_relative_abundVIP_final

########################################################################
########################################################################
# Heat Map
########################################################################
########################################################################
library(dplyr)
library(tidyr)
library(pheatmap)
library(RColorBrewer)

#aggregate at the SuperClass Level

dry_is.log.df1[1:5,1:5] #sampleID and origin here 
dry_is.log.df1[1:5,1865:1873] #sample order here
dim(dry_is.log.df1)

SampleID<-dry_is.log.df1$SampleID
dry_is.log.df3<-dry_is.log.df1[,3:1872]
dry_is.log.df4<-cbind(SampleID,dry_is.log.df3)
dim(dry_is.log.df4)

# long format (you already did this earlier)
vip_features <- rownames(dry_PLSDA_filtered) #go back up to dry_PLSDA_filtered and change from 0.0 to 1.0
dry_is.log.df4_vips <- dry_is.log.df4[, colnames(dry_is.log.df4) %in% vip_features]
dry_is.log.df5_vips<-cbind(SampleID,dry_is.log.df4_vips)

data_long <- dry_is.log.df5_vips %>% #change to dry_is.log.df4_vips or dry_is.log.df4
  pivot_longer(
    cols = -SampleID,
    names_to = "Feature",
    values_to = "Intensity"
  )

# join annotation
data_annotated <- data_long %>%
  left_join(
    info_feature_filtered_s2_psms %>%
      dplyr::select(Feature, Sirius_NPC.superclass),
    by = "Feature"
  )

rownames(info_feature_filtered_s2_psms)<-info_feature_filtered_s2_psms$Feature
info_PSM_vips<- info_feature_filtered_s2_psms[rownames(info_feature_filtered_s2_psms) %in% vip_features, ]
info_PSM_vips2<-info_PSM_vips %>% dplyr::select(Feature, mz, RT, Sirius_NPC.superclass)
has_rownames(info_PSM_vips2)
info_PSM_vips2<-remove_rownames(info_PSM_vips2)

# collapse to superclass per sample (sum or mean — sum is typical)
data_superclass <- data_annotated %>%
  group_by(SampleID, Sirius_NPC.superclass) %>%
  summarize(
    Intensity = sum(Intensity, na.rm = TRUE),
    .groups = "drop"
  )

#pivot to wide
heatmap_mat <- data_superclass %>%
  pivot_wider(
    names_from = Sirius_NPC.superclass,
    values_from = Intensity,
    values_fill = 0
  ) %>%
  column_to_rownames("SampleID") %>%
  as.matrix()

metadata_ordered <- metadata_metabolomics %>%
  dplyr::filter(SampleID %in% rownames(heatmap_mat)) %>%
  dplyr::arrange(match(SampleID, rownames(heatmap_mat)))

annotation_df <- data.frame(
  Origin = metadata_ordered$Origin
)

rownames(annotation_df) <- metadata_ordered$SampleID

my_colors <- colorRampPalette(c("navy", "white", "firebrick3"))(100)

pheatmap(
  heatmap_mat,
  annotation_row = annotation_df,
  clustering_method = "ward.D2",
  show_rownames = FALSE,
  fontsize_col = 8,
  color = my_colors
)

###### more complex heat map ######
library(microViz) #v. 0.13.1
library(ComplexHeatmap)

#build obj 
otu_mat <- dry_is.log.df5_vips %>%
  column_to_rownames("SampleID") %>%   # skip if SampleID already rownames
  as.matrix() %>%
  t()  # now features x samples

otu <- otu_table(otu_mat, taxa_are_rows = TRUE)

tax_mat <- info_PSM_vips2 %>%
  tibble::column_to_rownames("Feature") %>%
  dplyr::select(Sirius_NPC.superclass) %>%
  as.matrix()

stopifnot(all(rownames(otu_mat) %in% rownames(tax_mat)))
tax_tab <- tax_table(tax_mat)

samp <- metadata_ordered %>%
  column_to_rownames("SampleID") %>%
  sample_data()

ps_obj <- phyloseq(otu, tax_tab, samp)

ps_obj@sam_data$Origin[ps_obj@sam_data$Origin=="LH1"]<-"Little Hoko River"
ps_obj@sam_data$Origin[ps_obj@sam_data$Origin=="ELK"]<-"Pysht River"
ps_obj@sam_data$Origin<- factor(ps_obj@sam_data$Origin, levels=c("Little Hoko River","Pysht River"))

#aggregate and transform
ps_scaled <- ps_obj %>%
  tax_agg(rank = "Sirius_NPC.superclass") %>%
  tax_transform("identity") %>%
  tax_scale() 

ht <- ps_scaled %>%
  comp_heatmap(
    sample_seriation = "OLO_ward",
    sample_ser_dist = "euclidean",
    tax_seriation = "OLO_ward",
    tax_ser_dist = "euclidean",
    colors = heat_palette(palette = my_colors, sym = TRUE),
    sample_names_show = FALSE,
    name = "z-score",
    row_split = 4,
    column_split = 2
  ) #works

ht %>% ComplexHeatmap::draw()

args(microViz::anno_sample)
args(microViz::sampleAnnotation)

origin_cols <- c(
  "Little Hoko River" = "#C4CBCA",
  "Pysht River" = "black"
)

my_colors <- colorRampPalette(c("#8B4513", "white", "#2D5F2D"))(100)
my_colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, "BrBG")))(100)
my_colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, "PuOr")))(100)

ht <- ps_scaled %>%
  comp_heatmap(
    sample_anno = sampleAnnotation(
      Origin = anno_sample_cat("Origin", col = origin_cols, legend_title = "Origin")
    ),
    sample_seriation = "OLO_ward",
    sample_ser_dist = "euclidean",
    tax_seriation = "OLO_ward",
    tax_ser_dist = "euclidean",
    taxa_side = "right",
    sample_side = "top",
    colors = heat_palette(palette = my_colors, sym = TRUE, range = c(-2.5, 2.5)),
    sample_names_show = FALSE,
    name = "z-score",
    row_split = 3,
    column_split = 2,
    row_names_gp = grid::gpar(fontsize = 13),
    heatmap_legend_param = list(
      title_gp = grid::gpar(fontsize = 14, fontface = "bold"),
      labels_gp = grid::gpar(fontsize = 12)
    ))

hm.plot<- ht %>% ComplexHeatmap::draw(annotation_legend_list = attr(ht, "AnnoLegends"), merge_legends = TRUE)

#sample order
ht_drawn <- ComplexHeatmap::draw(
  ht,
  annotation_legend_list = attr(ht, "AnnoLegends"),
  merge_legends = TRUE
)

# column order for each split
col_order <- ComplexHeatmap::column_order(ht_drawn)

col_order

sample_names <- rownames(otu_table(ps_scaled))

meta <- data.frame(sample_data(ps_scaled))

cluster_samples <- lapply(col_order, function(x) sample_names[x])

lapply(cluster_samples, function(x) {
  table(meta[x, "Origin"])
})

lapply(col_order, function(x) sample_names[x])

for(i in seq_along(cluster_samples)) {
  cat("\nCluster", i, "\n")
  print(data.frame(
    Sample = cluster_samples[[i]],
    Origin = meta[cluster_samples[[i]], "Origin"],
    Tree_ID = meta[cluster_samples[[i]],"Tree_ID"]
  ))
}

#######################################################
#######################################################
#####Diarylheptanoid specific analyses
#######################################################
#######################################################
Candidate_diaryl_features <- info_feature_filtered_s2_psms %>%
  filter(Sirius_NPC.superclass %in% c("Diarylheptanoids")) %>%
  pull(Feature)

pls_XSacc <- dry_is.log.df1 %>% dplyr::select(
  # keep all non-numeric columns
  where(~ !is.numeric(.x)),
  # keep only numeric features that are in your feature list
  all_of(Candidate_diaryl_features)
) #16 x 61 ? what are the two extra cols ?

pls_XSacc[1:5,55:61]
pls_XSacc[1:5,1:5] #sampleID and Origin

nzv_idx_sugar <- caret::nearZeroVar(pls_XSacc %>% dplyr::select(-SampleID, -Origin), names = TRUE)

pls_X_filteredSacc <- pls_XSacc %>%
  dplyr::select(-SampleID, -Origin) %>%
  dplyr::select(-all_of(nzv_idx_sugar)) #16 x 59

pls_YSacc <- factor(dry_is.log.df1$Origin)

if (nrow(pls_X_filteredSacc) != length(pls_YSacc)) stop("Mismatch between X rows and Y length.")

dry_PLSDA_sugar <- mixOmics::plsda(
  X = pls_X_filteredSacc,
  Y = pls_YSacc,
  ncomp = 2,
  scale = TRUE
)

perf_PLSDA_diaryl <- perf(dry_PLSDA_sugar, 
                          validation = "Mfold",  
                          folds = 3,             # adjust to what you used
                          nrepeat = 100,   
                          progressBar = TRUE)


cumsum(dry_PLSDA_sugar$prop_expl_var$X) #cumsum comp2
perf_PLSDA_diaryl$error.rate$BER #max.dist comp2
perf_PLSDA_diaryl$error.rate$overall #max.dist comp2

#Extract PLS-DA scores and reattach metadata
diaryl_PLSDA_scores <- as.data.frame(dry_PLSDA_sugar$variates$X) %>%
  rownames_to_column("tmp_rowname") %>%    # mixOmics rownames map to original row order
  mutate(SampleID = pls_YSacc,
         Origin = as.character(pls_Y)) %>%
  left_join(metadata_metabolomics, by = "SampleID")

#Outlier detection (Mahalanobis) on comp1 & comp2 — identify, but do NOT drop
plsda_mahal <- mahalanobis(diaryl_PLSDA_scores[, c("comp1", "comp2")],
                           colMeans(diaryl_PLSDA_scores[, c("comp1", "comp2")]),
                           cov(diaryl_PLSDA_scores[, c("comp1", "comp2")]))

threshold <- qchisq(0.90, df = 4) #re-examine the degrees of freedom here. Should it be 2?
outlier_idx <- which(plsda_mahal > threshold)
message("Identified ", length(outlier_idx), " PLS-DA outlier(s) by Mahalanobis (90% cutoff).") #0 outlier

# custom site colors
site_colors <- c("black","#C4CBCA")
diaryl_PLSDA_scores$Origin.x[diaryl_PLSDA_scores$Origin.x=="LH1"]<-"Little Hoko River"
diaryl_PLSDA_scores$Origin.x[diaryl_PLSDA_scores$Origin.x=="ELK"]<-"Pysht River"
diaryl_PLSDA_scores$Origin.x<-factor(diaryl_PLSDA_scores$Origin.x,levels=c("Little Hoko River","Pysht River"))

diaryl_PLSDA_scores$Origin.x<-as.character(diaryl_PLSDA_scores$Origin.x)
diaryl_PLSDA_scores$Origin.x[diaryl_PLSDA_scores$Origin.x=="Little Hoko River"]<-"Non-local"
diaryl_PLSDA_scores$Origin.x[diaryl_PLSDA_scores$Origin.x=="Pysht River"]<-"Local"
diaryl_PLSDA_scores$Origin.x<-factor(diaryl_PLSDA_scores$Origin.x,levels=c("Local","Non-local"))

library(scales)

# ---- Build Plot ----
Diaryl_PLSDA.plot <- ggplot(diaryl_PLSDA_scores,
                            aes(x = comp1, y = comp2)) +
  scale_x_continuous(limits=c(-6,12),breaks=c(-6, -3, 0, 3, 6, 9, 12))+ 
  scale_y_continuous(limits=c(-6,6),breaks=c(-6, -3, 0, 3, 6))+
  # ELLIPSES (68% CI)
  stat_ellipse(
    aes(group = Origin.x, color = Origin.x),
    type = "t",
    level = 0.68,
    linewidth = 0.8,
    linetype = 1,
    show.legend = FALSE
  ) +
  
  # POINTS
  geom_point(
    aes(fill = Origin.x, shape = Origin.x),
    size = 3.25,
    alpha = 0.9,
    color = "#119822",
    stroke = 2
  ) +
  
  # MANUAL SCALES
  scale_fill_manual(name = "Origin River", values = site_colors) +
  scale_color_manual(name = "Origin River", values = site_colors) +
  scale_shape_manual(name = "Origin River", values = c("Non-local" = 21, "Local" = 21)
  ) +
  
  # THEMING
  theme_bw(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
  ) +
  
  # LABELS
  labs(
    x = paste0("PLS-DA Comp 1 (", round(dry_PLSDA_sugar$prop_expl_var$X[1] * 100, 1), "%)"),
    y = paste0("PLS-DA Comp 2 (", round(dry_PLSDA_sugar$prop_expl_var$X[2] * 100, 1), "%)")
  )

Diaryl_PLSDA.plot

#get the loadings
diaryl_PLSDA_loadings <- mixOmics::plotLoadings(
  dry_PLSDA_sugar, 
  plot = FALSE, 
  contrib = "max")

#convert into a df
diaryl_PLSDA_loadings.df <- as.data.frame(diaryl_PLSDA_loadings$X) %>%
  rownames_to_column("Feature") %>%
  dplyr::select(Feature, GroupContrib)

#generate VIP scores
diaryl_PLSDA_VIP <- as.data.frame(mixOmics::vip(dry_PLSDA_sugar))

#filter for scores > 1.2 (anything > 1 is contributing to group separation)
diaryl_PLSDA_filtered <- dplyr::filter(diaryl_PLSDA_VIP, diaryl_PLSDA_VIP$comp1 > 1.0) #set threshold to zero for main figure

diaryl_PLSDA_filtered$ID <- rownames(diaryl_PLSDA_filtered)

#select comp 1 (similar to PC axis for PCA) and comp 2
diaryl_PLSDA_select <- diaryl_PLSDA_filtered %>% dplyr::select(ID, comp1, comp2)

#create new VIP df
diaryl_PLSDA_VIP.df <- diaryl_PLSDA_select %>% 
  left_join(diaryl_PLSDA_loadings.df, by = c("ID" = "Feature")) %>%
  left_join(info_feature_filtered_s2_clean, by = c("ID" = "Feature")) %>% 
  arrange(desc(comp1))

dim(diaryl_PLSDA_VIP.df) #VIP > 1: 21 x 136

#filter for the top 20 VIP features
diaryl_PLSDA_20_VIP <- diaryl_PLSDA_VIP.df %>%
  head(20)

#filter NA / aggregate with annotations / arrange by median VIP
diaryl_PLSDA_aggregated_data <- diaryl_PLSDA_VIP.df %>%
  filter(!is.na(Sirius_NPC.superclass)) %>%
  dplyr::group_by(Sirius_NPC.superclass) %>%
  dplyr::summarise(median_VIP = median(comp1, na.rm = TRUE)) %>%
  arrange(desc(median_VIP))

############Diarylheptanoid VIP relative abundance models
#filter
dry_vips<-diaryl_PLSDA_VIP.df$ID
dry_vip_table <- pls_X_filteredSacc[,colnames(pls_X_filteredSacc) %in% dry_vips]
dim(dry_vip_table) #16 x 21

#rowsum intensities
diaryl_VIPintensities <- rowSums(dry_vip_table) 

#getting a weighted score instead
vip_features <- intersect(colnames(dry_vip_table),diaryl_PLSDA_VIP.df$ID)

length(vip_features)
vip_matrix <- dry_vip_table[, vip_features]

vip_df <- diaryl_PLSDA_VIP.df[match(vip_features, diaryl_PLSDA_VIP.df$ID), ]

vip_weights <- vip_df$comp1 / sum(vip_df$comp1) #weighted vip scores. 
vip_weights2 <- vip_df$comp2 / sum(vip_df$comp1) #weighted vip scores. 
length(vip_weights)

str(vip_matrix)
vip_matrix <- as.matrix(vip_matrix)
str(vip_weights)

#Per sample VIP-weighted Score
vip_score <- as.vector(vip_matrix %*% vip_weights) #A weighted sum of all VIP features, Essentially: “how much does this sample express the features that most drive group separation?”
vip_score_comp2 <- as.vector(vip_matrix %*% vip_weights2)
vip_total_weighted_score<-vip_score+vip_score_comp2

meta_y<-dry_is.log.df[,1:2]

meta_xy<-cbind(meta_y,diaryl_VIPintensities,vip_score,vip_score_comp2,vip_total_weighted_score)
str(meta_xy)

hist(meta_xy$diaryl_VIPintensities)
hist(meta_xy$vip_score)
hist(vip_score_comp2)
hist(vip_total_weighted_score)

meta_xy$Origin[meta_xy$Origin=="LH1"]<-"Little Hoko River"
meta_xy$Origin[meta_xy$Origin=="ELK"]<-"Pysht River"
meta_xy$Origin<-factor(meta_xy$Origin,levels=c("Little Hoko River","Pysht River"))

meta_xy$Origin<-as.character(meta_xy$Origin)
meta_xy$Origin[meta_xy$Origin=="Little Hoko River"]<-"Non-local"
meta_xy$Origin[meta_xy$Origin=="Pysht River"]<-"Local"
meta_xy$Origin<-factor(meta_xy$Origin,levels=c("Local","Non-local"))

#model
library(car)
library(multcompView)

Anova(lm(meta_xy$diaryl_VIPintensities~meta_xy$Origin)) #summed intensities of VIPs
Anova(lm(meta_xy$vip_score~meta_xy$Origin)) #weighted summed VIP scores.
Anova(lm(meta_xy$vip_score_comp2~meta_xy$Origin))
Anova(lm(meta_xy$vip_total_weighted_score~meta_xy$Origin))

vp1<-aov(meta_xy$diaryl_VIPintensities~meta_xy$Origin)
vp1b<-TukeyHSD(vp1)
vp1b

multcompLetters4(vp1, vp1b)
#$`meta_xy$Origin`
#Pysht River Little Hoko River 
#"a"               "b" 

diaryl_relative_abundVIP<-ggplot(data=meta_xy, aes(x=Origin,y=diaryl_VIPintensities,fill=Origin,shape=Origin))+
  geom_boxplot(aes(fill=Origin),alpha=0.5, outlier.size=0) +
  geom_point(aes(group=Origin, fill=Origin), color = "#119822", stroke = 1.5, size=3.25, position=position_jitterdodge(jitter.width=1.25)) +
  #facet_grid(~ Origin_River)+
  scale_shape_manual(values = c(21,21,21,21,21)) +
  scale_fill_manual(values = c(site_colors)) + # Boxplot fill color
  scale_y_continuous(label=scales::comma,limits=c(230,290),breaks=c(230,240,250,260,270,280,290))+
  theme_classic() +
  labs(y="Cumulative Feature Abundance\n of Diarylheptanoids",x="Origin River") +
  theme(legend.position = "right")

diaryl_relative_abundVIP_final<-diaryl_relative_abundVIP + theme(axis.text.x = element_text(size = 12)) + theme(text = element_text(size = 15),legend.text.align = 0)

diaryl_relative_abundVIP_final 

#######################################################
#######################################################
###ellagitannins, placing useful lines here so that i can read object names without scrolling.
length(ellagic_keeps1) #42
#######################################################
#######################################################

pls_X_ellagic <- dry_is.log.df1 %>% dplyr::select(
  all_of(ellagic_keeps1)
) #16 x 42

pls_YSellagic <- factor(dry_is.log.df1$Origin)

if (nrow(pls_X_ellagic) != length(pls_YSellagic)) stop("Mismatch between X rows and Y length.")

dry_PLSDA_ellagic <- mixOmics::plsda(
  X = pls_X_ellagic,
  Y = pls_YSellagic,
  ncomp = 2,
  scale = TRUE
)

perf_PLSDA_ellag <- perf(dry_PLSDA_ellagic, 
                         validation = "Mfold",  
                         folds = 3,             # adjust to what you used
                         nrepeat = 100,   
                         progressBar = TRUE)


cumsum(dry_PLSDA_ellagic$prop_expl_var$X) #cumsum comp2
perf_PLSDA_ellag$error.rate$BER #max.dist comp2
perf_PLSDA_ellag$error.rate$overall #max.dist comp2

#Extract PLS-DA scores and reattach metadata
ellagic_PLSDA_scores <- as.data.frame(dry_PLSDA_ellagic$variates$X) %>%
  rownames_to_column("tmp_rowname") %>%    # mixOmics rownames map to original row order
  mutate(SampleID = pls_YSellagic,
         Origin = as.character(pls_YSellagic)) %>%
  left_join(metadata_metabolomics, by = "SampleID")

#Outlier detection (Mahalanobis) on comp1 & comp2 — identify, but do NOT drop
plsda_mahal <- mahalanobis(ellagic_PLSDA_scores[, c("comp1", "comp2")],
                           colMeans(ellagic_PLSDA_scores[, c("comp1", "comp2")]),
                           cov(ellagic_PLSDA_scores[, c("comp1", "comp2")]))

threshold <- qchisq(0.90, df = 4) #re-examine the degrees of freedom here. Should it be 2?
outlier_idx <- which(plsda_mahal > threshold)
message("Identified ", length(outlier_idx), " PLS-DA outlier(s) by Mahalanobis (90% cutoff).") #0 outlier

# custom site colors
site_colors <- c("black","#C4CBCA")
ellagic_PLSDA_scores$Origin.x[ellagic_PLSDA_scores$Origin.x=="LH1"]<-"Little Hoko River"
ellagic_PLSDA_scores$Origin.x[ellagic_PLSDA_scores$Origin.x=="ELK"]<-"Pysht River"
ellagic_PLSDA_scores$Origin.x<-factor(ellagic_PLSDA_scores$Origin.x,levels=c("Little Hoko River","Pysht River"))

# ---- Build Plot ----
Ellagic_PLSDA.plot <- ggplot(ellagic_PLSDA_scores,
                             aes(x = comp1, y = comp2)) +
  scale_x_continuous(limits=c(-5,7.5),breaks=c(-5, -2.5, 0, 2.5, 5, 7.5))+
  scale_y_continuous(limits=c(-4,4),breaks=c(-4, -2, 0, 2, 4))+
  
  # ELLIPSES (68% CI)
  stat_ellipse(
    aes(group = Origin.x, color = Origin.x),
    type = "t",
    level = 0.68,
    linewidth = 0.8,
    linetype = 1,
    show.legend = FALSE
  ) +
  
  # POINTS
  geom_point(
    aes(fill = Origin.x, shape = Origin.x),
    size = 4,
    alpha = 0.9,
    color = "#119822",
    stroke = 2
  ) +
  
  # MANUAL SCALES
  scale_fill_manual(name = "Origin River", values = site_colors) +
  scale_color_manual(name = "Origin River", values = site_colors) +
  scale_shape_manual(name = "Origin River", values = c("Little Hoko River" = 21, "Pysht River" = 21)
  ) +
  
  # THEMING
  theme_bw(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
  ) +
  
  # LABELS
  labs(
    x = paste0("PLS-DA Comp 1 (", round(dry_PLSDA_ellagic$prop_expl_var$X[1] * 100, 1), "%)"),
    y = paste0("PLS-DA Comp 2 (", round(dry_PLSDA_ellagic$prop_expl_var$X[2] * 100, 1), "%)")
  ) 

Ellagic_PLSDA.plot

#get the loadings
ellagic_PLSDA_loadings <- mixOmics::plotLoadings(
  dry_PLSDA_ellagic, 
  plot = FALSE, 
  contrib = "max")

#convert into a df
ellagic_PLSDA_loadings.df <- as.data.frame(ellagic_PLSDA_loadings$X) %>%
  rownames_to_column("Feature") %>%
  dplyr::select(Feature, GroupContrib)

#generate VIP scores
ellagic_PLSDA_VIP <- as.data.frame(mixOmics::vip(dry_PLSDA_ellagic))

#filter for scores > 1.2 (anything > 1 is contributing to group separation)
ellagic_PLSDA_filtered <- dplyr::filter(ellagic_PLSDA_VIP, ellagic_PLSDA_VIP$comp1 > 1.0) 

ellagic_PLSDA_filtered$ID <- rownames(ellagic_PLSDA_filtered)

#select comp 1 (similar to PC axis for PCA) and comp 2
ellagic_PLSDA_select <- ellagic_PLSDA_filtered %>% dplyr::select(ID, comp1, comp2)

#create new VIP df
ellagic_PLSDA_VIP.df <- ellagic_PLSDA_select %>% 
  left_join(ellagic_PLSDA_loadings.df, by = c("ID" = "Feature")) %>%
  left_join(info_feature_filtered_s2_clean, by = c("ID" = "Feature")) %>% 
  arrange(desc(comp1))

dim(ellagic_PLSDA_VIP.df) #VIP > 1: 19 x 136

#filter for the top 20 VIP features
ellagic_PLSDA_20_VIP <- ellagic_PLSDA_VIP.df %>%
  head(20)

#filter NA / aggregate with annotations / arrange by median VIP
ellagic_PLSDA_aggregated_data <- ellagic_PLSDA_VIP.df %>%
  filter(!is.na(Sirius_NPC.superclass)) %>%
  dplyr::group_by(Sirius_NPC.superclass) %>%
  dplyr::summarise(median_VIP = median(comp1, na.rm = TRUE)) %>%
  arrange(desc(median_VIP))

#filter
ellag_vips<-ellagic_PLSDA_VIP.df$ID
ellag_vip_table <- pls_X_ellagic[,colnames(pls_X_ellagic) %in% ellag_vips]
dim(ellag_vip_table) #16 x 19

#rowsum intensities
ellag_VIPintensities <- rowSums(ellag_vip_table) 

#getting a weighted score instead
vip_features <- intersect(colnames(ellag_vip_table), ellagic_PLSDA_VIP.df$ID)

length(vip_features) #19
vip_matrix <- ellag_vip_table[, vip_features]

vip_df <- ellagic_PLSDA_VIP.df[match(vip_features, ellagic_PLSDA_VIP.df$ID), ]

vip_weights <- vip_df$comp1 / sum(vip_df$comp1) #weighted vip scores. 
vip_weights2 <- vip_df$comp2 / sum(vip_df$comp1) #weighted vip scores. 
length(vip_weights)

str(vip_matrix)
vip_matrix <- as.matrix(vip_matrix)
str(vip_weights)

#Per sample VIP-weighted Score
vip_score <- as.vector(vip_matrix %*% vip_weights) #A weighted sum of all VIP features, Essentially: “how much does this sample express the features that most drive group separation?”
vip_score_comp2 <- as.vector(vip_matrix %*% vip_weights2)
vip_total_weighted_score<-vip_score+vip_score_comp2

meta_y<-dry_is.log.df[,1:2]

meta_xy_ellag<-cbind(meta_y,ellag_VIPintensities,vip_score,vip_score_comp2,vip_total_weighted_score)
str(meta_xy_ellag)

hist(meta_xy_ellag$ellag_VIPintensities)
hist(meta_xy_ellag$vip_score)
hist(vip_score_comp2)
hist(vip_total_weighted_score)

meta_xy_ellag$Origin[meta_xy_ellag$Origin=="LH1"]<-"Little Hoko River"
meta_xy_ellag$Origin[meta_xy_ellag$Origin=="ELK"]<-"Pysht River"
meta_xy_ellag$Origin<-factor(meta_xy_ellag$Origin,levels=c("Little Hoko River","Pysht River"))

Anova(lm(meta_xy_ellag$ellag_VIPintensities~meta_xy_ellag$Origin)) #summed intensities of VIPs
Anova(lm(meta_xy_ellag$vip_score~meta_xy_ellag$Origin)) #weighted summed VIP scores.
Anova(lm(meta_xy_ellag$vip_score_comp2~meta_xy_ellag$Origin))
Anova(lm(meta_xy_ellag$vip_total_weighted_score~meta_xy_ellag$Origin))

vp1<-aov(meta_xy_ellag$ellag_VIPintensities~meta_xy_ellag$Origin)
vp1b<-TukeyHSD(vp1)
vp1b

multcompLetters4(vp1, vp1b)
#$`meta_xy$Origin`
#Pysht River Little Hoko River 
#"a"               "a" 

ellag_relative_abundVIP<-ggplot(data=meta_xy_ellag, aes(x=Origin,y=ellag_VIPintensities,fill=Origin,shape=Origin))+
  geom_boxplot(aes(fill=Origin),alpha=0.5, outlier.size=0) +
  geom_point(aes(group=Origin, fill=Origin), color = "#119822", stroke = 1.5, size=3.0, position=position_jitterdodge(jitter.width=1.25)) +
  #facet_grid(~ Origin_River)+
  scale_shape_manual(values = c(21,21,21,21,21)) +
  scale_fill_manual(values = c(site_colors)) + # Boxplot fill color
  scale_y_continuous(label=scales::comma,limits=c(170,230),breaks=c(170,180,190,200,210,220,230))+
  theme_classic() +
  labs(y="Cumulative Feature Abundance",x="Origin Site") +
  theme(legend.position = "none")

ellag_relative_abundVIP_final<-ellag_relative_abundVIP + theme(axis.text.x = element_text(size = 12)) + theme(text = element_text(size = 15),legend.text.align = 0)

ellag_relative_abundVIP_final 

#######################################################
#######################################################
#Flavonoids
#######################################################
#######################################################
Candidate_flav_features <- info_feature_filtered_s2_psms %>%
  filter(Sirius_NPC.superclass %in% c("Flavonoids")) %>%
  pull(Feature)

pls_XFlav <- dry_is.log.df1 %>% dplyr::select(
  all_of(Candidate_flav_features))

nzv_idxflav <- caret::nearZeroVar(pls_XFlav, names = TRUE)

pls_X_filteredFlav <- pls_XFlav %>%
  #dplyr::select(-SampleID, -Origin) %>%
  dplyr::select(-all_of(nzv_idxflav))

pls_YFlav <- factor(dry_is.log.df1$Origin)

if (nrow(pls_X_filteredFlav) != length(pls_YFlav)) stop("Mismatch between X rows and Y length.")

dry_PLSDA_flav <- mixOmics::plsda(
  X = pls_X_filteredFlav,
  Y = pls_YFlav,
  ncomp = 2,
  scale = TRUE
)

perf_PLSDA_flav <- perf(dry_PLSDA_flav, 
                        validation = "Mfold",  
                        folds = 3,             # adjust to what you used
                        nrepeat = 100,   
                        progressBar = TRUE)


cumsum(dry_PLSDA_flav$prop_expl_var$X) #cumsum comp2
perf_PLSDA_flav$error.rate$BER #max.dist comp2
perf_PLSDA_flav$error.rate$overall #max.dist comp2

#Extract PLS-DA scores and reattach metadata
flav_PLSDA_scores <- as.data.frame(dry_PLSDA_flav$variates$X) %>%
  rownames_to_column("tmp_rowname") %>%    # mixOmics rownames map to original row order
  mutate(SampleID = pls_YFlav,
         Origin = as.character(pls_YFlav)) %>%
  left_join(metadata_metabolomics, by = "SampleID")

#Outlier detection (Mahalanobis) on comp1 & comp2 — identify, but do NOT drop
plsda_mahal <- mahalanobis(flav_PLSDA_scores[, c("comp1", "comp2")],
                           colMeans(flav_PLSDA_scores[, c("comp1", "comp2")]),
                           cov(flav_PLSDA_scores[, c("comp1", "comp2")]))

threshold <- qchisq(0.90, df = 4) #re-examine the degrees of freedom here. Should it be 2?
outlier_idx <- which(plsda_mahal > threshold)
message("Identified ", length(outlier_idx), " PLS-DA outlier(s) by Mahalanobis (90% cutoff).") #0 outlier

# custom site colors
site_colors <- c("black","#C4CBCA")
flav_PLSDA_scores$Origin.x[flav_PLSDA_scores$Origin.x=="LH1"]<-"Little Hoko River"
flav_PLSDA_scores$Origin.x[flav_PLSDA_scores$Origin.x=="ELK"]<-"Pysht River"
flav_PLSDA_scores$Origin.x<-factor(flav_PLSDA_scores$Origin.x,levels=c("Little Hoko River","Pysht River"))

# ---- Build Plot ----
flav_PLSDA.plot <- ggplot(flav_PLSDA_scores,
                          aes(x = comp1, y = comp2)) +
  scale_x_continuous(limits=c(-6,6),breaks=c(-6, -3, 0, 3, 6))+
  scale_y_continuous(limits=c(-7.5,5.25),breaks=c(-7.5, -5, -2.5, 0, 2.5, 5))+
  
  # ELLIPSES (68% CI)
  stat_ellipse(
    aes(group = Origin.x, color = Origin.x),
    type = "t",
    level = 0.68,
    linewidth = 0.8,
    linetype = 1,
    show.legend = FALSE
  ) +
  
  # POINTS
  geom_point(
    aes(fill = Origin.x, shape = Origin.x),
    size = 4,
    alpha = 0.9,
    color = "#119822",
    stroke = 2
  ) +
  
  # MANUAL SCALES
  scale_fill_manual(name = "Origin River", values = site_colors) +
  scale_color_manual(name = "Origin River", values = site_colors) +
  scale_shape_manual(name = "Origin River", values = c("Little Hoko River" = 21, "Pysht River" = 21)
  ) +
  
  # THEMING
  theme_bw(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
  ) +
  
  # LABELS
  labs(
    x = paste0("PLS-DA Comp 1 (", round(dry_PLSDA_flav$prop_expl_var$X[2] * 100, 1), "%)"),
    y = paste0("PLS-DA Comp 2 (", round(dry_PLSDA_flav$prop_expl_var$X[1] * 100, 1), "%)")
  )

flav_PLSDA.plot

#get the loadings
flav_PLSDA_loadings <- mixOmics::plotLoadings(
  dry_PLSDA_flav, 
  plot = FALSE, 
  contrib = "max")

#convert into a df
flav_PLSDA_loadings.df <- as.data.frame(flav_PLSDA_loadings$X) %>%
  rownames_to_column("Feature") %>%
  dplyr::select(Feature, GroupContrib)

#generate VIP scores
flav_PLSDA_VIP <- as.data.frame(mixOmics::vip(dry_PLSDA_flav))

#filter for scores > 1.2 (anything > 1 is contributing to group separation)
flav_PLSDA_filtered <- dplyr::filter(flav_PLSDA_VIP, flav_PLSDA_VIP$comp1 > 1.0) #set threshold to zero for main figure

flav_PLSDA_filtered$ID <- rownames(flav_PLSDA_filtered)

#select comp 1 (similar to PC axis for PCA) and comp 2
flav_PLSDA_select <- flav_PLSDA_filtered %>% dplyr::select(ID, comp1, comp2)

#create new VIP df
flav_PLSDA_VIP.df <- flav_PLSDA_select %>% 
  left_join(flav_PLSDA_loadings.df, by = c("ID" = "Feature")) %>%
  left_join(info_feature_filtered_s2_clean, by = c("ID" = "Feature")) %>% 
  arrange(desc(comp1))

dim(flav_PLSDA_VIP.df) #VIP > 1: 45 x 136

#filter for the top 20 VIP features
flav_PLSDA_20_VIP <- flav_PLSDA_VIP.df %>%
  head(20)

#filter NA / aggregate with annotations / arrange by median VIP
flav_PLSDA_aggregated_data <- flav_PLSDA_VIP.df %>%
  filter(!is.na(Sirius_NPC.superclass)) %>%
  dplyr::group_by(Sirius_NPC.superclass) %>%
  dplyr::summarise(median_VIP = median(comp1, na.rm = TRUE)) %>%
  arrange(desc(median_VIP))

#filter
flav_vips<-flav_PLSDA_VIP.df$ID
flav_vip_table <- pls_X_filteredFlav[,colnames(pls_X_filteredFlav) %in% flav_vips]
dim(flav_vip_table) #16 x 45

#rowsum intensities
flav_VIPintensities <- rowSums(flav_vip_table) 

#getting a weighted score instead
vip_features <- intersect(colnames(flav_vip_table), flav_PLSDA_VIP.df$ID)

length(vip_features) #45
vip_matrix <- flav_vip_table[, vip_features]

vip_df <- flav_PLSDA_VIP.df[match(vip_features, flav_PLSDA_VIP.df$ID), ]

vip_weights <- vip_df$comp1 / sum(vip_df$comp1) #weighted vip scores. 
vip_weights2 <- vip_df$comp2 / sum(vip_df$comp1) #weighted vip scores. 
length(vip_weights)

str(vip_matrix)
vip_matrix <- as.matrix(vip_matrix)
str(vip_weights)

#Per sample VIP-weighted Score
vip_score <- as.vector(vip_matrix %*% vip_weights) #A weighted sum of all VIP features, Essentially: “how much does this sample express the features that most drive group separation?”
vip_score_comp2 <- as.vector(vip_matrix %*% vip_weights2)
vip_total_weighted_score<-vip_score+vip_score_comp2

meta_y<-dry_is.log.df[,1:2]

meta_xy_flav<-cbind(meta_y,flav_VIPintensities,vip_score,vip_score_comp2,vip_total_weighted_score)
str(meta_xy_flav)

hist(meta_xy_flav$flav_VIPintensities)
hist(meta_xy_flav$vip_score)
hist(vip_score_comp2)
hist(vip_total_weighted_score)

meta_xy_flav$Origin[meta_xy_flav$Origin=="LH1"]<-"Little Hoko River"
meta_xy_flav$Origin[meta_xy_flav$Origin=="ELK"]<-"Pysht River"
meta_xy_flav$Origin<-factor(meta_xy_flav$Origin,levels=c("Little Hoko River","Pysht River"))

Anova(lm(meta_xy_flav$flav_VIPintensities~meta_xy_flav$Origin)) #summed intensities of VIPs
Anova(lm(meta_xy_flav$vip_score~meta_xy_flav$Origin)) #weighted summed VIP scores.
Anova(lm(meta_xy_flav$vip_score_comp2~meta_xy_flav$Origin))
Anova(lm(meta_xy_flav$vip_total_weighted_score~meta_xy_flav$Origin))

vp1<-aov(meta_xy_flav$flav_VIPintensities~meta_xy_flav$Origin)
vp1b<-TukeyHSD(vp1)
vp1b

multcompLetters4(vp1, vp1b)
#$`meta_xy$Origin`
#Pysht River Little Hoko River 
#"a"               "a" 

flav_relative_abundVIP<-ggplot(data=meta_xy_flav, aes(x=Origin,y=flav_VIPintensities,fill=Origin,shape=Origin))+
  geom_boxplot(aes(fill=Origin),alpha=0.5, outlier.size=0) +
  geom_point(aes(group=Origin, fill=Origin), color = "#119822", stroke = 1.5, size=3.0, position=position_jitterdodge(jitter.width=1.25)) +
  #facet_grid(~ Origin_River)+
  scale_shape_manual(values = c(21,21,21,21,21)) +
  scale_fill_manual(values = c(site_colors)) + # Boxplot fill color
  scale_y_continuous(label=scales::comma,limits=c(400,500),breaks=c(400,420,440,460,480,500))+
  theme_classic() +
  labs(y="Cumulative Feature Abundance",x="Origin Site") +
  theme(legend.position = "none")

flav_relative_abundVIP_final<-flav_relative_abundVIP + theme(axis.text.x = element_text(size = 12)) + theme(text = element_text(size = 15),legend.text.align = 0)

flav_relative_abundVIP_final

plot_grid(Ellagic_PLSDA.plot,flav_PLSDA.plot, nrow=1, ncol =2, labels= LETTERS)

########now that i have all possible figures of interest: let's build a multipanel figure. 

library(cowplot)

hm_grob <- grid::grid.grabExpr(
  draw(hm.plot)
)

left_col <- plot_grid(
  leaf_PLSDA.plot, hm_grob,
  nrow = 2, ncol = 1,
  align = "hv",
  axis = "tl",
  labels = c("a", "c")
)

left_col_alt <- plot_grid(
  leaf_PLSDA.plot, vip_plot1,
  nrow = 1, ncol = 2,
  #align = "hv",
  #axis = "tl",
  labels = c("a", "b"))

diaryls <- plot_grid(Diaryl_PLSDA.plot, diaryl_relative_abundVIP_final, 
                     nrow = 1, ncol = 2,rel_widths = c(1,1), 
                     align = "hv",
                     axis = "tl",
                     labels = c("d", "e"))

same <- plot_grid(left_col_alt, diaryls, nrow=2, ncol=1, rel_heights = c(1,1), align = "hv", axis = "tlbr")

plot_grid(
  same, hm_grob,
  nrow = 1, ncol = 2,
  labels = c("", "c"),
  rel_widths = c(1, 0.6),
  label_size= 14) #13 x 22.5 pdf

plot_grid(
  same, hm_grob,
  nrow = 1, ncol = 2,
  labels = c("", "c"),
  rel_widths = c(1, 0.75),
  label_size = 14) #13 x 22.5 pdf

leftcol3 <- plot_grid(
  leaf_PLSDA.plot, Diaryl_PLSDA.plot, diaryl_relative_abundVIP_final,
  nrow = 3, ncol = 1,
  #align = "hv",
  #axis = "tl",
  rel_widths = c(1,1, 1),
  labels = c("", "", ""))

plot_grid(leftcol3, vip_plot1, hm_grob,
          nrow=1, ncol = 3,
          labels = c("","",""),
          rel_widths = c(0.75,1,1),
          label_size = 16)