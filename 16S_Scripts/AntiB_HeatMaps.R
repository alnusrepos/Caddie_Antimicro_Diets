library(phyloseq); packageVersion("phyloseq") #1.26.1
library(ggplot2); packageVersion("ggplot2") #3.4.1
#remotes::install_github("david-barnett/microViz")
library("microViz")
#install.packages("/Users/lab/rjson",repos=NULL,type="source")
library(rjson)

setwd("/Users/lab/Documents/Former_Researchers/Loomis_Caddies/Final_Qiime2_outputs")
#to load in the global environment and skip lines 26-248 run the line below.
load(file="CaddiesSeqAnalysis072825.rda")

CaddieD_phylo.rm2 #non rarefied data
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 2895 taxa and 37 samples ]
# sample_data() Sample Data:       [ 37 samples by 8 sample variables ]
# tax_table()   Taxonomy Table:    [ 2895 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 2895 tips and 2889 internal nodes ]

#examine and modify asv table, sample data, and write taxonomy table to make changes and read back in
colSums(CaddieD_phylo.rm2@otu_table) #no zeros
taxa_names(CaddieD_phylo.rm2)<-paste0("ASV",seq(ntaxa(CaddieD_phylo.rm2)))
pss.rm<-CaddieD_phylo.rm2 #rename
pss.rm@otu_table[1:10,1:10]

mean(rowSums(pss.rm@otu_table)) 
median(rowSums(pss.rm@otu_table)) 

###### phyloseq-class experiment-level object
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 2895 taxa and 37 samples ]
# sample_data() Sample Data:       [ 37 samples by 8 sample variables ]
# tax_table()   Taxonomy Table:    [ 2895 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 2895 tips and 2889 internal nodes ]

colnames(pss.rm@sam_data)
#"Sample.Name.", "Sample.Type.and.." ,"Treatment.Type","Case.Mass..g.","Body.Mass..g.","Case.Width..cm." ,"Abodomen.and.Thorax.Length..cm." ,"Entire.Length..cm." 
str(pss.rm@sam_data)
pss.rm.md<-pss.rm@sam_data
pss.rm.md$Treatment.Type[pss.rm.md$Treatment.Type=="treated"]<-"Antibiotic_Treated"
pss.rm.md$Treatment.Type[pss.rm.md$Treatment.Type=="untreated"]<-"Control"
pss.rm.md$Treatment.Type<-factor(pss.rm.md$Treatment.Type, levels=c("Antibiotic_Treated","Control"))
pss.rm.tab<-pss.rm@otu_table
pss.rm.tax<-pss.rm@tax_table

#write this to create concatenated column.
#write.csv(pss.rm.tax,file="pss.rm.tax.csv") 

#read back in
pss.rm.tax1<-read.csv(file="pss.rm.tax.concat.csv") #this csv file contains an extra column where taxonomic information were combined into a single column using TEXTJOIN() in excel. For example, for each ASV, the column would read "Phylum; Family; Genus; species; ASV#" or whatever hierarchical compoents you would like. 
pss.rm.phy<-pss.rm@phy_tree
colnames(pss.rm.tax)
colnames(pss.rm.tax1)
rownames(pss.rm.tax)
asv.id<-pss.rm.tax1$X
pss.rm.tax1<-pss.rm.tax1[,2:10]
colnames(pss.rm.tax1)
rownames(pss.rm.tax1)<-asv.id
class(pss.rm.tax1)
pss.rm.tax1m<-as.matrix(pss.rm.tax1)

#steps to put back into a phyloseq object for downstream heat map use
taxa_names(pss.rm.tab)
taxa_names(pss.rm.tax1m)<-taxa_names(pss.rm.tab)
taxa_names(pss.rm.md)<-taxa_names(pss.rm.tab)
taxa_names(pss.rm.phy)<-taxa_names(pss.rm.tab)

sn<-pss.rm.md$Sample.Name.
sample_names(pss.rm.tab)<-sn
sample_names(pss.rm.tax1m)<-sn
sample_names(pss.rm.md)<-sn
sample_names(pss.rm.phy)<-sn

pss.rm2<-phyloseq(otu_table(pss.rm.tab,taxa_are_rows=FALSE),sample_data(pss.rm.md),tax_table(pss.rm.tax1m),phy_tree(pss.rm.phy))
str(pss.rm2@sam_data)

#view the goodies
colSums(pss.rm2@otu_table) #minimum total abundance -- minimum total abund of a taxon, summed across all samples.
rowSums(pss.rm2@otu_table) 
colnames(pss.rm2@tax_table)

library(dplyr)
#when and if you need to pre-filter ASVs based on abundance here is where to do it. Running 82-85 as is does no filtering. 
psq <- pss.rm2 %>% #2895 taxa
  #tax_mutate(ConCat = NULL) %>%
  tax_filter(min_prevalence = 1) %>% #number or proportion of samples that a taxon must be present in
  tax_fix() 

#view the goodies, anything weird? Yeah maybe the Uncultured, Unknowns, Unassigned that need to be resolved in csv. Are they useful?
table(tax_table(psq)[, "Domain"], exclude = NULL)
table(tax_table(psq)[, "Phylum"], exclude = NULL) 
table(tax_table(psq)[, "Class"], exclude = NULL) 

##############################
##############################
#For all ASVs grouped by Family
##############################
##############################
htmp7 <- psq %>%
  #sort all samples by similarity
 ps_seriate(rank = "Family", tax_transform = "log10") %>% 
  # arrange 
  ps_arrange(Treatment.Type) %>% 
  tax_transform("log10", rank = "Family") %>%
  comp_heatmap(
    colors = heat_palette(palette = viridis::viridis(10)),
    row_names_side = "right", row_dend_side = "left", row_dend_reorder = TRUE, row_dend_width= grid::unit(5,"cm"), sample_side = "bottom",
    sample_anno = sampleAnnotation(
      Treatment.Type = anno_sample_cat(var = "Treatment.Type", col = c(Antibiotic_Treated = "#d81159", Control = "#b5e2fa"), 
        border_col = "black", box_col = NA, legend_title = "Treatment Type", size = grid::unit(4, "mm"))
    ),
    sample_seriation = "Identity" # suppress sample reordering
  )

htmp7 %>% ComplexHeatmap::draw(
  annotation_legend_list = attr(htmp7, "AnnoLegends"),
)


#An ASV must be present in at least two samples
psq1 <- pss.rm2 %>%
  #tax_mutate(ConCat = NULL) %>%
  tax_filter(min_prevalence = 0.05) %>% #number or proportion of samples that a taxon must be present in
  tax_fix()

htmp7b <- psq1 %>%
  #sort all samples by similarity
  ps_seriate(rank = "Family", tax_transform = "log10") %>% 
  # arrange 
  ps_arrange(Treatment.Type) %>% 
  tax_transform("log10", rank = "Family") %>%
  comp_heatmap(
    colors = heat_palette(palette = viridis::viridis(10)),
    row_names_side = "right", row_dend_side = "left", row_dend_reorder = TRUE, row_dend_width= grid::unit(5,"cm"), sample_side = "bottom",
    sample_anno = sampleAnnotation(
      Treatment.Type = anno_sample_cat(var = "Treatment.Type", col = c(Antibiotic_Treated = "#d81159", Control = "#b5e2fa"), 
                                       border_col = "black", box_col = NA, legend_title = "Treatment Type", size = grid::unit(4, "mm"))
    ),
    sample_seriation = "Identity" # suppress sample reordering
  )

htmp7b %>% ComplexHeatmap::draw(
  annotation_legend_list = attr(htmp7b, "AnnoLegends"),
)


#Set a threshold for low abundance. For example, you might want to remove taxa that appear in less than a certain percentage of your samples or those that have a total abundance below a certain count.

#Calculate Abundance: You can calculate the total abundance for each taxon across all samples.
total_abundance <- taxa_sums(psq)

# Set the abundance threshold
threshold <- 33 #the top 20% most abundant ASVs across all samples
thresh <- 21 #the top 25% most abundant ASVs across all samples
thresh2 <- 138 #the top 10% most abundant ASVs across all samples

# Filter taxa based on the threshold
physeq_filtered <- prune_taxa(total_abundance >= threshold, psq) #580
physeq_filtered1 <- prune_taxa(total_abundance >= thresh, psq) #722
physeq_filtered2 <- prune_taxa(total_abundance >= thresh2, psq) #291 taxa

#Top 25% most abundant ASVs across all samples, at the Family level
htmp7c <- physeq_filtered1 %>%
  #sort all samples by similarity
  ps_seriate(rank = "Family", tax_transform = "log10") %>% 
  # arrange 
  ps_arrange(Treatment.Type) %>% 
  tax_transform("log10", rank = "Family") %>%
  comp_heatmap(
    colors = heat_palette(palette = viridis::viridis(10)),
    row_names_side = "right", row_dend_side = "left", row_dend_reorder = TRUE, row_dend_width= grid::unit(5,"cm"), sample_side = "bottom",
    sample_anno = sampleAnnotation(
      Treatment.Type = anno_sample_cat(var = "Treatment.Type", col = c(Antibiotic_Treated = "#d81159", Control = "#b5e2fa"), 
                                       border_col = "black", box_col = NA, legend_title = "Treatment Type", size = grid::unit(4, "mm"))
    ),
    sample_seriation = "Identity" # suppress sample reordering
  )

htmp7c %>% ComplexHeatmap::draw(
  annotation_legend_list = attr(htmp7c, "AnnoLegends"),
)

#Top 20% most abundant ASVs across all samples, at the Family level
htmp7d <- physeq_filtered %>%
  #sort all samples by similarity
  ps_seriate(rank = "Family", tax_transform = "log10") %>% 
  # arrange 
  ps_arrange(Treatment.Type) %>% 
  tax_transform("log10", rank = "Family") %>%
  comp_heatmap(
    colors = heat_palette(palette = viridis::viridis(10)),
    row_names_side = "right", row_dend_side = "left", row_dend_reorder = TRUE, row_dend_width= grid::unit(5,"cm"), sample_side = "bottom",
    sample_anno = sampleAnnotation(
      Treatment.Type = anno_sample_cat(var = "Treatment.Type", col = c(Antibiotic_Treated = "#d81159", Control = "#b5e2fa"), 
                                       border_col = "black", box_col = NA, legend_title = "Treatment Type", size = grid::unit(4, "mm"))
    ),
    sample_seriation = "Identity" # suppress sample reordering
  )

htmp7d %>% ComplexHeatmap::draw(
  annotation_legend_list = attr(htmp7d, "AnnoLegends"),
)

#Top 10% most abundant ASVs across all samples, at the Famly level
htmp7e <- physeq_filtered2 %>%
  #sort all samples by similarity
  ps_seriate(rank = "Family", tax_transform = "log10") %>% 
  # arrange 
  ps_arrange(Treatment.Type) %>% 
  tax_transform("log10", rank = "Family") %>%
  comp_heatmap(
    colors = heat_palette(palette = viridis::viridis(10)),
    row_names_side = "right", row_dend_side = "left", row_dend_reorder = TRUE, row_dend_width= grid::unit(5,"cm"), sample_side = "bottom",
    sample_anno = sampleAnnotation(
      Treatment.Type = anno_sample_cat(var = "Treatment.Type", col = c(Antibiotic_Treated = "#d81159", Control = "#b5e2fa"), 
                                       border_col = "black", box_col = NA, legend_title = "Treatment Type", size = grid::unit(4, "mm"))
    ),
    sample_seriation = "Identity" # suppress sample reordering
  )

htmp7e %>% ComplexHeatmap::draw(
  annotation_legend_list = attr(htmp7e, "AnnoLegends"),
)

# fix
str(tax_table(physeq_filtered2))
tax_table_df <- as.data.frame(tax_table(physeq_filtered2))
tax_table(physeq_filtered2) <- as.matrix(tax_table_df)
rownames(physeq_filtered2@tax_table)==colnames(physeq_filtered2@otu_table)


#Top 10% most abundant ASVs across all samples, at the Famly level
htmp7e <- physeq_filtered2 %>%
  #sort all samples by similarity
  ps_seriate(rank = "PhyOrdFam", tax_transform = "log10") %>% 
  # arrange 
  ps_arrange(Treatment.Type) %>% 
  tax_transform("log10", rank = "PhyOrdFam") %>%
  comp_heatmap(
    colors = heat_palette(palette = viridis::viridis(10)),
    row_names_side = "right", row_dend_side = "left", row_dend_reorder = TRUE, row_dend_width= grid::unit(5,"cm"), sample_side = "bottom",
    sample_anno = sampleAnnotation(
      Treatment.Type = anno_sample_cat(var = "Treatment.Type", col = c(Antibiotic_Treated = "#d81159", Control = "#b5e2fa"), 
                                       border_col = "black", box_col = NA, legend_title = "Treatment Type", size = grid::unit(4, "mm"))
    ),
    sample_seriation = "Identity" # suppress sample reordering
  )

htmp7e %>% ComplexHeatmap::draw(
  annotation_legend_list = attr(htmp7e, "AnnoLegends"),
)

#Top 10% most abundant ASVs across all samples, at the Genus Level
htmp7f <- physeq_filtered2 %>%
  #sort all samples by similarity
  ps_seriate(rank = "Genus", tax_transform = "log10") %>% 
  # arrange 
  ps_arrange(Treatment.Type) %>% 
  tax_transform("log10", rank = "Genus") %>%
  comp_heatmap(
    colors = heat_palette(palette = viridis::viridis(10)),
    row_names_side = "right", row_dend_side = "left", row_dend_reorder = TRUE, row_dend_width= grid::unit(5,"cm"), sample_side = "bottom",
    sample_anno = sampleAnnotation(
      Treatment.Type = anno_sample_cat(var = "Treatment.Type", col = c(Antibiotic_Treated = "#d81159", Untreated = "#b5e2fa"), 
                                       border_col = "black", box_col = NA, legend_title = "Treatment Type", size = grid::unit(4, "mm"))
    ),
    sample_seriation = "Identity" # suppress sample reordering
  )

htmp7f %>% ComplexHeatmap::draw(
  annotation_legend_list = attr(htmp7f, "AnnoLegends"),
)

