#AntiB and Caddie guts analyses

#load packages and data
library(microbiome)
library(vegetarian) 
library(phyloseq); packageVersion("phyloseq") #1.42.0
library(ggplot2); packageVersion("ggplot2") #3.5.1
library(gdata)
library(ecodist)
library(vegan)
library("car") 
library(dplyr)
library(biomformat)
library("ape") #phylogenetic tools packages
library(phytools) #phylogenetic tools packages
library(castor)
library(doParallel) #used for UniFrac but meh not super necessary for small data sets.
library(lme4)

#set seed
set.seed(45)

setwd("/Users/lab/Documents/Former_Researchers/Loomis_Caddies/Final_Qiime2_outputs")
#to load in the global environment and skip lines 26-248 run the line below.
#load(file="CaddiesSeqAnalysis072825.rda")

#Read in biom table from working directory
setwd("/Users/lab/Documents/Former_Researchers/Loomis_Caddies/Final_Qiime2_outputs")
ASV_reads<-read_biom("AntiB_Caddie_feature-table.biom")
ASV_table<-as.data.frame(as.matrix(biom_data(ASV_reads)))
ASV_table[1:10,1:10]
otu_tab<-t(ASV_table)
dim(otu_tab) #40 x 3814
otu_tab[1:10, 1:10]
rownames(otu_tab) #ordered as DAL1, 10, 11, 12, etc.. 

#Read in metadata file
meta_data<-read.csv("caddisfly-measurements-metadata1.csv",header=TRUE) #40 obs x 8 vars
head(meta_data)
str(meta_data)
meta_data<-meta_data[order(meta_data$Sample.Name),] 
dim(meta_data) #40 x 8
meta_data$Sample.Name #order matches line 31. 
meta_data$Sample.Name==rownames(otu_tab) #bien

#Read in taxonomic information. #First go into Excel, open tsv, and resave as csv.
bac.taxa<-read.csv(file="AntiB_Caddie_taxonomy1.csv")
str(bac.taxa)
rownames(bac.taxa) #whole numbers
colnames(bac.taxa) # "Feature.ID" "Domain"     "Phylum"     "Class"      "Order"      "Family"     "Genus"      "species"
rownames(bac.taxa)<-bac.taxa[,1] #assign the feature IDs to the row names
bac.taxa<-bac.taxa[,2:8] #get rid of the first column. 
colnames(bac.taxa) # "Domain"  "Phylum"  "Class"   "Order"   "Family"  "Genus"   "species"

#Reading in phylogenetic tree
dated.16Stree<-read.tree(file="AntiB_Caddie_Tree.nwk") #reading in tree, uses ape package
is.rooted(dated.16Stree) #TRUE
sample_names(dated.16Stree) #NULL
dated.16Stree$tip.label #for example "41e44b590a0d2c9f9c4e3ea3993981ac" 

#Creating phyloseq object
str(bac.taxa) #data.frame 3814 x 7 
bac.taxa.m<-as.matrix(bac.taxa) #VERY NECESSARY TO DO, DON'T SKIP. 
str(bac.taxa.m)
colnames(bac.taxa.m)
rownames(bac.taxa.m)
str(otu_tab) #37 x 3814
df.bacterial.OTU<-as.data.frame(otu_tab)
str(meta_data) #data.frame 37 x 8

#Matching row names
rownames(df.bacterial.OTU)<-as.character(meta_data[,1])
colnames(df.bacterial.OTU) #accession numbers
rownames(meta_data)<-as.character(meta_data[,1]) #sample names
rownames(bac.taxa.m)<-as.character(rownames(bac.taxa.m)) #these the accession numbers
samp.names<-as.character(meta_data[,1]) 

#Matching sample names (originally marked NULL)
sample_names(df.bacterial.OTU)<-samp.names
sample_names(meta_data)<-samp.names
sample_names(bac.taxa.m)<-samp.names
sample_names(dated.16Stree)<-samp.names

#Matching taxa names (originally marked NULL)
taxa_names(df.bacterial.OTU)<-colnames(df.bacterial.OTU)
taxa_names(dated.16Stree)<-colnames(df.bacterial.OTU)
taxa_names(meta_data)<-colnames(df.bacterial.OTU)
taxa_names(bac.taxa.m)<-colnames(df.bacterial.OTU)

#Here is the actual phyloseq object
Bacterial_phylo<-phyloseq(otu_table(df.bacterial.OTU, taxa_are_rows=FALSE), sample_data(meta_data), tax_table(bac.taxa.m), phy_tree(dated.16Stree))
Bacterial_phylo@otu_table[1:10,1:10]
dim(Bacterial_phylo@otu_table) #40 x 3814
rowSums(Bacterial_phylo@otu_table) #
sort(colSums(Bacterial_phylo@otu_table),decreasing=FALSE) # a handful of singletons ##singletons: a sequence that cannot be matched to its counter strand 

#Examining taxonomic ranks to examine where chloroplasts and mitochondria are nested within

table(tax_table(Bacterial_phylo)[, "Domain"], exclude = NULL) #29 unassigned
table(tax_table(Bacterial_phylo)[, "Phylum"], exclude = NULL) #237 unassigned, potentially of interest to remove
table(tax_table(Bacterial_phylo)[, "Class"], exclude = NULL) #examine #244 unassigned
table(tax_table(Bacterial_phylo)[, "Order"], exclude = NULL) #317 reads assigned as chloroplast at this taxonomic rank
#15 uncultured and 284 unassigned
table(tax_table(Bacterial_phylo)[, "Family"], exclude = NULL) #297 reads assigned as mitochondria at this taxonomic rank
table(tax_table(Bacterial_phylo)[, "Genus"], exclude = NULL) #25 unknown family, #683 unassigned
table(tax_table(Bacterial_phylo)[, "Species"], exclude = NULL) #1 unidentified marine bacteria, 2392 unassigned

#Replace taxa names to short hand
#taxa_names(Bacterial_phylo) <- paste0("ASV", seq(ntaxa(Bacterial_phylo)))

#remove other stuff
p1<-subset_taxa(Bacterial_phylo,  !Domain %in% "Unassigned") #requires dplyr (%in% syntax)
#p2<-subset_taxa(p1,  !Domain %in% "Archaea") #be no more!
p3<-subset_taxa(p1,  !Domain %in% "Eukaryota") #sayonara euks!
p4<-subset_taxa(p3,  !Order %in% "Chloroplast") #Chlorplasts be gone! 
p5<-subset_taxa(p4,  !Family %in% "Mitochondria") #Adios amigos!

# #double check removal
# table(tax_table(p5)[, "Domain"], exclude = NULL)
# table(tax_table(p5)[, "Phylum"], exclude = NULL) #could be beneficial to blastN the unassigned ones
# table(tax_table(p5)[, "Family"], exclude = NULL)
table(tax_table(p5)[, "Domain"], exclude = NULL) #arachaea 23 and Bacteria 3140


#time to remove the blanks and all the ASVs found in them from all of our samples. 
rowSums(Bacterial_phylo@otu_table)
rowSums(p5@otu_table)
dim(p5@otu_table) #40 x 3163

PCRBlank1<-subset_samples(p5, Treatment.Type=="BLANK1")
PCRBlank1.rm<-prune_taxa(taxa_sums(PCRBlank1) > 0, PCRBlank1)
bad_taxa1<-colnames(PCRBlank1.rm@otu_table) #length = 9

pop_taxa = function(physeq, badTaxa){
  allTaxa = taxa_names(physeq)
  allTaxa <- allTaxa[!(allTaxa %in% badTaxa)] #requires dplyr (%in% syntax)
  return(prune_taxa(allTaxa, physeq))
}

dim(p5@otu_table) #40 x 3140
p6 = pop_taxa(p5, bad_taxa1)
dim(p6@otu_table) #40 x 3153
str(p6@tax_table)
str(p6@phy_tree) #works!
rowSums(p6@otu_table)

PCRBlank2<-subset_samples(p6, Treatment.Type=="BLANK2")
PCRBlank2.rm<-prune_taxa(taxa_sums(PCRBlank2) > 0, PCRBlank2)
bad_taxa2<-colnames(PCRBlank2.rm@otu_table) #length = 4

dim(p6@otu_table) #40 x 3135
p7 = pop_taxa(p6, bad_taxa2)
dim(p7@otu_table) #40 x 3149; #4 less
rowSums(p7@otu_table)

PCRBlank3<-subset_samples(p6, Treatment.Type=="BLANK3")
PCRBlank3.rm<-prune_taxa(taxa_sums(PCRBlank3) > 0, PCRBlank3)
bad_taxa3<-colnames(PCRBlank3.rm@otu_table) #length = 13

dim(p7@otu_table) #40 x 3149
p8 = pop_taxa(p7, bad_taxa3)
dim(p8@otu_table) #40 x 3138
rowSums(p8@otu_table)

#Remove samples with low read depth
CaddieD_phylo<-subset_samples(p8, Sample.Name. != "NEGATIVE_CONTROL_1" & Sample.Name. != "NEGATIVE_CONTROL_2" & Sample.Name. != "NEGATIVE_CONTROL_3")
CaddieD_phylo<-subset_taxa(CaddieD_phylo,  !Phylum %in% "Unassigned") #be no more!
CaddieD_phylo.rm<-prune_taxa(taxa_sums(CaddieD_phylo) > 0, CaddieD_phylo) #this is the object that we want to work with


dim(CaddieD_phylo.rm@sam_data) #37 x 8
CaddieD_meta_pruned<-CaddieD_phylo.rm@sam_data

dim(CaddieD_phylo.rm@otu_table) #37 x 2934

#ASV table with out zeros and singletons ## I dont think this part of the code is significant for me since my read depths were good
CaddieD_phylo.rm2<- prune_taxa(taxa_sums(CaddieD_phylo.rm) > 1, CaddieD_phylo.rm) #remove singletons
dim(CaddieD_phylo.rm2@otu_table) #37 x 2917
obj45<-rowSums(CaddieD_phylo.rm2@otu_table)

CaddieD_phylo.rm2.ASVtable<-CaddieD_phylo.rm2@otu_table

#Rarefying
##accounting for uneven sampling efforts and biases during the sampling process such as if one sample gets amplifyed a lot more for some reason compared to other samples
##will rarefy for the lowest number of reads in our case 
min(rowSums(CaddieD_phylo.rm2.ASVtable)) #28484 
mean(rowSums(CaddieD_phylo.rm2.ASVtable)) #56,389.62
median(rowSums(CaddieD_phylo.rm2.ASVtable)) #51,430
sort(x=rowSums(CaddieD_phylo.rm2.ASVtable), decreasing=FALSE)
bac.tab.df<-as.data.frame(CaddieD_phylo.rm2.ASVtable)
rdat<-rrarefy(bac.tab.df,28484) #rarefy! to minimum number of read depths. 

#Standardize abundances into proportions. 
rowSums(rdat) #Check if it worked
std.bac.tab<-decostand(rdat,"total") #replace object to rdat after rarefying. 
std.bac.tab[1:10,1:10] #Looks groovy! 

#Smush back together into single phyloseq object; this contains rarefied and standardized data
bacterial.phylo.4analysis<-phyloseq(otu_table(std.bac.tab, taxa_are_rows=FALSE), sample_data(CaddieD_phylo.rm2@sam_data), tax_table(CaddieD_phylo.rm2@tax_table), phy_tree(CaddieD_phylo.rm2@phy_tree))

#Create Quant Jaccard distance matrix
drdat<-vegdist(std.bac.tab,"bray")

#metadata
str(bacterial.phylo.4analysis@sam_data)
bacterial.phylo.4analysis@otu_table[1:10,1:10]
bacterial.phylo.4analysis@sam_data$Treatment.Type[bacterial.phylo.4analysis@sam_data$Treatment.Type=="treated"]<-"Antibiotic Treated"
bacterial.phylo.4analysis@sam_data$Treatment.Type[bacterial.phylo.4analysis@sam_data$Treatment.Type=="untreated"]<-"Control"
bacterial.phylo.4analysis@sam_data$Treatment.Type<-factor(bacterial.phylo.4analysis@sam_data$Treatment.Type, levels=c("Antibiotic Treated","Control"))
#bacterial.phylo.4analysis@sam_data$Filter.type<-factor(bacterial.phylo.4analysis@sam_data$Filter.type)

Caddie_meta<-bacterial.phylo.4analysis@sam_data
str(Caddie_meta)
str(Caddie_meta) #37 x 8
colnames(Caddie_meta)

basic.mod<-dbrda(drdat~Caddie_meta$Treatment.Type)
h<-how(nperm=10000)
anova(basic.mod, permutations = h, by="margin")
# Model: dbrda(formula = drdat ~ Caddie_meta$Treatment.Type)
#                              Df SumOfSqs      F    Pr(>F)    
# Caddie_meta$Treatment.Type   1   4.3778 58.643 9.999e-05 ***
#   Residual                   35  2.6128                     

wdistUnif<-UniFrac(physeq = bacterial.phylo.4analysis, weighted = TRUE, parallel = TRUE, fast = TRUE)
wUnimod<-dbrda(wdistUnif~Caddie_meta$Treatment.Type) #this model matches 177 on the right side of the ~
h<-how(nperm=10000)
anova(wUnimod,permutations=h,by="margin")
# Model: dbrda(formula = wdistUnif ~ Caddie_meta$Treatment.Type)
# Df SumOfSqs      F    Pr(>F)    
# Caddie_meta$Treatment.Type  1  0.45789 21.912 9.999e-05 ***
#   Residual                   35  0.73138                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

uwdistUnif<-UniFrac(physeq = bacterial.phylo.4analysis, weighted = FALSE, parallel = TRUE, fast = TRUE)
uwUnimod<-dbrda(uwdistUnif~Caddie_meta$Treatment.Type)
h<-how(nperm=10000)
anova(uwUnimod,permutations=h,by="margin")
# Model: dbrda(formula = uwdistUnif ~ Caddie_meta$Treatment.Type)
# Df SumOfSqs      F    Pr(>F)    
# Caddie_meta$Treatment.Type  1  0.79845 11.055 9.999e-05 ***
#   Residual                   35  2.52797                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
       
---
#Unconstrained Ordinations 
hfa_bac_pcoa <- ordinate(physeq = bacterial.phylo.4analysis, method = "PCoA", distance = "bray")
hfa_jac_pcoa <- ordinate(physeq = bacterial.phylo.4analysis, method = "PCoA", distance = "jaccard") 
wu_hfa_bac_pcoa <- ordinate(physeq = bacterial.phylo.4analysis, method = "PCoA", distance = "wunifrac") 
uwu_hfa_bac_pcoa <- ordinate(physeq = bacterial.phylo.4analysis, method = "PCoA", distance = "unifrac")

#Plot bray curtis pcoa 
bc1<-plot_ordination(physeq = bacterial.phylo.4analysis, ordination = hfa_bac_pcoa, shape= "Treatment.Type", color = "Treatment.Type", axes = c(1,2)) + 
  scale_fill_manual(values = c("#d81159","#b5e2fa")) +
  scale_color_manual(values =c("black","black"), guide = "none") +
  scale_shape_manual(values = c(22,21)) +
  scale_x_continuous(limits=c(-0.42,0.42),breaks=c(-0.40,-0.20,0.0,0.20,0.4))+
  scale_y_continuous(limits=c(-0.32,0.47),breaks=c(-0.3,-0.15,0.0,0.15,0.30,0.45))+
  geom_point(aes(fill = Treatment.Type), alpha = 1, size = 4.5, stroke=1)+
  theme_test()+
  labs(color="Treatment Type")+ labs(shape = "Treatment Type") + labs(fill = "Treatment Type") +
  labs(x="Axis 1 [63.3%] Bray-Curtis",y="Axis 2 [9.7%] Bray-Curtis")+
  theme(text = element_text(size = 15))+
  theme(legend.position = "right")
bc1$layers<-bc1$layers[-1]

bc1a <- bc1 + theme(
  axis.text.x = element_text(size = 12),
  axis.title.x = element_text(size = 15),
  axis.text.y = element_text(size = 12),
  axis.title.y = element_text(size = 15),
  legend.position = "right")


#plot jaccard pcoa
bc2<-plot_ordination(physeq = bacterial.phylo.4analysis, ordination = hfa_jac_pcoa, shape = "Treatment.Type", color = "Treatment.Type", axes = c(1,2)) + 
  scale_fill_manual(values = c("#d81159","#b5e2fa")) +
  scale_color_manual(values =c("black","black")) +
  scale_shape_manual(values = c(22,21)) +
  scale_x_continuous(limits=c(-0.42,0.42),breaks=c(-0.40,-0.20,0.0,0.20,0.4))+
  scale_y_continuous(limits=c(-0.27,0.52),breaks=c(-0.25,0.0,0.25,0.50))+
  geom_point(aes(fill = Treatment.Type), alpha = 1, size = 6, stroke=1)+
  theme_test()+
  labs(color="Treatment Type")+ labs(shape = "Treatment Type") + labs(fill = "Treatment Type") +
  labs(x="Axis 1 [33.6%]\n Quantitative Jaccard",y="Axis 2 [8.3%]\n Quantitative Jaccard")+
  theme(text = element_text(size = 15))+
  theme(legend.position = "right")
bc2$layers<-bc2$layers[-1]

# bc2a <- bc2 + theme(
#   axis.text.x = element_text(size = 12),
#   axis.title.x = element_text(size = 15),
#   axis.text.y = element_text(size = 12),
#   axis.title.y = element_text(size = 15),
#   legend.position = "none")

library(scales)
#weighted unifrac
bc3<-plot_ordination(physeq = bacterial.phylo.4analysis, ordination = wu_hfa_bac_pcoa, shape = "Treatment.Type", color = "Treatment.Type", axes = c(1,2)) + 
  scale_fill_manual(values = c("#d81159","#b5e2fa")) +
  scale_color_manual(values =c("black","black"), guide = "none") +
  scale_shape_manual(values = c(22,21)) +
  scale_x_continuous(limits=c(-0.4,0.4),breaks=c(-0.4,-0.2,0.0,0.2,0.4),labels=label_number(accuracy = 0.01))+
  scale_y_continuous(limits=c(-0.3,0.3),breaks=c(-0.3,-0.15,0.0,0.15,0.3))+
  geom_point(aes(fill = Treatment.Type), alpha = 1, size = 6, stroke=1)+
  theme_test()+
  labs(color="Treatment Type")+ labs(shape = "Treatment Type") + labs(fill = "Treatment Type") +
  labs(x="Axis 1 [48.6%] Weighted-UniFrac",y="Axis 2 [28.8%] Weighted-UniFrac")+
  theme(text = element_text(size = 15))+
  theme(legend.position = "right") 
bc3$layers<-bc3$layers[-1]

bc3a <- bc3 + theme(
  axis.text.x = element_text(size = 12),
  axis.title.x = element_text(size = 15),
  axis.text.y = element_text(size = 12),
  axis.title.y = element_text(size = 15),
  legend.position = "none")


#unweighted unifrac
bc4<-plot_ordination(physeq = bacterial.phylo.4analysis, ordination = uwu_hfa_bac_pcoa, shape = "Treatment.Type", color = "Treatment.Type", axes = c(1,2)) + 
  scale_fill_manual(values = c("#d81159","#b5e2fa")) +
  scale_color_manual(values =c("black","black"), guide = "none") +
  scale_shape_manual(values = c(22,21)) +
  scale_x_continuous(limits=c(-0.3,0.3),breaks=c(-0.3,-0.15,0.0,0.15,0.3))+
  scale_y_continuous(limits=c(-0.15,0.3),breaks=c(-0.15,0.0,0.15,0.3))+
  geom_point(aes(fill = Treatment.Type), alpha = 1, size = 6, stroke=1) +
  theme_test()+
  labs(color="Gut Treatment")+ labs(shape = "Gut Treatment") + labs(fill = "Gut Treatment") +
  labs(x="Axis 1 [24.9%] Unweighted-UniFrac",y="Axis 2 [7.9%] Unweighted-UniFrac")+
  theme(text = element_text(size = 15))+
  theme(legend.position = "right") 
bc4$layers<-bc4$layers[-1]

bc4a <- bc4 + theme(
  axis.text.x = element_text(size = 12),
  axis.title.x = element_text(size = 15),
  axis.text.y = element_text(size = 12),
  axis.title.y = element_text(size = 15),
  legend.position = "right",
  legend.key.size = unit(1,"cm"),
  legend.text = element_text(size = 12),
  legend.title = element_text(size = 15),
  legend.box.margin = margin(0,10,0,0))

library(cowplot)
plot_grid(bc1a,bc3a,bc4a,nrow=1,ncol=3,rel_widths=c(1,1,1.5),labels= c("A","B","C"))

plot_grid(bc3a,bc4a,nrow=1,ncol=2,rel_widths=c(1,1.5),labels= c("A","B"))

#plot unweighted unifrac pcoa
# bc3<-plot_ordination(physeq = bacterial.phylo.4analysis, ordination = uwu_hfa_bac_pcoa, color = "Water_Source", shape = "Timepoint", axes = c(1,2)) + 
#   scale_fill_manual(values = c("#060282","#5aafed","#ffc001","#ff0000")) +
#   scale_color_manual(values =c("black","black","black","black")) +
#   scale_shape_manual(values = c(21,22,23,24,25)) +
#   scale_x_continuous(limits=c(-0.3,0.3),breaks=c(-0.3,-0.15,0.0,0.15,0.3))+
#   #scale_y_continuous(limits=c(-0.3,0.4),breaks=c(-0.3,-0.1,0.0,0.1,0.4))+
#   geom_point(aes(fill = Water_Source), alpha = 1, size = 4, stroke=1)+
#   theme_test()+
#   labs(shape="Time Point",fill="Origin")+
#   labs(x="Axis 1 [16.3%]",y="Axis 2 [14.9%]")+
#   theme(text = element_text(size = 15))+
#   theme(legend.position = "none")
# bc3$layers<-bc3$layers[-1]

#save(list=ls(),file="AntiB_Caddie_16SAnalyses.rda")


### PERMANOVA TUTORIAL ###
# Calculate bray curtis distance matrix
caddie_bray <- phyloseq::distance(bacterial.phylo.4analysis, method = "bray")

# make a data frame from the sample_data
caddie_sampledf <- data.frame(sample_data(Caddie_meta))
str(caddie_sampledf)
adonis2(caddie_bray ~ Treatment.Type, data = caddie_sampledf)

#Permutation test for adonis under reduced model
#adonis2(formula = caddie_bray ~ Treatment.Type, data = caddie_sampledf)
#               Df SumOfSqs      R2      F Pr(>F)    
#Treatment.Type  1   4.3826 0.61914 56.897  0.001 ***
#Residual       35   2.6959 0.38086                  
#Total          36   7.0785 1.00000  

# Homogeneity of dispersion test
dist.hill<-function(dat,q=1){
  dmat<-array(dim=c(dim(dat)[1],dim(dat)[1]))
  for(i in 1:dim(dmat)[1]){
    for(j in 1:dim(dmat)[1]){
      dmat[i,j]<-d(rbind(dat[i,],dat[j,]),lev="beta",q=q)-1
    }
  }
  res<-as.dist(dmat)
  return(res)
}

caddie_sampledf$Treatment.Type<-as.character(caddie_sampledf$Treatment.Type)
#caddie_sampledf$Treatment.Type[caddie_sampledf$Treatment.Type=="Untreated"]<-"Control"
caddie_sampledf$Treatment.Type<-factor(caddie_sampledf$Treatment.Type, levels=c("Antibiotic Treated","Control"))
ddatq0<-dist.hill(bacterial.phylo.4analysis@otu_table,q=0) #pairwise beta-diversity
caddiebeta <- betadisper(ddatq0, caddie_sampledf$Treatment.Type)
mod.beta<-lm(caddiebeta$distances ~ caddiebeta$group)
Anova(mod.beta)

# Anova Table (Type II tests)
# 
# Response: caddiebeta$distances
# Sum Sq Df F value    Pr(>F)    
# caddiebeta$group 0.065060  1  46.772 6.171e-08 ***
#   Residuals        0.048686 35                      

### TRYING ALPHA DIVERSITY A DIFFERENT WAY ###
alpha0 <- apply(rdat, 1, vegetarian::d, q =0)
alpha1 <- apply(rdat, 1, vegetarian::d, q =1)
alpha2 <- apply(rdat, 1, vegetarian::d, q =2)

#attaching alpha diversity values to the meta data 
caddie_alpha_meta <- cbind(Caddie_meta, alpha0, alpha1, alpha2)

str(caddie_alpha_meta)
# caddie_alpha_meta$Treatment.Type<-as.character(caddie_alpha_meta$Treatment.Type)
# caddie_alpha_meta$Treatment.Type[caddie_alpha_meta$Treatment.Type=="Untreated"]<-"Control"

#creating a linear model where the effects of whether or not antibiotics were added
caddie_lm_alpha <- lm(alpha0 ~ Treatment.Type, dat = caddie_alpha_meta)
caddie_lm_alpha1<- lm(alpha1 ~ Treatment.Type, dat = caddie_alpha_meta)
caddie_lm_alpha2<- lm(alpha2 ~ Treatment.Type, dat = caddie_alpha_meta)

#running ANOVA for the model just created 
Anova(caddie_lm_alpha) #richness
# Anova Table (Type II tests)
# 
# Response: alpha0
# Sum Sq Df F value  Pr(>F)  
# Treatment.Type   6626  1  4.1845 0.04837 *
#   Residuals       55425 35                  

Anova(caddie_lm_alpha1) #Exponential Shannon
# Anova Table (Type II tests)
# 
# Response: alpha1
# Sum Sq Df F value    Pr(>F)    
# Treatment.Type 320.52  1  15.107 0.0004327 ***
#   Residuals      742.59 35                      

Anova(caddie_lm_alpha2) #Inverse simpsons index
# Anova Table (Type II tests)
# 
# Response: alpha2
# Sum Sq Df F value    Pr(>F)    
# Treatment.Type  55.552  1  13.674 0.0007409 ***
#   Residuals      142.191 35                      

#creating a boxplot figure 
a0.plot <- ggplot(data = caddie_alpha_meta, aes(x = Treatment.Type, y = alpha0, fill = Treatment.Type)) +
  geom_boxplot(alpha = 0.65, outlier.size = 0) +
  geom_point(aes(shape = Treatment.Type), size = 3.5, position = position_jitterdodge(), color = "black", stroke = 0.5) +
  scale_fill_manual(values = c("#d81159", "#b5e2fa")) +
  scale_shape_manual(values = c(22, 21)) +  # 22 = square filled, 21 = circle filled
  scale_y_continuous(limits = c(150, 350), breaks = c(150, 200, 250, 300, 350)) +
  theme_classic() +
  labs(y = "ASV Richness", x = "Gut Treatment") +
  theme(legend.position = "right")+
  guides(shape="none")

a0.plota <- a0.plot + labs(fill = "Gut Treatment")+ 
  theme(
  axis.title.y = element_text(size=15),
  axis.title.x = element_text(size=15),
  axis.text.x = element_text(size=12),
  axis.text.y = element_text(size=12),
  legend.key.size = unit(1,"cm"),
  legend.text = element_text(size = 12),
  legend.title = element_text(size = 15),
  legend.box.margin = margin(0,10,0,0)
)

a1.plot <- ggplot(data = caddie_alpha_meta, aes(x = Treatment.Type, y = alpha1, fill = Treatment.Type)) +
  geom_boxplot(aes(fill = Treatment.Type), alpha = 0.65, outlier.size = 0) +
  geom_point(aes(shape = Treatment.Type), size = 3.5, position = position_jitterdodge(), color = "black", stroke = 0.5) +
  scale_fill_manual(values = c("#d81159","#b5e2fa")) +
  scale_shape_manual(values = c(22, 21)) +  # 22 = square filled, 21 = circle filled
  scale_y_continuous(limits=c(0,30),breaks=c(0,10,20,30))+
  theme_classic()+
  labs(y = "Expontential Shannon's Diversity", x = "Gut Treatment") +
  theme(legend.position = "none")

a1.plota<-a1.plot + labs(fill = "Gut Treatment") +  labs(fill = "Gut Treatment")+ 
  theme(
    axis.title.y = element_text(size=15),
    axis.title.x = element_text(size=15),
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12),
    legend.key.size = unit(1,"cm"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 15),
    legend.box.margin = margin(0,10,0,0)
  )


a2.plot <- ggplot(data = caddie_alpha_meta, aes(x = Treatment.Type, y = alpha2, fill = Treatment.Type)) +
  geom_boxplot(aes(fill = Treatment.Type), alpha = 0.65, outlier.size = 0) +
  geom_point(aes(shape = Treatment.Type), size = 3.5, position = position_jitterdodge(), color = "black", stroke = 0.5) +
  scale_fill_manual(values = c("#d81159","#b5e2fa")) +
  scale_shape_manual(values = c(22, 21)) +  # 22 = square filled, 21 = circle filled
  scale_y_continuous(limits=c(0,15),breaks=c(0,5,10,15))+
  theme_classic()+
  labs(y = "Inverse Simpson's Diversity", x = "Gut Treatment") +
  theme(legend.position = "right")+
  guides(shape = "none")

a2.plota<-a2.plot + labs(fill = "Gut Treatment") +  labs(fill = "Gut Treatment")+ 
  theme(
    axis.title.y = element_text(size=15),
    axis.title.x = element_text(size=15),
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12),
    legend.key.size = unit(1,"cm"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 15),
    legend.box.margin = margin(0,10,0,0)
  )

library(cowplot)
plot_grid(a0.plota,a1.plota,a2.plota,nrow=1,ncol=3,rel_widths=c(1,1,1.5),labels= c("A","B","C"))
plot_grid(a1.plota,a2.plota,nrow=1,ncol=2,rel_widths=c(1,1.5),labels= c("A","B"))

#CaddieD_phylo.rm2
depth<-rowSums(CaddieD_phylo.rm2@otu_table)
mmmmme<-CaddieD_phylo.rm2@sam_data
klklklk<-cbind(mmmmme,depth)
str(klklklk)
klklklk$Treatment.Type[klklklk$Treatment.Type=="treated"]<-"Antibiotic Treated"
klklklk$Treatment.Type[klklklk$Treatment.Type=="untreated"]<-"Control"
klklklk$Treatment.Type<-factor(klklklk$Treatment.Type, levels=c("Antibiotic Treated","Control"))

rd.plot <- ggplot(data = klklklk, aes(x = Treatment.Type, y = depth, fill = Treatment.Type)) +
  geom_boxplot(aes(fill = Treatment.Type), alpha = 0.65, outlier.size = 0) +
  geom_point(aes(shape = Treatment.Type), size = 3.5, position = position_jitterdodge(), color = "black", stroke = 0.5) +
  scale_shape_manual(values = c(22, 21)) +
  scale_fill_manual(values = c("#d81159","#b5e2fa")) +
  scale_y_continuous(limits=c(25000,125000),breaks=c(25000,50000,75000,100000,125000))+
  theme_classic()+
  labs(y = "Raw Read Depth", x = "Gut Treatment") +
  theme(legend.position = "right")

rd.plota<-rd.plot  + labs(fill = "Gut Treatment") +  labs(fill = "Gut Treatment")+ 
  theme(
    axis.title.y = element_text(size=15),
    axis.title.x = element_text(size=15),
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12))

klklklk.mod <- lm(depth ~ Treatment.Type, dat = klklklk)
Anova(klklklk.mod)
# Anova Table (Type II tests)
# 
# Response: depth
# Sum Sq Df F value    Pr(>F)    
# Treatment.Type 4.5646e+09  1  14.977 0.0004539 ***
#   Residuals      1.0667e+10 35                      

klklklk.antiB<-klklklk[which(klklklk$Treatment.Type=="Antibiotic Treated"),]
klklklk.control<-klklklk[which(klklklk$Treatment.Type=="Control"),]

mean(klklklk.antiB$depth) #45,578.74
mean(klklklk.control$depth) #67,801.11

dim(klklklk.antiB) #19 x 9
dim(klklklk.control) #18 x 9

sd(klklklk.antiB$depth) #10,936.66
sd(klklklk.control$depth) #22,379.12

10936.66/(sqrt(19)) #2509.042
22379.12/(sqrt(18)) #5274.809

#venn diagram for treatments using rarefied phyloseq object , bacterial.phylo.4analysis
library(ggvenn)
#pull out your two treats
antiB_cadd<-subset_samples(bacterial.phylo.4analysis, Treatment.Type=="Antibiotic Treated")
antiB_cadd.rm<-prune_taxa(taxa_sums(antiB_cadd) > 0, antiB_cadd)

Cont_cadd<-subset_samples(bacterial.phylo.4analysis, Treatment.Type=="Control")
Cont_cadd.rm<-prune_taxa(taxa_sums(Cont_cadd) > 0, Cont_cadd)


#pull out just the colnames
lowASV<-colnames(antiB_cadd.rm@otu_table)
medASV<-colnames(Cont_cadd.rm@otu_table)

#make list and plot fig
new_list<-list("Antibiotic Treated" = lowASV, "Control" = medASV)
ggvenn(new_list, fill_color = c("#d81159","#b5e2fa"), fill_alpha = 0.8, stroke_size = 0.8, set_name_size = 6, text_size = 6, show_percentage = FALSE)

#save(list=ls(),file="CaddiesSeqAnalysis072825.rda")
