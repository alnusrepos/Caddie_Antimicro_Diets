#Load libraries
library("phyloseq"); packageVersion("phyloseq") #1.42.0
library("edgeR"); packageVersion("edgeR") #4.2.2
library("ggplot2"); packageVersion("ggplot2") #3.5.1
library(forcats)
library(dplyr)

#internal functions that are reqiured
phyloseq_to_edgeR = function(physeq, group, method="RLE", ...){
  require("edgeR")
  require("phyloseq")
  # Enforce orientation.
  if( !taxa_are_rows(physeq) ){ physeq <- t(physeq) }
  x = as(otu_table(physeq), "matrix")
  # Add one to protect against overflow, log(0) issues.
  x = x + 1
  # Check `group` argument
  if( identical(all.equal(length(group), 1), TRUE) & nsamples(physeq) > 1 ){
    # Assume that group was a sample variable name (must be categorical)
    group = get_variable(physeq, group)
  }
  # Define gene annotations (`genes`) as tax_table
  taxonomy = tax_table(physeq, errorIfNULL=FALSE)
  if( !is.null(taxonomy) ){
    taxonomy = data.frame(as(taxonomy, "matrix"))
  } 
  # Now turn into a DGEList
  y = DGEList(counts=x, group=group, genes=taxonomy, remove.zeros = TRUE, ...)
  # Calculate the normalization factors
  z = calcNormFactors(y, method=method)
  # Check for division by zero inside `calcNormFactors`
  if( !all(is.finite(z$samples$norm.factors)) ){
    stop("Something wrong with edgeR::calcNormFactors on this data,
         non-finite $norm.factors, consider changing `method` argument")
  }
  # Estimate dispersions
  return(estimateTagwiseDisp(estimateCommonDisp(z)))
}

scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}

#set seed
set.seed(45)

#load data
setwd("/Users/lab/Documents/Former_Researchers/Loomis_Caddies/Final_Qiime2_outputs")
#to load in the global environment and skip lines 26-248 run the line below.
load(file="CaddiesSeqAnalysis072825.rda")

#phyloseq object that contains a non-rarefied, non-standardized ASV table. Remove unassigned Phyla and Families
edgeR.phyloseq.1<-phyloseq(otu_table(CaddieD_phylo.rm2@otu_table, taxa_are_rows=FALSE), sample_data(CaddieD_phylo.rm2@sam_data), tax_table(CaddieD_phylo.rm2@tax_table), phy_tree(CaddieD_phylo.rm2@phy_tree)) #replaced rarefied data (rdat) with non rarefied. 
table(tax_table(edgeR.phyloseq.1)[, "Phylum"], exclude = NULL) #no unassigned phyla
table(tax_table(edgeR.phyloseq.1)[, "Family"], exclude = NULL)
p1.a<-subset_taxa(edgeR.phyloseq.1,  !Family %in% "Unassigned")
table(tax_table(p1.a)[, "Family"], exclude = NULL)
p1.b<-subset_taxa(p1.a,  !Family %in% "Unknown_Family")
table(tax_table(p1.b)[, "Family"], exclude = NULL)
p1.c<-subset_taxa(p1.b,  !Family %in% "uncultured")
table(tax_table(p1.c)[, "Family"], exclude = NULL)
p1.a
p1.b
p1.c
table(tax_table(p1.c)[, "Domain"], exclude = NULL) #21 taxa present as Archaea. *thumbs up emoji*

#fixing metadata
p1.c@sam_data$Treatment.Type
p1.c@sam_data$Treatment.Type[p1.c@sam_data$Treatment.Type=="treated"]<-"Antibiotic Treated"
p1.c@sam_data$Treatment.Type[p1.c@sam_data$Treatment.Type=="untreated"]<-"Untreated" #leaving as such because i do not want to mess up the order of the figure downstream. 
p1.c@sam_data$Treatment.Type<-factor(p1.c@sam_data$Treatment.Type, levels=c("Antibiotic Treated","Untreated"))


#breaking apart two treatments to understand dimensions
p1.c.antib<-subset_samples(p1.c, Treatment.Type !="Antibiotic Treated")
p1.c.antib.rm<-prune_taxa(taxa_sums(p1.c.antib) > 0, p1.c.antib) 
p1.c.untreated<-subset_samples(p1.c, Treatment.Type !="Untreated")
p1.c.untreated.rm<-prune_taxa(taxa_sums(p1.c.untreated) > 0, p1.c.untreated) 
dim(p1.c.antib.rm@otu_table) #18 x 1518
dim(p1.c.untreated.rm@otu_table) #19 x 1547

#independent filtering
p1.c #37 samples by 2720 taxa; nice (07/30/25)

#transform sample counts and view variance per ASV
edgeR.phyloseq.p = transform_sample_counts(p1.c, function(x){x/sum(x)})

#view
hist(log10(apply(otu_table(edgeR.phyloseq.p), 1, var)),
     xlab="log10(variance)", breaks=50,
     main="A large fraction of OTUs have very low variance")

#apply taxa filtering
AbundThreshold = 0.0001
f4Keeps= names(which(apply(otu_table(edgeR.phyloseq.p),2, sum) > AbundThreshold))
edgeR.phyloseq.f = prune_taxa(f4Keeps, p1.c)

#edgeR.phyloseq.f<-prune_taxa(taxa_sums(edgeR.phyloseq.p) > 0.001, edgeR.phyloseq.p)
edgeR.phyloseq.f #1673 taxa by 37 samples

#now lets use this newly defined function to transform phyloseq to edgeR's "DGE" type object
dge = phyloseq_to_edgeR(edgeR.phyloseq.f, group="Treatment.Type") #works!

# Perform binary test
et = exactTest(dge)

# Extract values from test results
tt = topTags(et, n=nrow(dge$table), adjust.method="BH", sort.by="PValue")
res = tt@.Data[[1]]
alpha_edger = 0.001 #0.001, 0.01, 0.05 
sigtab = res[(res$FDR < alpha_edger), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(edgeR.phyloseq.f)[rownames(sigtab), ], "matrix"))
dim(sigtab) #255 x 18
head(sigtab)

#figure time! 
theme_set(theme_bw())
sigtabgen = subset(sigtab, !is.na(Family))

# Phylum order
x = tapply(sigtabgen$logFC, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels = names(x))

# Family order
x = tapply(sigtabgen$logFC, sigtabgen$Family, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Family = factor(as.character(sigtabgen$Family), levels = names(x))

ggplot(sigtabgen, aes(x = logFC, y = Family, color = Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = , hjust = 0, vjust = 0.5))

#write.csv(sigtabgen, "Caddie_EdgeR_07302025.csv") #write csv to add shape column for negative values and positive log fold  change values

#more complex figure
custom_col42 = c("#781156","#A51876","#D21E96","#E43FAD","#EA6CC0","#F098D3","#114578","#185EA5","#1E78D2","#3F91E4","#6CABEA","#98C4F0","#117878","#18A5A5","#3FE4E4","#6CEAEA","#98F0F0", "#117845","#18A55E","#1ED278","#3FE491","#6CEAAB","#98F0C4","#787811","#A5A518","#D2D21E","#E4E43F","#EAEA6C","#F0F098","#F7F7C5","#784511","#A55E18","#D2781E","#E4913F","#EAAB6C","#F0C498","#781122","#A5182F","#D21E2C","#E43F5B","#EA6C81","#F098A7","black")

Phylum_cols<-c("#781156","#ff0000","#ff4d00","#ff7400","#ff9a00","#ffc100","#117845","#1ED278","#6CEAAB","#117878","#3FE4E4","#98F0F0","#0200b9","#007dff","#2b89ea") #15

Rain15<-c("#ff00cc","#ee34d2","#9c27b0","#50bfe6","#16d0cb","#aaf0d1","#66ff66","#ccff00","#ffff66","#ffcc33","#ff9933","#ff9966","#ff6037","#fd5b78","#ff355e")
Adapted15<-c("#4c00a4","#9c27b0","#ff00cc","#114578","#117878","#3FE4E4","#004242","#117845","#66ff66","#ffc100","#ff9005","#f9530b","#ff5577","#ff0000",	"#781156")
Adapted13b<-c("#781156","#ff0000","#ff5577","#f9530b","#ff9005","#ffc100","#004242","#117845","#66ff66","#114578","#117878","#3FE4E4","#4c00a4","#9c27b0","#ff00cc")
Adapted13c<-c("#781156","#ff0000","#ff5577","#f24b04","#ff9005","#ffc100","#004242","#117845","#0200b9","#007dff","#3FE4E4","#531cb3","#9c27b0","#ff00cc")
Adapted13<-c("#4c00a4","#9c27b0","#ff00cc","#114578","#3FE4E4","#004242","#117845","#66ff66","#f9530b","#ff9005","#ffc100","#781156","#ff5577")
Adapted10<-c("#4c00a4","#ff00cc","#114578","#3FE4E4","#004242","#66ff66","#ffc100","#f9530b","#ff0000","#781156")
Adapted13d<-c("#781156","#ff00cc","#ff5400","#ff9005","#fcf300","#004b23","#38b000","#072ac8","#1e96fc","#3FE4E4","#531cb3","#c86bfa", "#ffcbf2")


fig.data.edgR<-read.csv(file="Caddie_EdgeR_07302025.csv")
fig.data.edgR$Phylum<-factor(fig.data.edgR$Phylum,levels=sort(unique(fig.data.edgR$Phylum)))
fig.data.edgR$Order<-factor(fig.data.edgR$Order,levels=sort(unique(fig.data.edgR$Order)))
fig.data.edgR$Family<-factor(fig.data.edgR$Family,levels=sort(unique(fig.data.edgR$Family)))
fig.data.edgR$FAM_ASV<-factor(fig.data.edgR$FAM_ASV,levels=sort(unique(fig.data.edgR$FAM_ASV)))

fig.data.edgR$Shape<-as.factor(fig.data.edgR$Shape)
fig.data.edgR$Shape1<-as.factor(fig.data.edgR$Shape1)
fig.data.edgR$Treat<-factor(fig.data.edgR$Treat,levels=c("Antibiotic Treated", "Control"))

ggplot(fig.data.edgR, aes(x = logFC, y = Family, color = Phylum)) + geom_point(size=4) + 
  geom_vline(xintercept = c(0,0,0,0),linetype="dashed",colour="black",size=0.5)+
  theme(axis.text.x = element_text(angle = , hjust = 0, vjust = 0.5))+
  theme_bw() +
  scale_color_manual(values=Phylum_cols) +
  labs(y="Bacterial Family",x="Log Fold Change") +
  theme(legend.position = "right")+ labs(color="Order")+ theme(text = element_text(size = 15))

shape1<-c(22,21)



#perfect figure with color by phylum
# fig.data.edgR %>%
#   mutate(FAM_ASV = fct_inorder(FAM_ASV)) %>%
#   ggplot(aes(x = logFC, y = FAM_ASV)) + 
#   geom_point(aes(fill=Phylum,shape=Shape1,color=Phylum), color="black", stroke=1, size=3.50, alpha=1) + 
#   geom_vline(xintercept = c(0,0,0,0),linetype="dashed",colour="black",size=0.5)+
#   theme(axis.text.x = element_text(angle = , hjust = 0, vjust = 0.5))+
#   theme_bw() +
#   scale_shape_manual(values=shape1)+
#   scale_fill_manual(values=Adapted13c)+
#   scale_color_manual(values=Adapted13c)+
#   scale_y_discrete(limits=rev) +
#   scale_x_continuous(limits=c(-10,10),breaks=c(-10.0,-7.5,-5.0,-2.5,0.0,2.5,5.0,7.5,10.0))+
#   #geom_point(data = fig.data.edgR[fig.data.edgR$Shape == 21, ], aes(x = logFC, y = FAM_ASV,fill=Phylum),stroke=1, pch=21,size=3.50,color="black")+
#   labs(y="Bacterial Taxonomy",x="Log Fold Change") +
#   theme(legend.position = "right")+ labs(color="Phylum")+ theme(text = element_text(size = 14))
library(tidyverse)
edgE<- fig.data.edgR %>%
  mutate(FAM_ASV = fct_inorder(FAM_ASV)) %>%
  ggplot(aes(x = logFC, y = FAM_ASV)) + 
  geom_point(aes(fill = Phylum,  shape = Treat),
             stroke = 1, size = 3.5, alpha = 1) +
  
  geom_vline(xintercept = 0, linetype = "dashed", colour = "black", size = 0.5) +
  
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
    text = element_text(size = 15),
    legend.position = "right"
  ) +
  
  scale_shape_manual(values = shape1) +  # Ensure shape1 is keyed by Treat levels
  scale_fill_manual(values = Adapted13d) +
  #scale_color_manual(values = Adapted13d) +
  
  scale_y_discrete(limits = rev) +
  scale_x_continuous(limits = c(-10, 10),
                     breaks = c(-10.0, -7.5, -5.0, -2.5, 0.0, 2.5, 5.0, 7.5, 10.0)) +
  
  geom_point(data = fig.data.edgR[fig.data.edgR$Shape1 == 22, ], aes(x = logFC, y = FAM_ASV,fill=Phylum),stroke=1, pch=21,size=3.50,color="black")+
  geom_point(data = fig.data.edgR[fig.data.edgR$Shape1 == 21, ], aes(x = logFC, y = FAM_ASV,fill=Phylum),stroke=1, pch=22,size=3.50,color="black")+
  
  labs(
    y = "Bacterial Taxonomy",
    x = "Log Fold Change",
    fill = "Phylum",
    #color = "Phylum",
    shape = "Treatment"
  )

lefcol<-plot_grid(bc1a, a0.plota, nrow=2,ncol=1,rel_widths=c(1,1),labels= c("A","B"), align="tl")
right_col<-plot_grid(edgE,nrow=1,ncol=1,rel_widths=c(1),labels= c("C"))
FINAL_FIG1<-plot_grid(lefcol,right_col,ncol=2,nrow=1,rel_widths=c(1.5,2.5),align="tl")


#test figure for color by treatment
treat.col<-c("#d81159","#b5e2fa")
fig.data.edgR$PhyFam<-factor(fig.data.edgR$PhyFam,levels=sort(unique(fig.data.edgR$PhyFam)))

fig.data.edgR %>%
  mutate(PhyFam = fct_inorder(PhyFam)) %>%
  ggplot(aes(x = logFC, y = PhyFam)) + 
  geom_point(aes(fill=Shape,shape=Shape,color=Shape), stroke=1, size=3.50, alpha=1) + 
  geom_vline(xintercept = c(0,0,0,0),linetype="dashed",colour="black",size=0.5)+
  theme(axis.text.x = element_text(angle = , hjust = 0, vjust = 0.5))+
  theme_bw() +
  scale_shape_manual(values=shape1)+
  scale_fill_manual(values=treat.col)+
  scale_color_manual(values=treat.col)+
  scale_y_discrete(limits=rev) +
  scale_x_continuous(limits=c(-10,10),breaks=c(-10.0,-7.5,-5.0,-2.5,0.0,2.5,5.0,7.5,10.0))+
  geom_point(data = fig.data.edgR[fig.data.edgR$Shape == 21, ], aes(x = logFC, y = PhyFam,fill=Shape),stroke=1, pch=21,size=3.50,color="black")+
  labs(y="Bacterial Taxonomy",x="Log Fold Change")+
  theme(legend.position = "right")+ labs(color="Phylum")+ theme(text = element_text(size = 14))
  