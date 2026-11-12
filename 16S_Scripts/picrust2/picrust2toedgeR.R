#load packages
library("KEGGREST"); packageVersion("KEGGREST") #1.42.0
library(readr)
library(ggpicrust2)
library(tibble)
library(tidyverse)
library(ggprism)
library(patchwork)

#picrust2 output: KO to KEGG abundance conversion. Move into EdgeR. Make figures. 
#load(file="ggpicrust2.rda")
setwd("~/Documents/UCSD_EBE_PostDoc/Mentee_Materials/Dahlia_Loomis/Caddisfly_HFA/analyses")
load(file="picrust2_edgeR.rda")

#useful global environment objects
str(metadata) #tbl 37 x 9
dim(kegg_abundance) #259 x 37
dim(ko_abundance) #7194
#dim(metacyc_abundance) #413 x 38

###############################################################
###############################################################
###############################################################
#                        Begin EdgeR                          #
###############################################################
###############################################################
###############################################################
library("phyloseq"); packageVersion("phyloseq") #1.50.0
library("edgeR"); packageVersion("edgeR") #4.4.0
library("ggplot2"); packageVersion("ggplot2") #3.5.1
library(forcats); packageVersion("forcats") #1.0.0
library(dplyr); packageVersion("dplyr") #1.1.4

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

#phyloseq object that contains KEGG Abundances (kegg_abundance),
kegg_abundance.t<-base::t(kegg_abundance) #37 x 259
dd<-sort(colSums(kegg_abundance.t),decreasing=TRUE)
#write.csv(dd,file="KEGG_PathwayAbundSorted.csv")
KEGG.phyloseq<-phyloseq(otu_table(kegg_abundance.t, taxa_are_rows=FALSE), sample_data(CaddieD_phylo.rm2@sam_data))

#fixing metadata
KEGG.phyloseq@sam_data$Treatment.Type
rownames(KEGG.phyloseq@otu_table)==KEGG.phyloseq@sam_data$Sample.Name. #everything matches

#transform sample counts and view variance per ASV
KEGG.phyloseq.p = transform_sample_counts(KEGG.phyloseq, function(x){x/sum(x)})
KEGG.phyloseq@otu_table[1:10,1:10]
rowSums(KEGG.phyloseq@otu_table)
KEGG.phyloseq.p@otu_table[1:10,1:10]

#view
hist(log10(apply(otu_table(KEGG.phyloseq.p), 1, var)),
     xlab="log10(variance)", breaks=50) #log10(x)=-4.4 x=0.0000398107

#apply KO filtering
AbundThreshold = 0.00001
f4Keeps= names(which(apply(otu_table(KEGG.phyloseq.p),2, sum) > AbundThreshold))
KEGG.phyloseq.f = prune_taxa(f4Keeps, KEGG.phyloseq)

#now lets use this newly defined function to transform phyloseq to edgeR's "DGE" type object
KEGG.dge = phyloseq_to_edgeR(KEGG.phyloseq.f, group="Treatment.Type") #works!

# Perform binary test
KEGG.et = exactTest(KEGG.dge)

# Extract values from test results
KEGG.tt = topTags(KEGG.et, n=nrow(KEGG.dge$table), adjust.method="BH", sort.by="PValue")
KEGG.res = KEGG.tt@.Data[[1]]
alpha_edger = 0.05 #0.001, 0.01, 0.05 
KEGG.sigtab = KEGG.res[(KEGG.res$FDR < alpha_edger), ] #165 x 4
dim(KEGG.sigtab) #165 x 4
head(KEGG.sigtab)
#write.csv(KEGG.sigtab,file="KEGG.sigtab.csv")
kegg.edgeR4annotation<-read.csv(file="KEGG.sigtab.csv")
dim(kegg.edgeR4annotation) #165 x 9

kegg.edgeR4annotation1<-kegg.edgeR4annotation[1:50,]
dim(kegg.edgeR4annotation1) #50x9
kegg.edgeR4annotation2<-kegg.edgeR4annotation[51:102,]
dim(kegg.edgeR4annotation2) #52x9
kegg.edgeR4annotation3<-kegg.edgeR4annotation[103:165,]
dim(kegg.edgeR4annotation3) #63x9

Annotated_kegg1 <- pathway_annotation(pathway = "KO", daa_results_df = kegg.edgeR4annotation1, ko_to_kegg = TRUE)
Annotated_kegg2 <- pathway_annotation(pathway = "KO", daa_results_df = kegg.edgeR4annotation2, ko_to_kegg = TRUE)
Annotated_kegg3 <- pathway_annotation(pathway = "KO", daa_results_df = kegg.edgeR4annotation3, ko_to_kegg = TRUE)
Final_KOtoKEGGwAnnos<-rbind(Annotated_kegg1,Annotated_kegg2,Annotated_kegg3)
#write.csv(Final_KOtoKEGGwAnnos,file="Final_KOtoKEGGwAnnos.csv")

#figure time! 
theme_set(theme_bw())
Final_KOtoKEGGwAnnos2<-read.csv(file="Final_KOtoKEGGwAnnos.csv")
str(Final_KOtoKEGGwAnnos2)
sigtabgen<-Final_KOtoKEGGwAnnos2
#sigtabgen = subset(Final_KOtoKEGGwAnnos, !is.na("pathway_mechanism")) #use if you want to get rid of NAs

# pathway_class order
x = tapply(sigtabgen$logFC, sigtabgen$pathway_class, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$pathway_class = factor(as.character(sigtabgen$pathway_class), levels = names(x))

# pathway_mechanism order
x = tapply(sigtabgen$logFC, sigtabgen$pathway_mechanism, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$pathway_mechanism = factor(as.character(sigtabgen$pathway_mechanism), levels = names(x))

# pathway_map order
x = tapply(sigtabgen$logFC, sigtabgen$pathway_map, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$pathway_map = factor(as.character(sigtabgen$pathway_map), levels = names(x))

ggplot(sigtabgen, aes(x = logFC, y = pathway_mechanism, color = pathway_class)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = , hjust = 0, vjust = 0.5))

Phylum_cols<-c("#781156","#ff0000","#ff4d00","#ff7400","#ff9a00","#ffc100","#117845","#1ED278","#6CEAAB","#117878","#3FE4E4","#98F0F0","#0200b9","#007dff","#2b89ea") #15

#write.csv(sigtabgen,file="sigtabgen.csv")
sigtabgen<-read.csv(file="sigtabgen.csv")
str(sigtabgen)
sigtabgen$pathway_class<-factor(sigtabgen$pathway_class,levels=sort(unique(sigtabgen$pathway_class)))
sigtabgen$pathway_mechanism<-factor(sigtabgen$pathway_mechanism,levels=sort(unique(sigtabgen$pathway_mechanism)))
sigtabgen$pathway_map<-factor(sigtabgen$pathway_map,levels=sort(unique(sigtabgen$pathway_map)))
sigtabgen$Combined<-factor(sigtabgen$Combined,levels=sort(unique(sigtabgen$Combined)))
sigtabgen$shape<-as.factor(sigtabgen$shape)

ggplot(sigtabgen, aes(x = logFC, y = pathway_mechanism, color = pathway_class)) + geom_point(size=4) + 
  geom_vline(xintercept = c(0,0,0,0),linetype="dashed",colour="black",size=0.5)+
  theme(axis.text.x = element_text(angle = , hjust = 0, vjust = 0.5))+
  theme_bw() +
  scale_color_manual(values=Phylum_cols) +
  labs(y="Predicted Metabolic Mechanism",x="Log Fold Change") +
  theme(legend.position = "right")+ labs(color="Pathway Class")+ theme(text = element_text(size = 15))

shape1<-c(22,21)
Phylum_7<-c("#ff0000","#ff7400","#ffc100","#1ED278","#007dff","#0200b9","#7b2cbf") #781156

#perfect figure with color by pathway class
sigtabgen %>%
  mutate(pathway_mechanism = fct_inorder(pathway_mechanism)) %>%
  ggplot(aes(x = logFC, y = pathway_mechanism)) + 
  geom_point(aes(fill=pathway_class,shape=shape,color=pathway_class), stroke=1, size=3.50, alpha=1) + 
  geom_vline(xintercept = c(0,0,0,0),linetype="dashed",colour="black",size=0.5)+
  theme(axis.text.x = element_text(angle = , hjust = 0, vjust = 0.5))+
  theme_bw() +
  scale_shape_manual(values=shape1)+
  scale_fill_manual(values=Phylum_7)+
  scale_color_manual(values=Phylum_7)+
  scale_y_discrete(limits=rev,expand = expansion(0.02)) +
  scale_x_continuous(limits=c(-5.0,3.75),breaks=c(-5.0,-3.75,-2.5,-1.25,0.0,1.25,2.5,3.75),expand = expansion(0.01))+
  geom_point(data = sigtabgen[sigtabgen$shape == 21, ], aes(x = logFC, y = pathway_mechanism, fill=pathway_class),stroke=1, pch=21,size=3.50,color="black")+
  labs(y="Predicted Bacterial Subclass",x="Log Fold Change") +
  theme(legend.position = "none")+ labs(color="Pathway Class")+ theme(text = element_text(size = 14))

sigtabgen %>%
  mutate(pathway_mechanism = fct_inorder(pathway_mechanism)) %>%
  ggplot(aes(x = logFC, y = pathway_mechanism)) + 
  geom_point(aes(fill = pathway_class,  shape = shape), 
             stroke = 1, size = 3.5, alpha = 1) +
  geom_vline(xintercept = c(0,0,0,0),linetype="dashed",colour="black",size=0.5)+
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
    text = element_text(size = 14),
    legend.position = "right"
  ) +
  
  scale_shape_manual(values = shape1) +  # Ensure shape1 is keyed by Treat levels
  scale_fill_manual(values = Phylum_7) +
  
  scale_y_discrete(limits = rev) +
  scale_x_continuous(limits=c(-5.0,3.75),breaks=c(-5.0,-3.75,-2.5,-1.25,0.0,1.25,2.5,3.75),expand = expansion(0.01))+
  
  geom_point(data = sigtabgen[sigtabgen$Shape1 == 22, ], aes(x = logFC, y = pathway_mechanism,fill=pathway_class),stroke=1, pch=21,size=3.50,color="black")+
  geom_point(data = sigtabgen[sigtabgen$Shape1 == 21, ], aes(x = logFC, y = pathway_mechanism,fill=pathway_class),stroke=1, pch=22,size=3.50,color="black")+
  
  labs(
    y = "Predicted Bacterial Subclass",
    x = "Log Fold Change",
    fill = "Pathway Class",
    #color = "Phylum",
    shape = "Treatment"
  )

#in combined data remove class
sigtabgen %>%
  mutate(Combined = fct_inorder(Combined)) %>%
  ggplot(aes(x = logFC, y = Combined)) + 
  geom_point(aes(fill=pathway_class,shape=shape,color=pathway_class), stroke=1, size=3.50, alpha=1) + 
  geom_vline(xintercept = c(0,0,0,0),linetype="dashed",colour="black",size=0.5)+
  theme(axis.text.x = element_text(angle = , hjust = 0, vjust = 0.5))+
  theme_bw() +
  scale_shape_manual(values=shape1)+
  scale_fill_manual(values=Phylum_7)+
  scale_color_manual(values=Phylum_7)+
  scale_y_discrete(limits=rev) +
  scale_x_continuous(limits=c(-5,3.75),breaks=c(-5.0,-3.75,-2.5,-1.25,0.0,1.25,2.5,3.75))+
  geom_point(data = sigtabgen[sigtabgen$shape == 21, ], aes(x = logFC, y = Combined, fill=pathway_class),stroke=1, pch=21,size=3.50,color="black")+
  labs(y="Predicted Metabolic Mechanism",x="Log Fold Change") +
  theme(legend.position = "right")+ labs(color="Pathway Class")+ theme(text = element_text(size = 14))

##phyloseq object that contains KO Abundances (KO),
#write.csv(ko_abundance,"ko_abundance.csv")
ko_abundance.t<-read.csv("ko_abundance_t.csv")
rownames(ko_abundance.t)<-ko_abundance.t$X
ko_abundance.t<-ko_abundance.t[1:37,2:7195]
str(ko_abundance.t)
ko_abundance.t2<-mutate_all(ko_abundance.t, function(x) as.numeric(as.character(x)))
ko.phyloseq<-phyloseq(otu_table(ko_abundance.t2, taxa_are_rows=FALSE), sample_data(CaddieD_phylo.rm2@sam_data))

#fixing metadata
ko.phyloseq@sam_data$Treatment.Type
rownames(ko.phyloseq@otu_table)==ko.phyloseq@sam_data$Sample.Name. #everything matches

#transform sample counts and view variance per ASV
ko.phyloseq.p = transform_sample_counts(ko.phyloseq, function(x){x/sum(x)})
ko.phyloseq@otu_table[1:10,1:10]
rowSums(ko.phyloseq@otu_table)
ko.phyloseq.p@otu_table[1:10,1:10]

#view
hist(log10(apply(otu_table(ko.phyloseq.p), 1, var)),
     xlab="log10(variance)", breaks=50) #log10(x)=-7.15 x=0.0000000707946

#apply KO filtering
AbundThreshold = 0.0000000707946
f4Keeps= names(which(apply(otu_table(ko.phyloseq.p),2, sum) > AbundThreshold))
ko.phyloseq.f = prune_taxa(f4Keeps, ko.phyloseq)

#now lets use this newly defined function to transform phyloseq to edgeR's "DGE" type object
ko.dge = phyloseq_to_edgeR(ko.phyloseq.f, group="Treatment.Type") #works!

# Perform binary test
ko.et = exactTest(ko.dge)

# Extract values from test results
ko.tt = topTags(ko.et, n=nrow(ko.dge$table), adjust.method="BH", sort.by="PValue")
ko.res = ko.tt@.Data[[1]]
alpha_edger = 0.05 #0.001, 0.01, 0.05 
ko.sigtab = ko.res[(ko.res$FDR < alpha_edger), ] #5303 x 4
dim(ko.sigtab) #5303 x 4
head(ko.sigtab)
#write.csv(ko.sigtab,file="ko.sigtab.csv")
ko.edgeR4annotation<-read.csv(file="ko.sigtab.csv")
dim(ko.edgeR4annotation) #5303 x 9

ko.edgeR4annotation1<-ko.edgeR4annotation[1:100,]
ko.edgeR4annotation2<-ko.edgeR4annotation[101:200,]
ko.edgeR4annotation3<-ko.edgeR4annotation[201:300,]
ko.edgeR4annotation4<-ko.edgeR4annotation[301:400,]
ko.edgeR4annotation5<-ko.edgeR4annotation[401:500,]
ko.edgeR4annotation6<-ko.edgeR4annotation[501:600,]
ko.edgeR4annotation7<-ko.edgeR4annotation[601:700,]
ko.edgeR4annotation8<-ko.edgeR4annotation[701:800,]
ko.edgeR4annotation9<-ko.edgeR4annotation[801:900,]
ko.edgeR4annotation10<-ko.edgeR4annotation[901:1000,]
ko.edgeR4annotation11<-ko.edgeR4annotation[1001:1100,]
ko.edgeR4annotation12<-ko.edgeR4annotation[1101:1200,]
ko.edgeR4annotation13<-ko.edgeR4annotation[1201:1300,]
ko.edgeR4annotation14<-ko.edgeR4annotation[1301:1400,]
ko.edgeR4annotation15<-ko.edgeR4annotation[1401:1500,]
ko.edgeR4annotation16<-ko.edgeR4annotation[1501:1600,]
ko.edgeR4annotation17<-ko.edgeR4annotation[1601:1700,]
ko.edgeR4annotation18<-ko.edgeR4annotation[1701:1800,]
ko.edgeR4annotation19<-ko.edgeR4annotation[1801:1900,]
ko.edgeR4annotation20<-ko.edgeR4annotation[1901:2000,]
ko.edgeR4annotation21<-ko.edgeR4annotation[2001:2100,]
ko.edgeR4annotation22<-ko.edgeR4annotation[2101:2200,]
ko.edgeR4annotation23<-ko.edgeR4annotation[2201:2300,]
ko.edgeR4annotation24<-ko.edgeR4annotation[2301:2400,]
ko.edgeR4annotation25<-ko.edgeR4annotation[2401:2500,]
ko.edgeR4annotation26<-ko.edgeR4annotation[2501:2600,]
ko.edgeR4annotation27<-ko.edgeR4annotation[2601:2700,]
ko.edgeR4annotation28<-ko.edgeR4annotation[2701:2800,]
ko.edgeR4annotation29<-ko.edgeR4annotation[2801:2900,]
ko.edgeR4annotation30<-ko.edgeR4annotation[2901:3000,]
ko.edgeR4annotation31<-ko.edgeR4annotation[3001:3100,]
ko.edgeR4annotation32<-ko.edgeR4annotation[3101:3200,]
ko.edgeR4annotation33<-ko.edgeR4annotation[3201:3300,]
ko.edgeR4annotation34<-ko.edgeR4annotation[3301:3400,]
ko.edgeR4annotation35<-ko.edgeR4annotation[3401:3500,]
ko.edgeR4annotation36<-ko.edgeR4annotation[3501:3600,]
ko.edgeR4annotation37<-ko.edgeR4annotation[3601:3700,]
ko.edgeR4annotation38<-ko.edgeR4annotation[3701:3800,]
ko.edgeR4annotation39<-ko.edgeR4annotation[3801:3900,]
ko.edgeR4annotation40<-ko.edgeR4annotation[3901:4000,]
ko.edgeR4annotation41<-ko.edgeR4annotation[4001:4100,]
ko.edgeR4annotation42<-ko.edgeR4annotation[4101:4200,]
ko.edgeR4annotation43<-ko.edgeR4annotation[4201:4300,]
ko.edgeR4annotation44<-ko.edgeR4annotation[4301:4400,]
ko.edgeR4annotation45<-ko.edgeR4annotation[4401:4500,]
ko.edgeR4annotation46<-ko.edgeR4annotation[4501:4600,]
ko.edgeR4annotation47<-ko.edgeR4annotation[4601:4700,]
ko.edgeR4annotation48<-ko.edgeR4annotation[4701:4800,]
ko.edgeR4annotation49<-ko.edgeR4annotation[4801:4900,]
ko.edgeR4annotation50<-ko.edgeR4annotation[4901:5000,]
ko.edgeR4annotation51<-ko.edgeR4annotation[5001:5100,]
ko.edgeR4annotation52<-ko.edgeR4annotation[5101:5200,]
ko.edgeR4annotation53<-ko.edgeR4annotation[5201:5300,]
ko.edgeR4annotation54<-ko.edgeR4annotation[5301:5303,]

ko.Annotated_kegg1 <- pathway_annotation(pathway = "KO", daa_results_df = ko.edgeR4annotation1, ko_to_kegg = FALSE)
ko.Annotated_kegg2 <- pathway_annotation(pathway = "KO", daa_results_df = ko.edgeR4annotation2, ko_to_kegg = FALSE)
ko.Annotated_kegg3 <- pathway_annotation(pathway = "KO", daa_results_df = ko.edgeR4annotation3, ko_to_kegg = FALSE)
ko.Annotated_kegg4 <- pathway_annotation(pathway = "KO", daa_results_df = ko.edgeR4annotation4, ko_to_kegg = FALSE)
ko.Annotated_kegg5 <- pathway_annotation(pathway = "KO", daa_results_df = ko.edgeR4annotation5, ko_to_kegg = FALSE)
ko.Annotated_kegg6 <- pathway_annotation(pathway = "KO", daa_results_df = ko.edgeR4annotation6, ko_to_kegg = FALSE)
ko.Annotated_kegg7 <- pathway_annotation(pathway = "KO", daa_results_df = ko.edgeR4annotation7, ko_to_kegg = FALSE)
ko.Annotated_kegg8 <- pathway_annotation(pathway = "KO", daa_results_df = ko.edgeR4annotation8, ko_to_kegg = FALSE)
ko.Annotated_kegg9 <- pathway_annotation(pathway = "KO", daa_results_df = ko.edgeR4annotation9, ko_to_kegg = FALSE)
ko.Annotated_kegg10 <- pathway_annotation(pathway = "KO", daa_results_df = ko.edgeR4annotation10, ko_to_kegg = FALSE)
ko.Annotated_kegg11 <- pathway_annotation(pathway = "KO", daa_results_df = ko.edgeR4annotation11, ko_to_kegg = FALSE)
ko.Annotated_kegg12 <- pathway_annotation(pathway = "KO", daa_results_df = ko.edgeR4annotation12, ko_to_kegg = FALSE)
ko.Annotated_kegg13 <- pathway_annotation(pathway = "KO", daa_results_df = ko.edgeR4annotation13, ko_to_kegg = FALSE)
ko.Annotated_kegg14 <- pathway_annotation(pathway = "KO", daa_results_df = ko.edgeR4annotation14, ko_to_kegg = FALSE)
ko.Annotated_kegg15 <- pathway_annotation(pathway = "KO", daa_results_df = ko.edgeR4annotation15, ko_to_kegg = FALSE)
ko.Annotated_kegg16 <- pathway_annotation(pathway = "KO", daa_results_df = ko.edgeR4annotation16, ko_to_kegg = FALSE)
ko.Annotated_kegg17 <- pathway_annotation(pathway = "KO", daa_results_df = ko.edgeR4annotation17, ko_to_kegg = FALSE)
ko.Annotated_kegg18 <- pathway_annotation(pathway = "KO", daa_results_df = ko.edgeR4annotation18, ko_to_kegg = FALSE)
ko.Annotated_kegg19 <- pathway_annotation(pathway = "KO", daa_results_df = ko.edgeR4annotation19, ko_to_kegg = FALSE)
ko.Annotated_kegg20 <- pathway_annotation(pathway = "KO", daa_results_df = ko.edgeR4annotation20, ko_to_kegg = FALSE)
ko.Annotated_kegg21 <- pathway_annotation(pathway = "KO", daa_results_df = ko.edgeR4annotation21, ko_to_kegg = FALSE)
ko.Annotated_kegg22 <- pathway_annotation(pathway = "KO", daa_results_df = ko.edgeR4annotation22, ko_to_kegg = FALSE)
ko.Annotated_kegg23 <- pathway_annotation(pathway = "KO", daa_results_df = ko.edgeR4annotation23, ko_to_kegg = FALSE)
ko.Annotated_kegg24 <- pathway_annotation(pathway = "KO", daa_results_df = ko.edgeR4annotation24, ko_to_kegg = FALSE)
ko.Annotated_kegg25 <- pathway_annotation(pathway = "KO", daa_results_df = ko.edgeR4annotation25, ko_to_kegg = FALSE)
ko.Annotated_kegg26 <- pathway_annotation(pathway = "KO", daa_results_df = ko.edgeR4annotation26, ko_to_kegg = FALSE)
ko.Annotated_kegg27 <- pathway_annotation(pathway = "KO", daa_results_df = ko.edgeR4annotation27, ko_to_kegg = FALSE)
ko.Annotated_kegg28 <- pathway_annotation(pathway = "KO", daa_results_df = ko.edgeR4annotation28, ko_to_kegg = FALSE)
ko.Annotated_kegg29 <- pathway_annotation(pathway = "KO", daa_results_df = ko.edgeR4annotation29, ko_to_kegg = FALSE)
ko.Annotated_kegg30 <- pathway_annotation(pathway = "KO", daa_results_df = ko.edgeR4annotation30, ko_to_kegg = FALSE)
ko.Annotated_kegg31 <- pathway_annotation(pathway = "KO", daa_results_df = ko.edgeR4annotation31, ko_to_kegg = FALSE)
ko.Annotated_kegg32 <- pathway_annotation(pathway = "KO", daa_results_df = ko.edgeR4annotation32, ko_to_kegg = FALSE)
ko.Annotated_kegg33 <- pathway_annotation(pathway = "KO", daa_results_df = ko.edgeR4annotation33, ko_to_kegg = FALSE)
ko.Annotated_kegg34 <- pathway_annotation(pathway = "KO", daa_results_df = ko.edgeR4annotation34, ko_to_kegg = FALSE)
ko.Annotated_kegg35 <- pathway_annotation(pathway = "KO", daa_results_df = ko.edgeR4annotation35, ko_to_kegg = FALSE)
ko.Annotated_kegg36 <- pathway_annotation(pathway = "KO", daa_results_df = ko.edgeR4annotation36, ko_to_kegg = FALSE)
ko.Annotated_kegg37 <- pathway_annotation(pathway = "KO", daa_results_df = ko.edgeR4annotation37, ko_to_kegg = FALSE)
ko.Annotated_kegg38 <- pathway_annotation(pathway = "KO", daa_results_df = ko.edgeR4annotation38, ko_to_kegg = FALSE)
ko.Annotated_kegg39 <- pathway_annotation(pathway = "KO", daa_results_df = ko.edgeR4annotation39, ko_to_kegg = FALSE)
ko.Annotated_kegg40 <- pathway_annotation(pathway = "KO", daa_results_df = ko.edgeR4annotation40, ko_to_kegg = FALSE)
ko.Annotated_kegg41 <- pathway_annotation(pathway = "KO", daa_results_df = ko.edgeR4annotation41, ko_to_kegg = FALSE)
ko.Annotated_kegg42 <- pathway_annotation(pathway = "KO", daa_results_df = ko.edgeR4annotation42, ko_to_kegg = FALSE)
ko.Annotated_kegg43 <- pathway_annotation(pathway = "KO", daa_results_df = ko.edgeR4annotation43, ko_to_kegg = FALSE)
ko.Annotated_kegg44 <- pathway_annotation(pathway = "KO", daa_results_df = ko.edgeR4annotation44, ko_to_kegg = FALSE)
ko.Annotated_kegg45 <- pathway_annotation(pathway = "KO", daa_results_df = ko.edgeR4annotation45, ko_to_kegg = FALSE)
ko.Annotated_kegg46 <- pathway_annotation(pathway = "KO", daa_results_df = ko.edgeR4annotation46, ko_to_kegg = FALSE)
ko.Annotated_kegg47 <- pathway_annotation(pathway = "KO", daa_results_df = ko.edgeR4annotation47, ko_to_kegg = FALSE)
ko.Annotated_kegg48 <- pathway_annotation(pathway = "KO", daa_results_df = ko.edgeR4annotation48, ko_to_kegg = FALSE)
ko.Annotated_kegg49 <- pathway_annotation(pathway = "KO", daa_results_df = ko.edgeR4annotation49, ko_to_kegg = FALSE)
ko.Annotated_kegg50 <- pathway_annotation(pathway = "KO", daa_results_df = ko.edgeR4annotation50, ko_to_kegg = FALSE)
ko.Annotated_kegg51 <- pathway_annotation(pathway = "KO", daa_results_df = ko.edgeR4annotation51, ko_to_kegg = FALSE)
ko.Annotated_kegg52 <- pathway_annotation(pathway = "KO", daa_results_df = ko.edgeR4annotation52, ko_to_kegg = FALSE)
ko.Annotated_kegg53 <- pathway_annotation(pathway = "KO", daa_results_df = ko.edgeR4annotation53, ko_to_kegg = FALSE)
ko.Annotated_kegg54 <- pathway_annotation(pathway = "KO", daa_results_df = ko.edgeR4annotation54, ko_to_kegg = FALSE)

Final_KOAnnos<-rbind(ko.Annotated_kegg1,ko.Annotated_kegg2,ko.Annotated_kegg3,ko.Annotated_kegg4,ko.Annotated_kegg5,ko.Annotated_kegg6,ko.Annotated_kegg7,ko.Annotated_kegg8,ko.Annotated_kegg9,ko.Annotated_kegg10,ko.Annotated_kegg11,ko.Annotated_kegg12,ko.Annotated_kegg13,ko.Annotated_kegg14,ko.Annotated_kegg15,ko.Annotated_kegg16,ko.Annotated_kegg17,ko.Annotated_kegg18,ko.Annotated_kegg19,ko.Annotated_kegg20,ko.Annotated_kegg21,ko.Annotated_kegg22,ko.Annotated_kegg23,ko.Annotated_kegg24,ko.Annotated_kegg25,ko.Annotated_kegg26,ko.Annotated_kegg27,ko.Annotated_kegg28,ko.Annotated_kegg29,ko.Annotated_kegg30,ko.Annotated_kegg31,ko.Annotated_kegg32,ko.Annotated_kegg33,ko.Annotated_kegg34,ko.Annotated_kegg35,ko.Annotated_kegg36,ko.Annotated_kegg37,ko.Annotated_kegg38,ko.Annotated_kegg39,ko.Annotated_kegg40,ko.Annotated_kegg41,ko.Annotated_kegg42,ko.Annotated_kegg43,ko.Annotated_kegg44,ko.Annotated_kegg45,ko.Annotated_kegg46,ko.Annotated_kegg47,ko.Annotated_kegg48,ko.Annotated_kegg49,ko.Annotated_kegg50,ko.Annotated_kegg51,ko.Annotated_kegg52,ko.Annotated_kegg53,ko.Annotated_kegg54)
#write.csv(Final_KOAnnos,file="Final_KOAnnos.csv")
#save(list=ls(),file="picrust2_edgeR.rda")



# I need to filter and select only the top 10% most significant KO pathways to annotate. 
