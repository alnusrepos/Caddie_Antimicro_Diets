library(dplyr)
library(ggplot2)
library(forcats)
metab_keggs<-read.csv(file = "MetabolismKeggsEdgeR.csv")

str(metab_keggs)
metab_keggs$pathway_mechanism<-as.factor(metab_keggs$pathway_mechanism)
metab_keggs$pathway_map<-as.factor(metab_keggs$pathway_map)
metab_keggs$Shape<-as.factor(metab_keggs$Shape)

shape1<-c(22,21)
shape1<-as.integer(shape1)
Phylum_7<-c("#ff0000","#ff7400","#ffc100","#1ED278","#007dff","#0200b9","#7b2cbf") #781156
Adapted15<-c("#4c00a4","#9c27b0","#ff00cc","#114578","#117878","#3FE4E4","#004242","#117845","#66ff66","#ffc100","#ff9005","#f9530b","#ff5577","#ff0000",	"#781156")
Adapted13d<-c("#781156","#ff00cc","#ff5400","#fcf300","#004b23","#38b000","#072ac8","#1e96fc","#3FE4E4","#531cb3","#c86bfa","#ffcbf2")

metab_keggs %>%
  mutate(pathway_map = fct_inorder(pathway_map)) %>%
  ggplot(aes(x = logFC, y = pathway_map)) + 
  geom_point(aes(fill = pathway_mechanism,  shape = Shape), 
             stroke = 1, size = 3.5, alpha = 1) +
  geom_vline(xintercept = c(0,0,0,0),linetype="dashed",colour="black",size=0.5)+
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
    text = element_text(size = 14),
    legend.position = "right"
  ) +
  
  scale_shape_manual(values = shape1) +  # Ensure shape1 is keyed by Treat levels
  scale_fill_manual(values = Adapted13d) +
  #scale_color_manual(values = Adapted13d) +
  
  scale_y_discrete(limits = rev) +
  scale_x_continuous(limits = c(-5, 1),
                     breaks = c(-5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 1.0)) +
  
  geom_point(data = metab_keggs[metab_keggs$Shape1 == 22, ], aes(x = logFC, y = pathway_map,fill=pathway_mechanism),stroke=1, pch=21,size=3.50,color="black")+
  geom_point(data = metab_keggs[metab_keggs$Shape1 == 21, ], aes(x = logFC, y = pathway_map,fill=pathway_mechanism),stroke=1, pch=22,size=3.50,color="black")+
  
  labs(
    y = "Predicted Metabolic Pathways",
    x = "Log Fold Change",
    fill = "Pathway Class",
    #color = "Phylum",
    shape = "Treatment"
  )


