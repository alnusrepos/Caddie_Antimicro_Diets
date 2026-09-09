load(file="CaddiesSeqAnalysis072825.rda")
colnames(caddie_alpha_meta)

#case mass
cm.plot <- ggplot(data = caddie_alpha_meta, aes(x = Treatment.Type, y = Case.Mass..g., fill = Treatment.Type)) +
  geom_boxplot(aes(fill = Treatment.Type), alpha = 0.65, outlier.size = 0) +
  geom_point(aes(shape = Treatment.Type), size = 3.5, position = position_jitterdodge(), color = "black", stroke = 0.5) +
  scale_fill_manual(values = c("#d81159","#b5e2fa")) +
  scale_shape_manual(values = c(22, 21)) +  # 22 = square filled, 21 = circle filled
  scale_y_continuous(limits=c(0,1),breaks=c(0.0,0.25,0.50,0.75,1.0))+
  theme_classic()+
  labs(y = "Case Mass (g)", x = "Gut Treatment") +
  theme(legend.position = "none")

cm.plota <- cm.plot + labs(fill = "Gut Treatment")+ 
  theme(
    axis.title.y = element_text(size=15),
    axis.title.x = element_text(size=15),
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12))


caddie_alpha_meta$Treatment.Type<-as.factor(caddie_alpha_meta$Treatment.Type)
caddie_treated<-caddie_alpha_meta[which(caddie_alpha_meta$Treatment.Type=="Antibiotic Treated"),]
caddie_control<-caddie_alpha_meta[which(caddie_alpha_meta$Treatment.Type=="Control"),]

hist((caddie_alpha_meta$Case.Mass..g.))
fBasics::normalTest(caddie_alpha_meta$Case.Mass..g., method = "da", na.rm = FALSE) #appears to me normally distributed
caddie_lm_casemass <- lm(Case.Mass..g. ~ Treatment.Type, dat = caddie_alpha_meta)
Anova(caddie_lm_casemass) #Treatment.Type 0.1254  1  3.3486 0.07579 .

t.test(caddie_control$Case.Mass..g.,caddie_treated$Case.Mass..g.,alternative = "two.sided",paired=FALSE) #it is a very similar result

#body mass
bm.plot <- ggplot(data = caddie_alpha_meta, aes(x = Treatment.Type, y = Body.Mass..g., fill = Treatment.Type)) +
  geom_boxplot(aes(fill = Treatment.Type), alpha = 0.65, outlier.size = 0) +
  geom_point(aes(shape = Treatment.Type), size = 3.5, position = position_jitterdodge(), color = "black", stroke = 0.5) +
  scale_fill_manual(values = c("#d81159","#b5e2fa")) +
  scale_shape_manual(values = c(22, 21)) +  # 22 = square filled, 21 = circle filled
  scale_y_continuous(limits=c(0,0.5),breaks=c(0,0.1,0.2,0.3,0.4,0.5))+
  theme_classic()+
  labs(y = "Body Mass (g)", x = "Gut Treatment") +
  theme(legend.position = "none")

bm.plota <- bm.plot + labs(fill = "Gut Treatment") +  labs(fill = "Gut Treatment")+ 
  theme(
    axis.title.y = element_text(size=15),
    axis.title.x = element_text(size=15),
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12)
  )

hist((caddie_alpha_meta$Body.Mass..g.))
fBasics::normalTest(caddie_alpha_meta$Body.Mass..g., method = "da", na.rm = FALSE) #fails all three tests, consider transformation

hist(log(caddie_alpha_meta$Body.Mass..g.))
fBasics::normalTest(log(caddie_alpha_meta$Body.Mass..g.), method = "da", na.rm = FALSE) #appears to be normally distributed

caddie_lm_bodymass<- lm(log(caddie_alpha_meta$Body.Mass..g.) ~ Treatment.Type, dat = caddie_alpha_meta)
Anova(caddie_lm_bodymass) #Treatment.Type 0.0063  1  0.0436 0.8359

t.test(log(caddie_control$Body.Mass..g.),log(caddie_treated$Body.Mass..g.),alternative = "two.sided",paired=FALSE) #it is a very similar result


#full body length
fbl.plot <- ggplot(data = caddie_alpha_meta, aes(x = Treatment.Type, y = Entire.Length..cm., fill = Treatment.Type)) +
  geom_boxplot(aes(fill = Treatment.Type), alpha = 0.65, outlier.size = 0) +
  geom_point(aes(shape = Treatment.Type), size = 3.5, position = position_jitterdodge(), color = "black", stroke = 0.5) +
  scale_fill_manual(values = c("#d81159","#b5e2fa")) +
  scale_shape_manual(values = c(22, 21)) +  # 22 = square filled, 21 = circle filled
  scale_y_continuous(limits=c(1.5,3.0),breaks=c(1.5,2.0,2.5,3.0))+
  theme_classic()+
  labs(y = "Full Body Length (cm)", x = "Gut Treatment") +
  theme(legend.position = "none")

fbl.plota<-fbl.plot + labs(fill = "Gut Treatment") +  labs(fill = "Gut Treatment")+ 
  theme(
    axis.title.y = element_text(size=15),
    axis.title.x = element_text(size=15),
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12),
  )


hist((caddie_alpha_meta$Entire.Length..cm.))
fBasics::normalTest(caddie_alpha_meta$Entire.Length..cm., method = "da", na.rm = FALSE) #normal!

caddie_lm_fbl<- lm(caddie_alpha_meta$Entire.Length..cm. ~ Treatment.Type, dat = caddie_alpha_meta)
Anova(caddie_lm_fbl) #Treatment.Type 0.01082  1   0.131 0.7196

t.test(caddie_control$Entire.Length..cm.,caddie_treated$Entire.Length..cm.,alternative = "two.sided",paired=FALSE) #it is a very similar result

#case width
cw.plot <- ggplot(data = caddie_alpha_meta, aes(x = Treatment.Type, y = Case.Width..cm., fill = Treatment.Type)) +
  geom_boxplot(aes(fill = Treatment.Type), alpha = 0.65, outlier.size = 0) +
  geom_point(aes(shape = Treatment.Type), size = 3.5, position = position_jitterdodge(), color = "black", stroke = 0.5) +
  scale_fill_manual(values = c("#d81159","#b5e2fa")) +
  scale_shape_manual(values = c(22, 21)) +  # 22 = square filled, 21 = circle filled
  scale_y_continuous(limits=c(0.4,1.4),breaks=c(0.4,0.6,0.8,1.0,1.2,1.4))+
  theme_classic()+
  labs(y = "Case Width (cm)", x = "Gut Treatment") +
  theme(legend.position = "none")

cw.plota<-cw.plot + labs(fill = "Gut Treatment") +  labs(fill = "Gut Treatment")+ 
  theme(
    axis.title.y = element_text(size=15),
    axis.title.x = element_text(size=15),
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12),
  )

hist((caddie_alpha_meta$Case.Width..cm.))
fBasics::normalTest(caddie_alpha_meta$Case.Width..cm., method = "da", na.rm = FALSE) #normal!

caddie_lm_cw<- lm(caddie_alpha_meta$Case.Width..cm. ~ Treatment.Type, dat = caddie_alpha_meta)
Anova(caddie_lm_cw) #Treatment.Type 0.02823  1  1.0573 0.3109

t.test(caddie_control$Case.Width..cm.,caddie_treated$Case.Width..cm.,alternative = "two.sided",paired=FALSE) #it is a very similar result

#Abdomen and Thorax length
at.plot <- ggplot(data = caddie_alpha_meta, aes(x = Treatment.Type, y = Abodomen.and.Thorax.Length..cm., fill = Treatment.Type)) +
  geom_boxplot(aes(fill = Treatment.Type), alpha = 0.65, outlier.size = 0) +
  geom_point(aes(shape = Treatment.Type), size = 3.5, position = position_jitterdodge(), color = "black", stroke = 0.5) +
  scale_fill_manual(values = c("#d81159","#b5e2fa")) +
  scale_shape_manual(values = c(22, 21)) +  # 22 = square filled, 21 = circle filled
  scale_y_continuous(limits=c(1.0,2.5),breaks=c(1.0,1.5,2.0,2.5))+
  theme_classic()+
  labs(y = "Abdomen Length (cm)", x = "Gut Treatment") +
  theme(legend.position = "none")

at.plota<-at.plot + labs(fill = "Gut Treatment") +  labs(fill = "Gut Treatment")+ 
  theme(
    axis.title.y = element_text(size=15),
    axis.title.x = element_text(size=15),
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12),
  )

hist((caddie_alpha_meta$Abodomen.and.Thorax.Length..cm.))
fBasics::normalTest(caddie_alpha_meta$Abodomen.and.Thorax.Length..cm., method = "da", na.rm = FALSE) #normal!

caddie_lm_atl<- lm(caddie_alpha_meta$Abodomen.and.Thorax.Length..cm. ~ Treatment.Type, dat = caddie_alpha_meta)
Anova(caddie_lm_atl) #Treatment.Type 0.00086  1  0.0115 0.9153

t.test(caddie_control$Abodomen.and.Thorax.Length..cm.,caddie_treated$Abodomen.and.Thorax.Length..cm.,alternative = "two.sided",paired=FALSE) #it is a very similar result

library(cowplot)
plot_grid(at.plota,
          fbl.plota,
          bm.plota,
          cw.plota,
          cm.plota,
          nrow=2,ncol=3,
          rel_widths=c(1,1,1,1,2,0),
          labels= c("A","B","C","D","E"))

#go back up to cm.plot and turn legend on to right

cm.plot.legend <- cm.plot + labs(fill = "Gut Treatment")+ 
  theme(
    axis.title.y = element_text(size=15),
    axis.title.x = element_text(size=15),
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12),
    legend.key.size = unit(1,"cm"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 15),
    legend.box.margin = margin(0,10,0,0))

