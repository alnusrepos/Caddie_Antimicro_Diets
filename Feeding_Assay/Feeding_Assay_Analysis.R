#July 2026
library(ggplot2)
library(emmeans)

#setwd
setwd("/Users/jrdickey/Documents/UCSD_EBE_PostDoc/Mentee_Materials/Dahlia_Loomis/Caddisfly_HFA/Diet_analysis")

#read in data
newdat<-read.csv(file="GOAL_4_PLATE_WEIGHT_Nov25_2024.csv")
str(newdat)
newdat$Mass.Consumed_1<-as.numeric(newdat$Mass.Consumed_1)
newdat1<-na.omit(newdat) #removing NAs, 838 x 12
table(newdat1$Treatment) #420 x 418

#arc sin transform
mass.cons.t<-sqrt(newdat1$Mass.Consumed_1)
mass.cons.t.arcsin<-asin(mass.cons.t)
newdat2<-cbind(newdat1,mass.cons.t.arcsin)

#convert to factor
newdat2$SAMPLE<-as.factor(newdat2$SAMPLE)

#make composite variable with Treatment and Home.Away
metab.treat<-newdat2$Home.Away

metab.treat<-replace(metab.treat,metab.treat=="H","Local")
metab.treat<-replace(metab.treat,metab.treat=="A","Non-local")
metab.treat<-replace(metab.treat,metab.treat=="C","Control")

Metab2.treat<-newdat2$Home.Away

Metab2.treat<-replace(Metab2.treat,Metab2.treat=="H","A. rubra") #Just A. rubra inbued now
Metab2.treat<-replace(Metab2.treat,Metab2.treat=="A","A. rubra") #non-location specific
Metab2.treat<-replace(Metab2.treat,Metab2.treat=="C","Control")

newdat3<-cbind(newdat2,metab.treat,Metab2.treat)

newdat3$metab.treat<-factor(newdat3$metab.treat,c("Control","Local","Non-local"))
newdat3$Metab2.treat<-factor(newdat3$Metab2.treat,c("A. rubra","Control"))

newdat3$Treatment[newdat3$Treatment=="Untreated"] <- "Control"
newdat3$Treatment<-factor(newdat3$Treatment,c("Antibiotic Treated","Control"))

library(dplyr)
newdat3 <- newdat3 %>%
  mutate(Comp_Treat = paste(Treatment, metab.treat, sep = "_"))

newdat3$Comp_Treat
table(newdat3$Comp_Treat)
str(newdat3)
newdat3$Comp_Treat<-as.factor(newdat3$Comp_Treat)

#models
library(lmerTest)
null<-lmer(mass.cons.t.arcsin~(1|SAMPLE),data=newdat3)
full<-lmer(mass.cons.t.arcsin~Comp_Treat+(1|SAMPLE),data=newdat3)

anova(null,full) #similarly, Treatment is highly significant with lmer
#null: mass.cons.t.arcsin ~ (1 | SAMPLE)
#full: mass.cons.t.arcsin ~ Comp_Treat + (1 | SAMPLE)
#npar     AIC     BIC logLik -2*log(L)  Chisq Df Pr(>Chisq)   
#null    3 -1838.2 -1824.0 922.12   -1844.2                        
#full    8 -1845.0 -1807.2 930.52   -1861.0 16.815  5   0.004864 **

anova(full)
#            Sum Sq  Mean Sq NumDF  DenDF F value  Pr(>F)   
#Comp_Treat  0.12929 0.025858     5 27.441  4.1095 0.00651 **

levels(newdat3$Comp_Treat)
library(emmeans)
treatment_means<-emmeans(full, ~Comp_Treat)
contrast_weights<-list("Away"=c(-1,1,0,0,0,0),"Control"=c(0,0,-1,1,0,0),"Home"=c(0,0,0,0,-1,1))
results<-contrast(treatment_means,method=contrast_weights)

summary(results,adjust="none")
# contrast estimate     SE   df t.ratio p.value
#Away      -0.0543 0.0178 13.2  -3.042  0.0093
#Control    0.0728 0.0178 13.2   4.079  0.0013
#Home       0.0130 0.0101 41.4   1.287  0.2053


#Fig 5:
null<-lmer(mass.cons.t.arcsin~(1|SAMPLE),data=newdat3)
full<-lmer(mass.cons.t.arcsin~Comp_Treat+(1|SAMPLE),data=newdat3)
anova(null,full)#similarly, Treatment is highly significant with lmer and I think this is what we didn't make clear in the original manuscript
anova(full)
#Comp_Treat 0.12929 0.025858     5 27.441  4.1095 0.00651 **

levels(newdat3$Comp_Treat)
treatment_means<-emmeans(full,~Comp_Treat)
contrast_weights<-list("antibiotic caddies"=c(1,-0.5,-0.5,0,0,0),"untreated caddies"=c(0,0,0,1,-0.5,-0.5))
cw_rubraonly<-list("irrespective_gutcond"=c(0, 1, -1, 0, 1, -1)) #contrast irrespective of gut microbiome, local vs non local
results<-contrast(treatment_means,method=contrast_weights)
summary(results,adjust="none") 
# antibiotic caddies   0.0561 0.0171 12.3   3.279  0.006
# untreated caddies    0.0626 0.0171 12.3   3.660  0.0031

results<-contrast(treatment_means,method=cw_rubraonly)
summary(results,adjust="none")
# irrespective_gutcond -0.00938 0.0164 18  -0.571  0.5754

#Fig. 6
null<-lmer(mass.cons.t.arcsin~(1|SAMPLE),data=newdat3)
full<-lmer(mass.cons.t.arcsin~Comp_Treat+(1|SAMPLE),data=newdat3)
anova(null,full)#similarly, Treatment is highly significant with lmer
anova(full)
#Comp_Treat 0.12929 0.025858     5 27.441  4.1095 0.00651 **

levels(newdat3$Comp_Treat)
treatment_means<-emmeans(full,~Comp_Treat)
contrast_weights<-list("Non-local"=c(-1,1,0,0,0,0),"Control"=c(0,0,-1,1,0,0),"Local"=c(0,0,0,0,-1,1))
contrast_weights<-list("Non-local"=c(0,0,-1,0,0,1),"Control"=c(-1,0,0,1,0,0),"Local"=c(0,-1,0,0,1,0))
results<-contrast(treatment_means,method=contrast_weights)
summary(results,adjust="none")
# contrast  estimate     SE   df t.ratio p.value
# Non-local  1.67e-02 0.00839 814   1.987  0.0472
# Control    1.48e-02 0.01450 814   1.024  0.3060
# Local     -5.15e-05 0.00836 814  -0.006  0.9951


#Plotting
newdat3_means <- newdat3 %>%
  group_by(Metab2.treat) %>% #a.rubra vs control
  summarise(mean_y = mean(mass.cons.t.arcsin, na.rm = TRUE)) #mean

newdat3_meansb <- newdat3 %>%
  group_by(Metab2.treat, Treatment) %>%
  summarise(
    mean_y = mean(mass.cons.t.arcsin, na.rm = TRUE),
    .groups = "drop"
  )

#  Metab2.treat mean_y
# A. rubra      0.283
# Control       0.342

grpmean<-mean(newdat3$mass.cons.t.arcsin) #0.29

library(ggplot2)
ggplot(newdat3, aes(x = Metab2.treat, y = mass.cons.t.arcsin)) +
  
  geom_boxplot(
    aes(fill = Metab2.treat),
    color = "black", alpha = 0.4, linewidth = 0.75,
    outliers = FALSE, fatten = NULL) +
  
  geom_jitter(
    aes(fill = Metab2.treat, color = Metab2.treat, shape = Treatment),
    width = 0.25, height = 0, size = 2.5, alpha = 1, stroke = 1) +
  
  geom_hline(yintercept = grpmean, linetype = "dashed",
             color = "grey60", linewidth = 0.7) +
  geom_point(
    data = newdat3_means,
    aes(x = Metab2.treat, y = mean_y),
    shape = 23, size = 4, fill = "white", color = "black", stroke = 1) +
    scale_fill_manual(values = c("A. rubra" = "#538d22", "Control" = "#7f4f24")) +
    scale_color_manual(values = c("A. rubra" = "#538d22", "Control" = "#7f4f24")) +
    scale_shape_manual(values = c("Control" = 21, "Antibiotic Treated" = 1),
                     labels = c("Control" = "Control", "Antibiotic Treated" = "Antibiotic Treated")) +
    scale_x_discrete(labels = c("A. rubra" = expression(italic("A. rubra")),
                              "Control"  = "Control")) +
    guides(fill  = "none", color = "none", shape = guide_legend( title = "Gut",
                                                               override.aes = list(
                                                                 fill  = c(NA, "black"),
                                                                 color = c("black", "black"),
                                                                 size  = 4))) +
  scale_y_continuous(limits=c(0,0.6),breaks=c(0.0,0.2,0.4,0.6))+
  
  theme_classic() +
  
  labs(y = "Mass of Diet Consumed (g)\n (arcsine transformed)",
       x = "Diet Type") +
  theme(legend.position = "right",
        axis.text.x  = element_text(size = 12),
        text         = element_text(size = 15),
        panel.grid.major.y = element_line(color = "grey80", linewidth = 0.5),
        panel.grid.minor.y = element_line(color = "grey90", linewidth = 0.3))

#treatment cols
treatment_colors <- c("Antibiotic Treated" = "#d81159", "Control" = "#b5e2fa")

#diet type
diet_type_fill  <- c(Control = "#7f4f24", Local = "black",   "Non-local" = "#C4CBCA")
diet_type_color <- c(Control = "black",   Local = "#538d22", "Non-local" = "#538d22")

#rename some things
# newdat3$metab.treat <- as.character(newdat3$metab.treat)
# newdat3$metab.treat[newdat3$metab.treat == "Home"] <- "Local"
# newdat3$metab.treat[newdat3$metab.treat == "Away"] <- "Non-local"
# newdat3$Treatment <- as.character(newdat3$Treatment)
# 
# plot_dat <- newdat3 %>% mutate(
#   metab.treat = factor(metab.treat, levels = c("Control", "Local", "Non-local")),
#   Treatment = factor(Treatment,   levels = c("Antibiotic Treated", "Control")))

plot_dat <- newdat3 

#Per-Treatment x diet means -> white diamond black outline
group_means <- plot_dat %>%
  group_by(metab.treat, Treatment) %>%
  summarise(mean_y = mean(mass.cons.t.arcsin), .groups = "drop")

#Overall per-diet mean (collapsed across Treatment) -> dashed grey line from Sara's fig
diet_means <- plot_dat %>%
  group_by(metab.treat) %>%
  summarise(mean_y = mean(mass.cons.t.arcsin), .groups = "drop")

panel_titles <- c(Control = "Control Diets", Local = "Local A. rubra Diets", `Non-local` = "Non-local A. rubra Diets")

panel_titles <- c(Control = "\"Control Diets\"", Local = "\"Local \"*italic(\"A. rubra\")*\" Diets\"", `Non-local` = "\"Non-local \"*italic(\"A. rubra\")*\" Diets\"")

##Edit these based on planned-contrast p-values differ
sig_labels <- c(Control = "N.S.", Local = "N.S.", `Non-local` = "p < 0.05")

y_max <- 0.6   #Edit if need be, but max is 0.6

make_violin_panel <- function(diet_level) {
  
  dat_sub <- plot_dat %>%
    filter(metab.treat == diet_level) %>%
    mutate(panel_label = panel_titles[[diet_level]])
  
  means_sub <- group_means %>%
    filter(metab.treat == diet_level) %>%
    mutate(panel_label = panel_titles[[diet_level]])
  
  hline_y   <- diet_means$mean_y[diet_means$metab.treat == diet_level]
  bracket_y <- y_max * 0.93
  
  ggplot(dat_sub, aes(x = Treatment, y = mass.cons.t.arcsin)) +
    facet_wrap(~ panel_label, labeller = label_parsed) + 
    
    geom_boxplot(aes(fill = Treatment), color = "black",
                alpha = 0.4, linewidth = 0.75, outliers = FALSE, fatten = NULL) +
    
    geom_jitter(aes(color = Treatment), width = 0.25, height = 0,
                size = 2.5, alpha = 1) +
    
    geom_hline(yintercept = hline_y, linetype = "dashed",
               color = "grey60", linewidth = 0.7) +
    
    geom_point(data = means_sub, aes(x = Treatment, y = mean_y),
               shape = 23, size = 4, fill = "white", color = "black", stroke = 1) +
    
    annotate("segment", x = 1, xend = 2, y = bracket_y, yend = bracket_y,
             linewidth = 0.5) +
    
    annotate("text", x = 1.5, y = bracket_y * 1.04,
             label = sig_labels[[diet_level]], size = 5) +
    
    scale_fill_manual(values = treatment_colors, guide = "none") +
    scale_color_manual(values = treatment_colors, guide = "none") +
    
    scale_y_continuous(limits = c(0, 0.6),
                       breaks = c(0, 0.2, 0.4, 0.6)) +
    scale_x_discrete(labels = c("Antibiotic Treated" = "Antibiotic-treated",
                                "Control" = "Control")) +
    
    labs(x = expression(italic("Dicosmoecus") ~ "Larvae"),
         y = "Mass of Diet Consumed (g)\n(arcsine transformed)") +
    
    theme_classic(base_size = 13) +
    theme(
      strip.background    = element_rect(fill = "grey85", color = "black"),
      strip.text = element_text(size = 12),
      panel.grid.major.y = element_line(color = "grey80", linewidth = 0.5),
      panel.grid.minor.y = element_line(color = "grey90", linewidth = 0.3),
      axis.text.x  = element_text(size = 12),
      axis.text.y = element_text(size = 12)
    )
} # weird way to do this but whatever

panel_A <- make_violin_panel("Control")
panel_B <- make_violin_panel("Local")
panel_C <- make_violin_panel("Non-local")

panel_A
panel_B
panel_C

library(cowplot)
library(tidyverse)
plot_grid(panel_A,panel_B,panel_C, nrow=1, ncol =3)

#D
str(plot_dat) #838 x 17

sample_summary <- plot_dat %>%
  group_by(SAMPLE, metab.treat, Treatment) %>%
  summarise(
    mean_y = mean(mass.cons.t.arcsin),
    se_y   = sd(mass.cons.t.arcsin) / sqrt(n()),
    .groups = "drop"
  )

sample_wide <- sample_summary %>%
  pivot_wider(
    id_cols     = c(SAMPLE, metab.treat),
    names_from  = Treatment,
    values_from = c(mean_y, se_y)
  ) %>%
  
  rename(
    mean_control    = `mean_y_Control`,
    se_control      = `se_y_Control`,
    mean_antibiotic = `mean_y_Antibiotic Treated`,
    se_antibiotic   = `se_y_Antibiotic Treated`
  ) %>%
  
  mutate(
    DietType = recode(as.character(metab.treat),
                      "Control" = "Control", "Local" = "Local", "Non-local" = "Non-local"),
    DietType = factor(DietType, levels = c("Control", "Local", "Non-local"))
  )

panel_D <- ggplot(sample_wide, aes(x = mean_control, y = mean_antibiotic)) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.9) +
  geom_errorbar(aes(ymin = mean_antibiotic - se_antibiotic,
                    ymax = mean_antibiotic + se_antibiotic),
                width = 0, color = "grey60") +
  geom_errorbarh(aes(xmin = mean_control - se_control,
                     xmax = mean_control + se_control),
                 height = 0, color = "grey60") +
  geom_point(aes(fill = DietType, color = DietType),
             shape = 21, size = 3.5, stroke = 1.2) +
  scale_x_continuous(limits = c(0.23, 0.37),
                     breaks = c(0.23,0.25,0.27,0.29,0.31,0.33,0.35,0.37)) +
  scale_y_continuous(limits = c(0.23, 0.37),
                     breaks = c(0.23,0.25,0.27,0.29,0.31,0.33,0.35,0.37)) +
  scale_fill_manual(name = "Diet Type", values = diet_type_fill) +
  scale_color_manual(name = "Diet Type", values = diet_type_color) +
  labs(x = "Mass of Diet Consumed (g) by\nControl Larvae",
       y = "Mass of Diet Consumed (g) by\nAntibiotic-treated Larvae") +
  theme_classic(base_size = 13) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))

panel_D
plot_grid(panel_A,panel_B,panel_C, panel_D, nrow=1, ncol =4, labels= LETTERS)
plot_grid(panel_B,panel_C,panel_A, panel_D, nrow=1, ncol =4, labels= LETTERS)
plot_grid(panel_B,panel_C,panel_A, panel_D, nrow=2, ncol =2, labels= LETTERS)
plot_grid(panel_B,panel_C, panel_D, nrow=1, ncol =3, labels= LETTERS)

# --- Prep: Away trees only, ordered LH-1 to LH-10 ---
away_dat <- newdat3 %>%
  filter(metab.treat == "Non-local") %>%
  mutate(
    Treatment = factor(Treatment, levels = c("Antibiotic Treated", "Control")),
    SAMPLE    = factor(SAMPLE, levels = paste0("LH-", 1:10))
  )

home_dat <- newdat3 %>%
  filter(metab.treat == "Local") %>%
  mutate(
    Treatment = factor(Treatment, levels = c("Antibiotic Treated", "Control")),
    SAMPLE    = factor(SAMPLE, levels = paste0("ELK-", 1:10))
  )

# Per-Treatment x SAMPLE means (diamonds)
away_means <- away_dat %>%
  group_by(SAMPLE, Treatment) %>%
  summarise(mean_y = mean(mass.cons.t.arcsin), .groups = "drop")

home_means <- home_dat %>%
  group_by(SAMPLE, Treatment) %>%
  summarise(mean_y = mean(mass.cons.t.arcsin), .groups = "drop")

# Overall mean across all Away wells (dashed line)
overall_mean <- mean(away_dat$mass.cons.t.arcsin)
overall_mean_h <- mean(home_dat$mass.cons.t.arcsin)

# --- Colors (matching main figure) ---
treatment_colors <- c("Antibiotic Treated" = "#d81159", "Control" = "#b5e2fa")

# --- Plot ---
si_violin <- ggplot(away_dat,
                    aes(x = Treatment, y = mass.cons.t.arcsin)) +
  
  #geom_boxplot(aes(fill = Treatment), color = "black",
  #            alpha = 0.4, linewidth = 0.75, outliers = FALSE, fatten = NULL) +
  
  geom_jitter(aes(fill = Treatment, color = Treatment),
              width = 0.25, height = 0,
              size = 2.5, alpha = 1, shape = 21) +
  
  geom_hline(yintercept = overall_mean, linetype = "dashed",
             color = "grey50", linewidth = 0.6) +
  
  geom_point(data = away_means,
             aes(x = Treatment, y = mean_y),
             shape = 23, size = 3.5,
             fill = "white", color = "black", stroke = 1) +
  
  facet_wrap(~ SAMPLE, nrow = 1) +
  
  scale_fill_manual(values  = treatment_colors,
                    name    = "Gut",
                    labels  = c("Antibiotic Treated", "Control")) +
  scale_color_manual(values = treatment_colors,
                     name   = "Gut",
                     labels = c("Antibiotic Treated", "Control")) +
  scale_y_continuous(limits = c(0.0, 0.6),
                     breaks = c(0.0,0.1, 0.2,0.3, 0.4, 0.5, 0.6)) +
  scale_x_discrete(labels = NULL) +   # suppress x-axis text; Treatment shown by color/legend
  
  labs(y = "Mass of Diet Consumed (g)\n(arcsine transformed)",
       x = NULL) +
  
  theme_classic(base_size = 12) +
  theme(
    strip.background    = element_rect(fill = "grey85", color = "black"),
    strip.text          = element_text(size = 12),
    axis.ticks.x        = element_blank(),
    panel.grid.major.y = element_line(color = "grey80", linewidth = 0.5),
    #panel.grid.minor.y = element_line(color = "grey90", linewidth = 0.3),
    legend.position     = "right",
    legend.title        = element_text(size = 12),
    legend.text         = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  )

si_violin

si_violin_h <- ggplot(home_dat,
                      aes(x = Treatment, y = mass.cons.t.arcsin)) +

  #geom_boxplot(aes(fill = Treatment), color = "black",
  #            alpha = 0.4, linewidth = 0.75, outliers = FALSE, fatten = NULL) +
  
  geom_jitter(aes(fill = Treatment, color = Treatment),
              width = 0.25, height = 0,
              size = 2.5, alpha = 1, shape = 21) +
  
  geom_hline(yintercept = overall_mean_h, linetype = "dashed",
             color = "grey50", linewidth = 0.6) +
  
  geom_point(data = home_means,
             aes(x = Treatment, y = mean_y),
             shape = 23, size = 3.5,
             fill = "white", color = "black", stroke = 1) +
  
  facet_wrap(~ SAMPLE, nrow = 1) +
  
  scale_fill_manual(values  = treatment_colors,
                    name    = "Gut",
                    labels  = c("Antibiotic Treated", "Control")) +
  scale_color_manual(values = treatment_colors,
                     name   = "Gut",
                     labels = c("Antibiotic Treated", "Control")) +
  scale_y_continuous(limits = c(0.0, 0.6),
                     breaks = c(0.0,0.1, 0.2,0.3, 0.4, 0.5, 0.6)) +
  scale_x_discrete(labels = NULL) +   # suppress x-axis text; Treatment shown by color/legend
  
  labs(y = "Mass of Diet Consumed (g)\n(arcsine transformed)",
       x = NULL) +
  
  theme_classic(base_size = 12) +
  theme(
    strip.background    = element_rect(fill = "grey85", color = "black"),
    strip.text          = element_text(size = 12),
    axis.ticks.x        = element_blank(),
    panel.grid.major.y  = element_line(color = "grey80", linewidth = 0.3),
    legend.position     = "right",
    legend.title        = element_text(size = 12),
    legend.text         = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  )

si_violin_h


###
ar_wide <- sample_wide %>%
  filter(metab.treat %in% c("Local", "Non-local"))   # n = 20 trees

# Spearman rank correlation
cor.test(ar_wide$mean_control,
         ar_wide$mean_antibiotic,
         method = "spearman",
         exact  = FALSE) 

newdat3 <- newdat3 %>% #raw mass consumed in mg and percent consumed
  mutate(
    mass_mg      = Mass.Consumed_1 * 1000,        # grams -> mg
    pct_consumed = pmax(0, as.numeric(X..Consumed))
  )

ff<-newdat3 %>%
  group_by(SAMPLE) %>%
  summarise(
    mean_arcsin = mean(mass.cons.t.arcsin),
    mean_mg     = mean(mass_mg),
    count_l = length(mass_mg),
    std = sd(mass_mg),
    se = std/sqrt(count_l),
    .groups = "drop"
  )

as.data.frame(ff)

ss<-newdat3 %>% #gray dashed lines for panel abc
  group_by(metab.treat) %>%
  summarise(
    mean_arcsin = mean(mass.cons.t.arcsin),
    mean_mg     = mean(mass_mg),
    count_l = length(mass_mg),
    std = sd(mass.cons.t.arcsin),
    se = std/sqrt(count_l),
    .groups = "drop"
  )

as.data.frame(ss)

ff1<-newdat3 %>% #values of white daimonds
  group_by(metab.treat, Treatment) %>%
  summarise(
    mean_arcsin = mean(mass.cons.t.arcsin),
    mean_pct    = mean(pct_consumed),
    mean_mg     = mean(mass_mg),
    .groups = "drop"
  ) %>%
  arrange(metab.treat, Treatment)

as.data.frame(ff1)

newdat3 %>%
  group_by(Metab2.treat) %>%
  summarise(
    mean_arcsin = mean(mass.cons.t.arcsin),
    mean_pct    = mean(pct_consumed),
    mean_mg     = mean(mass_mg),
    .groups = "drop"
  ) %>%
  arrange(Metab2.treat)

away_summ <- newdat3 %>%
  filter(metab.treat == "Non-local") %>%
  group_by(Treatment) %>%
  summarise(mean_pct = mean(pct_consumed), mean_mg = mean(mass_mg), .groups = "drop")

print(away_summ)

away_diff_pct <- away_summ %>%
  pivot_wider(names_from = Treatment, values_from = c(mean_pct, mean_mg)) %>%
  mutate(
    diff_pct_consumed = `mean_pct_Antibiotic Treated` - `mean_pct_Control`,
    pct_less          = ((`mean_pct_Control` - `mean_pct_Antibiotic Treated`) /
                           `mean_pct_Control`) * 100
  )

print(away_diff_pct)

diet_summ <- newdat3 %>%
  group_by(Metab2.treat) %>%
  summarise(mean_pct = mean(mass.cons.t.arcsin),
            std = sd(mass.cons.t.arcsin),
            se = std/sqrt(length(mass.cons.t.arcsin)),
            mean(pct_consumed), 
            mean_mg = mean(mass_mg), .groups = "drop")

print(diet_summ)

# % less consumed on A. rubra vs control diet
arubra_mean  <- mean(newdat3$pct_consumed[newdat3$metab.treat != "Control"])
control_mean <- mean(newdat3$pct_consumed[newdat3$metab.treat == "Control"])
pct_less_arubra <- ((control_mean - arubra_mean) / control_mean) * 100
cat("A. rubra vs Control -- % less consumed:", round(pct_less_arubra, 2), "\n")
