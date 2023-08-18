library(tidyverse)

dependencies <- c("dplyr")

loadDependencies <- function(dependencies) {
  suppressAll <- function(expr) {
    suppressPackageStartupMessages(suppressWarnings(expr))
  }
  
  lapply(dependencies,
         function(x)
         {
           suppressAll(library(x, character.only = TRUE))
         }
  )
  invisible()
}

loadDependencies(dependencies)


# RESULTAT avec les données actuelle sans ajouter les nouveaux taxons
Results <- tibble(read.csv2("IBD_v1.2.4_resultats.csv", sep = ";", skip = 1)) %>%
  mutate(RESULTAT = as.numeric(RESULTAT),
         grp = 1,
         year = lubridate::year(DATE))

Results %>% dplyr::filter(LIB_PAR == "IndiceBioDiat") %>%
  ggplot(aes(x = factor(year), y = RESULTAT)) +
  geom_boxplot(aes(color = year)) +
  theme_minimal()

Results %>% dplyr::filter(LIB_PAR == "NbTaxonsIBDcontributifs") %>%
  ggplot(aes(x = factor(year), y = RESULTAT)) +
  geom_boxplot(aes(color = year)) +
  theme_minimal()

Results %>% group_by(CODE_STATION) %>%
  dplyr::filter(LIB_PAR == "IndiceBioDiat") %>%
  mutate(mean_result = mean(RESULTAT)) %>%
  ggplot(aes(y = mean_result)) +
  geom_boxplot()

# chargement des HER

load("R/Scripts_IBD/Stat_IBD.Rda")

HER_TAB <- Results %>% 
  mutate(DATE = as.Date(DATE)) %>%
  left_join(Stat_IBD %>% distinct(), 
            by = c("CODE_STATION" = "CODE_STATION", "DATE" = "DATE"),
            relationship = "many-to-many") %>%
  drop_na()

# visualisation HER/ années

HER_TAB %>% 
  dplyr::filter(LIB_PAR == "IndiceBioDiat") %>%
  ggplot(aes(y = RESULTAT, fill = NomHER1)) +
  geom_boxplot()+
  facet_wrap(~NomHER1)+
  theme_minimal() +
  theme(legend.position = "none")


# RESULTAT avec les données actuelle en ajoutant les nouveaux taxons
Results2 <- tibble(read.csv2("IBD_v1.2.4_resultats2.csv", sep = ";", skip = 1)) %>%
  mutate(grp = 2,
         year = lubridate::year(DATE),
         RESULTAT = as.numeric(RESULTAT))

Results2 %>% dplyr::filter(LIB_PAR == "IndiceBioDiat") %>%
  ggplot(aes(x = factor(year), y = RESULTAT)) +
  geom_boxplot(aes(color = year)) +
  theme_minimal()

Results2 %>% dplyr::filter(LIB_PAR == "NbTaxonsIBDcontributifs") %>%
  ggplot(aes(x = factor(year), y = RESULTAT)) +
  geom_boxplot(aes(color = year)) +
  theme_minimal()

Results2 %>% group_by(CODE_STATION) %>%
  dplyr::filter(LIB_PAR == "IndiceBioDiat") %>%
  mutate(mean_result = mean(RESULTAT)) %>%
  ggplot(aes(y = mean_result)) +
  geom_boxplot()

HER_TAB2 <- Results2 %>% 
  mutate(DATE = as.Date(DATE)) %>%
  left_join(Stat_IBD %>% distinct(), 
            by = c("CODE_STATION" = "CODE_STATION", "DATE" = "DATE"),
            relationship = "many-to-many") %>%
  drop_na()

# visualisation HER/ années

HER_TAB2 %>% 
  dplyr::filter(LIB_PAR == "IndiceBioDiat") %>%
  ggplot(aes(y = RESULTAT, fill = NomHER1)) +
  geom_boxplot()+
  facet_wrap(~NomHER1)+
  theme_minimal() +
  theme(legend.position = "none")



# Comparaison des deux 

tab <- rbind(Results, Results2)

tab %>%
  dplyr::filter(LIB_PAR == "NbTaxonsIBDcontributifs") %>%
  mutate(year = lubridate::year(DATE)) %>%
  ggplot(aes(x = factor(year), y = RESULTAT)) +
  geom_boxplot(aes(color = year, fill = grp)) +
  theme_minimal() +
  facet_wrap(~grp, nrow = 2)

tab %>%
  dplyr::filter(LIB_PAR == "IndiceBioDiat") %>%
  mutate(year = lubridate::year(DATE)) %>%
  ggplot(aes(x = factor(year), y = RESULTAT, fill = factor(grp))) +  # Mapping grp to fill
  geom_boxplot(color = "black") +  # Adding black border to boxplots
  scale_fill_discrete(name = "grp") +  # Customizing legend title
  theme_inrae()

tab %>%
  dplyr::filter(LIB_PAR == "IndiceBioDiat") %>%
  group_by(CODE_STATION, grp) %>%
  mutate(mean_result = mean(RESULTAT)) %>%
  ungroup() %>%
  ggplot(aes(y = mean_result)) +
  geom_boxplot() +
  facet_wrap(~grp)

tab_HER <- rbind(HER_TAB, HER_TAB2) 
 

tab_HER %>%
  dplyr::filter(LIB_PAR == "IndiceBioDiat") %>%
  mutate(year = lubridate::year(DATE)) %>%
  ggplot(aes(
    # x = factor(year), 
    y = RESULTAT)) +
  geom_boxplot(aes(
    # color = year, 
    fill = grp)) +
  theme_minimal() +
  facet_grid(grp~NomHER1)+
  theme(legend.position = "none")

tab_HER %>%
  dplyr::filter(LIB_PAR == "IndiceBioDiat") %>%
  mutate(year = lubridate::year(DATE)) %>%
  ggplot(aes(x = factor(NomHER1), y = RESULTAT, fill = factor(grp))) +  # Mapping grp to fill
  geom_boxplot(color = "black") +  # Adding black border to boxplots
  scale_fill_discrete(name = "grp") +  # Customizing legend title
  theme_inrae() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))
  


# Données Normales ?
library(ggpubr)
ggqqplot(tab %>% dplyr::filter(grp == 1, LIB_PAR == "IndiceBioDiat") %>% pull(RESULTAT))
ggqqplot(tab %>% dplyr::filter(grp == 2, LIB_PAR == "IndiceBioDiat") %>% pull(RESULTAT))

# NON, En plus, données dépendantes car échantillonnées sur les mêmes stations. Test de Wilcoxon


tab %>%
  dplyr::filter(LIB_PAR == "IndiceBioDiat") %>%
  mutate(year = lubridate::year(DATE)) %>%
  ggplot(aes(x = factor(grp), y = RESULTAT, fill = factor(grp))) +
  geom_boxplot(color = "black") +
  geom_hline(aes(yintercept = mean(RESULTAT)), color = "black", linetype = "dashed", size = 0.8) +  # Add horizontal lines for means
  scale_fill_discrete(name = "IBD") +
  labs(y = "Notes") +
  theme_inrae() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))

wilcox.test(tab %>% dplyr::filter(grp == 1, LIB_PAR == "IndiceBioDiat") %>% pull(RESULTAT),
            tab %>% dplyr::filter(grp == 2, LIB_PAR == "IndiceBioDiat") %>% pull(RESULTAT), paired = TRUE)

# Assuming your tibble is named 'your_tibble'
tab %>% dplyr::filter(LIB_PAR == "IndiceBioDiat") %>%
  dplyr::select(grp, RESULTAT, CODE_OPERATION) %>%
  pivot_wider(names_from = grp, values_from = RESULTAT) %>%
  unnest() %>%
  rename(grp1 = `1`, grp2 = `2`) %>%
  mutate(IBD1 = rank(grp1),
         IBD2 = rank(grp2)) %>%
  ggplot(aes(x = IBD2, y = IBD1)) +
  geom_point(size = 0.5) +
  geom_smooth(method = "lm", color = "red", se = FALSE) +
  theme_inrae()

glm_tab <- tab %>% dplyr::filter(LIB_PAR == "IndiceBioDiat") %>%
  dplyr::select(grp, RESULTAT, CODE_OPERATION) %>%
  pivot_wider(names_from = grp, values_from = RESULTAT) %>%
  unnest() %>%
  rename(grp1 = `1`, grp2 = `2`)

result <- cor.test(glm_tab$grp1, glm_tab$grp2, method = "spearman")
print(result)

mod=lm(grp1~grp2,data=glm_tab)
abline(mod,col='red',lw=2)
hist(mod$res)

plot(mod)


tab %>% dplyr::filter(LIB_PAR == "IndiceBioDiat") %>%
  dplyr::select(grp, RESULTAT, CODE_OPERATION) %>%
  pivot_wider(names_from = grp, values_from = RESULTAT) %>%
  unnest() %>%
  rename(grp1 = `1`, grp2 = `2`) %>%
  mutate(diff = grp1-grp2) %>%
  ggplot(aes(y = diff)) +
  geom_boxplot() +
  theme_inrae() +
  labs(y = "Différence de notes IBD1/IB2")

diff_cor_tab <- tab %>% dplyr::filter(LIB_PAR == "IndiceBioDiat") %>%
  dplyr::select(grp, RESULTAT, CODE_OPERATION) %>%
  pivot_wider(names_from = grp, values_from = RESULTAT) %>%
  unnest() %>%
  rename(grp1 = `1`, grp2 = `2`) %>%
  mutate(diff = grp1-grp2)

p1 <- diff_cor_tab %>%
  ggplot(aes(y = diff, x = grp1)) +
  geom_point(size = 0.8) +
  theme_inrae() +
  labs(y = "différence", x = "Notes IBD1")

result2 <- cor.test(diff_cor_tab$diff, diff_cor_tab$grp1, method = "spearman")
result2

grob1 <- grobTree(textGrob("A", x=0.02,  y=0.95, hjust=0,gp=gpar(col="red", fontsize=15, fontface="italic")))

p1 <- p1 + annotation_custom(grob1)


p2 <- diff_cor_tab %>%
  ggplot(aes(y = diff, x = grp2)) +
  geom_point(size = 0.8) +
  theme_inrae() +
  labs(y = "Différence IBD1/IBD2", x = "Notes IBD2")

result3 <- cor.test(diff_cor_tab$grp2, diff_cor_tab$diff, method = "spearman")
result3

grob2 <- grobTree(textGrob("B", x=0.02,  y=0.95, hjust=0,gp=gpar(col="red", fontsize=15, fontface="italic")))

p2 <- p2 + annotation_custom(grob2)

p1 + p2

# differences entre les indices

diff_ind <- tab %>% dplyr::filter(LIB_PAR == "IndiceBioDiat") %>%
  dplyr::select(grp, RESULTAT, CODE_OPERATION) %>%
  pivot_wider(names_from = grp, values_from = RESULTAT) %>%
  unnest() %>%
  rename(grp1 = `1`, grp2 = `2`) %>%
  mutate(diff = grp1-grp2) %>% pull(diff)

qqnorm(diff_ind, pch = 1, frame = FALSE)
boxplot(diff_ind)

p_val = NULL

for(i in 1:1000){
  print(i)
  p_value = shapiro.test(sample(diff_ind, 5000, replace = FALSE))$p.value
  p_val = c(p_val, p_value)
}


wilcox.test(diff_ind, mu = 0, alternative = "two.sided")

mean(diff_ind) # Difference significative mais à peine de -0.1 !
sd(diff_ind)



# Test simple entre les années des deux tableaux
Wilcox_tab_test <- tab %>% dplyr::filter(LIB_PAR == "IndiceBioDiat") %>% select(RESULTAT, year, grp)

# Step 1: Filter the data for each year and perform the Wilcoxon rank-sum test
test_results <- Wilcox_tab_test %>%
  group_by(year) %>%
  summarise(wilcox_p_value = wilcox.test(RESULTAT ~ grp, paired = TRUE)$p.value)

# Step 2: Print the test results
print(test_results)

# Nous avons regardé si il existe une différence d'indice entre les années des deux jeux de données

# Test entre les HER
HER_Wilcox_tab_test <- tab_HER %>% 
  dplyr::filter(LIB_PAR == "IndiceBioDiat") %>% select(RESULTAT, NomHER1, grp)

# Step 1: Filter the data for each year and perform the Wilcoxon rank-sum test
HER_test_results <- HER_Wilcox_tab_test %>%
  group_by(NomHER1) %>%
  summarise(wilcox_p_value = wilcox.test(RESULTAT ~ grp, paired = TRUE)$p.value)


