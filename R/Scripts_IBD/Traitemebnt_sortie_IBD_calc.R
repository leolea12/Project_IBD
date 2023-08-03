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


# RESULTAT avec les données actuelle en ajoutant les nouveaux taxons
Results2 <- tibble(read.csv2("IBD_v1.2.4_resultats.csv", sep = ";", skip = 1)) %>%
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
  ggplot(aes(x = factor(year), y = RESULTAT)) +
  geom_boxplot(aes(color = year, fill = grp)) +
  theme_minimal() +
  facet_wrap(~grp, nrow = 2)

tab %>%
  dplyr::filter(LIB_PAR == "IndiceBioDiat") %>%
  group_by(CODE_STATION, grp) %>%
  mutate(mean_result = mean(RESULTAT)) %>%
  ungroup() %>%
  ggplot(aes(y = mean_result)) +
  geom_boxplot() +
  facet_wrap(~grp)


# Données Normales ?
library(ggpubr)
ggqqplot(tab %>% dplyr::filter(grp == 1, LIB_PAR == "IndiceBioDiat") %>% pull(RESULTAT))
ggqqplot(tab %>% dplyr::filter(grp == 2, LIB_PAR == "IndiceBioDiat") %>% pull(RESULTAT))

# NON, En plus, données dépendantes car échantillonnées sur les mêmes stations. Test de Wilcoxon

boxplot(tab %>% dplyr::filter(grp == 1, LIB_PAR == "IndiceBioDiat") %>% pull(RESULTAT))
boxplot(tab %>% dplyr::filter(grp == 2, LIB_PAR == "IndiceBioDiat") %>% pull(RESULTAT))

wilcox.test(tab %>% dplyr::filter(grp == 1, LIB_PAR == "IndiceBioDiat") %>% pull(RESULTAT),
            tab %>% dplyr::filter(grp == 2, LIB_PAR == "IndiceBioDiat") %>% pull(RESULTAT))


Wilcox_tab_test <- tab %>% dplyr::filter(LIB_PAR == "IndiceBioDiat") %>% select(RESULTAT, year, grp)

# Step 1: Filter the data for each year and perform the Wilcoxon rank-sum test
test_results <- Wilcox_tab_test %>%
  group_by(year) %>%
  summarise(wilcox_p_value = wilcox.test(RESULTAT ~ grp)$p.value)

# Step 2: Print the test results
print(test_results)

# Nous avons regardé si il existe une différence d'indice entre les années des deux jeux de données








  



