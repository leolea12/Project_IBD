# Chargement des packages

library(data.table)
library(ade4)
library(vegan)
library(factoextra)
library(corrplot)
library(RVAideMemoire)
require(tidyverse)
conflicted::conflict_prefer(name = "filter", winner = "dplyr")
conflicted::conflict_prefer(name = "rename", winner = "dplyr")
conflicted::conflict_prefer(name = "mutate", winner = "dplyr")
conflicted::conflict_prefer(name = "arrange", winner = "dplyr")
require(InraeThemes)
require(readxl)
require(cooccur)
library(lubridate)
conflicted::conflict_prefer(name = "year", winner = "lubridate")
library(ggrepel)
library(devtools)
library(patchwork)
library(sf)
library(leaflet)
library(packcircles)
library(viridis)
library(htmlwidgets)
library(oceanis)
library(png)
library(leafpop)
library(RMySQL)
library(RPostgreSQL)
library(RODBC)
library(tmap)
library(igraph)
library(dplyr)
library(DescTools)
library(GGally)
library(network)
library(sna)
library(ggplot2)
library(intergraph)
library(ggraph)
library(FactoMineR)
library(ade4)
library(ggnetwork)
library(labdsv)
library(tictoc)
library(vegan)
library(fastcluster)
library(ggdendro)
library(cluster)
library (vegan3d)
library(geometry)
library(rgl)
library(shipunov)
library(RODBC)
library(RMySQL)
library(RSQLite)
library(sqldf)
library(svMisc)
library(labdsv)
library(hBayesDM)
library(overlapping)
# library(DataScienceR)


# # Chargement des données 
load("data/Donnees_utilisables/Donnees_completes.RData")

# Construction des fichiers et des objets qui serviront soit dans les calculs, soit dans les graphiques --------------------------------------------------------------------

Chim <- Chimie %>% filter(!is.na(Mediane)) %>% 
  mutate(site = paste0(CODE_STATION,"_", 
                       lubridate::year(DATE), "_0", 
                       lubridate::month(DATE), "_0", 
                       lubridate::day(DATE))) %>% 
  mutate(Nom_parametre = case_when(Code_parametre == 1304 ~ "Conduct_20",
                                   Code_parametre == 1302 ~ "PH",
                                   Code_parametre == 1433 ~ "PO4",
                                   Code_parametre == 1335 ~ "Ammonium",
                                   Code_parametre == 1314 ~ "DCO",
                                   Code_parametre == 1340 ~ "Nitrates",
                                   Code_parametre == 1311 ~ "Ox_dissous")) %>% 
  select(-DATE,-CODE_STATION,-Code_Unite_mesure)

`%notin%` <- Negate(`%in%`)

Communes <- c("Mana", "Roura", "Regina", "Maripasoula", "Papaichton", "Kourou",
              "Saint-Laurent-du-Maroni", "Iracoubo", "Montsinery-Tonnegrande", "Apatou")

code <- as_tibble(read.csv2("data/Donnees_utilisables/CODE_OMNIDIA_traites.csv", stringsAsFactors = FALSE)) %>% 
  select(-X) %>% 
  mutate(code_taxon = code_omnidia)

# On recupere les classes trophiques des taxons qui nous serviront a posteriori
Trophie <- read_excel(path = "data/Donnees_utilisables/Classes_trophie.xlsx", 
           sheet = "Numbers") %>% filter(parameter == "PO4") %>% 
  dplyr::rename("code_taxon" = "code") %>% 
  mutate(Class = as.character(Class)) %>% 
  select(code_taxon, Class)
# 

tab_NMDS <- tibble(Flore %>% ungroup() %>% distinct() %>% 
                     filter(year(DATE) >= 2015) %>% 
                     left_join(Trophie, by = "code_taxon") %>%
                     mutate(code_taxon = paste0(code_taxon,"_",Class)) %>%
                     # filter(pourcent_count_tot > 5) %>%
                     # filter(pourcent_count_year <= 60) %>% # On enleve les espèces rares et ubiquistes de l'annee
                     mutate(Location = paste0(CODE_STATION,"_", year(DATE), "_0", lubridate::month(DATE), "_0", lubridate::day(DATE))) %>%
                     dplyr::select(Location, code_taxon, RESULTAT) %>%
                     pivot_wider(names_from = code_taxon, values_from = RESULTAT, values_fill = 0))

tab_NMDS_CEUO <- tibble(Flore %>% ungroup() %>% distinct() %>% 
                          filter(year(DATE) <= 2009) %>% 
                          left_join(Trophie, by = "code_taxon") %>%
                          mutate(code_taxon = paste0(code_taxon,"_",Class)) %>%
                          # filter(pourcent_count_tot > 5) %>%
                          # filter(pourcent_count_year <= 60) %>% # On enleve les espèces rares et ubiquistes de l'annee
                          mutate(Location = paste0(CODE_STATION,"_", year(DATE), "_0", lubridate::month(DATE), "_0", lubridate::day(DATE))) %>%
                          dplyr::select(Location, code_taxon, RESULTAT) %>%
                          pivot_wider(names_from = code_taxon, values_from = RESULTAT, values_fill = 0))

NMDS_test <- data.frame(tab_NMDS) %>%
  column_to_rownames(var = "Location") %>% 
  filter(rowSums(across(where(is.numeric)))!=0)

colnames(NMDS_test) <- sub("_.*","", colnames(NMDS_test))

NMDS_test_CEUO <- data.frame(tab_NMDS_CEUO) %>%
  column_to_rownames(var = "Location") %>% 
  filter(rowSums(across(where(is.numeric)))!=0)

NMDS_test <- decostand(NMDS_test, method = "hellinger")
dist_bc <- parDist(as.matrix(NMDS_test), 
                   method="bray", 
                   threads = 20)

NMDS_test_CEUO <- decostand(NMDS_test_CEUO, method = "hellinger")
dist_bc_CEUO <- parDist(as.matrix(NMDS_test_CEUO), 
                        method="bray", 
                        threads = 20)

clus_ward <- hclust(dist_bc, method="ward.D2")
clus_ward_CEUO <- hclust(dist_bc_CEUO, method="ward.D2")

Dendo <- clus_ward %>% as.dendrogram()
Dendo_CEUO <- clus_ward_CEUO %>% as.dendrogram()


New_taxons <- as_tibble(read_excel("data/Donnees_utilisables/Fichier_DREAL_1.xlsx",
                                   sheet = "Non_contributifs")) %>% 
  dplyr::select(abre, Profil_assoc) %>%
  dplyr::rename(AFNOR = Profil_assoc) %>% 
  left_join(as_tibble(read.csv2("data/Donnees_utilisables/table_transcodage.csv", stringsAsFactors = FALSE)) %>% # On change l'ancien code 4 lettres
              dplyr::select(AFNOR = "abre", True_name = "CodeValid"), by = "AFNOR") %>% 
  dplyr::mutate(AFNOR = if_else(is.na(True_name == T), AFNOR, True_name)) %>% 
  dplyr::select(-True_name) %>% 
  left_join(as_tibble(read.csv2("data/Donnees_utilisables/table_transcodage.csv", stringsAsFactors = FALSE)) %>% # On change l'ancien code 4 lettres
              dplyr::select(abre = "abre", True_name = "CodeValid"), by = "abre") %>% 
  mutate(abre = if_else(is.na(True_name == T), abre, True_name)) %>% 
  dplyr::select(-True_name)


# On recupere pour chaque taxon ses sites de presence
site_taxon <- NMDS_test %>%
  rownames_to_column(var = "site") %>%
  pivot_longer(!site, names_to = "taxon", values_to = "Ab_relative") %>%
  mutate(taxon = sub("_.*", "", taxon)) %>% 
  filter(Ab_relative > 0)

site_taxon_CEUO <- NMDS_test_CEUO %>%
  rownames_to_column(var = "site") %>%
  pivot_longer(!site, names_to = "taxon", values_to = "Ab_relative") %>%
  mutate(taxon = sub("_.*", "", taxon)) %>% 
  filter(Ab_relative > 0)

save(dist_bc, file = "data/Donnees_utilisables/Dist_mat_2015_2021.RData")
save(dist_bc_CEUO, file = "data/Donnees_utilisables/Dist_mat_2015_2021_CEUO.RData")
save(clus_ward, file = "data/Donnees_utilisables/Cluster_2015_2021.RData")
save(clus_ward, file = "IBD_2022_shiny/Data/Cluster_2015_2021.RData")
save(clus_ward_CEUO, file = "data/Donnees_utilisables/Cluster_CEUO_2015_2021.RData")
save(clus_ward_CEUO, file = "IBD_2022_shiny/Data/Cluster_CEUO_2015_2021.RData")
save(NMDS_test, file = "data/Donnees_utilisables/Donnees_test_serveur.RData")
save(NMDS_test_CEUO, file = "data/Donnees_utilisables/Donnees_test_serveur_CEUO.RData")
save(NMDS_test, file = "IBD_2022_shiny/Data/Donnees_test_serveur.RData")
save(NMDS_test_CEUO, file = "IBD_2022_shiny/Data/Donnees_test_serveur_CEUO.RData")
save(New_taxons, file = "data/Donnees_utilisables/New_taxons.RData")
save(New_taxons, file = "IBD_2022_shiny/Data/New_taxons.RData")
save(Dendo_CEUO, file = "data/Donnees_utilisables/Dendo_CEUO.RData")
save(site_taxon_CEUO, file = "data/Donnees_utilisables/site_taxon_CEUO.RData")
save(Dendo, file = "data/Donnees_utilisables/Dendo.RData")
save(site_taxon, file = "data/Donnees_utilisables/site_taxon.RData")

load("data/Donnees_utilisables/Cluster_2015_2021.RData")
load("data/Donnees_utilisables/Donnees_test_serveur.RData")
load("data/Donnees_utilisables/Donnees_test_serveur_CEUO.RData")

# Traitement -----------------------------------------------------

load("data/Donnees_utilisables/Dendo_CEUO.RData")
load("data/Donnees_utilisables/site_taxon_CEUO.RData")
load("data/Donnees_utilisables/Dendo.RData")
load("data/Donnees_utilisables/site_taxon.RData")
load("data/Donnees_utilisables/site_taxon_CEUO.RData")
load("data/Donnees_utilisables/Donnees_test_serveur_CEUO.RData")
load("data/Donnees_utilisables/Donnees_test_serveur.RData")
# load("data/Donnees_utilisables/Mat_dist_CAH.RData")



get_density <- function(x, y) {
  dens <- MASS::kde2d(x, y)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

extract_density <- function(Dend, tax, tab){
  
  select_tax <- tab %>% filter(taxon == tax)
  
  no_select_tax <- tab %>% filter(site %notin% unique(select_tax %>% pull(site))) %>%
    distinct(site)
  
  label_select_tax <- data.frame(site = labels(Dend)) %>%
    left_join(select_tax, by = "site") %>%
    mutate(site = if_else(taxon == tax, site, ""))
  
  # Construire un pas entre chaque site pour rééchantillonner ensuite (sorte de probabilité)
  
  tax_density <- label_select_tax %>% mutate(site = seq(1:length(label_select_tax$site)),
                                             Ab_relative = if_else(is.na(Ab_relative) == TRUE, 0, Ab_relative)) %>% 
    drop_na() %>% 
    mutate(density = get_density(site, Ab_relative))
  
  return(tax_density)
}

# On ne va récupérer que les taxons vu au minimum 10 fois, et encore niveau stats
# c'est chaud donc on appliquera d'autres filtres par la suite 
# FULL_TAX <- data.frame(x = sub("_.*","", colnames(NMDS_test))) %>% 
#   filter(x %in% unique(site_taxon %>% select(-Ab_relative) %>% group_by(taxon) %>% 
#                          summarise(count = n()) %>% 
#                          filter(count>10) %>% pull(taxon))) %>% 
#   pull(x)

FULL_TAX <- data.frame(x = sub("_.*","", colnames(NMDS_test))) %>%
  filter(x %in% unique(site_taxon %>% select(-Ab_relative) %>% group_by(taxon) %>%
                         summarise(count = n()) %>%
                         filter(count>10) %>% pull(taxon))) %>%
  pull(x)

# tab_density = do.call(rbind, 
#                       map(FULL_TAX, 
#                           function(FULL_TAX) 
#                             extract_density(Dend = Dendo, tax = FULL_TAX, tab = site_taxon)))

tab_density = do.call(rbind,
                      map(FULL_TAX,
                          function(FULL_TAX)
                            extract_density(Dend = Dendo, tax = FULL_TAX, tab = site_taxon)))

save(tab_density, file = "data/Donnees_utilisables/tab_density.RData")

# tab_density %>% group_by(taxon) %>% 
#   summarise(somme = sum(density))

# Le tableau des densités récupéré permettra la comparaison statistique des densités 
# lorsque d'un taxon sera sans profil DREAL.

# Taxon sans profil: 
# Les taxons pour lesquels moins de 30 recensements ont eu lieu sont enlevés 
# de l'analayse. Pour le rest, les taxons vu 50 fois ou plus sont echantillonnés à 50,
# Ceux vu moins de 50 fois sont boostrappés pour atteindre la taille d'échantillon
# minimale nécessaire pour appliquer un tets de kolmogorov smirnov. Avoir des tailles
# d'échantillons similaire sera plus intéressant.
# On construit donc un jeu avec des echantillons tous de taille 50 et on applique des 
# pairwise ks test 1000 fois. 

# Récupérer les densités et comparer avec kolmogorov 


library(svMisc)
library(parallel)
library(doParallel)
library(foreach)
library(DataScienceR)
library(tidyverse)
library(progress)
library(igraph)
library(tictoc)
load("data/Donnees_utilisables/tab_density.RData")

#Setup backend to use many processors
totalCores = detectCores()

#Leave one core to avoid overload your computer
cluster <- makePSOCKcluster(20, outfile="") 
registerDoParallel(cluster)

N <- 1000
# Density_list <- vector("list", N)
Density_list<- vector("list", N)
# dt<-tibble(value=NA,group=NA,i=NA)

pb <- txtProgressBar(min=1, max=N, style=3)

res <- foreach(
  i = 1:N,
  .combine = 'rbind',
  .packages= c('tidyverse','DataScienceR', 'svMisc','base', 'dplyr','igraph')
) %dopar% {
  
  # Sleep for 0.1 seconds
  # Sys.sleep(0.01)
  # print(i)
  
  setTxtProgressBar(pb, i)
  
  ech1 <- tab_density %>% 
    dplyr::group_by(taxon) %>% 
    dplyr::mutate(size = n()) %>% 
    filter(size >= 30) %>% 
    dplyr::slice_sample(n=50, weight_by = density, replace = TRUE) %>% 
    dplyr::ungroup() %>% 
    dplyr::arrange(taxon)
  
  full_ech <- ech1
  
  value <- full_ech$density
  group <- as.factor(full_ech$taxon)
  
  sink("nul")
  Density_compare <- invisible(pairwise_ks_test(value, group, warning = -1))
  sink()
  # 
  as_data_frame(graph_from_adjacency_matrix(Density_compare, weighted = TRUE)) %>%
    rename(taxon = from, taxon2 = to, pval = weight) %>% 
    select(taxon, taxon2, pval) %>%
    filter(pval >= 0.05) %>%
    filter(taxon != taxon2)
}


Density_tab <- res %>% 
  dplyr::arrange(taxon, taxon2, desc(pval))

stopCluster(cluster)

# for(i in 1:N){
#
#   progress(i, N)
#
#   # On echantillonne 50 pour tout le monde
#   ech1 <- tab_density %>% # tab_density
#     group_by(taxon) %>%
#     # mutate(s_density = sum(density),
#     #        rapport_density = (density/s_density)*100) %>%
#     mutate(size = n()) %>%
#     filter(size >= 30) %>%
#     slice_sample(n=50, weight_by = density, replace = TRUE) %>%
#     ungroup() %>%
#     dplyr::arrange(taxon)
#
#   # ech2 <- tab_density %>% # tab_density
#   #   group_by(taxon) %>%
#   #   mutate(size = n()) %>%
#   #   filter(size %in% (30:49)) %>%
#   #   mutate(samp_size = 50 - size) %>%
#   #   dplyr::arrange(taxon)
#   #
#   # vec_ech2 <- ech2 %>%
#   #   distinct(samp_size) %>%
#   #   ungroup() %>%
#   #   pull(samp_size)
#   #
#   # ech2_sampled <- ech2 %>% group_split() %>%
#   #   map2_dfr(vec_ech2, ~ slice_sample(.x, n = .y, weight_by = density, replace = TRUE))
#   #
#   # ech2 %>% slice_sample(n=50, weight_by = density, replace = TRUE)
#
#   full_ech <- ech1
#     # add_row(ech2 %>%
#     #           ungroup() %>%
#     #           add_row(ech2_sampled) %>%
#     #           select(-samp_size)) %>%
#     # group_by(taxon) %>%
#     # dplyr::arrange(taxon, site) %>%
#     # ungroup()
#
#   value <- full_ech$density
#   group <- as.factor(full_ech$taxon)
#   Density_compare <- pairwise_ks_test(value, group, warning = -1)
#
#   # dt<-bind_rows(dt,tibble(value=value,group=group,i=i))
#   #
#   Density_extract <- data.frame(Density_compare) %>%
#     rownames_to_column(var = "taxon") %>%
#     pivot_longer(cols = 2:nrow(.), names_to = "taxon2", values_to = "pval") %>%
#     select(taxon, taxon2, pval) %>%
#     filter(pval >= 0.05) %>%
#     filter(taxon != taxon2)
#   #
#   Density_list[[i]] <- Density_extract
# }

# Vu avec seb, proposition
# tic()
# t <- data.frame(pairwise_ks_test(dt$value[-1], factor(dt$group[-1]), warning = -1)) %>%
#     rownames_to_column(var = "taxon") %>%
#     pivot_longer(cols = 2:nrow(.), names_to = "taxon2", values_to = "pval") %>%
#     select(taxon, taxon2, pval) %>%
#     filter(pval >= 0.05) %>%
#     filter(taxon != taxon2)

# Density_tab <- bind_rows(res) %>% dplyr::arrange(taxon, taxon2, desc(pval))

# toc()


# Après avoir récupéré et fusionné les tableaux des test contenant les pvalues significatives,
# on va regarder avec quel taxon chaque taxon est le plus proche (pour ceux qui n'ont pas de profil)

#Tableau contenant toutes les paires significatives des taxons
# Density_tab <- bind_rows(Density_list) %>% dplyr::arrange(taxon, taxon2, desc(pval)) 


# save(Density_tab, file = "data/Donnees_utilisables/Densité_comparées_bootsrap.RData")
save(Density_tab, file = "data/Donnees_utilisables/Densité_comparées_bootsrap.RData")

# load("data/Donnees_utilisables/Densité_comparées_bootsrap.RData")

# Density_tab
# Recuperation ds pourcentage pour 500 simulations
# Check_profile <- Density_tab %>% group_by(taxon, taxon2) %>% # Density_tab
#   mutate(mean_pval = mean(pval)) %>% 
#   mutate(number = n()) %>% 
#   ungroup() %>% 
#   distinct(taxon, taxon2, mean_pval, number) %>% 
#   dplyr::arrange(taxon, desc(number)) %>% 
#   mutate(percent_over_500 = (number*100)/500)

Check_profile <- Density_tab %>% group_by(taxon, taxon2) %>% # Density_tab_CEUO
  mutate(med_pval = median(pval)) %>% 
  mutate(number = n()) %>% 
  ungroup() %>% 
  distinct(taxon, taxon2, med_pval, number) %>% 
  dplyr::arrange(taxon, desc(number))
# mutate(percent_over_500 = (number*100)/500)

# save(Check_profile, file = "data/Donnees_utilisables/Fiche_identification.RData")
save(Check_profile, file = "data/Donnees_utilisables/Chech_profile.RData")

# Maintenant, il s'agit de récupérer le pourcentage de sites en commun dans la CAH entre
# deux taxons car en effet, deux taxons peuvent provenir de la même distribution mais
# rester complètement différents

# A partir du tab density intialement construit, on peut aller chercher les numéros de sites 
# qui se correspondent entre les taxons 

paired_site <- function(tax1, tax2, tab){
  pair1 <- tab %>% dplyr::arrange(site) %>% 
    filter(taxon == tax1) %>% 
    select(site, taxon) %>% 
    mutate(n_site1 = length(site))
  
  pair2 <- tab %>% dplyr::arrange(site) %>% 
    filter(taxon == tax2) %>% 
    select(site, taxon2 = taxon) %>% 
    mutate(n_site2 = length(site))
  
  joined <- pair1 %>% left_join(pair2, by = "site") %>% 
    mutate(n_site2 = pair2$n_site2[1]) %>% #On remet le nombre de site au cas ou ils n'en n'ont pas en commun
    drop_na() %>% mutate(nb_site_com = length(site)) %>% 
    distinct(taxon,taxon2,nb_site_com, n_site1, n_site2) %>% 
    mutate(part_tax1 = round((nb_site_com/length(pair1$site))*100,3),
           part_tax2 = round((nb_site_com/length(pair2$site))*100,3))
  
  return(joined)
}

# RECUPERER LES POURCENTAGE ET LES NOMBRES DE SITES COMMUN DANS LE TABLEAU 
# CHECK PROFILE

data_to_fill <- data.frame(taxon = "taxon", taxon2= "taxon2", nb_site_com= 0, part_tax1= 0, part_tax2= 0, n_site1 = 0, n_site2 = 0)

for(i in 1:length(Check_profile$taxon)){
  print(i)
  test <- paired_site(tax1 = Check_profile[i,]$taxon, tax2 = Check_profile[i,]$taxon2, tab = tab_density) # tab_density_CEUO
  data_to_fill <- data_to_fill %>% add_row(test)
}

# Tableau contenant la part des sites partagés entre deux taxons aux distributions similaires
# FINAL_PROFILES <- Check_profile %>% left_join(data_to_fill, by = c("taxon","taxon2"))
FINAL_PROFILES <- Check_profile %>% left_join(data_to_fill, by = c("taxon","taxon2"))

# save(FINAL_PROFILES, file = "data/Donnees_utilisables/FINAL_PROFILES.RData")
save(FINAL_PROFILES, file = "data/Donnees_utilisables/FINAL_PROFILES.RData")

# On va construire un dernière variable qui va évaluer la distance entre les sites de la CAH pour deux taxons. 
# Cela devrait nous permettre de mieux décider quel profil affilier et permettre de faciliter une affiliation
# lorsque les distributions sont similaires et que les sites ne sont pas partagés mais très proches

# data_overlap <- data.frame(taxon = "NULL", taxon2 = "NULL", overlap = 0)


# On récupère pour chaque paire de taxon l'overlapping des densités, de façon à 
# ce que deux taxons qui n'ont pas ou très de peu de sites communs puissent être comparables 
# sur le pourcentage de recouvrement de leur densité de sites.

overlaps <- function(x,y){
  l <- list(X1 = sort(tab_density %>% filter(taxon == x) %>% pull(site)),# tab_density
            X2 = sort(tab_density %>% filter(taxon == y) %>% pull(site)))
  out <- round(overlapping::overlap(l,plot=FALSE, type = "2")$OV,2)
}

# FINALE_PROFILES_DENSITY <- FINAL_PROFILES %>% mutate(overlapping = mapply(function(x, y) overlaps(x,y), 
# taxon, 
# taxon2))

FINALE_PROFILES_DENSITY <- FINAL_PROFILES %>% mutate(overlapping = mapply(function(x, y) overlaps(x,y), 
                                                                          taxon, 
                                                                          taxon2))

# save(FINALE_PROFILES_DENSITY, file = "data/Donnees_utilisables/FINALE_PROFILES_DENSITY.RData")
save(FINALE_PROFILES_DENSITY, file = "data/Donnees_utilisables/FINALE_PROFILES_DENSITY.RData")


# ASSOCIER UN RANG A CHAQUE PARAMETRE
# rank_table <- FINALE_PROFILES_DENSITY %>% 
#   replace(is.na(.), 0) %>% 
#   # On trouve les paires et on les ecrit dans le meme ordre pour supprimer les doublons
#   mutate(paire = map2_chr( taxon, taxon2, ~str_flatten(sort(c(.x,.y))) ) ) %>% 
#   distinct(paire, .keep_all = TRUE) %>% 
#   mutate_at(vars(mean_pval:overlapping), ~as.numeric(as.character(.))) %>% 
#   group_by(taxon) %>% 
#   mutate(across(c(mean_pval:overlapping), ~ rank(.x))) %>% 
#   ungroup() %>% 
#   group_by(taxon, taxon2) %>% 
#   rowwise() %>% 
#   mutate(mean_rank = base::mean(c(mean_pval,number,percent_over_500,
#                                     nb_site_com, part_tax1, part_tax2, overlapping))) %>% 
#   select(taxon, taxon2, mean_rank) %>% 
#   ungroup() %>% 
#   group_by(taxon) %>% 
#   dplyr::arrange(taxon, desc(mean_rank)) 

rank_table <- FINALE_PROFILES_DENSITY %>%
  group_by(taxon) |> 
  mutate(dif_n_site = abs(n_site1-n_site2)) |> 
  ungroup() |> 
  left_join(FINALE_PROFILES_DENSITY %>% 
              replace(is.na(.), 0) %>% 
              # On trouve les paires et on les ecrit dans le meme ordre pour supprimer les doublons
              mutate(paire = map2_chr( taxon, taxon2, ~str_flatten(sort(c(.x,.y))) ) ) %>% 
              distinct(paire, .keep_all = TRUE) %>% 
              mutate_at(vars(med_pval:overlapping), ~as.numeric(as.character(.))) %>% 
              mutate(dif_n_site = abs(n_site1-n_site2)) %>%
              group_by(taxon) %>% 
              mutate(across(c(med_pval:overlapping), ~ rank(.x))) %>% 
              mutate(across(c(dif_n_site), ~ rank(-.x))) %>%
              ungroup() %>% 
              # select(-percent_over_1000) %>%
              group_by(taxon, taxon2) %>% 
              rowwise() %>% 
              mutate(mean_rank = base::mean(c(med_pval,number,nb_site_com, part_tax1, part_tax2, overlapping, dif_n_site))) %>% 
              select(taxon, taxon2, mean_rank) %>% 
              ungroup(), by = c("taxon", "taxon2")) |> 
  group_by(taxon) %>% 
  dplyr::arrange(taxon, desc(mean_rank)) |> 
  ungroup() |> 
  select(taxon,taxon2,med_pval,number,nb_site_com,part_tax1, part_tax2, overlapping, dif_n_site, mean_rank, n_site1, n_site2)

save(rank_table, file = "data/Donnees_utilisables/rank_table.RData")

# On va remanier les rangs en ajoutant la pvalue du test de wilcoxon entre 
# les concentrations de PO4 sur les sites pour chaque paire de taxon
rank_table_PO4 <- tab_density %>% 
  select(site, taxon, Ab_relative, density) %>%
  left_join(Chim %>% distinct(site, Nom_parametre, Mediane) %>% # Jointure des donnees phosphore
              filter(Nom_parametre == "PO4"), by = "site")

# Tets de wilcoxon entre les paires
taxons <- as.factor(rank_table_PO4$taxon)
Wilcox_mat <- as.matrix(pairwise.wilcox.test(rank_table_PO4$Mediane, taxons, p.adj = "BH")$p.value)
paires_wilcox <- reshape2::melt(Wilcox_mat) %>%
  rename(taxon = Var1, taxon2 = Var2, pval_PO4 = value) %>% drop_na()

# Ajout des pvalues des tests de wilcoxon et filtration des pvalues non-significative (ce qu'on veut)
Phosphore_pval <- rank_table %>% 
  left_join(paires_wilcox %>% 
              rename(taxon2 = taxon, taxon = taxon2), by = c("taxon", "taxon2")) %>%
  filter(pval_PO4 >= 0.05)

# Dans Phosphore_pval, on a filtre parmis les candidats ceux qui presente une distribution
# similaire de phosphore avec le taxon cible. On conserve cependant les taxons n'ayant pas 
# de profil ecologique car ils permettent de se rendre compte tout de même pour certain 
# qu'ils ressemblent très fortement à des taxons sans profil.

Prof_eco <- tibble(read.csv("data/Donnees_utilisables/IBD_params.csv", header = TRUE, sep = ";")) %>%
  mutate(across(CL1:Val.Ind., ~as.numeric(gsub(",", ".", .)))) %>% 
  select(-DENOMINATION, -Origine, -SANDRE)


# 
# Transcoding <- tibble(read.csv("data/Donnees_utilisables/table_transcodage.csv", header = TRUE, sep = ";")) %>%
#               select(abre, CodeValid)

# On récupère les profils ecologiques des taxons candidats dans l'ordre de leur classement 

test1 <- Phosphore_pval %>% select(taxon, taxon2, pval_PO4) %>% rename(AFNOR = taxon2) %>% 
  filter(taxon == "CSMU")  %>% left_join( 
Prof_eco %>% filter(AFNOR %in% unique(Phosphore_pval %>% 
                                        filter(taxon == "CSMU") %>% 
                                        pull(taxon2))) %>%
  distinct(AFNOR, .keep_all = TRUE), by = "AFNOR") 

test1 <- rank_table %>% select(taxon, taxon2) %>% rename(AFNOR = taxon2) %>% 
  filter(taxon == "CSMU") %>% 
  left_join( 
    Prof_eco %>% filter(AFNOR %in% unique(rank_table %>% 
                                            filter(taxon == "CSMU") %>% 
                                            pull(taxon2))) %>%
      distinct(AFNOR, .keep_all = TRUE), by = "AFNOR") %>%
  slice(1:10)

# On analyse la variabilité des profils écologiques proposés

p1 <- test1 %>% 
  select(AFNOR, CL1,CL2,CL3,CL4,CL5,CL6,CL7) %>%
  drop_na() %>%
  pivot_longer(cols = -AFNOR, 
               names_to = "group_variable",
               values_to = "value") %>%
  ggplot(aes(x = group_variable, y = value)) +
  geom_point(aes(color = as.factor(AFNOR)))+
  geom_line(aes(group = AFNOR, color = as.factor(AFNOR)))

# SORTIR un pourcentage du nombre de profil moyen qui collent avec le profil vrai
p2 <- Prof_eco %>% 
  select(AFNOR, CL1,CL2,CL3,CL4,CL5,CL6,CL7) %>%
  filter(AFNOR == "CSMU") %>%
  distinct() %>%
  pivot_longer(cols = -AFNOR, 
               names_to = "group_variable",
               values_to = "value") %>%
  ggplot(aes(x = group_variable, y = value)) +
  geom_point(aes(color = as.factor(AFNOR)))+
  geom_line(aes(group = AFNOR, color = as.factor(AFNOR)))

p3 <- test1 %>% 
  select(AFNOR, CL1,CL2,CL3,CL4,CL5,CL6,CL7) %>%
  drop_na() %>%
  pivot_longer(cols = -AFNOR, 
               names_to = "group_variable",
               values_to = "value") %>%
  group_by(group_variable) %>%
  summarise(mean_prof = median(value)) %>%
  ggplot(aes(x = group_variable, y = mean_prof)) +
  geom_point()+
  geom_line(group = 1)

p1 + p2 + p3

# LE TABLEAU que je sauvegarde devra être joint au tableau rank_table. Ensuite
# il faudra que les taxons ayant une pvalue > 0.05 soient extrait pour voir 
# si parmis eux, on ne peut pas trouver un taxon qui ressemble à une proposiiton DREAL 
# etant donné qu'aucun d'entre eux ne sort pendant les analyses.


# Il faut donc récupérer le fichier IBD params 

save(paires_wilcox, file = "data/Donnees_utilisables/paires_Wilcox.RData")



# Travaux CAH cut 25% abondance, non conservé par la suite ---------------------

# Initialisation du fichier qui contiendra le resume de tout les cluster
# des taxons
tot_cluster <- data.frame()
tot_cluster_CEUO <- data.frame()

# Initiation dendextend plot
noLabel <- function(x) {
  if (stats::is.leaf(x)) {
    attr(x, "label") <- NULL }
  return(x)
}

get_density <- function(x, y) {
  dens <- MASS::kde2d(x, y)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

plot_dendo <- function(Cluster, Dendo, tax, site_taxon, prob){
  select_tax <- site_taxon %>% filter(taxon == tax)
  no_select_tax <- site_taxon %>% filter(site %notin% unique(select_tax %>% pull(site))) %>%
    distinct(site)
  
  label_select_tax <- data.frame(site = labels(Dendo)) %>%
    left_join(select_tax, by = "site") %>%
    mutate(site = if_else(taxon == tax, site, ""))
  
  plot_dend <- Dendo %>%
    # dendextend::set("labels_cex", 0.5) %>%
    # # set("nodes_pch", 19) %>%
    # # set("nodes_cex", 0.5) %>%
    # dendextend::set("by_labels_branches_col", value = select_tax$site) %>%
    # dendextend::set("by_labels_branches_lwd", value = select_tax$site, TF_values = c(1,0.25))  %>%
    # dendextend::set("by_labels_branches_lty", value = no_select_tax$site, TF_values = c(1,Inf)) %>%
    # dendextend::hang.dendrogram(hang_height = 0.1) %>%
    # dendextend::set("labels_cex", 0.2) %>%
    dendextend::set("labels", "")
  # 
  # 
  # 
  # # Extraction des quantiles en fonction de la densité
  X <- Categorical(label_select_tax %>%
                     mutate(site = seq(1:length(label_select_tax$site)),
                            Ab_relative = if_else(is.na(Ab_relative) == TRUE, 0, Ab_relative)) %>%
                     drop_na() %>% pull(site),
                   p = label_select_tax %>%
                     mutate(site = seq(1:length(label_select_tax$site)),
                            Ab_relative = if_else(is.na(Ab_relative) == TRUE, 0, Ab_relative)) %>%
                     drop_na() %>%
                     mutate(density = get_density(site, Ab_relative)) %>%
                     drop_na() %>% pull(density))
  # 
  Quantil1 <- quantile(X, probs = c(0.45, 0.55))
  Quantil2 <- quantile(X, probs = c(0.375, 0.625))
  Quantil3 <- quantile(X, probs = c(0.25, 0.75))
  Quantil4 <- quantile(X, probs = c(0.125, 0.875))
  
  # vlines <- data.frame(value = c(Quantil1, Quantil2, Quantil3, Quantil4),
  #                      Mean = c("10%","10%", "25%", "25%", "50%","50%", "75%", "75%"))
  
  # plot_density <- label_select_tax %>%
  #   mutate(site = seq(1:length(label_select_tax$site)),
  #          Ab_relative = if_else(is.na(Ab_relative) == TRUE, 0, Ab_relative)) %>%
  #   drop_na() %>%
  #   mutate(density = get_density(site, Ab_relative)) %>%
  #   filter(taxon == tax) %>%
  #   ggplot(aes(x=site)) +
  #   geom_density(fill="grey", color="grey", alpha=0.4)+
  #   geom_vline(data=vlines,
  #              aes(xintercept=value, colour=Mean),
  #              size=1, linetype="dashed", show.legend=TRUE)+
  #   theme_classic()+
  #   theme(axis.text.x=element_blank())
  
  # On créer des lignes quantiles entre lesquelles 50% des points seront présents
  
  # Plot_ab <- label_select_tax %>%
  #   mutate(site = seq(1:length(label_select_tax$site)),
  #          Ab_relative = if_else(is.na(Ab_relative) == TRUE, 0, Ab_relative)) %>%
  #   drop_na() %>%
  #   mutate(density = get_density(site, Ab_relative)) %>%
  #   ggplot()+
  #   geom_point(aes(x = site, y = Ab_relative, color = density), alpha = 0.8)+
  #   scale_color_viridis()+
  #   theme_classic()+
  #   theme(axis.text.x=element_blank())
  
  # # Combiner les deux images 
  # # library(cowplot)
  # # library(magick)
  # # ggsave(Plot_ab, filename = 'data/Plot_ab.png', device = 'png', bg = 'transparent')
  # # ggsave(plot_density, filename = 'data/plot_density.png', device = 'png', bg = 'transparent')
  # # cowplot::plot_grid(Plot_ab, plot_density, nrow = 1)
  # # plot1 <- image_read('data/Plot_ab.png')
  # # plot2 <- image_read('data/plot_density.png')
  # # img <- c(plot1, plot2)
  # # png("data/PLOT.png")
  # # image_mosaic(img)
  # # dev.off()
  # # image_write(image_mosaic(img), path = "data/final.png", format = "png")
  # 
  # Combined_plot <- Plot_ab/plot_density
  
  # DETERMINER A QUEL ENDROIT ON DOIT COUPER
  
  if(prob[1] == 0.45 & prob[2] == 0.55){
    data_tax <- data.frame(height = get_leaves_attr(plot_dend, "height")[Quantil1[1]:Quantil1[2]],
                           label = get_leaves_attr(plot_dend, "label")[Quantil1[1]:Quantil1[2]]) %>% 
      drop_na()}
  
  if(prob[1] == 0.375 & prob[2] == 0.625){
    data_tax <- data.frame(height = get_leaves_attr(plot_dend, "height")[Quantil2[1]:Quantil2[2]],
                           label = get_leaves_attr(plot_dend, "label")[Quantil2[1]:Quantil2[2]]) %>% 
      drop_na()}
  
  if(prob[1] == 0.25 & prob[2] == 0.75){
    data_tax <- data.frame(height = get_leaves_attr(plot_dend, "height")[Quantil3[1]:Quantil3[2]],
                           label = get_leaves_attr(plot_dend, "label")[Quantil3[1]:Quantil3[2]]) %>% 
      drop_na()}
  
  if(prob[1] == 0.125 & prob[2] == 0.875){
    data_tax <- data.frame(height = get_leaves_attr(plot_dend, "height")[Quantil4[1]:Quantil4[2]],
                           label = get_leaves_attr(plot_dend, "label")[Quantil4[1]:Quantil4[2]]) %>% 
      drop_na()}
  
  # RECUPERER TOUS LES NOEUDS COMPRIS DANS L'INTERVAL DE 25%
  Nodes <- data.frame(get_nodes_xy(plot_dend))
  # Si jamais le premier noeud est compris dans l'interval il faut l'enlever
  First_node <- max(Nodes$X2)
  
  # On coupera juste au dessus le noeud le plus haut 
  
  if(prob[1] == 0.45 & prob[2] == 0.55){
    Point_coupe <- round(Nodes %>% filter(X1 <= Quantil1[2] & X1 >= Quantil1[1]) %>% 
                           arrange(X2) %>% 
                           filter(X2 != First_node) %>% 
                           slice(which.max(X2)) %>% 
                           pull(X2), 2) + 0.01}
  
  if(prob[1] == 0.375 & prob[2] == 0.625){
    Point_coupe <- round(Nodes %>% filter(X1 <= Quantil2[2] & X1 >= Quantil2[1]) %>% 
                           arrange(X2) %>% 
                           filter(X2 != First_node) %>% 
                           slice(which.max(X2)) %>% 
                           pull(X2), 2) + 0.01}
  
  if(prob[1] == 0.25 & prob[2] == 0.75){
    Point_coupe <- round(Nodes %>% filter(X1 <= Quantil3[2] & X1 >= Quantil3[1]) %>% 
                           arrange(X2) %>% 
                           filter(X2 != First_node) %>% 
                           slice(which.max(X2)) %>% 
                           pull(X2), 2) + 0.01}
  
  if(prob[1] == 0.125 & prob[2] == 0.875){
    Point_coupe <- round(Nodes %>% filter(X1 <= Quantil4[2] & X1 >= Quantil4[1]) %>% 
                           arrange(X2) %>% 
                           filter(X2 != First_node) %>% 
                           slice(which.max(X2)) %>% 
                           pull(X2), 2) + 0.01}
  
  
  cluster <- data.frame(cutree(Dendo, h = Point_coupe)) %>% 
    rownames_to_column(var = "site_coupe") %>% 
    dplyr::rename(cluster = "cutree.Dendo..h...Point_coupe.") %>% 
    filter(site_coupe %in% unique(data_tax %>% pull(label))) %>% 
    arrange(desc(cluster)) %>% 
    rename_with(~ paste0("site_coupe_",tax), 'site_coupe')
  
  tot_cluster <- rbind(Cluster, cluster)
  # 
  # # Pourcentage de sites en commun entre un taxon et son profil proposé
  # 
  # 
  # label_point_coupe <- data.frame(height = get_leaves_attr(plot_dend, "height")[Quantil4[1]:Quantil4[2]],
  #                                 label = get_leaves_attr(plot_dend, "label")[Quantil4[1]:Quantil4[2]]) %>%
  #   drop_na()
  # # # 
  # # # # Parmis tous les site ou le taxon est présent, on prend celui qui se situe le plus haut
  # # # # ce sera le point de coupe de base !
  # point_coupe_base <- data.frame(height = get_leaves_attr(plot_dend, "height")[Quantil4[1]:Quantil4[2]],
  #                                label = get_leaves_attr(plot_dend, "label")[Quantil4[1]:Quantil4[2]]) %>%
  #   drop_na() %>%
  #   slice_max(height) %>%
  #   mutate(height = round(height,2))
  # 
  # # save plot 
  # # grDevices::pdf(file = paste0("data/CAH_",tax,".pdf"))
  # grDevices::jpeg(file = paste0("data/density_",tax,".jpeg"),width = 800, height = 800)
  # # read_image(...)
  # print(Combined_plot)
  # grDevices::dev.off()
  # 
  # grDevices::jpeg(file = paste0("data/CAH_",tax,".jpeg"), width = 800, height = 800)
  # plot(plot_dend, main = paste0("CAH_",tax))
  # # Ajout des informations sur le taxon sur la CAH (plus haute feuille + bornes)
  # 
  # abline(v = Quantil1[1], col="red", lwd=1.5, lty=2)
  # abline(v = Quantil1[2], col="red", lwd=1.5, lty=2)
  # abline(v = Quantil2[1], col="green", lwd=1.5, lty=2)
  # abline(v = Quantil2[2], col="green", lwd=1.5, lty=2)
  # abline(v = Quantil3[1], col="blue", lwd=1.5, lty=2)
  # abline(v = Quantil3[2], col="blue", lwd=1.5, lty=2)
  # abline(v = Quantil4[1], col="purple", lwd=1.5, lty=2)
  # abline(v = Quantil4[2], col="purple", lwd=1.5, lty=2)
  # abline(h = Point_coupe, col="black", lwd=1.5, lty=1)
  # legend(20,30.5, legend=c("10%","25%", "50%", "75%", paste0("COUPE_", as.character((1 - (prob[2]-prob[1]))*100))),
  #        col=c("red", "green", "blue", "purple", "black"),title = "partition", lty=2, cex=0.8)
  # 
  
  # grDevices::dev.off()
  # 
  # 
  # tbl <- create_table(cluster) %>% 
  #   titles(paste0("Result_",tax,"cut_", as.character(1 - ((prob[1]+prob[2])*100))))
  # rpt <- create_report(paste0("data/Result_",tax,"cut_", as.character((1 - (prob[2]-prob[1]))*100), ".pdf"), output_type = "PDF") %>% 
  #   add_content(tbl)
  # write_report(rpt)
  
  return(tot_cluster)
  
}

tot_cluster = data.frame()

plot_dendo(Cluster = tot_cluster,
           Dendo = Dendo, 
           tax = "ADMO", 
           site_taxon = site_taxon,
           prob = c(0.45, 0.55))

# test
# plot_dendo(site_taxon, Dendo, tax)

# No interest tax

No_interest_tax <- c("ADMI", "AFBA", "CBPY", "CLTL", "ENEE", "PLRC", "PSPO", "CROX", "GPPY", "CPLA", "COPL")

# map
# On recupere les taxons qui nous interesse, on enlève ceux vu avec les DREAL
# Pour chaque taxon avec profil proposé, on va comparer les sorties de la fonction
all_tax <- c(New_taxons %>% 
               filter(abre %notin% No_interest_tax) %>% 
               pull(abre), 
             New_taxons %>% 
               filter(AFNOR %notin% No_interest_tax) %>% 
               drop_na() %>% 
               pull(AFNOR)) %>% 
  unique() %>% 
  sort()
tax_CEUO <- "CEUO"

# Permettra de recuperer la sortie return de la fonction
all_cut <- map(all_tax, function(all_tax) plot_dendo(Cluster = tot_cluster,
                                                     Dendo = Dendo, 
                                                     tax = all_tax, 
                                                     site_taxon = site_taxon,
                                                     prob = c(0.45, 0.55)))
cut_CEUO <- plot_dendo(cluster = tot_cluster_CEUO,
                       Dendo = Dendo_CEUO, 
                       tax = tax_CEUO, 
                       site_taxon = site_taxon_CEUO,
                       prob = c(0.25, 0.75))

save(all_cut, file = "data/Donnees_utilisables/All_cut.RData")
save(cut_CEUO, file = "data/Donnees_utilisables/cut_CEUO.RData")


# Essais pour le shiny, ks tests, khi²... -------------------------------------

load("data/Donnees_utilisables/All_cut.RData")
load("data/Donnees_utilisables/cut_CEUO.RData")
load("data/Donnees_utilisables/Donnees_test_serveur.RData")
load("data/Donnees_utilisables/Donnees_completes.RData")
load("data/Donnees_utilisables/New_taxons.RData")
load("data/Donnees_utilisables/Cluster_2015_2021.RData")

`%notin%` <- Negate(`%in%`)

No_interest_tax <- c("ADMI", "AFBA", "CBPY", "CLTL", "ENEE", "PLRC", "PSPO", "CROX", "GPPY", "CPLA", "COPL")

# On recupere la liste des variablkes de chaque tibble pour nommer indépendemment 
# chaque element de la liste retournée par la fonction
flatten <- function(x){  
  islist <- sapply(x, class) %in% "list"
  r <- c(x[!islist], unlist(x[islist],recursive = F))
  if(!sum(islist))return(r)
  flatten(r)
}
out <- Map(colnames,flatten(all_cut))
# 
tax_order <- data.frame(out) %>% slice(1) %>% flatten_chr()

# Chaque element de la liste est extrait sous forme d'un dataframe unique
for (i in 1:length(all_cut)) {
  assign(sub(".*_", "",tax_order[i]), as.data.frame(all_cut[[i]]))
}

# Maintenant qu'on a extrait sous forme de tableau individuels, on peut combiner
# les paires de taxons (proposition DREAL et taxon nouveau)

# creation des paires 
paires <- New_taxons %>% 
  filter(abre %notin% No_interest_tax) %>% 
  filter(AFNOR %notin% c("ADMI", "CPLA", "COPL")) %>% 
  drop_na()

# Recuperation des tableaux recap 
paire1 <- c(table(ADMC$cluster),table(ADSA$cluster))
paire2 <- c(table(ADSK$cluster),table(ADSH$cluster))
paire3 <- c(table(ADTC$cluster),table(ADCT$cluster))
paire4 <- c(table(AZHA$cluster),table(APFI$cluster))
paire5 <- c(table(CEXF$cluster),table(CAFF$cluster))
paire6 <- c(table(CSBH$cluster),table(CHEL$cluster))
paire7 <- c(table(EARB$cluster),table(EARC$cluster))
paire8 <- c(table(ECAL$cluster),table(ECPM$cluster))
paire9 <- c(table(ECTO$cluster),table(EBLU$cluster))
paire10 <- c(table(EJUE$cluster),table(EBLU$cluster))
paire11 <- c(table(EPSG$cluster),table(EMIN$cluster))
paire12 <- c(table(FPDE$cluster),table(UDEL$cluster))
paire13 <- c(table(GMIS$cluster),table(GELG$cluster))
paire14 <- c(table(HPDA$cluster),table(NGRE$cluster))
paire15 <- c(table(LHLU$cluster),table(LGOP$cluster))
paire16 <- c(table(NSTS$cluster),table(NINC$cluster))
paire17 <- c(table(NFSO$cluster),table(SSVE$cluster))
paire18 <- c(table(PSXO$cluster),table(PHEL$cluster))
paire19 <- c(table(SBOS$cluster),table(SDIF$cluster))
paire20 <- c(table(SCRA$cluster),table(SNIG$cluster))
paire21 <- c(table(SECA$cluster),table(SPUP$cluster))
# Récupérer la trophie des sites des clusters de chaque paire

# On recupere les classes trophiques des taxons qui nous serviront a posteriori
Trophie <- read_excel(path = "data/Donnees_utilisables/Classes_trophie.xlsx", 
                      sheet = "Numbers") %>% filter(parameter == "PO4") %>% 
  dplyr::rename("code_taxon" = "code") %>% 
  mutate(Class = as.character(Class)) %>% 
  dplyr::select(code_taxon, Class)

# Recuperation des donnees de bases ayant servies au travaux
tab_data_base <- tibble(Flore %>% ungroup() %>% distinct() %>% 
                          filter(year(DATE) >= 2015)) %>% 
  left_join(Trophie, by = "code_taxon") %>% 
  mutate(Location = paste0(CODE_STATION,"_", year(DATE), "_0", lubridate::month(DATE), "_0", lubridate::day(DATE)))

# map(all_cut, ~ (pivot_wider(.x, names_from=1,values_from= 2)))

# Fonction pour changer le nom de la premiere colonne destableau pour pouvoir les 
# fusionner, recuperer les sites, realiser les boxplots.
tax_sibling <- function(tab1, tab2){
  new_name1 = stringr::str_sub(names(tab1[1]),start = 1,end = 7)
  title = stringr::str_sub(names(tab1[1]),start = 9,end = 13)
  tab1 = tab1 %>% 
    mutate(taxon = stringr::str_sub(names(tab1[1]),start = 9,end = 13)) %>% 
    rename_with(~c(new_name1), c(names(tab1[1])))
  nsite1 = length(unique(tab1$site_coupe))
  # Extraction des tout les taxons presents sur les memes sites
  tab1_site = tab_data_base %>% filter(Location %in% unique(tab1 %>% 
                                                              pull(site_coupe))) %>% 
    dplyr::select(Location, CODE_STATION, code_taxon, RESULTAT, Class)
  
  plot1 = tab1_site %>% 
    mutate_if(is.character, ~replace(., is.na(.), 0)) %>% 
    group_by(Class) %>% 
    dplyr::select(Class, RESULTAT) %>% 
    ggplot(aes(x = as.numeric(Class), y = log(RESULTAT), fill = Class))+
    geom_boxplot()+
    theme_bw()+
    labs(title = paste0(title, " Nombre de sites: ", nsite1),
         x = "", y = "Abondance (log)")+
    theme(legend.position = "none")
  
  new_name2 = stringr::str_sub(names(tab2[1]),start = 1,end = 7)
  title2 = stringr::str_sub(names(tab2[1]),start = 9,end = 13)
  tab2 = tab2 %>% 
    mutate(taxon = stringr::str_sub(names(tab2[1]),start = 9,end = 13)) %>% 
    rename_with(~c(new_name2), c(names(tab2[1])))
  nsite2 = length(unique(tab2$site_coupe))
  tab2_site = tab_data_base %>% filter(Location %in% unique(tab2 %>% 
                                                              pull(site_coupe))) %>% 
    dplyr::select(Location, CODE_STATION, code_taxon,RESULTAT, Class)
  
  plot2 = tab2_site %>% 
    mutate_if(is.character, ~replace(., is.na(.), 0)) %>% 
    group_by(Class) %>% 
    dplyr::select(Class, RESULTAT) %>% 
    ggplot(aes(x = as.numeric(Class), y = log(RESULTAT), fill = Class))+
    geom_boxplot()+
    theme_bw()+
    labs(title = paste0(title2, " Nombre de sites: ", nsite2),
         x = "", y = "Abondance (log)")+
    theme(legend.position = "none")
  
  plot3 = tab1_site %>% 
    mutate_if(is.character, ~replace(., is.na(.), 0)) %>% 
    group_by(Class) %>% 
    dplyr::select(Class) %>% 
    ggplot(aes(x = as.numeric(Class), fill = Class))+
    geom_bar(color = "black")+
    theme_bw()+
    labs(x = "Trophie", y = "Nombre de taxon")+
    theme(legend.position = "bottom")
  
  plot4 = tab2_site %>% 
    mutate_if(is.character, ~replace(., is.na(.), 0)) %>% 
    group_by(Class) %>% 
    dplyr::select(Class) %>% 
    ggplot(aes(x = as.numeric(Class), fill = Class))+
    geom_bar(color = "black")+
    theme_bw()+
    labs(x = "Trophie", y = "Nombre de taxon")+
    theme(legend.position = "none")
  
  pval_chi = NULL
  
  if(length(unique(tab1_site$Location)) > length(unique(tab2_site$Location))){
    
    for(i in 1:500){
      ech1 = tab1_site %>% 
        group_by(Location) %>% 
        nest() %>%
        ungroup() %>% 
        sample_n(length(unique(tab2_site$Location))) %>% 
        unnest()
      
      tab_ech1 = data.frame(table(ech1$Class))
      tab_ech2 = data.frame(table(tab2_site$Class))
      
      if(length(tab_ech1) > length(tab_ech2)){
        combine <- tab_ech1 %>% left_join(tab_ech2,by = "Var1")
      }else{combine <- tab_ech2 %>% left_join(tab_ech1,by = "Var1")}
      combine[is.na(combine)] <- 0
      
      Matrix <- matrix(c(combine$Freq.x,combine$Freq.y), 
                       nrow = 2, 
                       ncol = length(combine$Var1),
                       byrow = TRUE)
      
      X2 = chisq.test(Matrix)$p.value
      pval_chi = c(pval_chi,X2)}
  }else{
    for(i in 1:500){
      ech1 = tab2_site %>% 
        group_by(Location) %>% 
        nest() %>%
        ungroup() %>% 
        sample_n(length(unique(tab1_site$Location))) %>% 
        unnest()
      
      tab_ech1 = data.frame(table(ech1$Class))
      tab_ech2 = data.frame(table(tab1_site$Class))
      
      if(length(tab_ech1) > length(tab_ech2)){
        combine <- tab_ech1 %>% left_join(tab_ech2,by = "Var1")
      }else{combine <- tab_ech2 %>% left_join(tab_ech1,by = "Var1")}
      combine[is.na(combine)] <- 0
      
      Matrix <- matrix(c(combine$Freq.x,combine$Freq.y), 
                       nrow = 2, 
                       ncol = length(combine$Var1),
                       byrow = TRUE)
      
      X2 = chisq.test(Matrix)$p.value
      pval_chi = c(pval_chi,X2)}
  }
  
  conf_chi<-hBayesDM::HDIofMCMC(pval_chi)
  
  plot_chi = data.frame(pval_chi) %>% 
    ggplot(aes(x = pval_chi))+
    geom_histogram(color = "black", bins = 50)+
    geom_segment(x=conf_chi[1],xend=conf_chi[2],y=0,yend=0,
                 color="red",size=1,
                 lineend="round")+
    labs(title= "Test du X²", y = "", x = "p_value")+
    theme_bw()
  
  
  # performing the K-S
  # Test on x and x2
  x = data.frame(res = tab1_site$RESULTAT)
  x2 = data.frame(res =tab2_site$RESULTAT)
  
  Ks = ks.test(x$res, x2$res, alternative = "l")$p.value
  
  plot_Ks = ggplot(x, aes(res)) + stat_ecdf(geom = "step", col = "red")+
    stat_ecdf(data = x2, geom = "step", col = "blue")+
    theme_bw()+
    labs(title = paste0("Kolmogorov test: ", "p = ", round(Ks,5)),
         x = "", y = "")
  
  
  Chimie_tab <- Chimie %>% filter(!is.na(Mediane)) %>% 
    mutate(Location = paste0(CODE_STATION,"_", 
                             year(DATE), "_0", 
                             lubridate::month(DATE), "_0", 
                             lubridate::day(DATE))) %>% 
    mutate(Nom_parametre = case_when(Code_parametre == 1433 ~ "PO4")) %>%
    select(Location, Nom_parametre, Code_parametre, Mediane) %>% 
    filter(Nom_parametre == "PO4")
  
  plot5_tab <- Chimie_tab %>% 
    filter(Location %in% unique(tab1_site %>% pull(Location))) %>% 
    mutate(group = "1", taxon = tab1$taxon[1])
  
  plot5_tab2 <- Chimie_tab %>% 
    filter(Location %in% unique(tab2_site %>% pull(Location))) %>% 
    mutate(group = "2", taxon = tab2$taxon[1])
  
  if(length(plot5_tab2$Location)>5000){
    plot5_tab2 = sample_n(plot5_tab2, 5000)
  }
  
  # On construit les graphs et on applique des tests de comparaison de moyenne
  
  if(length(plot5_tab$Location) < 30 || length(plot5_tab2$Location) < 30){
    
    l1 = length(plot5_tab$Location)
    l2 = length(plot5_tab2$Location)
    
    if(l1<30){
      posy = "SW"
    }else{posy = "SE"}
    
    if(l1<30 & l2 < 30){
      posy = "N"}
    
    plot5 = bind_rows(plot5_tab, plot5_tab2) %>% 
      ggplot(aes(x = log(Mediane), fill = group))+
      geom_boxplot()+
      coord_flip()+
      theme_bw()+
      theme(legend.position = "none")+
      labs(x = "log Concentration", y = "", title = paste0("Concentration en PO4 dans les sites")) + 
      annotation_compass('Nb sites insuffisant',posy)+
      annotation_compass(plot5_tab$taxon[1],'NW')+
      annotation_compass(plot5_tab2$taxon[1],'NE')+
      geom_hline(yintercept = 0, 
                 color = "black", size=0.75)
    
    (plot1+plot2+plot_Ks+plot_layout(tag_level = 'new', ncol = 3))/
      (plot3+plot4+plot_chi+plot_layout(tag_level = 'new', ncol = 3))/plot5 +
      plot_annotation(tag_levels = c('A', '1'), 
                      title ='Comparaison des abondances taxonomiques (en log) (A), 
                  du nombre de taxons (B) par classe de trophie et de la quantité moyenne en ortophosphates (C)
                  sur les sites extraits de la CAH pour chacun des taxons',
                      theme = theme(plot.title = element_text(hjust = 0.5)),
                      caption = 'made by Leonard Heinry')
    
  }else{
    var_test = bartlett.test(Mediane ~ group, data = bind_rows(plot5_tab, plot5_tab2))$p.value
    Shpiro1 = shapiro.test(plot5_tab$Mediane)$p.value
    Shpiro2 = shapiro.test(plot5_tab2$Mediane)$p.value
    
    if(var_test & Shpiro1 & Shpiro2 >= 0.05){
      all_pvalue = NULL
      if(length(plot5_tab$Mediane) > length(plot5_tab2$Mediane)){
        for(i in 1:1000){
          samp1 = sample(plot5_tab$Mediane, length(plot5_tab2$Mediane), replace = FALSE)
          test = t.test(samp1, plot5_tab2$Mediane)$p.value
          all_pvalue = c(all_pvalue, test)}
      }
      else{
        for(i in 1:1000){
          samp1 = sample(plot5_tab2$Mediane, length(plot5_tab$Mediane), replace = FALSE)
          test = t.test(samp1, plot5_tab$Mediane)$p.value
          all_pvalue = c(all_pvalue, test)
        }
      }
      plot5 = bind_rows(plot5_tab, plot5_tab2) %>% 
        ggplot(aes(x = log(Mediane), fill = group))+
        geom_boxplot()+
        coord_flip()+
        theme_bw()+
        theme(legend.position = "none")+
        labs(x = "log Concentration", y = "", title = "Concentration de PO4 dans les sites") +
        annotation_compass(plot5_tab$taxon[1],'NW')+
        annotation_compass(plot5_tab2$taxon[1],'NE')+
        geom_hline(yintercept = 0,
                   color = "black", size=0.75)
      
      # Interval de confiance 95% automatique
      conf<-hBayesDM::HDIofMCMC(all_pvalue)
      plot6 = data.frame(all_pvalue) %>% 
        ggplot(aes(x = all_pvalue))+
        geom_histogram(color = "black", fill = "grey", bins = 100)+
        geom_segment(x=conf[1],xend=conf[2],y=0,yend=0,
                     color="red",size=1,
                     lineend="round")+
        labs(title = "Boostrap de la p_value pour 1000 répétitions avec interval de confiance à 95% (ligne rouge)")+
        theme_bw()
    }
    
    else{
      all_pvalue = NULL
      if(length(plot5_tab$Mediane) > length(plot5_tab2$Mediane)){
        for(i in 1:1000){
          samp1 = sample(plot5_tab$Mediane, length(plot5_tab2$Mediane), replace = FALSE)
          test = wilcox.test(samp1, plot5_tab2$Mediane)$p.value
          all_pvalue = c(all_pvalue, test)
        }
      }
      else{
        for(i in 1:1000){
          samp1 = sample(plot5_tab2$Mediane, length(plot5_tab$Mediane), replace = FALSE)
          test = wilcox.test(samp1, plot5_tab$Mediane)$p.value
          all_pvalue = c(all_pvalue, test)
        }
      }
      
      plot5 = bind_rows(plot5_tab, plot5_tab2) %>% 
        ggplot(aes(x = log(Mediane), fill = group))+
        geom_boxplot()+
        coord_flip()+
        theme_bw()+
        theme(legend.position = "none")+
        labs(x = "log Concentration", y = "", title = "Concentration de PO4 dans les sites") +
        annotation_compass(plot5_tab$taxon[1],'NW')+
        annotation_compass(plot5_tab2$taxon[1],'NE')+
        geom_hline(yintercept = 0, 
                   color = "black", size=0.75)
      
      # Interval de confiance 95% automatique
      conf<-hBayesDM::HDIofMCMC(all_pvalue)
      plot6 = data.frame(all_pvalue) %>% 
        ggplot(aes(x = all_pvalue))+
        geom_histogram(color = "black", fill = "grey", bins = 1000)+
        geom_segment(x=conf[1],xend=conf[2],y=0,yend=0,
                     color="red",size=1,
                     lineend="round")+
        ggtitle(paste0(round(conf,8)[1], " < ", "pvalue", " < ",round(conf,8)[2]))+
        theme_bw()
    }
    (plot1+plot2+plot_Ks+plot_layout(tag_level = 'new'))/
      (plot3+plot4+plot_chi+plot_layout(tag_level = 'new'))/(plot5 + plot6 + plot_layout(tag_level = 'new')) +
      plot_annotation(tag_levels = c('A', '1'), 
                      title ='Comparaison des abondances taxonomiques (en log) (A), 
                  du nombre de taxons (B) par classe de trophie et de la quantité moyenne 
                  en ortophosphates (C) sur les sites extraits de la CAH pour chacun des taxons',
                      theme = theme(plot.title = element_text(hjust = 0.5)),
                      caption = 'made by Leonard Heinry')
  }
  
}
tax_sibling(ADSK, EJUE)

tax_solo <- function(tab){
  new_name = stringr::str_sub(names(tab[1]),start = 1,end = 7)
  tab = tab %>% 
    mutate(taxon = stringr::str_sub(names(tab[1]),start = 9,end = 13)) %>% 
    rename_with(~c(new_name), c(names(tab[1])))
  return(tab)
}


# On recupere le tableau global de tous les taxons et leur 25% d'abondance et on
# le sauvegarde pour construire le tableau recap et y inclure les graphiques,
# on met en face d'un taxon son profil proposé
tab_recap <- bind_rows(lapply(all_cut, FUN = tax_solo)) %>% 
  filter(taxon %in% unique(New_taxons %>% 
                             filter(abre %notin% No_interest_tax) %>% 
                             pull(abre))) %>% 
  group_by(taxon) %>% 
  dplyr::summarise(nsites_25 = n()) %>% 
  left_join(New_taxons %>% rename(taxon = abre), by = "taxon") %>% 
  ungroup() %>% left_join(bind_rows(lapply(all_cut, FUN = tax_solo)) %>% 
                            filter(taxon %in% unique(New_taxons %>% 
                                                       filter(abre %notin% No_interest_tax) %>% 
                                                       pull(AFNOR))) %>% 
                            group_by(taxon) %>% 
                            dplyr::summarise(nsites_25_paired = n()) %>% 
                            rename(AFNOR = taxon), by = "AFNOR")

write_xlsx(tab_recap, "data/essaie.xlsx")

# On fait tourner la fonction pour chaque paire et on sauvegarde le resultat 
# pour le tableau recap --> AJOUTER UNE COLONNE COMPARAISON TROPHIE AVEC LES PLOTS

png_saver <- function(tax1,tax2){
  png(file = paste("data/Donnees_utilisables/",stringr::str_sub(names(tax1[1]),start = 9,end = 13),"_", 
                   stringr::str_sub(names(tax2[1]),start = 9,end = 13),".jpeg"),
      width = 600, height = 600)
  plot(tax_sibling(tax1,tax2))
  dev.off()
}

png_saver(ADMC,ADSA)
png_saver(ADSK,ADSH)
png_saver(ADTC,ADCT)
png_saver(AZHA,APFI)
png_saver(CEXF,CAFF)
png_saver(CSBH,CHEL)
png_saver(EARB,EARC)
png_saver(ECAL,ECPM)
png_saver(ECTO,EBLU)
png_saver(EJUE,EBLU)
png_saver(EPSG,EMIN)
png_saver(FPDE,UDEL)
png_saver(GMIS,GELG)
png_saver(HPDA,NGRE)
png_saver(LHLU,LGOP)
png_saver(NSTS,NINC)
png_saver(NFSO,SSVE)
png_saver(PSXO,PHEL)
png_saver(SBOS,SDIF)
png_saver(SCRA,SNIG)
png_saver(SECA,SPUP)

# Traitement CEUO

cut_CEUO

# A partir des images et du tableau excel créé, construire le tableau recap. 
# y ajouter toutes les informations necessaires, y inclure CPLA traité à part 

# Traiter CEUO à part !
#Recuperation des autres taxons presents sur les sites de CEUO entre 2007 et 2009

tab_CEUO_base <- tibble(Flore %>% ungroup() %>% distinct() %>% 
                          filter(year(DATE) <= 2009)) %>% 
  left_join(Trophie, by = "code_taxon") %>% 
  mutate(Location = paste0(CODE_STATION,"_", year(DATE), "_0", lubridate::month(DATE), "_0", lubridate::day(DATE)))

CEUO_site <- tab_CEUO_base %>% filter(Location %in% unique(cut_CEUO %>% 
                                                             pull(site_coupe_CEUO))) %>% 
  dplyr::select(Location, CODE_STATION, code_taxon, RESULTAT, Class)

nsite_CEUO <- length(unique(cut_CEUO$site_coupe_CEUO))

plot_CEUO = CEUO_site %>% 
  mutate_if(is.character, ~replace(., is.na(.), 0)) %>% 
  group_by(Class) %>% 
  dplyr::select(Class, RESULTAT) %>% 
  ggplot(aes(x = as.numeric(Class), y = log(RESULTAT), fill = Class))+
  geom_boxplot()+
  theme_bw()+
  labs(title = paste0("CEUO", " Nombre de sites: ", nsite_CEUO),
       x = "Trophie", y = "Abondance (log)")

plot_CEUO2 = CEUO_site %>% 
  mutate_if(is.character, ~replace(., is.na(.), 0)) %>% 
  group_by(Class) %>% 
  dplyr::select(Class) %>% 
  ggplot(aes(x = as.numeric(Class), fill = Class))+
  geom_bar(color = "black")+
  theme_bw()+
  labs(x = "Trophie", y = "Nombre de taxon")

png(file = paste("data/Donnees_utilisables/CEUO.jpeg"))
plot(plot_CEUO+plot_CEUO2)
dev.off()

Chim_plot <- Chimie %>% filter(!is.na(Mediane)) %>% 
  mutate(Location = paste0(CODE_STATION,"_", 
                           year(DATE), "_0", 
                           lubridate::month(DATE), "_0", 
                           lubridate::day(DATE))) %>% 
  mutate(Nom_parametre = case_when(Code_parametre == 1433 ~ "PO4")) %>%
  select(Location, Nom_parametre, Code_parametre, Mediane) %>% 
  filter(Nom_parametre == "PO4") %>% 
  filter(Location %in% unique(CEUO_site %>% pull(Location))) %>% 
  ggplot(aes(x = log(Mediane)))+
  geom_boxplot()+
  coord_flip()+
  theme_bw()+
  labs(x = "Concentration", title = "Concentration de PO4 dans les sites")




# INDVAL COMPUTATION ------------------------------------------------------

library(tidyverse)
library(labdsv)
library(vegan)
library(indicspecies)
library(parallel)
library(doParallel)
library(parallelDist)
library(dendextend)

load("data/Donnees_utilisables/site_taxon.RData")
load("data/Donnees_utilisables/New_taxons.RData")
`%notin%` <- Negate(`%in%`)


cluster <- parallel::makePSOCKcluster(4, outfile="") 
doParallel::registerDoParallel(cluster)
clusterExport(cluster, "site_taxon")

Giga_plot <- function(tax){
  
  select_tax <- site_taxon %>% filter(taxon == tax)
  
  Tax_mat <- data.frame(site_taxon %>% 
                          filter(site %in% select_tax$site) %>%
                          dplyr::select(site, taxon, Ab_relative) %>%
                          pivot_wider(names_from = taxon, values_from = Ab_relative, values_fill = 0)) %>%
    column_to_rownames(var = "site")
  
  Stand_Tax_mat <- decostand(Tax_mat, method = "hellinger")
  dist_bc <- parDist(as.matrix(Tax_mat), 
                     method="bray", 
                     threads = 20)
  
  clus_ward <- hclust(dist_bc, method="ward.D2")
  
  Dendo <- clus_ward %>% as.dendrogram()
  
  nodes <- sort(get_nodes_attr(Dendo, "height",
                               include_leaves = FALSE,
                               include_branches = TRUE,
                               na.rm = TRUE), decreasing = TRUE)
  
  Dendo_cut <- data.frame(cutree(Dendo, h = nodes))
  
  vecteur <- 1:ncol(Dendo_cut)
  
  N <- vecteur[seq(2, length(vecteur), by = 1)]
  
  pb <- txtProgressBar(min=0, max=N, style=3)
  
  list_indval <- foreach(
    i = seq_along(N),
    .packages = c("tidyverse", "indicspecies", "labdsv")
  ) %dopar% {
    
    gc()
    
    valeur = N[i]
    setTxtProgressBar(pb, i)
    
    indicator_species <- multipatt(Tax_mat, as.vector(Dendo_cut[,valeur]), duleg = TRUE,
                                   control = how(nperm=999))
    fidg <- indicator_species$sign %>%
      rownames_to_column(var = "esp") %>%
      rename(group = index, indval = stat) %>%
      mutate(group = as.character(group)) %>%
      dplyr::select(esp,group,indval,p.value) %>% filter(p.value <= 0.05) %>%
      left_join(data.frame(indicator_species$A) %>%
                  rownames_to_column(var = "esp") %>%
                  pivot_longer(!esp, values_to = "A", names_to = "group") %>%
                  mutate(group = gsub("^X","",group)), by = c("esp", "group")) %>%
      left_join(data.frame(indicator_species$B) %>%
                  rownames_to_column(var = "esp") %>%
                  pivot_longer(!esp, values_to = "B", names_to = "group") %>%
                  mutate(group = gsub("^X","",group)), by = c("esp", "group")) %>%
      left_join(site_taxon %>%
                  filter(site %in% select_tax$site) %>%
                  dplyr::select(site, taxon, Ab_relative) %>%
                  group_by(taxon) %>%
                  summarize(freq = n()) %>%
                  rename(esp = taxon), by = "esp")
    
    # indicator_species <- indval(Tax_mat, as.vector(Dendo_cut[,i]))
    # gr <- indicator_species$maxcls[indicator_species$pval <= 0.05]
    # iv <- indicator_species$indcls[indicator_species$pval <= 0.05]
    # pv <- indicator_species$pval[indicator_species$pval <= 0.05]
    # fr <- apply(Tax_mat > 0, 2, sum)[indicator_species$pval <= 0.05]
    # fidg2 <- data.frame(group = gr, indval = iv, pvalue = pv, freq = fr)
    # fidg2 <- fidg2[order(fidg2$group,-fidg2$freq, -fidg2$indval),]
    # fidg2 <- fidg2 %>% rownames_to_column(var = "esp")
    # 
    # data.frame(fidg2)
    
  }
  
  # Récupération des coupes dans lesquelles le taxon est INDVAL
  interest_list <- lapply(list_indval, function(x) {
    if (tax %in% x$esp) {
      return(x)
    }
  })
  
  # enlever les éléments nuls de la liste
  interest_list <- interest_list[sapply(interest_list, function(x) !is.null(x))] 
  
  # Ordre des groupes dans lesquels le taxon est INDVAL
  grp_order <- lapply(interest_list, function(x) filter(x, esp == tax) %>% pull(group)) %>%
    unlist()
  
  
  if(length(grp_order) == 0){
    ma_liste_finale <- NULL}else{
      
      # Récupération pour des groupes dans lesquels le taxon est vu en INDVAL
      
      list_group_interest <- list()
      
      for(i in 1:length(grp_order)){
        list_group_interest[[i]] <- interest_list[[i]] %>% filter(group == grp_order[i])
      }
      
      # On retire les doublons de colonnes esp pou éviter une surcharge de plot à la fin
      duplicates <- duplicated(lapply(list_group_interest, "[[", "esp"))
      # Obtenir les index des doublons
      duplicates_index <- which(duplicates)
      
      if(length(duplicates_index) == 0){
        list_group_interest <- list_group_interest
      }else{list_group_interest <- list_group_interest[-duplicates_index]}# Supprimer les doublons de la liste
      
      # Spécification des colonnes pour la jointure
      col_join <- "esp"
      
      # Fonction de jointure des profils écologique
      ma_jointure <- function(x) {
        left_join(x, Transcoded_profiles, by = col_join)
      }
      
      # Appliquer la fonction de jointure à chaque élément de la liste
      ma_liste_jointe <- lapply(list_group_interest, ma_jointure)
      
      # On enleve les endroits ou le taxon d'interet est seul indval
      ma_liste_finale <- Filter(function(df) nrow(df) > 1, ma_liste_jointe)
    }
  # save(ma_liste_finale, file = paste0("data/Donnees_utilisables/indval",tax,".RData"))
  return(ma_liste_finale)
}

taxons <- site_taxon %>% filter(taxon %in% unique(New_taxons %>% pull(abre))) %>%
  group_by(taxon) %>% summarise(count = n()) %>% arrange(count) %>% filter(count >= 30) %>%
  pull(taxon)

length(taxons)
sort(taxons[1:39])

tic()
Giga_plot("ADTC")
toc()


load("data/Donnees_utilisables/site_taxon_CEUO.RData")
view(site_taxon_CEUO %>% distinct(site))



# Traitement Sorties INDVAL du serveur de calcul  ----------------------------------------------------------------
IBD_profiles <- tibble(read.csv(file = "data/Donnees_utilisables/IBD_params.csv", 
                                sep = ";", dec = ",")) %>%
  dplyr::select(-DENOMINATION, -SANDRE, -Origine) %>%
  rename(abre = AFNOR)

transcodage <- read.csv(file = "data/Donnees_utilisables/table_transcodage.csv", 
                        sep = ";", dec = ",") %>% dplyr::select(-name, -name_valid)


Transcoded_profiles <- IBD_profiles %>% left_join(transcodage, by = "abre") %>% 
  left_join(IBD_profiles %>% rename(CodeValid = abre), by = "CodeValid") %>%
  dplyr::select(10:18) %>%
  rename(esp = CodeValid, CL1 = CL1.y,    CL2 = CL2.y,   CL3 = CL3.y,   CL4 = CL4.y,  
         CL5 = CL5.y,   CL6 = CL6.y,    CL7 = CL7.y, Val_ind = Val.Ind..y) %>%
  distinct(esp, .keep_all = TRUE) %>%
  drop_na()

Transcoded_profiles %>% arrange(esp) %>%
  anti_join(profiles %>% rename(esp = AFNOR), by = "esp")

all.equal(data.frame(Transcoded_profiles %>% 
            filter(esp %notin% c("CCTL", "CMLF", "GMUT", "HAMP", "SSMU")) %>%
                     rename(AFNOR = esp) %>% arrange(AFNOR)),
          profiles %>% select(-true_profile))

load("data/Donnees_utilisables/site_taxon.RData")
load("data/Donnees_utilisables/New_taxons.RData")

profiles <- read.csv(file = "data/Donnees_utilisables/IBD_params.csv", 
                     sep = ";", dec = ",") %>%
  select(-SANDRE, -DENOMINATION, -Origine) %>% 
  left_join(as_tibble(read.csv2("data/Donnees_utilisables/table_transcodage.csv", stringsAsFactors = FALSE)) %>%
                                                           select(AFNOR = "abre", True_name = "CodeValid"), by = "AFNOR") %>%
  
  mutate(true_profile = if_else(AFNOR == True_name, 1,0)) %>% 
  filter(true_profile == 1) %>%
  mutate(AFNOR = if_else(is.na(True_name) == T, AFNOR, True_name)) %>%
  select(-True_name) %>% filter(!is.na(AFNOR))


Candidats_Evaluation <- function(taxon){
  
  gc() 
  
  load(paste0("data/Donnees_utilisables/indval",taxon,".RData"))
  
  prof_prop <- ifelse(is.na(New_taxons %>% filter(abre == taxon) %>% pull(abre)) == TRUE, NA, New_taxons %>% filter(abre == taxon) %>% pull(AFNOR))
  
  if(is.na(prof_prop) == TRUE){
    prof_DREAL <- NA}else{prof_DREAL <- unique(profiles %>% filter(AFNOR == prof_prop)) %>% select(-AFNOR, -true_profile, -Val.Ind.)}
  
  if(length(prof_DREAL) == 1){
    
    candidats <- do.call("rbind", ma_liste_finale) %>% 
      filter(esp != taxon) %>%
      group_by(esp) %>%
      mutate(mean_indval = mean(indval, na.rm = TRUE), # Moyenne indval de l'espèce
             
             mean_A = mean(A, na.rm = TRUE), # Moyenne des abondances dans les relevés par rapport a tous les groupes (spécificité)
             
             mean_B = mean(B, na.rm = TRUE), # Moyenne fidélité de l'espèce
             med_freq = median(freq, na.rm = TRUE), # Nombre de coupes effectuées
             n_occ = n(), # Nombre de coupes effectuées
             n_coupes = length(ma_liste_finale)) %>% # Nombre de coupes effectuées
      distinct(esp, .keep_all = TRUE) %>%
      select(esp, mean_indval, mean_A, mean_B, med_freq, n_occ) %>%
      arrange(desc(med_freq)) %>%
      left_join(profiles %>% rename(esp = AFNOR) %>% distinct(), by = "esp") %>%
      drop_na()
    
    # save(candidats, file = paste0("data/Donnees_utilisables/Candidats_",taxon,".RData"))
    
  }else{
    candidats <- do.call("rbind", ma_liste_finale) %>% 
      filter(esp != taxon) %>%
      group_by(esp) %>%
      mutate(mean_indval = mean(indval, na.rm = TRUE), # Moyenne indval de l'espèce
             
             mean_A = mean(A, na.rm = TRUE), # Moyenne des abondances dans les relevés par rapport a tous les groupes (spécificité)
             
             mean_B = mean(B, na.rm = TRUE), # Moyenne fidélité de l'espèce
             med_freq = median(freq, na.rm = TRUE), # Nombre de coupes effectuées
             n_occ = n(), # Nombre de coupes effectuées
             n_coupes = length(ma_liste_finale)) %>% # Nombre de coupes effectuées
      distinct(esp, .keep_all = TRUE) %>%
      select(esp, mean_indval, mean_A, mean_B, med_freq, n_occ) %>%
      arrange(desc(med_freq)) %>%
      left_join(profiles %>% rename(esp = AFNOR) %>% distinct(), by = "esp") %>%
      drop_na() %>% 
      mutate(prof_dif = sum(abs(c(CL1,CL2,CL3,CL4,CL5,CL6,CL7) - prof_DREAL))) %>%
      select(-true_profile) %>% arrange(prof_dif)

   
  }
  
  return(candidats)
  
}

data.frame(Candidats_Evaluation("LHLU")) %>% arrange(esp)

load(paste0("data/Donnees_utilisables/indval","SSBG",".RData"))
ma_liste_finale

profiles %>% filter(AFNOR %in% unique(New_taxons %>% filter(abre == "LHLU") %>% pull(AFNOR)))

profiles %>% filter(AFNOR == "GADC")
IBD_profiles %>% filter(abre == "AHOF")
Prof_eco %>% filter(AFNOR == "FGRA")

site_taxon %>% filter(taxon == "CPLA")

taxons <- site_taxon %>% filter(taxon %in% unique(New_taxons %>% pull(abre))) %>%
  group_by(taxon) %>% summarise(count = n()) %>% arrange(count) %>% filter(count >= 30) %>%
  pull(taxon)

c(taxons[46:53],taxons[56:57])

profiles %>% filter(AFNOR == "SEAT")

New_taxons
