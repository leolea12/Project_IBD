
library(parallel)
library(parallelDist)
library(doSNOW)
library(NbClust)
library(tidyverse)
library(progress)
library(vegan)
library(vegan)
library(tictoc)
library(dendextend)
library(patchwork)
library(readxl)
library(distributions3)
library(gridExtra)
library(reporter)
library(writexl)
library(lubridate)
library(xlsx)
library(viridis)
library(labdsv)
library(hBayesDM)

load("data/Donnees_utilisables/Donnees_test_serveur_CEUO.RData")
load("data/Donnees_utilisables/Donnees_test_serveur.RData")
load("data/Donnees_utilisables/Donnees_completes.RData")
load("data/Donnees_utilisables/Cluster_2015_2021.RData")
load("data/Donnees_utilisables/New_taxons.RData")

# Fonction d'annotation de graphs (positionnement des légendes etc) 

annotation_compass <- function(label,
                               position = c('N','NE','E','SE','S','SW','W','NW'),
                               padding = grid::unit(c(0.5,0.5),"line"), ...){
  position <- match.arg(position)
  x <- switch (position,
               N = 0.5,
               NE = 1,
               E = 1,
               SE = 1,
               S = 0.5, 
               SW = 0,
               W = 0, 
               NW = 0
  )
  y <- switch (position,
               N = 1,
               NE = 1,
               E = 0.5,
               SE = 0,
               S = 0, 
               SW = 0,
               W = 0.5, 
               NW = 1
  )
  hjust <- switch (position,
                   N = 0.5,
                   NE = 1,
                   E = 1,
                   SE = 1,
                   S = 0.5, 
                   SW = 0,
                   W = 0, 
                   NW = 0
  )
  vjust <- switch (position,
                   N = 1,
                   NE = 1,
                   E = 0.5,
                   SE = 0,
                   S = 0, 
                   SW = 0,
                   W = 0.5, 
                   NW = 1
  )
  f1 <- switch (position,
                N = 0,
                NE = -1,
                E = -1,
                SE = -1,
                S = 0, 
                SW = 1,
                W = 1, 
                NW = 1
  )
  f2 <- switch (position,
                N = -1,
                NE = -1,
                E = 0,
                SE = 1,
                S = 1, 
                SW = 1,
                W = 0, 
                NW = -1
  )
  annotation_custom(grid::textGrob(label, 
                                   x=grid::unit(x,"npc") + f1*padding[1] , 
                                   y=grid::unit(y,"npc") + f2*padding[2],
                                   hjust=hjust,vjust=vjust, ...))
}

# Suppression des labels
noLabel <- function(x) {
  if (stats::is.leaf(x)) {
    attr(x, "label") <- NULL }
  return(x)
}

# Densité des points sur un plot
get_density <- function(x, y) {
  dens <- MASS::kde2d(x, y)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

# Rendre indépendant chaque element d'une liste 

flatten <- function(x){  
  islist <- sapply(x, class) %in% "list"
  r <- c(x[!islist], unlist(x[islist],recursive = F))
  if(!sum(islist))return(r)
  flatten(r)
}

# Creation d'objet nécessaires avant traitement et standardisation des données

colnames(NMDS_test) <- sub("_.*","", colnames(NMDS_test))

`%notin%` <- Negate(`%in%`)

# TRAVAIL SUR CAH DONNEES 2015-2021 (CEUO séparé) ---------------------

Dendo <- clus_ward %>% as.dendrogram()
save(Dendo, file = "data/Donnees_utilisables/Dendogramme.RData")
save(Dendo, file = "IBD_2022_shiny/Data/Dendogramme.RData")
Dendo_CEUO <- clus_ward_CEUO %>% as.dendrogram()
save(Dendo, file = "data/Donnees_utilisables/Dendogramme_CEUO.RData")
save(Dendo, file = "IBD_2022_shiny/Data/Dendogramme_CEUO.RData")

# On recupere pour chaque taxon ses sites de presence
site_taxon <- NMDS_test %>%
  rownames_to_column(var = "site") %>%
  pivot_longer(!site, names_to = "taxon", values_to = "Ab_relative") %>%
  mutate(taxon = sub("_.*", "", taxon)) %>% 
  filter(Ab_relative > 0)

save(site_taxon, file = "data/Donnees_utilisables/site_taxon.RData")
save(site_taxon, file = "IBD_2022_shiny/Data/site_taxon.RData")


site_taxon_CEUO <- NMDS_test_CEUO %>%
  rownames_to_column(var = "site") %>%
  pivot_longer(!site, names_to = "taxon", values_to = "Ab_relative") %>%
  mutate(taxon = sub("_.*", "", taxon)) %>% 
  filter(Ab_relative > 0)

save(site_taxon, file = "data/Donnees_utilisables/site_taxon_CEUO.RData")
save(site_taxon, file = "IBD_2022_shiny/Data/site_taxon_CEUO.RData")

# Initialisation du fichier qui contiendra le resume de tout les cluster
# des taxons

tot_cluster <- data.frame()
tot_cluster_CEUO <- data.frame()

# Fonction qui permettrae le plot des arbres par taxon, la recuperation des quantiles 
# la representation graphique et statisitque des resultats

plot_dendo <- function(Cluster, Dendo, tax, site_taxon, prob){
  select_tax <- site_taxon %>% filter(taxon == tax)
  no_select_tax <- site_taxon %>% filter(site %notin% unique(select_tax %>% pull(site))) %>% 
    distinct(site)
  
  label_select_tax <- data.frame(site = labels(Dendo)) %>% 
    left_join(select_tax, by = "site") %>% 
    mutate(site = if_else(taxon == tax, site, ""))
  
  plot_dend <- Dendo %>% 
    dendextend::set("labels_cex", 0.5) %>% 
    # set("nodes_pch", 19) %>% 
    # set("nodes_cex", 0.5) %>% 
    dendextend::set("by_labels_branches_col", value = select_tax$site) %>% 
    dendextend::set("by_labels_branches_lwd", value = select_tax$site, TF_values = c(1,0.25))  %>% 
    dendextend::set("by_labels_branches_lty", value = no_select_tax$site, TF_values = c(1,Inf)) %>%
    dendextend::hang.dendrogram(hang_height = 0.1) %>% 
    dendextend::set("labels_cex", 0.2) %>% 
    dendextend::set("labels", label_select_tax$site)
  
  # Extraction des quantiles en fonction de la densité
  X <- Categorical(label_select_tax %>% 
                     mutate(site = seq(1:length(label_select_tax$site)),
                            Ab_relative = if_else(is.na(Ab_relative) == TRUE, 0, Ab_relative)) %>% 
                     drop_na() %>% 
                     mutate(density = get_density(site, Ab_relative)) %>% 
                     drop_na() %>% pull(site),
                   p = label_select_tax %>% 
                     mutate(site = seq(1:length(label_select_tax$site)),
                            Ab_relative = if_else(is.na(Ab_relative) == TRUE, 0, Ab_relative)) %>% 
                     drop_na() %>% 
                     mutate(density = get_density(site, Ab_relative)) %>% 
                     drop_na() %>% pull(density))
  
  Quantil1 <- quantile(X, probs = c(0.45, 0.55))
  Quantil2 <- quantile(X, probs = c(0.375, 0.625))
  Quantil3 <- quantile(X, probs = c(0.25, 0.75))
  Quantil4 <- quantile(X, probs = c(0.125, 0.875))
  
  vlines <- data.frame(value = c(Quantil1, Quantil2, Quantil3, Quantil4),
                       Mean = c("10%","10%", "25%", "25%", "50%","50%", "75%", "75%"))
  
  plot_density <- label_select_tax %>% 
    mutate(site = seq(1:length(label_select_tax$site)),
           Ab_relative = if_else(is.na(Ab_relative) == TRUE, 0, Ab_relative)) %>% 
    drop_na() %>% 
    mutate(density = get_density(site, Ab_relative)) %>% 
    filter(taxon == tax) %>% 
    ggplot(aes(x=site)) +
    geom_density(fill="grey", color="grey", alpha=0.4)+
    geom_vline(data=vlines,
               aes(xintercept=value, colour=Mean),
               size=1, linetype="dashed", show.legend=TRUE)+
    theme_classic()+
    theme(axis.text.x=element_blank())
  
  # On créer des lignes quantiles entre lesquelles 50% des points seront présents
  
  Plot_ab <- label_select_tax %>% 
    mutate(site = seq(1:length(label_select_tax$site)),
           Ab_relative = if_else(is.na(Ab_relative) == TRUE, 0, Ab_relative)) %>% 
    drop_na() %>% 
    mutate(density = get_density(site, Ab_relative)) %>% 
    ggplot()+
    geom_point(aes(x = site, y = Ab_relative, color = density), alpha = 0.8)+
    scale_color_viridis()+
    theme_classic()+
    theme(axis.text.x=element_blank())
  
  # Combiner les deux images 
  # library(cowplot)
  # library(magick)
  # ggsave(Plot_ab, filename = 'data/Plot_ab.png', device = 'png', bg = 'transparent')
  # ggsave(plot_density, filename = 'data/plot_density.png', device = 'png', bg = 'transparent')
  # cowplot::plot_grid(Plot_ab, plot_density, nrow = 1)
  # plot1 <- image_read('data/Plot_ab.png')
  # plot2 <- image_read('data/plot_density.png')
  # img <- c(plot1, plot2)
  # png("data/PLOT.png")
  # image_mosaic(img)
  # dev.off()
  # image_write(image_mosaic(img), path = "data/final.png", format = "png")
  
  Combined_plot <- Plot_ab/plot_density
  
  # DETERMINER A QUEL ENDROIT ON DOIT COUPER
  
  if(prob[1] == 0.45 & prob[2] == 0.55){
    data <- data.frame(height = get_leaves_attr(plot_dend, "height")[Quantil1[1]:Quantil1[2]],
                       label = get_leaves_attr(plot_dend, "label")[Quantil1[1]:Quantil1[2]]) %>% 
      drop_na()}
  
  if(prob[1] == 0.375 & prob[2] == 0.625){
    data <- data.frame(height = get_leaves_attr(plot_dend, "height")[Quantil2[1]:Quantil2[2]],
                       label = get_leaves_attr(plot_dend, "label")[Quantil2[1]:Quantil2[2]]) %>% 
      drop_na()}
  
  if(prob[1] == 0.25 & prob[2] == 0.75){
    data <- data.frame(height = get_leaves_attr(plot_dend, "height")[Quantil3[1]:Quantil3[2]],
                       label = get_leaves_attr(plot_dend, "label")[Quantil3[1]:Quantil3[2]]) %>% 
      drop_na()}
  
  if(prob[1] == 0.125 & prob[2] == 0.875){
    data <- data.frame(height = get_leaves_attr(plot_dend, "height")[Quantil4[1]:Quantil4[2]],
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
    rename(cluster = "cutree.Dendo..h...Point_coupe.") %>% 
    filter(site_coupe %in% unique(data %>% pull(label))) %>% 
    arrange(desc(cluster)) %>% 
    rename_with(~ paste0("site_coupe_",tax), 'site_coupe')
  
  tot_cluster <- rbind(Cluster, cluster)
  
  # Pourcentage de sites en commun entre un taxon et son profil proposé
  
  
  # label_point_coupe <- data.frame(height = get_leaves_attr(plot_dend, "height")[Quantil4[1]:Quantil4[2]],
  #                                 label = get_leaves_attr(plot_dend, "label")[Quantil4[1]:Quantil4[2]]) %>% 
  #   drop_na()
  # 
  # # Parmis tous les site ou le taxon est présent, on prend celui qui se situe le plus haut
  # # ce sera le point de coupe de base !
  # point_coupe_base <- data.frame(height = get_leaves_attr(plot_dend, "height")[Quantil4[1]:Quantil4[2]],
  #                                label = get_leaves_attr(plot_dend, "label")[Quantil4[1]:Quantil4[2]]) %>% 
  #   drop_na() %>% 
  #   slice_max(height) %>% 
  #   mutate(height = round(height,2))
  
  # save plot 
  # grDevices::pdf(file = paste0("data/CAH_",tax,".pdf"))
  grDevices::jpeg(file = paste0("data/density_",tax,".jpeg"),width = 800, height = 800)
  # read_image(...)
  print(Combined_plot)
  grDevices::dev.off()
  
  grDevices::jpeg(file = paste0("data/CAH_",tax,".jpeg"), width = 800, height = 800)
  plot(plot_dend, main = paste0("CAH_",tax))
  # Ajout des informations sur le taxon sur la CAH (plus haute feuille + bornes)
  
  abline(v = Quantil1[1], col="red", lwd=1.5, lty=2)
  abline(v = Quantil1[2], col="red", lwd=1.5, lty=2)
  abline(v = Quantil2[1], col="green", lwd=1.5, lty=2)
  abline(v = Quantil2[2], col="green", lwd=1.5, lty=2)
  abline(v = Quantil3[1], col="blue", lwd=1.5, lty=2)
  abline(v = Quantil3[2], col="blue", lwd=1.5, lty=2)
  abline(v = Quantil4[1], col="purple", lwd=1.5, lty=2)
  abline(v = Quantil4[2], col="purple", lwd=1.5, lty=2)
  abline(h = Point_coupe, col="black", lwd=1.5, lty=1)
  legend(20,30.5, legend=c("10%","25%", "50%", "75%", paste0("COUPE_", as.character((1 - (prob[2]-prob[1]))*100))),
         col=c("red", "green", "blue", "purple", "black"),title = "partition", lty=2, cex=0.8)
  
  # Boucle qui test tous les points de coupe pour avoir 80% des sites dans le même cluster
  # A VOIR AVEC SEBASTIEN
  # for(i in seq(point_coupe_base$height,round(max(get_branches_heights(plot_dend),1)),1)){
  #   n_clust <- length(unique(cutree(plot_dend, h = i)))
  #   val_group <- length(table(as.vector(cutree(plot_dend, h = i)[label_point_coupe$label])))
  #   
  #   if(val_group == 1){
  #     seuil = i
  #     break
  #   }
  #   print(seuil)
  # }
  
  grDevices::dev.off()
  
  
  tbl <- create_table(cluster) %>% 
    titles(paste0("Result_",tax,"cut_", as.character(1 - ((prob[1]+prob[2])*100))))
  rpt <- create_report(paste0("data/Result_",tax,"cut_", as.character((1 - (prob[2]-prob[1]))*100), ".pdf"), output_type = "PDF") %>% 
    add_content(tbl)
  write_report(rpt)
  
  return(tot_cluster)
  
}

# Taxons non retenus pour le traitement des données

No_interest_tax <- c("ADMI", "AFBA", "CBPY", "CLTL", "ENEE", "PLRC", "PSPO", "CROX", "GPPY", "CPLA", "COPL")


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

# Permettra de recuperer la sortie return de la fonction plot_dendo

all_cut <- map(all_tax, function(all_tax) plot_dendo(Cluster = tot_cluster,
                                                     Dendo = Dendo, 
                                                     tax = all_tax, 
                                                     site_taxon = site_taxon,
                                                     prob = c(0.25, 0.75)))
cut_CEUO <- plot_dendo(cluster = tot_cluster_CEUO,
                       Dendo = Dendo_CEUO, 
                       tax = tax_CEUO, 
                       site_taxon = site_taxon_CEUO,
                       prob = c(0.25, 0.75))

save(all_cut, file = "data/Donnees_utilisables/All_cut.RData")
save(cut_CEUO, file = "data/Donnees_utilisables/cut_CEUO.RData")

# Construire un tableau recapitulatif -------------------------------------

load("data/Donnees_utilisables/All_cut.RData")
load("data/Donnees_utilisables/cut_CEUO.RData")
load("data/Donnees_utilisables/Donnees_test_serveur.RData")
load("data/Donnees_utilisables/Donnees_completes.RData")
load("data/Donnees_utilisables/New_taxons.RData")
load("data/Donnees_utilisables/Cluster_2015_2021.RData")

No_interest_tax <- c("ADMI", "AFBA", "CBPY", "CLTL", "ENEE", "PLRC", "PSPO", "CROX", "GPPY", "CPLA", "COPL")

out <- Map(colnames,flatten(all_cut))
 
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

# Fonction pour changer le nom de la premiere colonne des tableaux pour pouvoir les 
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
    
    (plot1+plot2+plot_layout(tag_level = 'new'))/(plot3+plot4+plot_layout(tag_level = 'new'))/plot5 +
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
    (plot1+plot2+plot_layout(tag_level = 'new'))/(plot3+plot4+plot_layout(tag_level = 'new'))/(plot5 + plot6 + plot_layout(tag_level = 'new')) +
      plot_annotation(tag_levels = c('A', '1'), 
                      title ='Comparaison des abondances taxonomiques (en log) (A), 
                  du nombre de taxons (B) par classe de trophie et de la quantité moyenne 
                  en ortophosphates (C) sur les sites extraits de la CAH pour chacun des taxons',
                      theme = theme(plot.title = element_text(hjust = 0.5)),
                      caption = 'made by Leonard Heinry')
  }
  
  # compact = rbind(tab1, tab2)
  # return(compact)
}

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

# Comparer les profils -----------------------------------------------------

# L'objectif ici est de récupérer les sites compris dans le 25% de l'interval 
# sur la CAH et de comparer les groupes auxquels les sites sont affiliés entre le
# profil proposé par les DREAL et le taxon concerné.


