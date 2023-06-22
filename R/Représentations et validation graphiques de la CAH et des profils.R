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
library(overlapping)
library(RODBC)
library(ggpubr)



load("data/Donnees_utilisables/Donnees_test_serveur_CEUO.RData")
load("data/Donnees_utilisables/Donnees_test_serveur.RData")
load("data/Donnees_utilisables/Donnees_completes.RData")
load("data/Donnees_utilisables/New_taxons.RData")
load("data/Donnees_utilisables/site_taxon.RData")
load("data/Donnees_utilisables/site_taxon_CEUO.RData")
load("data/Donnees_utilisables/Dendo.RData")
load("data/Donnees_utilisables/Dendo_CEUO.RData")

get_density <- function(x, y) {
  dens <- MASS::kde2d(x, y)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}


# Valider la CAH  ---------------------------------------------------------

load("data/Donnees_utilisables/Cluster_2015_2021.RData")
# load("data/Donnees_utilisables/Dist_mat_2015_2021.RData")
load("data/Donnees_utilisables/Donnees_test_serveur.RData")
load("data/Donnees_utilisables/Donnees_completes.RData")
load("data/Donnees_utilisables/sites_CAH.RData")

# Nombre de clusters souhaité
nb_clusters <- 15

# Du fait du volume de données que represente la matrice de distance, 
# on stock ses noms de sites pour y avoir accès facilement
sites <- rownames(as.matrix(dist_bc))
save(sites, file = "data/Donnees_utilisables/sites_CAH.RData")

# Boucle pour couper la CAH avec le nombre de clusters souhaité
liste_groupes <- list()

# Boucle pour couper la CAH avec le nombre de clusters souhaité
for (i in 1:(nb_clusters)) {
  clusters <- cutree(clus_ward, k = i)
  groupes <- data.frame(site = sites, groupe = clusters)
  liste_groupes[[i]] <- groupes
}

# Compilation des tableaux de groupes dans un data frame avec une colonne supplémentaire pour la valeur de i
resultat <- do.call(rbind, lapply(1:(nb_clusters), function(i) {cbind(liste_groupes[[i]], i)}))

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


to_plot = merge(resultat, Chim, by = "site", all.x = TRUE)

# Liste des groupes
groupes <- sort(unique(to_plot$i))[-1]

# Initialisation d'une liste pour stocker les boxplots de chaque groupe
boxplots <- list()

# Boucle pour créer un boxplot pour chaque groupe
for (g in groupes) {
  # Sous-ensemble des données pour le groupe actuel
  # Création du boxplot pour le groupe actuel
  dt <- as.data.frame(subset(to_plot, i == g)) |> 
    filter(Nom_parametre == "Nitrates") |> 
    mutate(Mediane = log(Mediane))
  
  n_grp <- unique(dt$groupe)
  
  # Générer toutes les paires possibles de noms de groupe
  paires <- combn(n_grp, 2, simplify = FALSE)
  
  bp <- # Affichage du boxplot avec les p-values
    ggboxplot(dt, x = "groupe", y = "Mediane", 
              color = "groupe", palette = "jco",
              add = "mean", legend = "none") +
    stat_compare_means(
      comparisons = paires,
      method = "wilcox.test",
      label = "p.format",
      paired = FALSE)
  
  # Ajout du boxplot à la liste
  boxplots[[as.character(g)]] <- bp
}

boxplots[[4]]

save(boxplots, file = "data/Donnees_utilisables/boxplots_significativite_CAH.RData")

# Creation d'un tableau recap pour le rapport -----------------------------

install.packages("webshot")
library(webshot)
library(knitr)
library(kableExtra)
library(tidyverse)

load("data/Donnees_utilisables/New_taxons.RData")
load("data/Donnees_utilisables/site_taxon.RData")
# Construction de GPPY absent de tous les relevés

Methode_tab1 <- site_taxon %>% 
  mutate(DATE = (sub(".*?_","", site)),
         year = sub("_.*", "", DATE)) %>% 
  filter(taxon %in% unique(New_taxons %>% pull(abre))) %>% 
  group_by(taxon, year) %>% summarise(count = as.character(n()),
                                      mean_ab = round(mean(Ab_relative)/10,1)) %>% 
  # mutate(count_ab = paste0(count,"/",mean_ab)) %>% 
  ungroup() %>% 
  select(-mean_ab) %>% 
  pivot_wider(names_from = year, values_from = count) %>% 
  replace(is.na(.), "absent") %>% left_join(
    site_taxon %>% 
      mutate(DATE = (sub(".*?_","", site)),
             year = sub("_.*", "", DATE)) %>% 
      filter(taxon %in% unique(New_taxons %>% pull(abre))) %>% 
      group_by(taxon, year) %>% summarise(count = as.character(n()),
                                          mean_ab = as.character(round(mean(Ab_relative)/10,1))) %>% 
      # mutate(count_ab = paste0(count,"/",mean_ab)) %>% 
      ungroup() %>% 
      select(-count) %>% 
      pivot_wider(names_from = year, values_from = mean_ab) %>% 
      replace(is.na(.), "absent"), by = 'taxon') %>% 
  add_row(tibble(taxon = "GPPY", 
                 `2015.x` = "absent",
                 `2016.x` = "absent",
                 `2017.x` = "absent",
                 `2018.x` = "absent",
                 `2019.x` = "absent",
                 `2020.x` = "absent",
                 `2021.x` = "absent",
                 `2015.y` = "absent",
                 `2016.y` = "absent",
                 `2017.y` = "absent",
                 `2018.y` = "absent",
                 `2019.y` = "absent",
                 `2020.y` = "absent",
                 `2021.y`= "absent")) %>% 
  dplyr::arrange(taxon)

Methode_tab2 <- site_taxon %>% 
  mutate(DATE = (sub(".*?_","", site)),
         year = sub("_.*", "", DATE)) %>% 
  filter(taxon %in% unique(New_taxons %>% pull(abre))) %>% 
  group_by(taxon, year) %>% summarise(count = n(),
                                      mean_ab = round(mean(Ab_relative)/10,1)) %>% 
  mutate(count_ab = paste0(count,"/",mean_ab)) %>% 
  ungroup() %>% group_by(taxon) %>% 
  mutate("Total d'apparitions" = sum(count),
         "Mediane d'abondance" = median(mean_ab),
         "Nombre d'années" = length(year)) %>% 
  select(-count, -mean_ab) %>% 
  ungroup() %>% 
  distinct(taxon, .keep_all = TRUE) %>% 
  add_row(tibble(taxon = "GPPY", year = c("none"),
                 `Total d'apparitions` = 0, `Mediane d'abondance` = 0,
                 `Nombre d'années` = 0
  )) %>% 
  dplyr::arrange(taxon)

Methode_tab <- Methode_tab1 %>% left_join(Methode_tab2, by= "taxon") %>% 
  as.data.frame() %>% 
  select(Taxon = taxon, "2015" = '2015.x', "2016" = '2016.x',"2017" = '2017.x',"2018" = '2018.x',"2019" = '2019.x',"2020" = '2020.x',"2021" = '2021.x',
         " 2015" = '2015.y'," 2016" = '2016.y'," 2017" = '2017.y'," 2018" = '2018.y'," 2019" = '2019.y'," 2020" = '2020.y'," 2021" = '2021.y',
         "Total d'apparitions", "Nombre d'années", "Mediane d'abondance") %>% 
  dplyr::arrange(Taxon)

save_kable(Methode_tab %>% 
             mutate(`Total d'apparitions` = as.numeric(`Total d'apparitions`)) %>% 
             mutate(`Total d'apparitions` = cell_spec(`Total d'apparitions`, color = ifelse(`Total d'apparitions` >= 30, "green","red"))) %>% 
             mutate(
               across(2:15, 
                      ~ cell_spec(.x, 
                                  color = ifelse(.x == "absent", "#990000", "black")))
             ) %>% 
             kable(booktabs = T, linesep = "", format = "html", escape = F, table.attr = "style='width:30%;'") %>% 
             kable_classic(full_width = F, html_font = "Cambria") %>% 
             add_header_above(c(" " = 1, "Nombre d'apparition dans les relevés" = 7, 
                                "Abondance relative moyenne par échantillon" = 7, 
                                "Résumé" = 3), font_size = 20) %>% 
             kable_styling(fixed_thead = T, font_size = 15) %>% 
             row_spec(0, font_size= 25,
                      bold = T) %>% 
             column_spec(2:8, background = "#CCCCCC") %>% 
             column_spec(9:15, background = "gainsboro") %>% 
             column_spec(16:17, background = "#CCCCCC") %>% 
             column_spec(18, color = "white", background = spec_color(Methode_tab$`Mediane d'abondance`, option = "D")) %>% 
             row_spec(0:69, align = "c") 
           ,file = "Table_methodo.html")

webshot("Table_methodo.html", "Table_methodo.pdf")
Sys.setenv("OPENSSL_CONF"="/dev/null")

load("data/Donnees_utilisables/clus_ward.RData")

plot(clus_ward, labels = FALSE, main = "", 
     xlab = "", ylab = "", sub = "", axes = FALSE, hang = -1)


# Plot des paires de taxons candidats, Density plot en 2D et 3D -----------------

#2D 

load("data/Donnees_utilisables/tab_density.RData")



tab_density %>% 
  filter(taxon %in% New_taxons[46,]$abre) %>%
  ggplot(aes(x = as.character(log(pos)), y = taxon)) +
  geom_tile(aes(fill = density)) +
  scale_fill_gradient(low = "white", high = "red") +
  theme_minimal() +
  theme(axis.text.x = element_blank()) +
  labs(x = "Site", y = "Taxon", fill = "density")

tab_density %>%
  filter(taxon == "AAMB")%>%
  ggplot(aes(x = log(pos), y = density))+
  geom_point()


ggsave("plot.png", plot = test)




plot_abondance <- function(Dendo, tax, site_taxon){
  select_tax <- site_taxon %>% filter(taxon == tax)
  label_select_tax <- data.frame(site = labels(Dendo)) %>% 
    left_join(select_tax, by = "site") %>% 
    mutate(site = if_else(taxon == tax, site, ""))
  
  label_select_tax %>% 
    mutate(site = seq(1:length(label_select_tax$site)),
           Ab_relative = if_else(is.na(Ab_relative) == TRUE, 0, Ab_relative)) %>% 
    drop_na() %>% 
    mutate(density = get_density(site, Ab_relative)) %>% 
    ggplot(aes(x = site, y = Ab_relative))+
    # stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white", contour = TRUE)+
    geom_point(aes(color = density), alpha = 0.8,
               size = 1)+
    scale_color_viridis()+
    theme_classic()+
    theme(axis.text.x=element_blank())
}


get_density <- function(x, y) {
  dens <- MASS::kde2d(x, y)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

?MASS::kde2d()
plot_abondance(Dend = Dendo, tax = "CEXF", site_taxon = site_taxon)




# -------------------------------------------------------------------

load("data/Donnees_utilisables/Mat_dist_CAH.RData")

select_tax <- tab %>% filter(taxon == tax)

label_select_tax <- data.frame(site = labels(Dend)) %>%
  left_join(select_tax, by = "site") %>%
  mutate(site = if_else(taxon == tax, site, ""))

all_sites <- label_select_tax %>% mutate(site = seq(1:length(label_select_tax$site)),
                                         Ab_relative = if_else(is.na(Ab_relative) == TRUE, 0, Ab_relative)) %>% 
  drop_na()

# On recupere dans la matrice de distance de branches de la CAH les sites 
# de présence du taxon et la matrice des distance qui va avec
essaie <- Mat_dist_CAH[all_sites$site,all_sites$site]

# Appliquer l'algorithme MDS pour projeter les sites dans un espace à une dimension
mds <- cmdscale(essaie, k = 1)

# Extraire les positions projetées
positions <- as.vector(mds[,1])
min <- abs(min(positions))
positions <- positions + min
# # Créer un jeu de données avec une seule observation pour chaque site
df <- data.frame(x = positions, y = all_sites$Ab_relative)
#
# # Tracer un graphique avec un axe x gradué en respectant les distances entre les sites
# p <- ggplot(df, aes(x = x, y = y)) +
#   geom_point(size = 0.1, position=position_jitter(h=0.5,w=0.5)) +
#   scale_x_continuous(breaks = positions, labels = rownames(essaie))
#   # stat_density_2d(
#   #   geom = "raster",
#   #   aes(fill = after_stat(density)),
#   #   contour = FALSE
#   # ) + scale_fill_viridis_c()
# 
# # Afficher le graphique
# print(p)

# On construit des colonnes qui permettront de visualiser les vraies distances entre les sites 
tax_density <- all_sites %>% 
  mutate(dist_site = positions,
         density = get_density(positions, all_sites$Ab_relative))

# créer un vecteur de données simulées
library(cluster)

sil_scores <- sapply(2:(length(unique(tax_density$dist_site))-1), function(k) {
  km <- kmeans(tax_density$dist_site, k)
  ss <- silhouette(km$cluster, dist(tax_density$dist_site))
  mean(ss[,3])
})

best_k <- which.max(sil_scores) + 1
print(paste0("Le nombre optimal de classes est ", best_k))

classes <- cut(tax_density$dist_site, best_k, labels = FALSE)
tax_density$group <- factor(classes)

library(ggplot2)

tax_density %>% 
  ggplot(aes(x = dist_site, y = Ab_relative, color = group)) +
  geom_point() +
  scale_fill_viridis_d() +
  theme_classic()


# -------------------------------------------------------------------


# 3D

library(plotly)
plot_density <- function(Dendo, tax, site_taxon){
  select_tax <- site_taxon %>% filter(taxon == tax)
  no_select_tax <- site_taxon %>% filter(site %notin% unique(select_tax %>% pull(site))) %>% 
    distinct(site)
  
  label_select_tax <- data.frame(site = labels(Dendo)) %>% 
    left_join(select_tax, by = "site") %>% 
    mutate(site = if_else(taxon == tax, site, ""))
  
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
  
  label_select_tax %>% 
    mutate(site = seq(1:length(label_select_tax$site)),
           Ab_relative = if_else(is.na(Ab_relative) == TRUE, 0, Ab_relative)) %>% 
    drop_na() %>% 
    mutate(density = get_density(site, Ab_relative)) %>% 
    filter(taxon == tax) %>% 
    ggplot(aes(x=site, y = Ab_relative)) +
    geom_point()+
    # geom_density(fill="grey", color="grey", alpha=0.4)+
    # geom_vline(data=vlines,
    #            aes(xintercept=value, colour=Mean),
    #            size=1, linetype="dashed", show.legend=TRUE)+
    theme_classic()+
    theme(axis.text.x=element_blank())
}

Density_3D <- function(tax,tax2, Dendo, tab){
  
  p1 = plot_density(Dendo, tax = tax, site_taxon)
  p2 = plot_density(Dendo, tax = tax2, site_taxon)
  
  
  p3D1 <- MASS::kde2d(p1$data$Ab_relative, p1$data$site)
  p3D2 <- MASS::kde2d(p2$data$Ab_relative, p2$data$site)
  
  df1 <- data.frame(x = data.frame(p3D1)$x, y = data.frame(p3D1)$y, z = data.frame(p3D1) %>% 
                      slice(1) %>% 
                      dplyr::select(-x,-y) %>% 
                      unlist() %>% 
                      unname())
  
  df2 <- data.frame(x = data.frame(p3D2)$x, y = data.frame(p3D2)$y, z = data.frame(p3D2) %>% 
                      slice(1) %>% 
                      dplyr::select(-x,-y) %>% 
                      unlist() %>% 
                      unname())
  
  
  fig1 <- plot_ly(x = p3D1$x, y = p3D1$y, z = p3D1$z, scene = 'scene1', type = "surface") 
  fig1 <- fig1 %>% add_surface(contours = list(z = list(show = TRUE, 
                                                        usecolormap = TRUE, 
                                                        project = list(z = TRUE)))) %>% 
    add_trace(data = df1, x = ~x, y = ~y, z = ~z, mode = "markers", type = "scatter3d", 
              marker = list(size = 4, color = "red", symbol = 104))
  
  fig2 <- plot_ly(x = p3D2$x, y = p3D2$y, z = p3D2$z, scene = 'scene2') 
  fig2 <- fig2 %>% add_surface(contours = list(z = list(show = TRUE, 
                                                        usecolormap = TRUE, 
                                                        project = list(z = TRUE)))) %>% 
    add_trace(data = df2, x = ~x, y = ~y, z = ~z, mode = "markers", type = "scatter3d", 
              marker = list(size = 4, color = "red", symbol = 104))
  
  fig3 <- plot_ly(x = p3D1$x, y = p3D1$y, z = p3D1$z, scene = 'scene3') 
  fig3 <- fig3 %>% add_surface(contours = list(z = list(show = TRUE, 
                                                        usecolormap = TRUE, 
                                                        project = list(z = TRUE)))) %>% 
    add_trace(data = df1, x = ~x, y = ~y, z = ~z, mode = "markers", type = "scatter3d", 
              marker = list(size = 4, color = "red", symbol = 104))
  
  
  fig4 <- plot_ly(x = p3D2$x, y = p3D2$y, z = p3D2$z, scene = 'scene4') 
  fig4 <- fig4 %>% add_surface(contours = list(z = list(show = TRUE, 
                                                        usecolormap = TRUE, 
                                                        project = list(z = TRUE)))) %>% 
    add_trace(data = df2, x = ~x, y = ~y, z = ~z, mode = "markers", type = "scatter3d", 
              marker = list(size = 4, color = "red", symbol = 104))
  
  
  fig <- plotly::subplot(fig1, fig2, fig3, fig4)
  
  axx <- list(
    gridcolor='rgb(255, 255, 255)',
    zerolinecolor='rgb(255, 255, 255)',
    showbackground=TRUE,
    backgroundcolor='rgb(230, 230,230)'
  )
  
  fig <- fig %>% plotly::layout(title = "3D Subplots",
                                
                                scene = list(domain=list(x=c(0,0.5),y=c(0.5,1)),
                                             camera = list(eye = list(x = 1.87, 
                                                                      y = 0.88, 
                                                                      z = -0.64)),
                                             xaxis=axx, yaxis=axx, zaxis=axx),
                                
                                scene2 = list(domain=list(x=c(0.5,1),y=c(0.5,1)),
                                              camera = list(eye = list(x = 1.87, 
                                                                       y = 0.88, 
                                                                       z = -0.64)),
                                              xaxis=axx, yaxis=axx, zaxis=axx),
                                
                                scene3 = list(domain=list(x=c(0,0.5),y=c(0,0.5)),
                                              camera = list(eye = list(x = 1.5, 
                                                                       y = 0.88, 
                                                                       z = 0.64)),
                                              xaxis=axx, yaxis=axx, zaxis=axx),
                                
                                scene4 = list(domain=list(x=c(0.5,1),y=c(0,0.5)),
                                              camera = list(eye = list(x = 1.5, 
                                                                       y = 0.88, 
                                                                       z = 0.64)),
                                              xaxis=axx, yaxis=axx, zaxis=axx))
  
  
  fig
}

plot_abondance(Dendo, "ADTC", site_taxon)
plot_abondance(Dendo, "NIGR", site_taxon)

Density_3D(tax = "AAMB", tax2 = "NIGR", Dendo = Dendo, tab = site_taxon)



