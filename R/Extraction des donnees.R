# Chargement des packages 

require(tidyverse)
require(InraeThemes)
require(readxl)
require(fuzzyjoin)
library(readr)
library(leaflet)
library(sf)
require(tictoc)
library(lubridate)
library(purrr)
library(data.table)
library(progress)
library(rlist)
library(openxlsx)
library(purrr)
library(cooccur)
library(magrittr)
library(RODBC)
library(RMySQL)
library(RSQLite)
library(sqldf)


# require(RSQLite)
# 
# drv <- dbDriver("MySQL")
# 
# conR <- dbConnect(drv,
#                   dbname = "pandore",
#                   host = "10.69.192.179",
#                   port = 3306,
#                   user = "",
#                   password = "")
# 
# sort(dbListTables(conR))
# 
# sqlutils::execQuery(connection = conR, query = "Requete_1_seb_le_meilleur")

# Recuperation des codes 4 lettres des taxons 
code <- as.tibble(read.csv2("data/Donnees_utilisables/Fichiers/CODE_OMNIDIA_traites.csv", stringsAsFactors = FALSE)) %>% 
  dplyr::select(-X) %>% 
  mutate(code_taxon = code_omnidia) %>% 
  dplyr::select(-code_omnidia)

# Base de donnees NAIADES ------------------------------------------------
# Importation de la base de donnees biologiques de NAIADES (Penser à télécharger 
# la dernière version en ammont)

Diatom <- as.tibble(data.table::fread("data/Donnees_de_base/fauneflore.csv")) %>%
  dplyr::select("CODE_STATION" = CdStationMesureEauxSurface,
                "Nom_groupe_taxo" = LbSupport,
                "DATE" = DateDebutOperationPrelBio,
                "SANDRE" = CdAppelTaxon,
                "Nom_latin_taxon" = NomLatinAppelTaxon,
                "RESULTAT" = RsTaxRep,
                "Code_groupe_taxo" = CdSupport) %>% 
  filter(Code_groupe_taxo == 10) %>% 
  dplyr::select(-Code_groupe_taxo) %>% 
  distinct(CODE_STATION, Nom_groupe_taxo, DATE, SANDRE, Nom_latin_taxon, RESULTAT) %>%
  arrange(DATE,
          CODE_STATION,
          Nom_latin_taxon, RESULTAT) %>% 
  filter(RESULTAT != 0)

# Recupération données pandore --------------------------------------------


# Requêtage de la base de données

# Lien SQL / R
channel=odbcConnect("pandore")


# JOINTURE PANDORE / DONNEES NAIADES +
# On va enlever les prélèvement qui ne vont pas jusqu'à 400 valves et ceux qui
# vont bien au dessus. En effet pour qu'un prélèvement soit viable, il
# faut qu'environ 400 valves aient été comptées, sinon on donne trop de poids
# au prélèvement faibles et forts 

Diatom2 <- Diatom %>% mutate(DATE = as.Date(DATE)) %>% dplyr::select(-Nom_groupe_taxo) %>% 
  bind_rows(as.tibble(sqlQuery(channel, "SELECT cd_opecont, cd_site, date_opecont, cd_taxon, comptage, nom_taxon FROM listes_diatomee 
                   left join pandore.compil using (cd_opecont) 
                   left join pandore.taxo_diatomee_sandre using (cd_taxon)
                               WHERE date_opecont BETWEEN '2007-01-01' AND '2022-09-28'") %>% #Recuperation des donnees
                        dplyr::select(CODE_STATION = "cd_site",
                                      DATE = "date_opecont", SANDRE = "cd_taxon", Nom_latin_taxon = "nom_taxon", RESULTAT = "comptage")) %>% 
              distinct(CODE_STATION, DATE, SANDRE, RESULTAT, .keep_all = T) %>% 
              filter(RESULTAT != 0)) %>%
  arrange(DATE) %>% 
  group_by(CODE_STATION, DATE) %>% 
  dplyr::mutate(tot=sum(RESULTAT)) %>% 
  ungroup() %>% 
  dplyr::mutate(RESULTAT = (RESULTAT/tot)*1000) %>% # Passage des abondances en relative pour 1000
  dplyr::select(-tot)

# Etape de transcodification du fichier FLORE pour recuperer les bons noms
# de taxons et ceux manquants

# Remplacement anciennes taxonommie par nouvelle
Transcoded_Flore <- Diatom2 %>% 
  left_join(code %>% dplyr::select(-name), by = "SANDRE") %>%
  filter(SANDRE != 0) %>% 
  # Fusion de la nouvelle taxonomie et remplaecement
  left_join(as.tibble(read.csv2("data/Donnees_utilisables/Fichiers/table_transcodage.csv", stringsAsFactors = FALSE)) %>% 
              dplyr::select(code_taxon = "abre", True_name = "CodeValid"), by = "code_taxon") %>% 
  mutate(code_taxon = if_else(is.na(True_name) == T, code_taxon, True_name)) %>% 
  dplyr::select(-True_name) %>% filter(!is.na(code_taxon)) %>% 
  group_by(CODE_STATION, DATE, code_taxon) %>% 
  dplyr::mutate(RESULTAT = sum(RESULTAT)) %>% 
  ungroup() %>% 
  distinct(CODE_STATION, DATE, code_taxon, .keep_all = TRUE)


channel=odbcConnect("pandore")
stations_pandore <- as.tibble(sqlQuery(channel, "SELECT cd_site, commune, x, y FROM site")) %>% 
  select(CODE_STATION = cd_site, commune, longitude = x, latitude = y) %>% 
  filter(!is.na(longitude)) %>% 
  filter(longitude > 0)

stations_NAIADES <- as.tibble(
  read.csv2("data/Donnees_utilisables/Fichiers/stations.csv", stringsAsFactors = FALSE,
            quote = "", na.strings=c("","NA"))
  ) %>% select(CODE_STATION = CdStationMesureEauxSurface, commune = LbCommune,
               longitude = CoordXStationMesureEauxSurface, latitude = CoordYStationMesureEauxSurface) %>% 
  mutate(longitude = as.numeric(longitude),
         latitude = as.numeric(latitude)) %>% 
  filter(!is.na(longitude)) %>% 
  filter(longitude > 0)

full_stations <- bind_rows(stations_NAIADES, stations_pandore) %>% 
  distinct(CODE_STATION, .keep_all = TRUE)

# load("data/Donnees_utilisables/Diatomees.RData")

# Recuperation des données chimiques de Pandore
tic()
channel=odbcConnect("pandore")

Chimie_pandore <- as.tibble(sqlQuery(channel, "SELECT cd_opecont, cd_site, date_opecont, cd_obsPhyChi, cd_support, resultat, cd_param, cd_unite, nom_param FROM sandre_parametre
    join pandore.listes_phychi using(cd_param)
    join pandore.opecont_phychi using(cd_obsPhyChi)
    join pandore.compil using(cd_opecont)
    WHERE cd_param = 1302 OR cd_param = 1304 OR cd_param = 1311 OR 
                        cd_param = 1314 OR cd_param = 1335 OR cd_param = 1433 OR cd_param = 1340 AND cd_support = 3") %>% 
                              dplyr::rename(CODE_STATION = cd_site,
                                            Code_Unite_mesure = cd_unite,
                                            Nom_parametre = nom_param,
                                            Code_parametre = cd_param,
                                            DATE = date_opecont,
                                            Concentration = resultat,
                                            Code_Prelevement = cd_obsPhyChi,
                                            CdSupport = cd_support) %>% 
                              arrange(DATE) %>% 
                              subset(CODE_STATION %in% c(unique(Transcoded_Flore$CODE_STATION))) %>% 
                              mutate(Code_Unite_mesure = as.character(Code_Unite_mesure),
                                     Code_Prelevement = as.character(Code_Prelevement)) %>% 
                              select(CODE_STATION, Code_Prelevement, DATE,       
                                     Code_parametre, Nom_parametre,
                                     Concentration, Code_Unite_mesure))



toc()

# save(Chimie_pandore, file = "data/Donnees_utilisables/Chimie_pandore.RData")
# save(Transcoded_Flore, file = "data/Donnees_utilisables/Diatomees.RData")
# save(Transcoded_Flore, file = "IBD_2022_shiny/Data/Diatomees.RData")

# # Regroupement de toutes les années et enregistrement des données ---------

# On utilise la fonction Pick_date pour lancer automatiquement l'extraction
# des donnees et la fusion chimie_biologie

# ATTENTION, besoin de 3h00 minimum pour tout récupérer

# chargement des fonctions pour filtration des donnees 

source("R/Fonction_fusion_chimie_biologie.R")
source("R/Fonction_pick_date.R")

# load("data/Donnees_utilisables/Diatomees.RData")
# load("data/Donnees_utilisables/Chimie_pandore.RData")

Pick_date(2007, Chimie = Chimie_pandore, Diatomees = Transcoded_Flore)
Pick_date(2008, Chimie = Chimie_pandore, Diatomees = Transcoded_Flore)
Pick_date(2009, Chimie = Chimie_pandore, Diatomees = Transcoded_Flore)
Pick_date(2010, Chimie = Chimie_pandore, Diatomees = Transcoded_Flore)
Pick_date(2011, Chimie = Chimie_pandore, Diatomees = Transcoded_Flore)
Pick_date(2012, Chimie = Chimie_pandore, Diatomees = Transcoded_Flore)
Pick_date(2013, Chimie = Chimie_pandore, Diatomees = Transcoded_Flore)
Pick_date(2014, Chimie = Chimie_pandore, Diatomees = Transcoded_Flore)
Pick_date(2015, Chimie = Chimie_pandore, Diatomees = Transcoded_Flore)
Pick_date(2016, Chimie = Chimie_pandore, Diatomees = Transcoded_Flore)
Pick_date(2017, Chimie = Chimie_pandore, Diatomees = Transcoded_Flore)
Pick_date(2018, Chimie = Chimie_pandore, Diatomees = Transcoded_Flore)
Pick_date(2019, Chimie = Chimie_pandore, Diatomees = Transcoded_Flore)
Pick_date(2020, Chimie = Chimie_pandore, Diatomees = Transcoded_Flore)
Pick_date(2021, Chimie = Chimie_pandore, Diatomees = Transcoded_Flore)

x <- as.tibble(list(file = c("data/Donnees_utilisables/Format donnees traitees/2007_traitee.xlsx", "data/Donnees_utilisables/Format donnees traitees/2007_traitee.xlsx",
                             "data/Donnees_utilisables/Format donnees traitees/2008_traitee.xlsx", "data/Donnees_utilisables/Format donnees traitees/2008_traitee.xlsx",
                             "data/Donnees_utilisables/Format donnees traitees/2009_traitee.xlsx", "data/Donnees_utilisables/Format donnees traitees/2009_traitee.xlsx",
                             "data/Donnees_utilisables/Format donnees traitees/2010_traitee.xlsx", "data/Donnees_utilisables/Format donnees traitees/2010_traitee.xlsx",
                             "data/Donnees_utilisables/Format donnees traitees/2011_traitee.xlsx", "data/Donnees_utilisables/Format donnees traitees/2011_traitee.xlsx",
                             "data/Donnees_utilisables/Format donnees traitees/2012_traitee.xlsx", "data/Donnees_utilisables/Format donnees traitees/2012_traitee.xlsx",
                             "data/Donnees_utilisables/Format donnees traitees/2013_traitee.xlsx", "data/Donnees_utilisables/Format donnees traitees/2013_traitee.xlsx",
                             "data/Donnees_utilisables/Format donnees traitees/2014_traitee.xlsx", "data/Donnees_utilisables/Format donnees traitees/2014_traitee.xlsx",
                             "data/Donnees_utilisables/Format donnees traitees/2015_traitee.xlsx", "data/Donnees_utilisables/Format donnees traitees/2015_traitee.xlsx",
                             "data/Donnees_utilisables/Format donnees traitees/2016_traitee.xlsx", "data/Donnees_utilisables/Format donnees traitees/2016_traitee.xlsx",
                             "data/Donnees_utilisables/Format donnees traitees/2017_traitee.xlsx", "data/Donnees_utilisables/Format donnees traitees/2017_traitee.xlsx",
                             "data/Donnees_utilisables/Format donnees traitees/2018_traitee.xlsx", "data/Donnees_utilisables/Format donnees traitees/2018_traitee.xlsx",
                             "data/Donnees_utilisables/Format donnees traitees/2019_traitee.xlsx", "data/Donnees_utilisables/Format donnees traitees/2019_traitee.xlsx",
                             "data/Donnees_utilisables/Format donnees traitees/2020_traitee.xlsx", "data/Donnees_utilisables/Format donnees traitees/2020_traitee.xlsx",
                             "data/Donnees_utilisables/Format donnees traitees/2021_traitee.xlsx", "data/Donnees_utilisables/Format donnees traitees/2021_traitee.xlsx"),
                    sheet = c("chimie", "flore", "chimie", "flore",
                              "chimie", "flore", "chimie", "flore",
                              "chimie", "flore", "chimie", "flore",
                              "chimie", "flore", "chimie", "flore",
                              "chimie", "flore", "chimie", "flore",
                              "chimie", "flore", "chimie", "flore",
                              "chimie", "flore", "chimie", "flore",
                              "chimie", "flore")))


# We read flora and chemical sheets independently
Files <- map2(x$file, x$sheet, ~ read_xlsx(path = .x, sheet = .y))

# We seperatly store chemical and flora sheets, chemical = odd, flore = even
chimie_sheet = seq(1,30,2)
flore_sheet = seq(2,30,2)

# Merge all chemical sheets together
Chimie <- bind_rows(Files[chimie_sheet]) %>%
  select(CODE_STATION,DATE,
         Code_parametre,Nom_parametre,Code_Unite_mesure,Mediane)

# Merge all flora sheets together
Flore <- bind_rows(Files[flore_sheet]) %>%
  select(CODE_STATION, DATE, code_taxon, SANDRE, RESULTAT, Commune, x, y) %>%
  group_by(CODE_STATION, DATE) %>% distinct()


save(Chimie, Flore, file = "data/Donnees_utilisables/Fichiers/Donnees_completes.RData")
write.csv2(Flore, "data/Donnees_utilisables/Fichiers/Flore.csv")


