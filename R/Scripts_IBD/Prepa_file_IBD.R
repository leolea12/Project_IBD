# Direct depuis NAIADES

library(downloader)

library(tidyverse)

library(utils)

library(sf)

library(data.table)

`%notin%` <- Negate(`%in%`)

options(timeout = 3600)

download(
  "https://naiades.eaufrance.fr/reports/reportsperyear/HB/Naiades_Export_France_Entiere_HB.zip",
  dest = "R/Scripts_IBD/dataset.zip",
  mode = "wb"
)
unzip("R/Scripts_IBD/dataset.zip",
      exdir = "./R/Scripts_IBD")
file.remove("R/Scripts_IBD/cep.csv")
file.remove("R/Scripts_IBD/resultat.csv")
file.remove("R/Scripts_IBD/DescriptionDonneesHB.pdf")

france_metropolitaine <- st_read("R/Scripts_IBD/FRA_adm0.shp")

Diatom <- as_tibble(fread("R/Scripts_IBD/fauneflore.csv")) %>%
  dplyr::select(
    "CODE_STATION" = CdStationMesureEauxSurface,
    "Nom_groupe_taxo" = LbSupport,
    "DATE" = DateDebutOperationPrelBio,
    "SANDRE" = CdAppelTaxon,
    "Nom_latin_taxon" = NomLatinAppelTaxon,
    "RESULTAT" = RsTaxRep,
    "Code_groupe_taxo" = CdSupport
  ) %>%
  filter(Code_groupe_taxo == 10) %>%
  dplyr::select(-Code_groupe_taxo) %>%
  distinct(CODE_STATION,
           Nom_groupe_taxo,
           DATE,
           SANDRE,
           Nom_latin_taxon,
           RESULTAT) %>%
  dplyr::select(-Nom_groupe_taxo) %>%
  arrange(DATE,
          CODE_STATION,
          Nom_latin_taxon, RESULTAT) %>%
  filter(RESULTAT != 0) %>%
  mutate(DATE = as.Date(DATE)) %>%
  arrange(DATE) %>%
  # group_by(CODE_STATION, DATE) %>%
  # dplyr::mutate(tot = sum(RESULTAT)) %>%
  # ungroup() %>%
  # dplyr::mutate(RESULTAT = (RESULTAT / tot) * 1000) %>%
  # dplyr::mutate(RESULTAT = round(RESULTAT, 2)) %>% # Passage des abondances en relative pour 1000
  # dplyr::select(-tot) %>%
  
  left_join(tibble(read.csv2("R/Scripts_IBD/operation.csv", sep = ";")) %>% 
              dplyr::select(CODE_STATION = CdStationMesureEauxSurface,
                            DATE =  DateDebutOperationPrelBio, 
                            CODE_OPERATION = RefOperationPrelBio,
                            CdSupport) %>% 
              dplyr::filter(lubridate::year(DATE) > 2006,
                            CdSupport == 10) %>% 
              dplyr::mutate(DATE = as.Date(DATE)) %>%
              dplyr::arrange(DATE) %>%
              rename(CODE_STATION_bis = CODE_STATION, DATE_bis = DATE) %>%
  dplyr::select(-CdSupport), by = c("CODE_STATION" = "CODE_STATION_bis", 
                                    "DATE" = "DATE_bis"), relationship = "many-to-many") %>%

  left_join(read.csv2("R/Scripts_IBD/CODE_OMNIDIA_traites.csv", stringsAsFactors = FALSE),
            by = "SANDRE") %>%
  filter(SANDRE != 0) %>%
  rename(taxon = code_omnidia) %>%
  mutate(CODE_STATION = str_remove(CODE_STATION, "^0+")) %>%
  left_join(
    as_tibble(
      read.csv2(
        "R/Scripts_IBD/stations.csv",
        stringsAsFactors = FALSE,
        quote = "",
        na.strings = c("", "NA")
      )
    ) %>% select(
      CODE_STATION = CdStationMesureEauxSurface,
      commune = LbCommune,
      longitude = CoordXStationMesureEauxSurface,
      latitude = CoordYStationMesureEauxSurface
    ) %>%
      mutate(
        longitude = as.numeric(longitude),
        latitude = as.numeric(latitude)
      ) %>%
      filter(!is.na(longitude)) %>%
      filter(longitude > 0) %>%
      mutate(CODE_STATION = str_remove(CODE_STATION, "^0+")),
    by = "CODE_STATION"
  ) %>%
  drop_na() %>%
  sf::st_as_sf(coords = c("longitude", "latitude"), crs = 2154) %>%
  st_transform(geometry, crs = 4326) %>%
  st_intersection(france_metropolitaine) %>%
  st_jitter(factor = 0.00001) %>%
  dplyr::select(CODE_OPERATION,
                CODE_STATION,
                DATE,
                SANDRE,
                Nom_latin_taxon,
                RESULTAT,
                taxon,
                commune) %>%
  tidyr::extract(geometry, c("long", "lat"), "\\((.*), (.*)\\)", convert = TRUE) %>%
  left_join(as_tibble(
    read.csv2("R/Scripts_IBD/table_transcodage.csv", stringsAsFactors = FALSE)
  ) %>%
    dplyr::select(taxon = "abre", True_name = "CodeValid"),
  by = "taxon") %>%
  mutate(taxon = if_else(is.na(True_name) == T, taxon, True_name)) %>%
  dplyr::select(-True_name) %>%
  filter(!is.na(taxon)) %>%
  left_join(
    read.csv2("R/Scripts_IBD/table_transcodage.csv", stringsAsFactors = FALSE) %>%
      select(taxon = CodeValid, full_name = name_valid) %>% distinct() %>%
      mutate(full_name = sub("\\_g.*", "", full_name)),
    by = "taxon"
  ) %>%
  mutate(full_name = str_replace_all(full_name, "[^[:alnum:]]", " ")) %>%
  mutate(full_name = paste0(full_name, " ", "(", taxon, ")")) %>%
  mutate(lon = round(long, 10), lat = round(lat, 10)) %>%
  left_join(
    read.csv2("R/Scripts_IBD/table_transcodage.csv", stringsAsFactors = FALSE) %>%
      select(abre, name, taxon = CodeValid) %>% unique() %>%
      group_by(taxon) %>% filter(abre %notin% taxon) %>% mutate(list = paste0(abre, " ", sub("\\_g.*", "", name))) %>%
      mutate(taxons_apparies = paste(list, collapse = " / ")) %>%
      select(-abre,-name,-list) %>% distinct(),
    by = "taxon"
  ) %>% 
  mutate(CodeValid = full_name) %>%
  separate_rows(taxons_apparies, sep = " / ") %>%
  group_by(CodeValid) %>%
  mutate(taxons_apparies = ifelse(is.na(taxons_apparies) == TRUE, "Aucun", paste0(str_sub(taxons_apparies,  start = 6)," (", str_sub(taxons_apparies,  start = 1, end = 4),")"))) %>%
  group_by(full_name) %>%
  mutate(grp = cur_group_id()) %>%
  mutate(taxons_apparies = ifelse(taxons_apparies == "Aucun", taxons_apparies, paste0(full_name, " / ", taxons_apparies))) %>%
  separate_rows(taxons_apparies, sep = " / ") %>%
  distinct(taxons_apparies, CODE_STATION, DATE, .keep_all = TRUE) %>%
  mutate(full_name = taxons_apparies) %>%
  ungroup() %>%
  group_by(grp, DATE, CODE_STATION) %>%
  mutate(taxons_apparies = map_chr(row_number(), ~paste(unique(taxons_apparies[-.x]), collapse = " / "))) %>%
  ungroup() %>%
  dplyr::mutate(taxons_apparies = ifelse(taxons_apparies == "", "Aucun", taxons_apparies)) %>%
  mutate(full_name = ifelse(full_name == "Aucun", paste0(Nom_latin_taxon, " (", taxon,")"), 
                            full_name)) %>% 
  dplyr::select(CODE_OPERATION, CODE_STATION, DATE, CODE_TAXON = taxon, RESULTAT) %>%
  distinct(CODE_OPERATION, CODE_TAXON, .keep_all = TRUE)


write.table(Diatom,
            file = "R/Scripts_IBD/Test_IBD.txt",
            sep = "\t", 
            col.names = TRUE, 
            row.names = FALSE)

file.remove("R/Scripts_IBD/stations.csv")
file.remove("R/Scripts_IBD/fauneflore.csv")
file.remove("R/Scripts_IBD/operation.csv")
file.remove("R/Scripts_IBD/dataset.zip")


