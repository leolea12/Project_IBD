
library(leaflet)
library(tidyverse)
library(RMySQL)
library(RPostgreSQL)
library(RODBC)
library(stars)
library(lubridate)

# Recuperation des stations

`%notin%` <- Negate(`%in%`)
Communes <- c(
  "Mana", "Roura", "Regina", "Maripasoula", "Papaichton", "Kourou",
  "Saint-Laurent-du-Maroni", "Iracoubo", "Montsinery-Tonnegrande", "Apatou",
  "Massangis"
)

load("data/Donnees_utilisables/Donnees_completes.RData")

channel=odbcConnect("pandore")
stations_pandore <- as.tibble(sqlQuery(channel, "SELECT cd_site, commune, x, y FROM site")) %>% 
  select(CODE_STATION = cd_site, commune, longitude = x, latitude = y) %>% 
  filter(!is.na(longitude)) %>% 
  filter(longitude > 0)

stations_NAIADES <- as.tibble(
  read.csv2("data/Donnees_utilisables/Station.csv", stringsAsFactors = FALSE,
            quote = "", na.strings=c("","NA"))
) %>% select(CODE_STATION = CdStationMesureEauxSurface, commune = LbCommune,
             longitude = CoordXStationMesureEauxSurface, latitude = CoordYStationMesureEauxSurface) %>% 
  mutate(longitude = as.numeric(longitude),
         latitude = as.numeric(latitude)) %>% 
  filter(!is.na(longitude)) %>% 
  filter(longitude > 0)

full_stations <- bind_rows(stations_NAIADES, stations_pandore) %>% 
  distinct(CODE_STATION, .keep_all = TRUE) %>% 
  mutate(CODE_STATION =  str_remove(CODE_STATION, "^0+"))

Station_Flore <- Flore %>% ungroup() %>% 
  filter(year(DATE) >= 2015) %>% 
  distinct(CODE_STATION) %>% 
  mutate(CODE_STATION =  str_remove(CODE_STATION, "^0+")) %>% 
  left_join(full_stations, by = "CODE_STATION") %>% 
  filter(!is.na(latitude)) %>% 
  mutate(
    longitude = if_else(commune == "Massangis", 4.00, longitude),
    latitude = if_else(commune == "Massangis", 47.6, latitude)
  ) %>%
  filter(commune %notin% Communes[-11]) %>% 
  select(CODE_STATION, longitude, latitude) %>% 
  filter(!is.na(longitude)) %>% 
  sf::st_as_sf(coords = c("longitude", "latitude"), crs = 2154) %>% 
  st_transform(geometry, crs = 4326)

# Carte des stations 

map_base <- leaflet() %>%
  addProviderTiles(providers$Esri.WorldGrayCanvas,
                   group = "Fond clair"
  ) %>%
  addProviderTiles(providers$CartoDB.DarkMatter,
                   group = "Fond noir"
  ) %>%
  addProviderTiles(providers$GeoportailFrance.orthos,
                   group = "Fond satellite"
  )

map_base %>%
  addCircleMarkers(
    data = Station_Flore,
    color = "blue",
    radius = 0.2,) %>%
  addLayersControl(
    position = "topleft",
    baseGroups = c(
      "Fond satellite",
      "Fond clair",
      "Fond noir"
    ))


sum_tax <- site_taxon %>% filter(taxon %in% unique(New_taxons %>% pull(abre))) %>% 
  group_by(taxon) %>% mutate(count = n(),
                             mean_ab = mean(Ab_relative)/10) %>% 
  ungroup() %>% 
  select(-site) %>% 
  distinct(taxon, count, mean_ab) %>% 
  dplyr::arrange(taxon)










