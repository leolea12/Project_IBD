
# Fonction pick_date pour selectionner la base de donnees chimique souhaitee
#' Pick_date, une fonction qui applique f_fusion pour une année donnée
#'
#' @param Date
#'
#' @return Renvoie un fichier xlsx contenant la liste floristique d'une année
#' donnée sur une feuille et la chimie associée sur l'autre
#' @export
#'
#' @examples
Pick_date <- function(Date, Chimie, Diatomees) {
  
  channel=odbcConnect("pandore")
  
  filename <- paste0("data/Donnees_de_base/analyses_", Date, ".csv")
  
  # Fusion NAIADES/PANDORE
  
  Chimie_tab <- Chimie %>%
    filter(year(DATE) == Date) %>%
    bind_rows((as.tibble(fread(filename)) %>%
                 select(
                   "CODE_STATION" = CdStationMesureEauxSurface, "Code_Prelevement" = CdPrelevement,
                   "DATE" = DatePrel, "Code_parametre" = CdParametre, "Nom_parametre" = LbLongParamètre,
                   "Concentration" = RsAna, "Code_Unite_mesure" = CdUniteMesure, "Nom_unite_mesure" = SymUniteMesure, CdSupport
                 ) %>%
                 subset(CODE_STATION %in% c(unique(Diatomees$CODE_STATION))) %>%
                 filter(Code_parametre == 1302 | Code_parametre == 1304 | Code_parametre == 1311 |
                          Code_parametre == 1314 | Code_parametre == 1335 | Code_parametre == 1433 |
                          Code_parametre == 1340, CdSupport == 3) %>%
                 mutate(DATE = as.Date(DATE, tryFormats = c("%Y-%m-%d", "%Y/%m/%d"))) %>%
                 mutate(Concentration = as.numeric(Concentration)))) %>%
    distinct() %>%
    select(-CdSupport) %>%
    group_nest(CODE_STATION, DATE)
  
  Diatomees %>%  
    group_nest(CODE_STATION, DATE) %>% 
    filter(year(DATE) == Date, between(month(DATE), 5, 9)) -> Biol
  
  Tab <- f_fusion(data_chimie = Chimie_tab, data_biologie = Biol)
  
  df1 <- Tab %>%
    unnest(chimie) %>%
    select(-flore)
  
  df2 <- Tab %>%
    unnest(flore) %>%
    select(-chimie)
  
  dataset_names <- list("chimie" = df1, "flore" = df2)
  
  write.xlsx(dataset_names, file = paste0("data/Donnees_utilisables/", Date, "_traitee.csv"))
  
  return(Tab)
}








