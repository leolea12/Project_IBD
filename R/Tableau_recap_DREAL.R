

# Ouverture du même tableau que celui pris pour le shiny

dt_Chorologie <- as.tibble(read.csv2("data/Donnees_utilisables/Chorologie.csv", stringsAsFactors = FALSE)) %>%
  select(-X)

tab_recap <- dt_Chorologie %>% group_by(code_taxon) %>% 
  summarise(
    Periode_recensement = paste(min(annee), "-", max(annee)),
    Nombre_annees_vu = length(unique(annee)),
    Nombre_stations = length(unique(CODE_STATION)),
    Tranche_occurence = paste(min(RESULTAT_1000), "-", max(RESULTAT_1000)),
    Importance_max = max(RESULTAT_1000),
    Priorité = case_when(max(RESULTAT_1000) < 100 ~ "++"))




library(xlsx)
write.xlsx(tab_recap, file = "data/Donnees_utilisables/Tab_Recap.xlsx")
