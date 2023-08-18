library(tidyverse)


IBD_old <- read.csv2("data/Donnees_utilisables/Fichiers/IBD_params.csv", header = TRUE, sep = ";")

trancoded_IBD <- as_tibble(read.csv2("data/Donnees_utilisables/Fichiers/table_transcodage.csv", stringsAsFactors = FALSE)) %>%
  rename(AFNOR = CodeValid) %>%
  select(AFNOR) %>% left_join(IBD_old, by = "AFNOR") %>%
  drop_na() %>% distinct(AFNOR, .keep_all = TRUE)

write.csv2(trancoded_IBD, "data/Donnees_utilisables/Fichiers/IBD_params.csv", row.names = FALSE)
# save(trancoded_IBD, file = "R/Scripts_IBD/IBD_profiles.RData")


trancoded_IBD2 <- trancoded_IBD %>%
  bind_rows(trancoded_IBD %>% filter(AFNOR == "ADEU") %>% mutate(AFNOR = "ADAM")) %>%
  bind_rows(trancoded_IBD %>% filter(AFNOR == "ACLI") %>% mutate(AFNOR = "ADCV")) %>%
  bind_rows(trancoded_IBD %>% filter(AFNOR == "ADMI") %>% mutate(AFNOR = "ADMC")) %>%
  bind_rows(trancoded_IBD %>% filter(AFNOR == "ADMI") %>% mutate(AFNOR = "ADMO")) %>%
  bind_rows(trancoded_IBD %>% filter(AFNOR == "ADSB") %>% mutate(AFNOR = "ADRU")) %>%
  bind_rows(trancoded_IBD %>% filter(AFNOR == "ADSH") %>% mutate(AFNOR = "ADSK")) %>%
  bind_rows(trancoded_IBD %>% filter(AFNOR == "ADCT") %>% mutate(AFNOR = "ADTC")) %>%
  bind_rows(trancoded_IBD %>% filter(AFNOR == "ADMI") %>% mutate(AFNOR = "AHOF")) %>%
  bind_rows(trancoded_IBD %>% filter(AFNOR == "ADMS") %>% mutate(AFNOR = "ALBL")) %>%
  bind_rows(trancoded_IBD %>% filter(AFNOR == "ACOP") %>% mutate(AFNOR = "AMCD")) %>%
  bind_rows(trancoded_IBD %>% filter(AFNOR == "APFI") %>% mutate(AFNOR = "AZHA")) %>%
  bind_rows(trancoded_IBD %>% filter(AFNOR == "CHEL") %>% mutate(AFNOR = "CEXF")) %>%
  bind_rows(trancoded_IBD %>% filter(AFNOR == "COPL") %>% mutate(AFNOR = "CROX")) %>%
  bind_rows(trancoded_IBD %>% filter(AFNOR == "CHEL") %>% mutate(AFNOR = "CSBH")) %>%
  bind_rows(trancoded_IBD %>% filter(AFNOR == "CHEL") %>% mutate(AFNOR = "CSUT")) %>%
  bind_rows(trancoded_IBD %>% filter(AFNOR == "CAFF") %>% mutate(AFNOR = "CTDE")) %>%
  bind_rows(trancoded_IBD %>% filter(AFNOR == "EARC") %>% mutate(AFNOR = "EARB")) %>%
  bind_rows(trancoded_IBD %>% filter(AFNOR == "ESUM") %>% mutate(AFNOR = "EBNA")) %>%
  bind_rows(trancoded_IBD %>% filter(AFNOR == "ECPM") %>% mutate(AFNOR = "ECAL")) %>%
  bind_rows(trancoded_IBD %>% filter(AFNOR == "EBLU") %>% mutate(AFNOR = "ECTO")) %>%
  bind_rows(trancoded_IBD %>% filter(AFNOR == "EBLU") %>% mutate(AFNOR = "EJUE")) %>%
  bind_rows(trancoded_IBD %>% filter(AFNOR == "HLAC") %>% mutate(AFNOR = "PUBL")) %>%
  bind_rows(trancoded_IBD %>% filter(AFNOR == "ESUM") %>% mutate(AFNOR = "ENEE")) %>%
  bind_rows(trancoded_IBD %>% filter(AFNOR == "ESLE") %>% mutate(AFNOR = "ENSI")) %>%
  bind_rows(trancoded_IBD %>% filter(AFNOR == "FVAU") %>% mutate(AFNOR = "FCAD")) %>%
  bind_rows(trancoded_IBD %>% filter(AFNOR == "FVAU") %>% mutate(AFNOR = "FMIV")) %>%
  bind_rows(trancoded_IBD %>% filter(AFNOR == "FSAX") %>% mutate(AFNOR = "FNEV")) %>%
  bind_rows(trancoded_IBD %>% filter(AFNOR == "FMES") %>% mutate(AFNOR = "FNIN")) %>%
  bind_rows(trancoded_IBD %>% filter(AFNOR == "UDEL") %>% mutate(AFNOR = "FPDE")) %>%
  bind_rows(trancoded_IBD %>% filter(AFNOR == "FBID") %>% mutate(AFNOR = "FPRU")) %>%
  bind_rows(trancoded_IBD %>% filter(AFNOR == "GEXL") %>% mutate(AFNOR = "GAGV")) %>%
  bind_rows(trancoded_IBD %>% filter(AFNOR == "GCOR") %>% mutate(AFNOR = "GAUR")) %>%
  bind_rows(trancoded_IBD %>% filter(AFNOR == "GPRI") %>% mutate(AFNOR = "GCUN")) %>%
  bind_rows(trancoded_IBD %>% filter(AFNOR == "GELG") %>% mutate(AFNOR = "GMIS")) %>%
  bind_rows(trancoded_IBD %>% filter(AFNOR == "GLAT") %>% mutate(AFNOR = "GTNO")) %>%
  bind_rows(trancoded_IBD %>% filter(AFNOR == "GRHB") %>% mutate(AFNOR = "GVRD")) %>%
  bind_rows(trancoded_IBD %>% filter(AFNOR == "NGRE") %>% mutate(AFNOR = "HPDA")) %>%
  bind_rows(trancoded_IBD %>% filter(AFNOR == "HLMO") %>% mutate(AFNOR = "HTHU")) %>%
  bind_rows(trancoded_IBD %>% filter(AFNOR == "LHLU") %>% mutate(AFNOR = "LGOP")) %>%
  bind_rows(trancoded_IBD %>% filter(AFNOR == "NGER") %>% mutate(AFNOR = "NDDF")) %>%
  bind_rows(trancoded_IBD %>% filter(AFNOR == "SSVE") %>% mutate(AFNOR = "NFSO")) %>%
  bind_rows(trancoded_IBD %>% filter(AFNOR == "NFON") %>% mutate(AFNOR = "NSTS")) %>%
  bind_rows(trancoded_IBD %>% filter(AFNOR == "NFON") %>% mutate(AFNOR = "NYCO")) %>%
  bind_rows(trancoded_IBD %>% filter(AFNOR == "PLFR") %>% mutate(AFNOR = "PLRC")) %>%
  bind_rows(trancoded_IBD %>% filter(AFNOR == "PTLA") %>% mutate(AFNOR = "PMNT")) %>%
  bind_rows(trancoded_IBD %>% filter(AFNOR == "PULA") %>% mutate(AFNOR = "POVA")) %>%
  bind_rows(trancoded_IBD %>% filter(AFNOR == "FCRO") %>% mutate(AFNOR = "SRBU")) %>%
  bind_rows(trancoded_IBD %>% filter(AFNOR == "POBL") %>% mutate(AFNOR = "PSXO")) %>%
  bind_rows(trancoded_IBD %>% filter(AFNOR == "SDIF") %>% mutate(AFNOR = "SBOS")) %>%
  bind_rows(trancoded_IBD %>% filter(AFNOR == "NAMP") %>% mutate(AFNOR = "SCRA")) %>%
  bind_rows(trancoded_IBD %>% filter(AFNOR == "SEBA") %>% mutate(AFNOR = "SEAT")) %>%
  bind_rows(trancoded_IBD %>% filter(AFNOR == "SPUP") %>% mutate(AFNOR = "SECA")) %>%
  bind_rows(trancoded_IBD %>% filter(AFNOR == "SPUP") %>% mutate(AFNOR = "SESP")) %>%
  bind_rows(trancoded_IBD %>% filter(AFNOR == "SRHE") %>% mutate(AFNOR = "SPDV")) %>%
  bind_rows(trancoded_IBD %>% filter(AFNOR == "NEXI") %>% mutate(AFNOR = "SRAE")) %>%
  bind_rows(trancoded_IBD %>% filter(AFNOR == "SSRT") %>% mutate(AFNOR = "SSBG")) %>%
  dplyr::select(-SANDRE,-DENOMINATION)


write.csv2(trancoded_IBD2, "data/Donnees_utilisables/Fichiers/IBD_params2.csv", row.names = FALSE)


