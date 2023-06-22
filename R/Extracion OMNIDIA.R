
# RECUPERATION ET FILTRATION DES CODES OMNIDIA
library(readr)
library(tidyverse)
library(readxl)

# DANS CODE1 selectionner avec start with
code1 <- as.tibble(fread("data/Donnees_utilisables/BaseOmnidia.csv", header = TRUE )) %>% 
  select(CdAppelTaxon, CdAlternatif1,CdAlternatif2,
         CdAlternatif3,CdAlternatif4,
         CdAlternatif5) %>% mutate(CdAppelTaxon = noquote(CdAppelTaxon), CdAlternatif1 = noquote(CdAlternatif1),CdAlternatif2 = noquote(CdAlternatif2),
                                   CdAlternatif3 = noquote(CdAlternatif3),CdAlternatif4 = noquote(CdAlternatif4),
                                   CdAlternatif5 = noquote(CdAlternatif5)) %>% 
  slice(-1) %>% 
  gather(key = "key", value = "value", CdAlternatif1:CdAlternatif5) %>% 
  group_by(CdAppelTaxon,key) %>% mutate(ID = cur_group_id()) %>% 
  distinct()

code1[code1 == ""] <- NA

code2 <- as.tibble(fread("data/Donnees_utilisables/BaseOmnidia.csv", header = TRUE )) %>% 
  select(CdAppelTaxon, OrgCdAlternatif1, 
         OrgCdAlternatif2, OrgCdAlternatif3, 
         OrgCdAlternatif4, 
         OrgCdAlternatif5) %>% mutate(OrgCdAlternatif1 = noquote(OrgCdAlternatif1), 
                                      OrgCdAlternatif2 = noquote(OrgCdAlternatif2), OrgCdAlternatif3 = noquote(OrgCdAlternatif3), 
                                      OrgCdAlternatif4 = noquote(OrgCdAlternatif4), 
                                      OrgCdAlternatif5 = noquote(OrgCdAlternatif5)) %>% 
  slice(-1) %>% 
  gather(key = "key", value = "value", OrgCdAlternatif1:OrgCdAlternatif5) %>% 
  group_by(CdAppelTaxon,key) %>% mutate(ID = cur_group_id()) %>% 
  distinct()

code2[code2 == ""] <- NA

code <- code1 %>% left_join(code2, by = "ID") %>% 
  arrange(key.x, key.y) %>% 
  select(-CdAppelTaxon.y, -key.x, -key.y, 
         SANDRE = CdAppelTaxon.x, code = value.x, 
         origine = value.y, ID) %>% mutate_all(na_if,"") %>% 
  drop_na() %>% 
  filter(origine %in% c("REF_OMNIDIA","TERA_OMNIDIA", "SYN_OMNIDIA")) %>% # AJOUTER SP_OMNIDIA
  arrange(SANDRE, origine) %>% 
  select(-ID) %>% 
  pivot_wider(names_from = origine, values_from = code) %>% 
  pivot_longer(cols=REF_OMNIDIA:TERA_OMNIDIA) %>% filter(!is.na(value)) %>% 
  rename(code_omnidia = value) %>% 
  mutate(SANDRE = as.numeric(SANDRE),
         code_omnidia = as.character(code_omnidia))

true_code = bind_rows(code,
                 as.tibble(read.csv2("data/Donnees_utilisables/Codes 4 lettres taxons.csv", stringsAsFactors = FALSE) %>% 
                             mutate(SANDRE = as.numeric(SANDRE))) %>% 
                   select("SANDRE", name = "auteur_apa", code_omnidia = "cd_apa"), 
                 
                 as.tibble(read.csv2("data/Donnees_utilisables/Complement_code4L.csv", stringsAsFactors = FALSE)) %>% 
                   select("SANDRE", code_omnidia = "TAXON") %>% 
                   drop_na() %>% 
                   mutate(name = "REF_OMNIDIA", SANDRE = as.numeric(SANDRE)) %>% 
                   select(SANDRE, name, code_omnidia),
                 
                 as.tibble(data.frame(SANDRE = c(4767, 38468, 20336, 45516, 66512, 31344, 31603), 
                                      name = "REF_OMNIDIA", 
                            code_omnidia = c("ERMS","FGNO","DIAS","CMFO","LSAP", "CROU", "ADMC")))) %>% 
  distinct(SANDRE, .keep_all = T)


write.csv2(true_code, "data/Donnees_utilisables/CODE_OMNIDIA_traites.csv")

