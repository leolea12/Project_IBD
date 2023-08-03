library(tidyverse)


IBD_old <- read.csv2("data/Donnees_utilisables/Fichiers/IBD_params.csv", header = TRUE, sep = ";")

trancoded_IBD <- as_tibble(read.csv2("data/Donnees_utilisables/Fichiers/table_transcodage.csv", stringsAsFactors = FALSE)) %>%
  rename(AFNOR = CodeValid) %>%
  select(AFNOR) %>% left_join(IBD_old %>% select(-SANDRE, -DENOMINATION, -Origine), by = "AFNOR") %>%
  drop_na() %>% distinct(AFNOR, .keep_all = TRUE)

trancoded_IBD
