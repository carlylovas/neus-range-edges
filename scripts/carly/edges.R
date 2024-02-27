## Linear relationship of range edges

# poleward 
poleward.edge.lm <- poldat.stats %>%
  dplyr::select(comname, est_year, spp.lat95) %>% 
  group_by(comname) %>% 
  nest() %>% 
  arrange(comname) %>% 
  mutate(model = purrr::map(data, ~lm(spp.lat95 ~ est_year, data = .x)), 
         tidymodel = purrr::map(model, tidy)) %>%
  unnest(tidymodel)  %>% 
  filter(!term=="(Intercept)") %>%
  dplyr::select(-data) %>%
  ungroup()%>% 
  mutate(model  = as.character(model),
         signif = ifelse(p.value < 0.05, "yes", "no"))

# signif poleward movement 
poleward.edge.lm %>% 
  filter(signif == "yes")

# equatorial 

