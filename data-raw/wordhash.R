wordhash <- lexicon::hash_grady_pos %>%
  as_tibble() %>%
  mutate(word = str_replace_all(word, "[:punct:]|[:space:]", "_")) %>%
  group_by(word) %>%
  summarise() %>%
  mutate(index = row_number(word)) %>%
  relocate(index)

save(wordhash, file = "data/wordhash.rda")
