# load("R/sysdata.rda")
save(list = c("wordhash",
              "zfin_granges",
              "Drerio",
              "meme",
              "hg38_granges",
              "Hsapiens",
              "hg38_full_model_gr",
              "dr11_full_model_gr"),
     file = "R/sysdata.rda", compress = TRUE)
