load("R/sysdata.rda")
save(list = c("wordhash",
              "zfin_granges",
              "Drerio",
              "meme",
              "hg38_granges",
              "Hsapiens"),
     file = "R/sysdata_temp.rda", compress = TRUE)
