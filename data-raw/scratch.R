library(blaseRtools)
library(blaseRdata)
library(tidyverse)

identical(vignette_cds, filter_cds(vignette_cds))
identical(monocle3::exprs(vignette_cds), monocle3::exprs(filter_cds(vignette_cds)))
filter_cds(vignette_cds,
           cells = bb_cellmeta(vignette_cds) |> filter(equipment == "chromium"),
           genes = bb_rowmeta(vignette_cds) |> filter(module == "1"))
