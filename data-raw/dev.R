
devtools::document()

# make the documents for the website
purrr::walk(.x = list.files(path = "vignettes", pattern = "*.Rmd", full.names = FALSE),
     .f = function(x) {
       rmarkdown::render(input = file.path("vignettes", x), output_dir = "docs/pages", output_format = "html_document")
       rfile <- str_replace(x, ".Rmd", ".R")
       knitr::purl(input = file.path("vignettes", x), output = file.path("docs", "pages",rfile))
     })

# commit and push
gert::git_add("*")
gert::git_commit("version 0.0.0.9090")
gert::git_push()

