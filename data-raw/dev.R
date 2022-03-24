
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
gert::git_commit("version 0.0.0.9064")
gert::git_push()

# build and insert into repo
pkg_build <- devtools::build()

drat::insertPackage(file = pkg_build,
                    repodir = "/home/OSUMC.EDU/blas02/network/X/Labs/Blaser/share/data/R/drat/",
                    action = "archive")

unlink(pkg_build)
