#'@rdname Image-class
Image_ <- setClass(
  Class = "Image_",
  slots = list(
    file_path = "character",
    species = "character",
    stage = "character",
    genetics = "character",
    treatment = "character",
    microscope = "character",
    mag = "character",
    filter = "character",
    use = "character",
    note = "character",
    md5sum = "character"
  ),
  prototype = list(
    file_path = NA_character_,
    species = NA_character_,
    stage = NA_character_,
    genetics = NA_character_,
    treatment = NA_character_,
    microscope = NA_character_,
    mag = NA_character_,
    filter = NA_character_,
    use = NA_character_,
    note = NA_character_,
    md5sum = NA_character_
  )
)

setValidity("Image_",
            function(object) {
              if (!fs::file_exists(object@file_path))
                return("The file does not exist!")
              # if (!"fs_path" %in% class(object))
              #   return("The file_path slot must be of type fs::path.")
            })

#' @title An S4 Class for Holding Image Metadata
#' @description
#' This is an S4 Class for holding the relevant metadata we need for key images.  The idea is to generate this object for each of the images we will use in a grant, paper or other important document.  That way when we want to reuse these images we know where they are.
#'
#' One can construct a single image object using the Image() constructor method.  When called as such, it will open a file chooser to identify the file from the network drive.  Then it will provide an interactive menu to add the needed metadata (see below).
#'
#' The workflow is to use this to create a new Image object when you are using a new key image in a grant or other document.  Then you will add the image to the image catalog using ImageCatalog.add.
#'
#' The Image constructor provides several validation checks, including that the image file must be accessible.
#'
#' Each slot has it's own getter and setter methods which are identical to the name of the slot.
#'
#' **It is critical that your images are stored in a common network drive.  Ideally this is X/Labs/Blaser/staff/keyence imaging data**
#'
#' **Avoid Duplication!! Please keep your all of your raw imaging data in X/Labs/Blaser/staff and subdirectories.**
#'
#' @slot file_path Path to file.  Should start with ~/network/X/Labs/Blaser...
#' @slot species The species being imaged.
#' @slot stage The stage of the sample.  Options include various times in hpf plus other for other timepoints.
#' @slot genetics Any genetic modifications.
#' @slot treatment Any treatments performed.
#' @slot microscope The microscope used.
#' @slot mag Magnification
#' @slot filter The filter or camera setup used.
#' @slot use The document the image is being used in.
#' @slot note Any additional notes.
#'
#' @export
Image <- setClass(
  Class = "Image",
  slots = list(
    file_path = "character",
    species = "character",
    stage = "character",
    genetics = "character",
    treatment = "character",
    microscope = "character",
    mag = "character",
    filter = "character",
    use = "character",
    note = "character",
    md5sum = "character"
  ),
  prototype = list(
    file_path = NA_character_,
    species = NA_character_,
    stage = NA_character_,
    genetics = NA_character_,
    treatment = NA_character_,
    microscope = NA_character_,
    mag = NA_character_,
    filter = NA_character_,
    use = NA_character_,
    note = NA_character_,
    md5sum = NA_character_
  ),
  contains = "Image_"
)



# image class methods -------------------------------
#' @importFrom cli cli_text
#' @importFrom stringr str_detect str_remove
#' @importFrom fs path path_home
#' @importFrom digest digest
setMethod("initialize",
          "Image",
          # function(.Object, ...) {
          function(.Object, file_path = NULL, ...) {
            .Object <- callNextMethod()
            if (is.null(file_path)) {
              cli::cli_text("Select an image file.  Note:  this file must be in a common network directory.")
              file_path <- file.choose()

            }

            if (stringr::str_detect(file_path, "~", negate = TRUE)) {
              file_path <-
                paste0("~", stringr::str_remove(file_path, fs::path_home()))

            }
            file_path <- fs::path(file_path)

            choices <- c("Zebrafish", "Human", "Mouse", "Other")
            species <-
              menu(choices, title = "Select the species of the specimen.")
            if (species == 4) {
              species <-
                readline("Enter a free text description of the species. ")
            } else {
              species <- choices[species]
            }

            choices <-
              c("24 hpf", "48 hpf", "72 hpf", "96 hpf", "Other")
            stage <-
              menu(choices, title = "Select the developmental stage.")
            if (stage == 5) {
              stage <-
                readline("Enter a free text description of the stage of the specimen if applicable. ")
            } else {
              stage <- choices[stage]
            }

            genetics <-
              readline("Enter a free-text description including all of the specimen's relevant genetics. ")
            treatment <-
              readline("Enter a free-text description including any of the specimen's chemical treatments. ")

            choices <-
              c("Keyence", "inverted", "dissecting", "Other")
            microscope <-
              menu(choices, title = "Select the microscope used for imaging.")
            if (microscope == 4) {
              microscope <- readline("Enter the name of the microscope. ")
            } else {
              microscope <- choices[microscope]
            }

            mag <-
              menu(c("1x", "2X", "4X", "10X", "20X", "40X", "Other"), title = "Enter the magnification used.")
            if (mag == 7) {
              mag <- readline("Enter the magnification.")
            } else {
              mag <- choices[mag]
            }

            filter <-
              menu(
                c(
                  "Green",
                  "Red",
                  "Blue",
                  "Brightfield - Greyscale",
                  "Brightfield - Full Color",
                  "Other"
                ),
                title  = "Enter the filter/camera options."
              )
            if (filter == 6) {
              filter <-
                readline("Enter a free-text description of the filter/camera used. ")
            } else {
              filter <- choices[filter]
            }

            use <-
              readline("Describe the document this image is being used in:  ")
            note <-
              readline("Add any other notes to this image. ")

            .Object@file_path <- file_path
            .Object@species <- species
            .Object@stage <- stage
            .Object@genetics <- genetics
            .Object@treatment <- treatment
            .Object@microscope <- microscope
            .Object@mag <- mag
            .Object@filter <- filter
            .Object@use <- use
            .Object@note <- note
            .Object@md5sum <-
              digest::digest(file = file_path, algo = "md5")



            .Object
          })


#' @rdname Image-class
#' @importFrom cli cli_text
#' @export
setMethod("show",
          "Image_",
          function(object) {
            cli::cli_text("An image file at {object@file_path} with hash {object@md5sum}.")
            cli::cli_text("Species:  {object@species}.")
            cli::cli_text("Stage:  {object@stage}.")
            cli::cli_text("Genetics:  {object@genetics}.")
            cli::cli_text("Treatment:  {object@treatment}.")
            cli::cli_text("Microscope:  {object@microscope}.")
            cli::cli_text("Magnification:  {object@mag}.")
            cli::cli_text("Filter:  {object@filter}.")
            cli::cli_text("Used in:  {object@use}.")
            cli::cli_text("Note:  {object@note}.")

          })


# setter and getter methods -----------------------
#' @rdname Image-class
#' @export
setGeneric("file_path", function(x)
  standardGeneric("file_path"))

#' @rdname Image-class
#' @export
setGeneric("file_path<-", function(x, value)
  standardGeneric("file_path<-"))
setMethod("file_path", "Image_", function(x)
  x@file_path)
setMethod("file_path<-", "Image_", function(x, value) {
  x@file_path <- value
  x
})

#' @rdname Image-class
#' @export
setGeneric("species", function(x)
  standardGeneric("species"))

#' @rdname Image-class
#' @export
setGeneric("species<-", function(x, value)
  standardGeneric("species<-"))
setMethod("species", "Image_", function(x)
  x@species)
setMethod("species<-", "Image_", function(x, value) {
  x@species <- value
  x
})

#' @rdname Image-class
#' @export
setGeneric("stage", function(x)
  standardGeneric("stage"))
#' @rdname Image-class
#' @export
setGeneric("stage<-", function(x, value)
  standardGeneric("stage<-"))
setMethod("stage", "Image_", function(x)
  x@stage)
setMethod("stage<-", "Image_", function(x, value) {
  x@stage <- value
  x
})

#' @rdname Image-class
#' @export
setGeneric("genetics", function(x)
  standardGeneric("genetics"))

#' @rdname Image-class
#' @export
setGeneric("genetics<-", function(x, value)
  standardGeneric("genetics<-"))
setMethod("genetics", "Image_", function(x)
  x@genetics)
setMethod("genetics<-", "Image_", function(x, value) {
  x@genetics <- value
  x
})

#' @rdname Image-class
#' @export
setGeneric("treatment", function(x)
  standardGeneric("treatment"))
#' @rdname Image-class
#' @export
setGeneric("treatment<-", function(x, value)
  standardGeneric("treatment<-"))
setMethod("treatment", "Image_", function(x)
  x@treatment)
setMethod("treatment<-", "Image_", function(x, value) {
  x@treatment <- value
  x
})

#' @rdname Image-class
#' @export
setGeneric("microscope", function(x)
  standardGeneric("microscope"))

#' @rdname Image-class
#' @export
setGeneric("microscope<-", function(x, value)
  standardGeneric("microscope<-"))
setMethod("microscope", "Image_", function(x)
  x@microscope)
setMethod("microscope<-", "Image_", function(x, value) {
  x@microscope <- value
  x
})

#' @rdname Image-class
#' @export
setGeneric("mag", function(x)
  standardGeneric("mag"))

#' @rdname Image-class
#' @export
setGeneric("mag<-", function(x, value)
  standardGeneric("mag<-"))
setMethod("mag", "Image_", function(x)
  x@mag)
setMethod("mag<-", "Image_", function(x, value) {
  x@mag <- value
  x
})

#' @rdname Image-class
#' @export
setGeneric("filter", function(x)
  standardGeneric("filter"))

#' @rdname Image-class
#' @export
setGeneric("filter<-", function(x, value)
  standardGeneric("filter<-"))

setMethod("filter", "Image_", function(x)
  x@filter)
setMethod("filter<-", "Image_", function(x, value) {
  x@filter <- value
  x
})

#' @rdname Image-class
#' @export
setGeneric("use", function(x)
  standardGeneric("use"))

#' @rdname Image-class
#' @export
setGeneric("use<-", function(x, value)
  standardGeneric("use<-"))

setMethod("use", "Image_", function(x)
  x@use)
setMethod("use<-", "Image_", function(x, value) {
  x@use <- value
  x
})


#' @rdname Image-class
#' @export
setGeneric("note", function(x)
  standardGeneric("note"))

#' @rdname Image-class
#' @export
setGeneric("note<-", function(x, value)
  standardGeneric("note<-"))

setMethod("note", "Image_", function(x)
  x@note)
setMethod("note<-", "Image_", function(x, value) {
  x@note <- value
  x
})

# imageCatalog -----------------------------------------
#' @title An S4 Image Catalog
#' @description
#' This object holds all of the individual Image objects.
#'
#' Methods are provided for
#'
#' * viewing as a tibble:  ImageCatalog.as_tibble
#' * writing to .json format:  ImageCatalog.write
#' * adding an image:  ImageCatalog.add
#' * deleting an image:  ImageCatalog.delete
#' * subsetting and extracting:  brackets and double brackets
#'
#' An ImageCatalog object can be made from a list of Image objects using ImageCatalog(list = <list of images>).
#'
#' More commonly, you will generate the ImageCatalog from a .json file. To make an image catalog this way, run ImageCatalog(json_path = <path>).
#'
#'
#' @export
ImageCatalog <- setClass(Class = "ImageCatalog",
                         slots = list(images = "list",
                                      json_path = "character"))

#' @importFrom purrr map list_rbind
#' @importFrom tibble tibble
#' @importFrom digest digest
#' @importFrom dplyr mutate filter count
setValidity("ImageCatalog",
            function(object) {
              if (!all(mapply(class, object@images) %in% c("Image", "Image_")))
                return("The catalog must only contain Images.")

              md5sum_tbl <-
                purrr::map(.x = object@images, .f = \(x) {
                  tibble::tibble(
                    file_path = x@file_path,
                    check = digest::digest(file = x@file_path, algo = "md5"),
                    md5sum = x@md5sum
                  )
                }) |>
                purrr::list_rbind() |>
                dplyr::mutate(match = ifelse(check == md5sum, "match", "mismatch"))
              md5sum_check <- md5sum_tbl |>
                dplyr::filter(match == "mismatch")
              dupe_check <- md5sum_tbl |>
                dplyr::count(md5sum)

              if (nrow(md5sum_check) > 0) {
                show(md5sum_check)
                return("There is a digest mismatch in the file(s) above.")
              }
              if (max(dupe_check$n > 1)) {
                show(dupe_check |> arrange(desc(n)))
                return("There are duplicates in the image catalog.")
              }

            })

# ImageCAtalog methods -------------------------------
#' @importFrom purrr map_chr map
#' @importFrom stringr str_sub
#' @importFrom fs path
setMethod("initialize",
          "ImageCatalog",
          function(.Object,
                   images = list(),
                   json_path = character(0),
                   ...) {
            # .Object <- callNextMethod()

            if (length(images) > 0) {
              names(images) <- purrr::map_chr(images, \(x) {
                stringr::str_sub(x@md5sum, start = 1L, end = 6L)
              })

              images <- purrr::map(images, \(x) {
                image_ <- new("Image_")
                image_@file_path <- fs::path(x@file_path)
                image_@species <- x@species
                image_@stage <- x@stage
                image_@genetics <- x@genetics
                image_@treatment <- x@treatment
                image_@microscope <- x@microscope
                image_@mag <- x@mag
                image_@filter <- x@filter
                image_@use <- x@use
                image_@note <- x@note
                image_@md5sum <- x@md5sum
                image_

              })
            } else {
              catalog_tibble <- rjson::fromJSON(file = json_path) |>
                tibble::as_tibble()
              images <- purrr::pmap(catalog_tibble,
                                    .f = \(
                                      file_path,
                                      species,
                                      stage,
                                      genetics,
                                      treatment,
                                      microscope,
                                      mag,
                                      filter,
                                      use,
                                      note,
                                      md5sum
                                    ) {
                                      image_ <- new("Image_")
                                      image_@file_path <- fs::path(file_path)
                                      image_@species <- species
                                      image_@stage <- stage
                                      image_@genetics <- genetics
                                      image_@treatment <- treatment
                                      image_@microscope <- microscope
                                      image_@mag <- mag
                                      image_@filter <- filter
                                      image_@use <- use
                                      image_@note <- note
                                      image_@md5sum <- md5sum
                                      image_
                                    })
              names(images) <- purrr::map_chr(images, \(x) {
                stringr::str_sub(x@md5sum, start = 1L, end = 6L)
              })
            }
              .Object@images <- images
              .Object@json_path <- json_path
              callNextMethod()
              .Object

          })

#' @rdname ImageCatalog-class
#' @importFrom cli cli_text
#' @export
setMethod("show",
          "ImageCatalog",
          function(object) {
            cli::cli_text("An image catalog with {length(object@images)} images.")

          })

#' @rdname ImageCatalog-class
#' @export
setMethod("[", c("ImageCatalog"),
          function(x, i, j, ..., drop = TRUE)
          {
            initialize(x, images = x@images[i])
          })

#' @rdname ImageCatalog-class
#' @export
setMethod("[[", c("ImageCatalog"),
          function(x, i, j, ..., drop = TRUE)
          {
            x@images[i]
          })

#' @rdname ImageCatalog-class
#' @export
setMethod("$", "ImageCatalog",
          function(x, name)
          {
            ## 'name' is a character(1)
            x@images[[name]]
          })


#' @rdname ImageCatalog-class
#' @importFrom purrr map list_rbind
#' @importFrom tibble tibble
#' @export
setGeneric("ImageCatalog.as_tibble", function(image_catalog)
  standardGeneric("ImageCatalog.as_tibble"))
setMethod("ImageCatalog.as_tibble", signature("ImageCatalog"), function(image_catalog) {
  purrr::map(.x = image_catalog@images, .f = \(x) {
    tibble::tibble(
      file_path = x@file_path,
      species = x@species,
      stage = x@stage,
      genetics = x@genetics,
      treatment = x@treatment,
      microscope = x@microscope,
      mag = x@mag,
      filter = x@filter,
      use = x@use,
      note = x@note,
      md5sum = x@md5sum
    )
  }) |>
    purrr::list_rbind() #|>
  # tidyr::nest(.by = c(use, note))
})


#' @rdname ImageCatalog-class
#' @importFrom rjson toJSON
#' @export
setGeneric("ImageCatalog.write", function(image_catalog, out)
  standardGeneric("ImageCatalog.write"))
setMethod("ImageCatalog.write", "ImageCatalog", function(image_catalog, out) {
  ImageCatalog.as_tibble(image_catalog) |>
    rjson::toJSON() |>
    cat(file = out)
})

#' @rdname ImageCatalog-class
#' @importFrom fs path
#' @importFrom dplyr filter
#' @importFrom purrr pmap
#' @import cli
#' @export
setGeneric("ImageCatalog.add", function(image_catalog, image)
  standardGeneric("ImageCatalog.add"))
setMethod("ImageCatalog.add", "ImageCatalog", function(image_catalog, image) {
  # check if the image is already in the catalog
  cli::cli_inform("Checking to see if this image is already catalogued...")

  image_hash <- image@md5sum

  cat <- ImageCatalog.as_tibble(image_catalog) |>
    dplyr::filter(md5sum == image_hash)

  if (nrow(cat) == 0) {
    cli::cli_inform("No it isn't.  Adding the image to the image catalog.")
    image_catalog@images <- c(image_catalog@images, image)
  } else if (nrow(cat) == 1) {
    cli::cli_alert_warning("This image is already catalogued!")
    cli::cli_div(theme = list(span.epnh = list(color = "orange")))
    cli::cli_inform(
      "Adding the {.emph use} and {.emph note} slots to the fields for this entry in the image catalog."
    )
    new_use <- image@use
    updated_image <-
      Image_(
        file_path = cat$file_path,
        species = cat$species,
        stage = cat$stage,
        genetics = cat$genetics,
        treatment = cat$treatment,
        microscope = cat$microscope,
        mag = cat$mag,
        filter = cat$filter,
        use = c(cat$use, image@use),
        note = c(cat$note, image@note),
        md5sum = cat$md5sum
      )

      image_list <- purrr::pmap(ImageCatalog.as_tibble(image_catalog) |>
                                dplyr::filter(md5sum != image_hash),
                              .f = \(
                                file_path,
                                species,
                                stage,
                                genetics,
                                treatment,
                                microscope,
                                mag,
                                filter,
                                use,
                                note,
                                md5sum
                              ) {
                                Image_(
                                  file_path = fs::path(file_path),
                                  species = species,
                                  stage = stage,
                                  genetics = genetics,
                                  treatment = treatment,
                                  microscope = microscope,
                                  mag = mag,
                                  filter = filter,
                                  use = use,
                                  note = note,
                                  md5sum = md5sum
                                )
                              })
    image_catalog <-
      ImageCatalog(images = c(image_list, updated_image))

  } else {
    cli::cli_abort("The image catalog has multiple duplicates!")
  }


  validObject(image_catalog)
  image_catalog
})

#' @rdname ImageCatalog-class
#' @importFrom fs path
#' @importFrom dplyr mutate filter select
#' @importFrom stringr str_sub
#' @importFrom purrr pmap
#' @export
setGeneric("ImageCatalog.delete", function(image_catalog, hash)
  standardGeneric("ImageCatalog.delete"))
setMethod("ImageCatalog.delete", "ImageCatalog", function(image_catalog, hash) {
  hash <- stringr::str_sub(hash, start = 1L, end = 6L)
  cat <- ImageCatalog.as_tibble(image_catalog) |>
    dplyr::mutate(hash_short = stringr::str_sub(md5sum, start = 1L, end = 6L)) |>
    dplyr::filter(hash_short != hash) |>
    dplyr::select(-hash_short)
  image_list <- purrr::pmap(cat,
                            .f = \(
                              file_path,
                              species,
                              stage,
                              genetics,
                              treatment,
                              microscope,
                              mag,
                              filter,
                              use,
                              note,
                              md5sum
                            ) {
                              image_ <- new("Image_")
                              image_@file_path <-
                                fs::path(file_path)
                              image_@species <- species
                              image_@stage <- stage
                              image_@genetics <- genetics
                              image_@treatment <- treatment
                              image_@microscope <- microscope
                              image_@mag <- mag
                              image_@filter <- filter
                              image_@use <- use
                              image_@note <- note
                              image_@md5sum <- md5sum
                              image_
                            })
  image_catalog <- ImageCatalog(images = image_list)
  validObject(image_catalog)
  image_catalog
})


