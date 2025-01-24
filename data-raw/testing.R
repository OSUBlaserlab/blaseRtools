devtools::load_all()
test1 <- Image(file = "/home/OSUMC.EDU/blas02/network/X/Labs/Blaser/Brad/test_img/1.png")
test2 <- Image(file = "/home/OSUMC.EDU/blas02/network/X/Labs/Blaser/Brad/test_img/2.png")
test3 <- Image(file = "/home/OSUMC.EDU/blas02/network/X/Labs/Blaser/Brad/test_img/3.png")
test4 <- Image(file = "~/network/X/Labs/Blaser/Brad/test_img/4.png")
test4.1 <- test4
test5 <- Image()
test1_ <- Image_()
ImageCatalog()
test_cat1 <- ImageCatalog(images = list(test4, test5))
test_cat1@json_path
test_cat1[[1]]
test_cat1["f912f7"]
test_cat1[["f912f7"]]
test_cat1$f912f7
as_tibble(test_cat1)
ImageCatalog.write(test_cat1, out = "test.json")

test_cat2 <- ImageCatalog(json_path = "test.json")
test_cat2@json_path

test_cat2[[1]]

waldo::compare(test_cat1, test_cat2)
validObject(test_cat2)
test_cat3 <- ImageCatalog.add(test_cat1, test3)
as_tibble(test_cat3)

test_cat3 <- ImageCatalog.delete(test_cat3, hash = "c0296b")
as_tibble(test_cat3)

test_cat <- ImageCatalog(images = list(test1, test2))
test1.1 <- test1
use(test1) <- "foo"
note(test1) <- "bar"
test_cat <- ImageCatalog.add(test_cat, test1)
as_tibble(test_cat)
as_tibble(test_cat) |> nest(.by = c(-use, -note), .key = "uses_and_notes")
as_tibble(test_cat) |> nest(.by = c(-use, -note), .key = "uses_and_notes") |> unnest(cols = c(`uses and notes`))
as_tibble(test_cat) |> nest(.by = c(-use, -note), .key = "uses_and_notes") |> dplyr::filter(str_detect(md5sum, "1dc7b"))
