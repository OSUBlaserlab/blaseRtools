devtools::load_all()

annika_analysis_file <- read_csv("~/network/X/Labs/Blaser/Annika/fluorescence scope/20250519_hgfa_drug_analysis-file.csv") |>
  mutate(filepaths = fs::path(bb_fix_file_path(filepaths)))

bb_blind_images(analysis_file = annika_analysis_file,
                output_dir = fs::path("./blinding_test"),
                file_column = "filepaths")

unblinded_result <- bb_unblind_images(
  directory = "~/network/X/Labs/Blaser/Annika/Fluorescence scope/20250513_20250513112841739153/",
  keyfile = "blinding_key.csv",
  scorefile = "scoresheet.csv",
  analysis_file = "~/network/X/Labs/Blaser/Annika/fluorescence scope/20250512_hgfa_drug_analysis-file.csv",
  file_column = "filepaths"
)




