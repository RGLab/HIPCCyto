## code to prepare `DATA` dataset goes here
DATA <- list(
  "SDY820" = list(
    map = data.frame(
      alias = c("APC-eFluor 780-A", "eFluor 450-A"),
      channels = c("APC-eFluor780-A", "eFluor450-A")
    ),
    markers = c(
      "APC-A" = "CD14",
      "Alexa Fluor 700-A" = "CD3E",
      "APC-eFluor 780-A" = "CD4",
      "V500-A" = "LD",
      "FITC-A" = "IL2RA",
      "PE-A" = "NCAM1",
      "PerCP-Cy5-5-A" = "CCR7",
      "PE-Cy7-A" = "CD19",
      "eFluor 450-A" = "PTPRC",
      "BV650-A" = "CD8A"
    ),
    batch = "CST SETUP DATE"
  )
)

usethis::use_data(DATA, internal = TRUE, overwrite = TRUE)
