## code to prepare `DATA` dataset goes here
DATA <- list(
  "SDY820" = list(
    map = data.frame(
      alias = c("APC-eFluor 780-A", "eFluor 450-A"),
      channels = c("APC-eFluor780-A", "eFluor450-A")
    ),
    # acquired date (to examine batches)
    # there are two batches, but don't correspond with either cohort or
    # study_time_collected
    batch = "CST SETUP DATE"
  ),
  "SDY514" = list(
    batch = "CST SETUP DATE"
  )
)


headers <- Rlabkey::labkey.selectRows(
  baseUrl = "https://www.immunespace.org",
  folderPath = "/Studies",
  schemaName = "immport",
  queryName = "fcs_header_marker",
  colFilter = Rlabkey::makeFilter(c("pns_preferred", "NOT_MISSING", "")),
  colNameOpt = "fieldname"
)
headers <- unique(headers[, c("pns_reported", "pns_preferred")])
MARKERS <- headers$pns_preferred
names(MARKERS) <- headers$pns_reported
MARKERS[c("TCRgd", "TCR GAMMA DELTA", "TCRGD", "gdTCR")] <- "TCRgd"
MARKERS[c("STAT1", "STAT 1", "STAT-1")] <- "STAT1"
MARKERS[c("STAT3", "STAT 3", "STAT-3")] <- "STAT3"
MARKERS[c("STAT5", "STAT 5", "STAT-5")] <- "STAT5A"
MARKERS[c("CD4/CD20", "CD4 / CD20", "CD4 /CD20", "CD4/ CD20", "CD20/CD4")] <- "CD4/CD20"
MARKERS[c("Cd85j")] <- "LILRB1"
warning(which(table(headers$pns_reported) > 1))


usethis::use_data(DATA, MARKERS, overwrite = TRUE, internal = TRUE)
