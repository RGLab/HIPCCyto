## code to prepare `DATA` dataset goes here
DATA <- list(
  "SDY820" = list(
    map = data.frame(
      alias = c("APC-eFluor 780-A", "eFluor 450-A"),
      channels = c("APC-eFluor780-A", "eFluor450-A")
    ),
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
warning(which(table(headers$pns_reported) > 1))


usethis::use_data(DATA, MARKERS, overwrite = TRUE, internal = TRUE)
