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
    batch = "CST SETUP DATE",
    live_method = "tailgate",
    live_args = "num_peaks = 3, ref_peak = 3"
  ),
  "SDY514" = list(
    batch = "CST SETUP DATE"
  ),
  "SDY144" = list(
    map = data.frame(
      alias = c("PE-Texas Red-A", "Qdot 565-A", "Qdot 565-A"),
      channels = c("PE-TxR-A", "QDOT 565-A", "QDot 565-A")
    ),
    target = c(60000, 5000)
  ),
  "SDY212" = list(
    batch = "EXPERIMENT NAME"
  ),
  "SDY296" = list(
    batch = "EXPERIMENT NAME"
  ),
  "SDY312" = list(
    batch = "EXPERIMENT NAME"
  ),
  "SDY887" = list(
    toRemove = 0.04
  ),
  "SDY702" = list(
    toRemove = 0.04,
    target = c(100000, 30000)
  ),
  "SDY1097" = list(
    toRemove = 0.04,
    target = c(100000, 30000)
  ),
  "SDY301" = list(
    target = c(40000, 5000)
  ),
  "SDY387" = list(
    target = c(60000, 6000),
    batch = "EXPERIMENT NAME"
  ),
  "SDY368" = list(
    target = c(60000, 6000)
  ),
  "SDY180" = list(
    target = c(45000, 30000),
    batch = "EXPERIMENT NAME"
  )
)

MARKERS <- character()
MARKERS[c("TCRgd", "TCR GAMMA DELTA", "TCRGD", "gdTCR")] <- "TCRgd" # SDY212 SDY514
# MARKERS[c("STAT1", "STAT 1", "STAT-1")] <- "STAT1"
# MARKERS[c("STAT3", "STAT 3", "STAT-3")] <- "STAT3"
# MARKERS[c("STAT5", "STAT 5", "STAT-5")] <- "STAT5A"
MARKERS[c("CD4/CD20", "CD4 / CD20", "CD4 /CD20", "CD4/ CD20", "CD20/CD4")] <- "CD4/CD20" # SDY212 SDY514
MARKERS[c("Cd85j")] <- "CD85j" # SDY514
MARKERS[c("Lin-1", "LIN1", "lin-1")] <- "LIN-1" # SDY144 SDY296
MARKERS[c("SLAM")] <- "SLAN" # SDY144

usethis::use_data(DATA, MARKERS, overwrite = TRUE, internal = TRUE)
