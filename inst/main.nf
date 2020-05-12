// Check that the input parameters are set
assert params.study, 'Please specify a study with --study'
assert params.inputDir, 'Please specify an input folder with --inputDir'
assert params.outputDir, 'Please specify an output folder with --outputDir'
assert params.username, 'Please specify ImmPort Username with --username'
assert params.password, 'Please specify ImmPort Password with --password'

process processStudy {
  container = 'rglab/hipccyto:latest'

  echo true

  input:
  val study from params.study
  file inputDir from Channel.fromPath(params.inputDir)
  val outputDir from params.outputDir
  val username from params.username
  val password from params.password

  output:
  file "v*"

  publishDir "${outputDir}", mode: 'copy', overwrite: true

  """
  #!/usr/bin/env Rscript
  library(HIPCCyto)

  Sys.setenv(ImmPortUsername="${username}")
  Sys.setenv(ImmPortPassword="${password}")

  ver <- paste0("v", packageVersion("HIPCCyto"))
  dir.create(ver)

  gsl <- process_study("${study}", "${inputDir}")

  save_gating_sets(gsl, ver)
  """
}
