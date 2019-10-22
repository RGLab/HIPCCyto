// Check that the input parameters are set
assert params.study, 'Please specify a study with --study'
assert params.inputDir, 'Please specify an input folder with --inputDir'
assert params.username, 'Please specify ImmPort Username with --username'
assert params.password, 'Please specify ImmPort Password with --password'

params.outputDir = 'gs'

process hipcCyto {
  container = 'hipccyto:latest'

  echo true

  input:
  val study from params.study
  file inputDir from Channel.fromPath(params.inputDir)
  val inputPath from params.inputDir
  val username from params.username
  val password from params.password
  val outputDir from params.outputDir

  output:
  file "${outputDir}"

  publishDir "${inputPath}", mode: 'copy', overwrite: true

  """
  #!/usr/local/bin/Rscript
  Sys.setenv(ImmPortUsername = "${username}")
  Sys.setenv(ImmPortPassword = "${password}")
  HIPCCyto:::process_study("${study}", "${inputDir}", "${outputDir}")
  """
}
