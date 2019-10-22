study = 'SDY820'
inputDir = Channel.fromPath('/home/jkim2345/ImmPort/SDY820')

process basicExample {

  echo true

  input:
  val study
  file inputDir

  """
  Rscript -e "HIPCCyto:::process_study('$study', '$inputDir')"
  """

}
