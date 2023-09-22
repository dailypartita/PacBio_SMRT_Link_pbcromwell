version 1.0

import "tasks/pbcoretools.wdl" as pbcoretools

workflow prepare_input {
  input {
    File? dataset_xml
    File? reads
    String dataset_ext = ".consensusreadset.xml"
    Int filter_min_qv = 20
    String dataset_filters = ""
    Int? downsample_factor

    Int nproc = 1
    String log_level = "INFO"
    Int max_nchunks = 40
    String tmp_dir = "/tmp"
  }

  if (defined(dataset_xml)) {
    call pbcoretools.dataset_filter {
      input:
        dataset = select_first([dataset_xml]),
        dataset_ext = dataset_ext,
        filters = dataset_filters,
        min_qv = filter_min_qv,
        downsample_factor = downsample_factor,
        nproc = 1,
        log_level = log_level
    }
  }
  Array[File?] all_input_reads = [dataset_filter.filtered, reads]

  output {
    File reads_file = select_first(all_input_reads)
  }
}
