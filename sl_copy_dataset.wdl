# Workflow for running auto-merge

version 1.0

import "tasks/sl_tasks.wdl" as sl_tasks
import "sl_dataset_reports.wdl"

task reheader_bams {
  input {
    File dataset
    String? biosample_name
    String? library_name

    Int nproc = 1
    String log_level = "DEBUG"
  }
  # note that although the library and biosample names are optional, Cromwell
  # has a bug that breaks quoting them in interpolations, so we have to pass
  # them always and the task itself treats empty strings as equivalent to None
  command <<<
    EXT=$(python3 -c "print('.'.join('~{dataset}'.split('.')[-2:]))")
    python3 \
      -m pbcoretools.tasks.reheader_bams \
      --log-level ~{log_level} \
      ~{dataset} \
      updated.${EXT} \
      --library-name "~{library_name}" \
      --biosample-name "~{biosample_name}"
  >>>
  runtime {
    cpu: 1
  }
  output {
    File updated = glob("*set.xml")[0]
  }
}

task dataset_update_counts {
  input {
    File dataset

    Int nproc = 1
    String log_level = "DEBUG"
  }

  command {
    dataset \
      --log-level ${log_level} \
      newuuid --updateCounts \
      ${dataset}
  }
  runtime {
    cpu: 1
  }
  output {
    File updated = "${dataset}"
  }
}

workflow sl_copy_dataset {
  input {
    File dataset
    Boolean make_reports = false
    String? biosample_name
    String? library_name

    Int nproc = 1
    String log_level = "DEBUG"
    String? tmp_dir
  }

  if (defined(biosample_name) || defined(library_name)) {
    call reheader_bams {
      input:
        dataset = dataset,
        biosample_name = biosample_name,
        library_name = library_name,
        nproc = 1,
        log_level = log_level
    }
  }

  if (!defined(reheader_bams.updated)) {
    call dataset_update_counts {
      input:
        dataset = dataset,
        nproc = 1,
        log_level = log_level
    }
  }

  Array[File?] output_xmls = [reheader_bams.updated, dataset_update_counts.updated]
  File output_xml = select_first(output_xmls)

  if (make_reports) {
    call sl_dataset_reports.import_dataset_reports {
      input:
        dataset_xml = output_xml,
        nproc = 1,
        log_level = log_level
    }
  }

  output {
    File updated = output_xml
    File? report_raw_data = import_dataset_reports.report_raw_data
    File? report_adapters = import_dataset_reports.report_adapter
    File? report_loading = import_dataset_reports.report_loading
    File? report_control = import_dataset_reports.report_control
    File? report_subread_stats = import_dataset_reports.report_subread_stats
    File? report_ccs2 = import_dataset_reports.report_ccs2
  }
}
