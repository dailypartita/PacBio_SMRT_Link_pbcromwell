# This is the middle layer of the minor variants workflows, starting from a
# CCS dataset that contains multiple barcoded samples, each of which needs
# to be analyzed separately.

version 1.0

import "wf_julietflow.wdl"
import "tasks/chunking.wdl" as chunking

task gather_html_zip {
  input {
    Array[File] html_files

    Int nproc = 1
    String log_level = "DEBUG"
  }
  command {
    cp ${sep=" " html_files} .
    tar -cvzf minor_variants_reports.zip *.html
  }
  runtime {
    cpu: 1
    backend: "Local"
    memory: "100MB"
  }
  output {
    File gathered = "minor_variants_reports.zip"
  }
}

# utility to wrap the dictionary output by juliet with the barcoded sample ID,
# to allow gathering multiple samples into a single file
# adapted from pysiv2.tasks.minor_variants
task export_juliet_outputs {
  input {
    File dataset
    File json
    File html
    String log_level = "INFO"
  }
  command <<<
    python3 <<EOF
    from pbcore.io import openDataSet
    import logging
    import shutil
    import json
    import os.path
    logging.basicConfig()
    ds_file = os.path.realpath("~{dataset}")
    json_file = "~{json}"
    html_file = "~{html}"
    with openDataSet(ds_file, strict=True) as ds:
      if ds.isBarcoded:
        bcs = list(set(zip(ds.index.bcForward, ds.index.bcReverse)))
        assert len(bcs) == 1
        label = "{f}--{r}".format(f=bcs[0][0], r=bcs[0][1])
        dict_key = label
        logging.info("Sample is '{s}'".format(s=label))
      else:
        label = "all"
        dict_key = "NA"
      json_output_file = "mv.{s}.json".format(s=label)
      html_output_file = "mv.{s}.html".format(s=label)
      with open(json_file) as json_in:
        d = json.loads(json_in.read())
        with open(json_output_file, "w") as json_out:
          json_out.write(json.dumps({dict_key: d}))
      shutil.copyfile(html_file, html_output_file)
    EOF
  >>>
  runtime {
    cpu: 1
    memory: "2GB"
  }
  output {
    File juliet_json = glob("mv.*.json")[0]
    File juliet_html = glob("mv.*.html")[0]
  }
}

task pbreports_minor_variants {
  input {
    File juliet_json

    Int nproc = 1
    String log_level = "INFO"
  }
  command {
    python3 \
      -m pbreports.report.minor_variants \
      --log-level ${log_level} \
      ${juliet_json} \
      minor_variants.report.json \
      minor_variants.report.csv
  }
  runtime {
    cpu: 1
    memory: "1GB"
  }
  output {
    File report = "minor_variants.report.json"
    File csv = "minor_variants.report.csv"
  }
}

workflow minorvariants {
  input {
    File ccs
    File reference
    String juliet_target_config = "none"
    Boolean juliet_mode_phasing = true
    String juliet_genomic_region = ""
    Boolean juliet_only_known_drms = false
    Float juliet_minimal_percentage = 0.1
    Float juliet_maximal_percentage = 100.0
    Float juliet_substitution_rate = 0.0
    Float juliet_deletion_rate = 0.0
    Boolean juliet_debug = false

    Int nproc
    String log_level = "DEBUG"
    String? tmp_dir
  }

  call chunking.split_dataset {
    input:
      ds_in = ccs,
      target_size = 1,
      split_mode = "barcodes",
      log_level = log_level
  }

  scatter (barcoded_sample in split_dataset.chunks) {
    call wf_julietflow.julietflow {
      input:
        ccs = barcoded_sample,
        reference = reference,
        target_config = juliet_target_config,
        mode_phasing = juliet_mode_phasing,
        genomic_region = juliet_genomic_region,
        only_known_drms = juliet_only_known_drms,
        minimal_percentage = juliet_minimal_percentage,
        maximal_percentage = juliet_maximal_percentage,
        substitution_rate = juliet_substitution_rate,
        deletion_rate = juliet_deletion_rate,
        juliet_debug = juliet_debug,
        nproc = nproc,
        log_level = log_level,
        tmp_dir = tmp_dir
    }

    call export_juliet_outputs {
      input:
        dataset = barcoded_sample,
        json = julietflow.juliet_json,
        html = julietflow.juliet_html,
        log_level = log_level
    }
  }

  call gather_html_zip {
    input:
      html_files = export_juliet_outputs.juliet_html
  }

  call chunking.gather_generic as gather_json {
    input:
      chunks = export_juliet_outputs.juliet_json,
      output_file_name = "minor_variants_data.json",
      log_level = log_level
  }

  call chunking.gather_datasets_lite as gather_alignments {
    input:
      chunks = julietflow.mapped,
      dataset_ext = ".consensusalignmentset.xml",
      dataset_name = "\"Minor Variants Alignments\"",
      log_level = log_level
  }

  call pbreports_minor_variants {
    input:
      juliet_json = gather_json.gathered,
      log_level = log_level
  }

  output {
    File juliet_html = gather_html_zip.gathered
    File report_minor_variants = pbreports_minor_variants.report
    File report_csv = pbreports_minor_variants.csv
    File mapped = gather_alignments.gathered
  }
}
