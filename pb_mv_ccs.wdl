# Version of the Minor Variants workflow that takes (barcoded) CCS as input.

version 1.0

import "tasks/pbcoretools.wdl" as pbcoretools
import "wf_minorvariants.wdl"

workflow pb_mv_ccs {
  input {
    File eid_ccs
    File eid_ref_dataset
    String dataset_filters = ""
    Int filter_min_qv = 20
    String juliet_target_config = "none"
    String? juliet_target_config_override
    Boolean juliet_mode_phasing = true
    String juliet_genomic_region = ""
    Boolean juliet_only_known_drms = false
    Float juliet_minimal_percentage = 0.1
    Float juliet_maximal_percentage = 100.0
    Float juliet_substitution_rate = 0.0
    Float juliet_deletion_rate = 0.0
    Boolean juliet_debug = false

    Int nproc
    # this is actually ignored internally due to the way barcoded samples
    # are handled
    Int max_nchunks = 100
    String log_level = "DEBUG"
    String? tmp_dir
  }

  call pbcoretools.dataset_filter {
    input:
      dataset = eid_ccs,
      dataset_ext = ".consensusreadset.xml",
      filters = dataset_filters,
      min_qv = filter_min_qv,
      nproc = 1,
      log_level = log_level
  }

  if ((juliet_target_config == "none") &&
      (defined(juliet_target_config_override)) &&
      (juliet_target_config_override != "")) {
    String? target_config_override = juliet_target_config_override
  }
  String target_config_final = select_first([target_config_override,
                                             juliet_target_config])

  call wf_minorvariants.minorvariants {
    input:
      ccs = dataset_filter.filtered,
      reference = eid_ref_dataset,
      juliet_target_config = target_config_final,
      juliet_mode_phasing = juliet_mode_phasing,
      juliet_genomic_region = juliet_genomic_region,
      juliet_only_known_drms = juliet_only_known_drms,
      juliet_minimal_percentage = juliet_minimal_percentage,
      juliet_maximal_percentage = juliet_maximal_percentage,
      juliet_substitution_rate = juliet_substitution_rate,
      juliet_deletion_rate = juliet_deletion_rate,
      juliet_debug = juliet_debug,
      nproc = nproc,
      log_level = log_level,
      tmp_dir = tmp_dir
  }

  output {
    File juliet_html = minorvariants.juliet_html
    File report_minor_variants = minorvariants.report_minor_variants
    File report_csv = minorvariants.report_csv
    File mapped = minorvariants.mapped
  }
}
