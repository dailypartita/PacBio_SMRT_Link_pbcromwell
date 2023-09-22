# Mock pseudo-mapping workflow.  Unlike the real mapping workflow, it takes
# an AlignmentSet as input, no reference, but we only use this to test the
# igv-files endpoint in SMRT Link services.

version 1.0

import "tasks/pbcoretools.wdl" as pbcoretools

task mock_mapping {
  input {
    File alignments

    String log_level = "INFO"
    Int nproc = 1
  }
  command {
    dataset \
      --log-level ${log_level} \
      create \
      --type AlignmentSet \
      --name "Mock Mapped Reads" \
      mapped.alignmentset.xml \
      "`readlink -f '${alignments}'`"
    dataset newuuid mapped.alignmentset.xml --random
  }
  runtime {
    cpu: 1
    backend: "Local"
  }
  output {
    File mapped = "mapped.alignmentset.xml"
  }
}

workflow dev_mock_mapping {
  input {
    File eid_alignment
    Boolean consolidate_aligned_bam = false

    Int nproc = 1
    String log_level = "INFO"
    String? tmp_dir
    Int max_nchunks = 1
  }

  call mock_mapping {
    input:
      alignments = eid_alignment,
      log_level = log_level
  }

  if (consolidate_aligned_bam) {
    call pbcoretools.auto_consolidate_alignments {
      input:
        mapped = mock_mapping.mapped,
        log_level = log_level
    }
  }

  output {
    File mapped = mock_mapping.mapped
    File? mapped_bam_datastore = auto_consolidate_alignments.datastore
  }
}
