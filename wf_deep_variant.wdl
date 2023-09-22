# Single-pass DeepVariant 1.4 workflow

version 1.0

task pb_collect_inputs {
  input {
    File eid_ref_dataset
    Int max_nchunks = 1
    Int genome_chunk_size = 10000000
  }
  command {
    python3 <<EOF
    import shutil
    import os.path as op
    from pbcore.io import openDataSet
    # copy the entire reference genome over so it can be hard-linked later
    ref = openDataSet(op.realpath("${eid_ref_dataset}"))
    fasta_name = ref.externalResources[0].resourceId
    fasta_fai_name = ref.externalResources[0].resourceId + ".fai"
    nbytes = op.getsize(fasta_name)
    nchunks_genome = max(1, nbytes // ${genome_chunk_size})
    with open("nchunks_actual.txt", "wt") as nchunks_out:
      # this allows us to test parallel execution on tiny genomes
      nchunks_out.write(str(min(${max_nchunks}, nchunks_genome)))
    shutil.copyfile(fasta_name, op.basename(fasta_name))
    shutil.copyfile(fasta_fai_name, op.basename(fasta_fai_name))
    EOF
  }
  runtime {
    cpu: 1
    memory: "500MB"
  }
  output {
    File ref_fasta_name = select_first(glob("*.fasta"))
    File ref_fasta_index = select_first(glob("*.fasta.fai"))
    Int nchunks_actual = read_int("nchunks_actual.txt")
  }
}

task deepvariant_make_examples {
  input {
    Array[File] mapped_bams
    Array[File] mapped_bam_bais
    File ref_fasta
    File ref_fasta_index
    Int i_chunk
    Int nchunks
    String docker_image
    Int base_memory_mb = 1024

    String log_level = "INFO"
    Int nproc = 1
  }
  Int dv_partition_size = 25000
  Int dv_max_reads_per_partition = 600
  # this is a property of the model, not mutable
  Int dv_pileup_image_width = 199
  Float dv_vsc_min_fraction_indels = 0.12
  command <<<
    set -e
    mkdir -p inputs
    for bam in ~{sep=" " mapped_bams}; do
      ln -s $bam inputs/
    done
    for bai in ~{sep=" " mapped_bam_bais}; do
      ln -s $bai inputs/
    done
    ln -s ~{ref_fasta} reference.fasta
    ln -s ~{ref_fasta_index} reference.fasta.fai
    READS=$(find inputs -name "*.bam" | tr '\n' ',' | sed 's/,$//;')
    echo "Input reads: $READS"
    /opt/deepvariant/bin/make_examples \
      --add_hp_channel \
      --alt_aligned_pileup=diff_channels \
      --min_mapping_quality=1 \
      --parse_sam_aux_fields \
      --partition_size=~{dv_partition_size} \
      --max_reads_per_partition=~{dv_max_reads_per_partition} \
      --phase_reads \
      --pileup_image_width ~{dv_pileup_image_width} \
      --norealign_reads \
      --sort_by_haplotypes \
      --track_ref_reads \
      --vsc_min_fraction_indels ~{dv_vsc_min_fraction_indels} \
      --mode calling \
      --ref reference.fasta \
      --reads $READS \
      --examples tfrecord@~{nchunks}.gz \
      --gvcf gvcf.tfrecord@~{nchunks}.gz \
      --task ~{i_chunk} > make_examples.out 2>&1
  >>>
  Int total_memory_mb = base_memory_mb + 8192
  runtime {
    cpu: 1
    docker: "${docker_image}"
    memory: "${total_memory_mb}MB"
  }
  output {
    # this gets flattened later
    Array[File] examples_gz = glob("tfrecord*.gz")
    Array[File] gvcf_gz = glob("gvcf.tfrecord*.gz")
  }
}

task deepvariant_call_variants {
  input {
    Array[File] examples_gz
    Int nchunks
    String docker_image
    Int base_memory_mb = 1024

    # Tensorflow behavior is weird: we can't use more than one GPU, but we also
    # can't stop the CPU version from taking every CPU regardless of how many
    # were allocated to the task
    Int nproc = 1
    Int ngpu = 0
  }
  String outfile = "call_variants_output.tfrecord.gz"
  command <<<
    set -e
    mkdir -p inputs
    # FIXME is there a better way to handle these?
    for gz_file in ~{sep=" " examples_gz}; do
      ln -s $gz_file inputs/
    done
    /opt/deepvariant/bin/call_variants \
      --outfile "~{outfile}" \
      --examples "inputs/tfrecord@~{nchunks}.gz" \
      --checkpoint /opt/models/pacbio/model.ckpt > call_variants.out 2>&1
  >>>
  Int total_memory_mb = base_memory_mb + 8192
  runtime {
    cpu: nproc
    gpu: ngpu
    docker: "${docker_image}"
    memory: "${total_memory_mb}MB"
  }
  output {
    File tfrecord = "${outfile}"
  }
}

task deepvariant_postprocess_variants {
  input {
    File tfrecord
    File ref_fasta
    File ref_fasta_index
    Array[File] gvcf_chunks
    String docker_image
    Int base_memory_mb = 1024
    Int nchunks
    Int nproc = 1
  }
  command <<<
    set -e
    mkdir -p inputs
    ln -s ~{ref_fasta} inputs/reference.fasta
    ln -s ~{ref_fasta_index} inputs/reference.fasta.fai
    for gvcf_file in ~{sep=" " gvcf_chunks}; do
      ln -s $gvcf_file inputs/
    done
    /opt/deepvariant/bin/postprocess_variants \
      --ref inputs/reference.fasta \
      --infile ~{tfrecord} \
      --nonvariant_site_tfrecord_path inputs/gvcf.tfrecord@~{nchunks}.gz \
      --outfile deepvariant.vcf.gz \
      --gvcf_outfile deepvariant.g.vcf.gz > postprocess_variants.out 2>&1
  >>>
  Int total_memory_mb = base_memory_mb + (32 * 1024)
  runtime {
    cpu: 1
    docker: "${docker_image}"
    memory: "${total_memory_mb}MB"
  }
  output {
    File variants_vcf = "deepvariant.vcf.gz"
    File variants_vcf_index = "deepvariant.vcf.gz.tbi"
    File variants_gvcf = "deepvariant.g.vcf.gz"
  }
}

task bcftools_stats {
  input {
    File variants_vcf
    Int base_memory_mb = 1024

    Int nproc = 1
    String log_level = "INFO"
  }
  command {
    bcftools stats \
      --threads ${nproc} \
      --apply-filters PASS \
      ${variants_vcf} > deepvariant.vcf.stats.txt
  }
  runtime {
    cpu: nproc
    memory: "${base_memory_mb}MB"
  }
  output {
    File stats_txt = "deepvariant.vcf.stats.txt"
  }
}

workflow wf_deep_variant {
  input {
    File eid_ref_dataset
    Array[File] mapped_bams
    Array[File] mapped_bams_bai
    String docker_image = "google/deepvariant:1.5.0"
    Boolean enable_gpu = false

    Int nproc = 1
    Int max_nchunks = 2
    Int base_memory_mb = 1024
    String log_level = "INFO"
    String? tmp_dir
  }

  call pb_collect_inputs {
    input:
      eid_ref_dataset = eid_ref_dataset,
      max_nchunks = max_nchunks
  }

  scatter (i_chunk in range(pb_collect_inputs.nchunks_actual)) {
    call deepvariant_make_examples {
      input:
        mapped_bams = mapped_bams,
        mapped_bam_bais = mapped_bams_bai,
        ref_fasta = pb_collect_inputs.ref_fasta_name,
        ref_fasta_index = pb_collect_inputs.ref_fasta_index,
        i_chunk = i_chunk,
        nchunks = pb_collect_inputs.nchunks_actual,
        docker_image = docker_image,
        base_memory_mb = base_memory_mb,
        nproc = nproc
    }
  }

  # relying on Google to be consistent here:
  # https://hub.docker.com/r/google/deepvariant/tags
  String docker_image_call = if (enable_gpu) then "${docker_image}-gpu" else docker_image
  Int number_of_gpus = if (enable_gpu) then 1 else 0
  Int nproc_call_variants = if (enable_gpu) then 2 else nproc

  call deepvariant_call_variants {
    input:
      examples_gz = flatten(deepvariant_make_examples.examples_gz),
      nchunks = pb_collect_inputs.nchunks_actual,
      docker_image = docker_image_call,
      base_memory_mb = base_memory_mb,
      nproc = nproc_call_variants,
      ngpu = number_of_gpus
  }

  call deepvariant_postprocess_variants {
    input:
      tfrecord = deepvariant_call_variants.tfrecord,
      ref_fasta = pb_collect_inputs.ref_fasta_name,
      ref_fasta_index = pb_collect_inputs.ref_fasta_index,
      gvcf_chunks = flatten(deepvariant_make_examples.gvcf_gz),
      nchunks = pb_collect_inputs.nchunks_actual,
      docker_image = docker_image,
      base_memory_mb = base_memory_mb,
      nproc = nproc,
  }

  call bcftools_stats {
    input:
      variants_vcf = deepvariant_postprocess_variants.variants_vcf,
      base_memory_mb = base_memory_mb,
      nproc = nproc
  }

  output {
    File variants_vcf = deepvariant_postprocess_variants.variants_vcf
    File variants_vcf_index = deepvariant_postprocess_variants.variants_vcf_index
    File variants_gvcf = deepvariant_postprocess_variants.variants_vcf
    File bcftools_stats_txt = bcftools_stats.stats_txt
    File reference_fasta = pb_collect_inputs.ref_fasta_name
    File reference_fasta_fai = pb_collect_inputs.ref_fasta_index
  }
}
