# whatshap phasing following DeepVariant

# TODO decide on deployment strategy

version 1.0

# this is tailored to work with the hg38 reference, where we are usually only
# interested in the first 25 contigs in the FASTA (and that's where the reads
# tend to map), but should handle other highly contiguous assemblies too.
# very fragmented assemblies will run as a single chunk because grouping
# chromosomes is messier
task split_by_chromosome {
  input {
    File variants_vcf
    File reference_fasta
    File reference_fasta_fai
    Int add_memory_mb = 0

    Int nproc = 1
    String log_level = "INFO"
  }
  Int max_nchunks = 50
  command {
    set -e
    ln -s ${reference_fasta} reference.fasta
    ln -s ${reference_fasta_fai} reference.fasta.fai
    python3 <<EOF
    import vcf
    with open("${variants_vcf}", "rb") as vcf_in:
      variants_vcf = vcf.Reader(vcf_in)
      variant_chroms = set((rec.CHROM for rec in variants_vcf))
    max_nchunks = int("${max_nchunks}")
    with open("chromosome_args.txt", "wt") as chunks_txt:
      if len(variant_chroms) > 1 and len(variant_chroms) < max_nchunks:
        chunks_txt.write("\n".join(list(variant_chroms)))
      else:
        chunks_txt.write("ALL")
    EOF
  }
  Int total_mem_mb = add_memory_mb + 512
  runtime {
    cpu: 1
    memory: "${total_mem_mb}MB"
  }
  output {
    Array[String] chromosome_list = read_lines("chromosome_args.txt")
  }
}

task tabix_extract {
  input {
    File variants_vcf
    File variants_vcf_tbi
    String chromosome
    Int add_memory_mb = 0

    Int nproc = 1
    String log_level = "INFO"
  }
  String output_vcf = "region.deepvariant.vcf"
  command {
    if [ "${chromosome}" = "ALL" ] || [ "${chromosome}" = "" ]; then
      ln -s ${variants_vcf} ${output_vcf}.gz
      ln -s ${variants_vcf_tbi} ${output_vcf}.gz.tbi
      echo "" > chromosome_args.txt
    else
      ln -s ${variants_vcf} in.vcf.gz
      ln -s ${variants_vcf_tbi} in.vcf.gz.tbi
      tabix -h in.vcf.gz ${chromosome} > ${output_vcf}
      bgzip ${output_vcf}
      tabix ${output_vcf}.gz
      echo "--chromosome ${chromosome}" > chromosome_args.txt
    fi
  }
  Int total_mem_mb = add_memory_mb + 1024
  runtime {
    cpu: 1
    memory: "${total_mem_mb}MB"
  }
  output {
    File region_variants_vcf = "${output_vcf}.gz"
    File region_variants_vcf_tbi = "${output_vcf}.gz.tbi"
    String whatshap_args = read_string("chromosome_args.txt")
  }
}

task whatshap_phase {
  input {
    File reference_fasta
    File reference_fasta_fai
    File variants_vcf
    File variants_vcf_tbi
    Array[File] mapped_bams
    Array[File] mapped_bams_bai
    String docker_image
    String whatshap_bin_path
    String? whatshap_overrides
    String? chromosome_args
    Int add_memory_mb = 0

    Int nproc = 1
    String log_level = "INFO"
  }
  command <<<
    set -e
    ln -s ~{reference_fasta} reference.fasta
    ln -s ~{reference_fasta_fai} reference.fasta.fai
    ln -s ~{variants_vcf} deepvariant.vcf.gz
    ln -s ~{variants_vcf_tbi} deepvariant.vcf.gz.tbi
    mkdir -p inputs
    for bam_file in "~{sep=" " mapped_bams}"; do
      ln -s ${bam_file} inputs/
    done
    for bai_file in "~{sep=" " mapped_bams_bai}"; do
      ln -s ${bai_file} inputs/
    done
    ~{whatshap_bin_path}/whatshap phase \
      --indels \
      ~{whatshap_overrides} \
      ~{chromosome_args} \
      --output deepvariant.phased.vcf \
      --reference reference.fasta \
      deepvariant.vcf.gz \
      inputs/mapped*.bam > whatshap-phase.log
  >>>
  Int total_mem_mb = add_memory_mb + 4096
  runtime {
    cpu: 1
    memory: "${total_mem_mb}MB"
    docker: "${docker_image}"
  }
  output {
    File phased_variants_vcf = "deepvariant.phased.vcf"
  }
}

task bcftools_concat {
  input {
    Array[File] variants_vcf
    String bcftools_params = "-a -Oz"
    Int add_memory_mb = 0

    Int nproc = 1
    String log_level = "INFO"
  }
  String vcf_out = "deepvariant.phased.vcf.gz"
  command <<<
    set -e
    VCF_FILES=$(for fn in ~{sep=" " variants_vcf}; do readlink -f $fn; done)
    echo "true VCF paths: $VCF_FILES"
    bcftools concat ~{bcftools_params} -o ~{vcf_out} \
      ${VCF_FILES}
    tabix ~{vcf_out}
  >>>
  Int total_mem_mb = add_memory_mb + 128
  runtime {
    cpu: 1
    memory: "${total_mem_mb}MB"
  }
  output {
    File phased_variants_vcf = "${vcf_out}"
    File phased_variants_vcf_tbi = "${vcf_out}.tbi"
  }
}

task whatshap_haplotag {
  input {
    File reference_fasta
    File reference_fasta_fai
    File phased_variants_vcf
    File phased_variants_vcf_tbi
    File mapped_bam
    File mapped_bam_bai
    String docker_image
    String whatshap_bin_path
    String? whatshap_overrides
    Int add_memory_mb = 0

    Int nproc = 1
    String log_level = "INFO"
  }
  command {
    set -e
    ln -s ${reference_fasta} reference.fasta
    ln -s ${reference_fasta_fai} reference.fasta.fai
    ln -s ${phased_variants_vcf} deepvariant.phased.vcf.gz
    ln -s ${phased_variants_vcf_tbi} deepvariant.phased.vcf.gz.tbi
    ln -s ${mapped_bam} mapped.bam
    ln -s ${mapped_bam_bai} mapped.bam.bai
    ${whatshap_bin_path}/whatshap haplotag \
      ${whatshap_overrides} \
      --output-threads ${nproc} \
      --output deepvariant.happlotagged.bam \
      --reference reference.fasta \
      deepvariant.phased.vcf.gz mapped.bam > whatshap-haplotag.log
    ${whatshap_bin_path}/python3 -c "import pysam; pysam.index(\"deepvariant.happlotagged.bam\")"
  }
  Int total_mem_mb = add_memory_mb + 4096
  runtime {
    cpu: 1
    memory: "${total_mem_mb}MB"
    docker: "${docker_image}"
  }
  output {
    File haplotagged_bam = "deepvariant.happlotagged.bam"
    File haplotagged_bam_bai = "deepvariant.happlotagged.bam.bai"
  }
}

task whatshap_stats {
  input {
    File phased_variants_vcf
    File phased_variants_vcf_tbi
    String docker_image
    String whatshap_bin_path
    Int add_memory_mb = 0

    Int nproc = 1
    String log_level = "INFO"
  }
  String output_gtf = "deepvariant.phased.gtf"
  String output_tsv = "deepvariant.phased.tsv"
  String output_blocklist = "deepvariant.phased.blocklist"
  command {
    set -e
    ln -s ${phased_variants_vcf} variants.vcf.gz
    ln -s ${phased_variants_vcf_tbi} variants.vcf.gz.tbi
    ${whatshap_bin_path}/whatshap stats \
      --gtf ${output_gtf} \
      --tsv ${output_tsv} \
      --block-list ${output_blocklist} \
      variants.vcf.gz > whatshap-stats.log
  }
  Int total_mem_mb = add_memory_mb + 4096
  runtime {
    cpu: 1
    memory: "${total_mem_mb}MB"
    docker: "${docker_image}"
  }
  output {
    File stats_gtf = "${output_gtf}"
    File stats_tsv = "${output_tsv}"
    File blocklist_txt = "${output_blocklist}"
  }
}

task tabix_index {
  input {
    File phased_variants_vcf
    Int nproc = 1
    Int add_memory_mb = 0
  }
  command {
    ln ${phased_variants_vcf} deepvariant.phased.vcf
    bgzip deepvariant.phased.vcf
    tabix deepvariant.phased.vcf.gz
  }
  Int total_mem_mb = add_memory_mb + 128
  runtime {
    cpu: 1
    memory: "${total_mem_mb}MB"
  }
  output {
    File indexed_variants_vcf = "deepvariant.phased.vcf.gz"
    File indexed_variants_vcf_tbi = "deepvariant.phased.vcf.gz.tbi"
  }
}

workflow wf_whatshap {
  input {
    Array[File] mapped_bams
    Array[File] mapped_bams_bai
    File variants_vcf
    File variants_vcf_tbi
    File reference_fasta
    File reference_fasta_fai
    Int add_memory_mb = 0

    String docker_image = "quay.io/biocontainers/whatshap:1.4--py39hc16433a_1"
    String whatshap_bin_path = "/usr/local/bin"

    Int nproc = 1
    String log_level = "INFO"
  }


  call split_by_chromosome {
    input:
      reference_fasta = reference_fasta,
      reference_fasta_fai = reference_fasta_fai,
      variants_vcf = variants_vcf,
      add_memory_mb = add_memory_mb
  }

  scatter (chromosome in split_by_chromosome.chromosome_list) {
    call tabix_extract {
      input:
        variants_vcf = variants_vcf,
        variants_vcf_tbi = variants_vcf_tbi,
        chromosome = chromosome,
        nproc = nproc,
        add_memory_mb = add_memory_mb
    }

    call whatshap_phase {
      input:
        reference_fasta = reference_fasta,
        reference_fasta_fai = reference_fasta_fai,
        variants_vcf = tabix_extract.region_variants_vcf,
        variants_vcf_tbi = tabix_extract.region_variants_vcf_tbi,
        mapped_bams = mapped_bams,
        mapped_bams_bai = mapped_bams_bai,
        chromosome_args = tabix_extract.whatshap_args,
        docker_image = docker_image,
        whatshap_bin_path = whatshap_bin_path,
        add_memory_mb = add_memory_mb,
        nproc = nproc
    }

    call tabix_index {
      input:
        phased_variants_vcf = whatshap_phase.phased_variants_vcf,
        add_memory_mb = add_memory_mb,
        nproc = 1
    }
  }

  call bcftools_concat {
    input:
      variants_vcf = tabix_index.indexed_variants_vcf,
      add_memory_mb = add_memory_mb,
      nproc = nproc
  }

  # more than this is unproductive
  Int NPROC_MAX = 4
  scatter (i_chunk in range(length(mapped_bams))) {
    call whatshap_haplotag {
      input:
        reference_fasta = reference_fasta,
        reference_fasta_fai = reference_fasta_fai,
        phased_variants_vcf = bcftools_concat.phased_variants_vcf,
        phased_variants_vcf_tbi = bcftools_concat.phased_variants_vcf_tbi,
        mapped_bam = mapped_bams[i_chunk],
        mapped_bam_bai = mapped_bams_bai[i_chunk],
        docker_image = docker_image,
        whatshap_bin_path = whatshap_bin_path,
        add_memory_mb = add_memory_mb,
        nproc = if (nproc > NPROC_MAX) then NPROC_MAX else nproc
    }
  }

  call whatshap_stats {
    input:
      phased_variants_vcf = bcftools_concat.phased_variants_vcf,
      phased_variants_vcf_tbi = bcftools_concat.phased_variants_vcf_tbi,
      docker_image = docker_image,
      whatshap_bin_path = whatshap_bin_path,
      add_memory_mb = add_memory_mb,
      nproc = nproc
  }

  output {
    File phased_variants_vcf = bcftools_concat.phased_variants_vcf
    File phased_variants_vcf_tbi = bcftools_concat.phased_variants_vcf_tbi
    Array[File] haplotagged_bams = whatshap_haplotag.haplotagged_bam
    Array[File] haplotagged_bams_bai = whatshap_haplotag.haplotagged_bam_bai
    File phasing_stats_gtf = whatshap_stats.stats_gtf
    File phasing_stats_tsv = whatshap_stats.stats_tsv
  }
}
