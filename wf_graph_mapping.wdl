version 1.0

import "tasks/pbreports.wdl" as pbreports
import "wf_coverage_reports.wdl"

task graph_mapping {
    input {
        File reference
        File gfa
        File reads
        String out_fn
        String param_overrides
        Int nproc
        Float min_idt
        Int min_length
    }
    command <<<
        graph_param="-g ~{gfa}"
        if [ -z "~{gfa}" ]; then
            graph_param=""
        fi

        # Downstream Pbcore tools cannot work with unmapped alignments, so we need to filter them out here.
        raptor -t ~{nproc} --out-fmt bam -r ~{reference} ${graph_param} -q ~{reads} --min-idt ~{min_idt} --min-qlen ~{min_length} ~{param_overrides} -o - | samtools view -bh -F 0x104 > ~{out_fn}.bam
        samtools sort --threads ~{nproc} -o ~{out_fn}.sorted.bam ~{out_fn}.bam
        samtools index ~{out_fn}.sorted.bam
        rm -f ~{out_fn}.bam
        pbindex ~{out_fn}.sorted.bam
        dataset create ~{out_fn} --type AlignmentSet ~{out_fn}.sorted.bam
    >>>
    Int total_mem_gb = 4 + (2 * nproc)
    runtime {
        cpu: nproc
        memory: "${total_mem_gb}GB"
    }
    output {
        File mapped = "${out_fn}"
    }
}

workflow RaptorGraphMapping {
    input {
        File reads_xml
        # Raptor can't read referencesets at the moment, just plain sequence files.
        File reference
        # The referenceset is needed for reports.
        File referenceset_xml
        File gfa
        String param_overrides
        String alignment_ext = ".alignmentset.xml"
        String mapping_stats_module = "mapping_stats"
        String coverage_report_module = "coverage"
        Int min_concordance = 70
        Int min_length = 50
        Int index_size_gb

        Int nproc = 16
        Int? max_nchunks = 1
        Int? target_size

        String log_level = "INFO"
    }

    call graph_mapping {
        input:
            reference = reference,
            gfa = gfa,
            reads = reads_xml,
            out_fn = "mapped.alignmentset.xml",
            param_overrides = param_overrides,
            nproc = nproc,
            min_idt = min_concordance,
            min_length = min_length,
    }

    call pbreports.mapping_stats as mapping_stats {
        input:
            mapped = graph_mapping.mapped,
            all_reads = reads_xml,
            report_module = mapping_stats_module,
            min_concordance = min_concordance,
            index_memory_gb = index_size_gb,
            nproc = nproc,
            log_level = log_level
    }

    call wf_coverage_reports.coverage_reports {
        input:
            mapped = graph_mapping.mapped,
            reference = referenceset_xml,
            nproc = nproc,
            log_level = log_level
    }

    output {
        File mapped = graph_mapping.mapped
        File report_mapping_stats = mapping_stats.report
        File coverage_gff = coverage_reports.coverage_gff
        File report_coverage = coverage_reports.report_coverage
    }
}
