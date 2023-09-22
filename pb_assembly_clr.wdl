version 1.0

import "tasks/assembly.wdl" as ipa

workflow pb_assembly_clr {
    input {
        File eid_subread
        String ctg_prefix = ""
        File? blacklist_fofn
        File config_sh_fn

        Int nproc = 1
        String log_level = "INFO"
        Int max_nchunks = 1
        String tmp_dir = "/tmp"
    }

    String raw_db_prefix = "rawreads"
    String erc_db_prefix = "ercreads"

    call ipa.build_raw_rdb {
        input:
            reads_fn=eid_subread,
            db_prefix=raw_db_prefix,
            blacklist_fofn=blacklist_fofn,
            config_sh_fn=config_sh_fn,
    }

    call ipa.ovlcns_raw_prepare {
        input:
            rdb=build_raw_rdb.out_rdb,
    }

    Array[String] blocks = read_lines(ovlcns_raw_prepare.blocks_fofn)

    scatter (block_id in blocks) {
        call ipa.ovlcns_raw_run {
            input:
                rdb=build_raw_rdb.out_rdb,
                db_prefix=raw_db_prefix,
                length_cutoff_fn=build_raw_rdb.out_cutoff,
                block_id=block_id,
                nproc=nproc,
                config_sh_fn=config_sh_fn,
        }
    }

    call ipa.ovlcns_raw_merge {
        input:
            in_fns = ovlcns_raw_run.out_ercreads_fasta,
    }

    call ipa.build_erc_rdb {
        input:
            reads_fn=ovlcns_raw_merge.gathered_fofn,
            db_prefix=erc_db_prefix,
            config_sh_fn=config_sh_fn,
    }

    call ipa.ovl_erc_prepare {
        input:
            rdb = build_erc_rdb.out_rdb,
    }

    Array[String] blocks2 = read_lines(ovl_erc_prepare.blocks_file)

    scatter (block_id in blocks2) {
        call ipa.ovl_erc_run {
            input:
                rdb=build_erc_rdb.out_rdb,
                block_id=block_id,
                nproc=nproc,
                config_sh_fn=config_sh_fn,
        }
    }

    call ipa.ovl_erc_merge {
        input:
            in_fns=ovl_erc_run.out_filtered_m4,
            db_prefix=erc_db_prefix,
            nproc=nproc,
            config_sh_fn=config_sh_fn,
    }

    call ipa.assemble {
        input:
            in_m4=ovl_erc_merge.gathered_m4,
            in_ercreads=ovlcns_raw_merge.gathered_fasta,
            ctg_prefix=ctg_prefix,
            config_sh_fn=config_sh_fn,
    }

    output {
            File p_ctg_fa = assemble.p_ctg_fa
            File a_ctg_fa = assemble.a_ctg_fa
            File circular_contigs = assemble.circular_contigs
            File raw_rdb_full = build_raw_rdb.out_rdb_full
            File raw_rdb = build_raw_rdb.out_rdb
            File raw_length_cutoff_fn = build_raw_rdb.out_cutoff
            File erc_rdb = build_erc_rdb.out_rdb
            File contig_gfa2 = assemble.contig_gfa2
    }
}
