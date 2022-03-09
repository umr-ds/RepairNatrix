rule constraint_repair:
    input:
        # todo
        # rules json
        "results/assembly/{sample}_{unit}/{sample}_{unit}_assembled.fastq"
    output:
        "repaired_seqs.fastq" # "repaired" bases will have a lower quality score and will be APPENDED!
        "only_repaired_seqs.fastq"  # this fastq contains only the repaired bases

    params:
        sequence_length=config["constraint_filtering"]["sequence_length"],
        min_gc_content=config["constraint_repair"]["min_gc_content"],
        max_gc_content=config["constraint_repair"]["max_gc_content"],
        max_homopolymer_length=config["constraint_repair"]["max_homopolymer_length"],
        illegal_sequences=config["constraint_repair"]["illegal_sequences"],

    conda:
        "../envs/constraint_repair.yaml"
    script:
        "../scripts/constraints/check_constraints.py"
