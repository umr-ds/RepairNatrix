rule constraint_repair:
    input:
        # "results/assembly/{sample}_{unit}/{sample}_{unit}_assembled.fastq"
        #derep="results/assembly/{sample}_{unit}/{sample}_{unit}_derep.fasta",
        derep="results/assembly/{sample}_{unit}/{sample}_{unit}.dereplicated.fasta",
        cluster="results/assembly/{sample}_{unit}/{sample}_{unit}_cluster.fasta"
    output:
        "results/assembly/{sample}_{unit}/{sample}_{unit}_repaired_cluster.fasta" # "repaired" bases will have a lower quality score (and will be APPENDED - or replace the original ones?)
        "results/assembly/{sample}_{unit}/{sample}_{unit}_repair_mapping.json"  # mapping repaired reads to original reads

    params:
        sequence_length=config["constraint_filtering"]["sequence_length"],
        min_gc_content=config["constraint_filtering"]["gc_content"]["gc_min"],
        max_gc_content=config["constraint_filtering"]["gc_content"]["gc_max"],
        max_homopolymer_length=config["constraint_filtering"]["homopolymer"]["count"],
        illegal_sequences=config["constraint_filtering"]["undesired_subsequences"]["file"],
        #todo add kmer checks to constraint_filtering code:
        kmer_counting_k=config["constraint_filtering"]["kmer_counting"]["k"],
        kmer_counting_upper_bound=config["constraint_filtering"]["kmer_counting"]["upper_bound"]

    conda:
        "../envs/constraint_repair.yaml"
    script:
        "../scripts/constraints/check_constraints.py"
