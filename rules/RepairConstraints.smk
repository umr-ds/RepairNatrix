rule repairseqs:
    output:
        temp(expand("demultiplexed/{unit.sample}_{unit.unit}_R{read}.fastq.gz", unit=units.reset_index().itertuples(), read=reads))
        "results/assembly/{sample}_{unit}/{sample}_{unit}_repaired_cluster.fasta"
        "results/assembly/{sample}_{unit}/{sample}_{unit}_repair_mapping.json" # mapping repaired reads to original reads
    input:
        "results/assembly/{sample}_{unit}/{sample}_{unit}_derep.fasta"
        "results/assembly/{sample}_{unit}/{sample}_{unit}_cluster.fasta"
    params:
        filename = config["general"]["filename"],
        primertable = config["general"]["primertable"],
        demultiplexing = config["general"]["demultiplexing"],
        read_sorting = config['general']['read_sorting'],
        assembled = config['general']['already_assembled'],
        name_ext = config['general']['name_ext']
    conda:
        "../envs/constraint_repair.yaml"
    script:
        "../scripts/constraints/check_constraints.py"
