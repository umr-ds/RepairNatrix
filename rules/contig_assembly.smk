rule data_assembly:
    input:
        derep = "results/assembly/{sample}_{unit}/{sample}_{unit}.dereplicated.fasta"
    output:
        blast_hits = "results/data_assembly/{sample}_{unit}/blast_hits",
        data_hits = "results/data_assembly/{sample}_{unit}/data_reads.fasta",
        contigs = "results/data_assembly/{sample}_{unit}/spades/contigs.fasta",
    params:
        spades_dir = "results/data_assembly/{sample}_{unit}/spades",
        blast_db = config["contig_assembly"]["blast_db"],
        sample = "{sample}"
    conda:
        "../envs/contig_assembly.yaml"
    script:
        "../scripts/contig_assembly.py"

rule rename_data_assembly:
    input:
        "results/data_assembly/{sample}_{unit}/spades/contigs.fasta"
    output:
        "results/data_assembly/assemblies/{sample}_{unit}.fasta",
    shell:
        "cp {input} {output}"