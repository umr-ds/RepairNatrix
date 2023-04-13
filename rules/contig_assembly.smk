
rule wrap_fasta:
    input:
        #"results/assembly/{sample}_{unit}/{sample}_{unit}.fastq"#, repaired = CONSTRAINT_REPAIRED, filtered = CONSTRAINT_FILTER)
        "results/assembly/{sample}_{unit}/{sample}_{unit}_derep.fasta" #if config["derep"]["centroid_selection"] == "frequency" else "results/assembly/{sample}_{unit}/{sample}_{unit}_sortedf.fasta"  #"results/assembly/{sample}_{unit}/{sample}_{unit}.dereplicated.fasta"
    output:
        "results/assembly/{sample}_{unit}/{sample}_{unit}.dereplicated.fasta"
    conda:
        "../envs/seqtk.yaml"
    shell:
        "seqtk seq {input} > {output}"


rule data_assembly:
    input:
        derep = "results/assembly/{sample}_{unit}/{sample}_{unit}_repaired.fasta" if config['constraint_filtering']['repair_after_assembly'] else "results/assembly/{sample}_{unit}/{sample}_{unit}.dereplicated.fasta"
    output:
        blast_hits = "results/contig_assembly/{sample}_{unit}/blast_hits",
        data_hits = "results/contig_assembly/{sample}_{unit}/data_reads.fasta",
        contigs = "results/contig_assembly/{sample}_{unit}/spades/contigs.fasta",
    params:
        spades_dir = "results/contig_assembly/{sample}_{unit}/spades",
        blast_db = config["contig_assembly"]["blast_db"],
        sample = "{sample}"
    conda:
        "../envs/contig_assembly.yaml"
    script:
        "../scripts/contig_assembly.py"

rule rename_data_assembly:
    input:
        "results/contig_assembly/{sample}_{unit}/spades/contigs.fasta"
    output:
        "results/contig_assembly/assemblies/{sample}_{unit}_assembled.fasta"
    shell:
        "cp {input} {output}"