def get_fastq(wildcards):
    if not is_single_end(wildcards.sample, wildcards.unit):
        return expand("demultiplexed/{sample}_{unit}_{group}.fastq",
                        group=[1,2], **wildcards)
    return "demultiplexed/{sample}_{unit}_1.fastq".format(**wildcards)


#def get_fastq(wildcards):
#    if not is_single_end(wildcards.sample, wildcards.unit):
#        return expand("demultiplexed/{sample}_{unit}_{group}.fastq",
#                        group=[1,2], **wildcards)
#    return expand("demultiplexed/{sample}_{unit}_1.fastq", **wildcards)


rule define_primer:
    input:
        primer_table=config["general"]["primertable"]
    output:
        "primer_table.csv"
    params:
       paired_end=config["general"]["paired_End"],
       offset=config["qc"]["primer_offset"],
       bar_removed=config["qc"]["barcode_removed"],
       all_removed=config["qc"]["all_primer"]
    conda:
        "../envs/define_primer.yaml"
    script:
        "../scripts/define_primer.py"

rule prinseq:
    input:
        sample=get_fastq
    output:
        expand(
        "results/assembly/{{sample}}_{{unit}}/{{sample}}_{{unit}}_{read}.fastq",
        read=reads)
    params:
        mq = config["qc"]["mq"],
        trim_to_length = config["qc"]["trim_to_length"]
    log:
        "results/logs/{sample}_{unit}/prinseq.log"
    conda:
        "../envs/prinseq.yaml"
    script:
        "../scripts/prinseq.py"

rule assembly:
    input:
        expand(
        "results/assembly/{{sample}}_{{unit}}/{{sample}}_{{unit}}_{read}.fastq",read=reads),
        primer_t="primer_table.csv"
    output:
        "results/assembly/{sample}_{unit}/{sample}_{unit}_assembled.fastq"
    threads: 20
    params:
        paired_end=config["general"]["paired_End"],
        threshold=config["qc"]["threshold"],
        minoverlap=config["qc"]["minoverlap"],
        minlen=config["qc"]["minlen"],
        maxlen=config["qc"]["maxlen"],
        minqual=config["qc"]["minqual"],
        prim_rm=config["qc"]["all_primer"]
    conda:
        "../envs/assembly.yaml"
    log:
        "results/logs/{sample}_{unit}/read_assembly.log"
    script:
        "../scripts/assembly.py"

if config["derep"]["centroid_selection"] == "frequency" or not config["general"]["clustering"]:
    rule copy_to_fasta:
        input:
            "results/assembly/{sample}_{unit}/{sample}_{unit}_assembled.fastq" if not config['constraint_filtering']['after_assembly'] else "results/assembly/{sample}_{unit}/{sample}_{unit}_filtered.fastq"
            #expand("results/assembly/{{sample}}_{{unit}}/{{sample}}_{{unit}}_{read}{repaired}{filtered}.fastq",read=reads ,repaired=CONSTRAINT_REPAIRED_2,filtered=CONSTRAINT_FILTER_2)
        output:
            "results/assembly/{sample}_{unit}/{sample}_{unit}_freq.fasta"
        conda:
            "../envs/seqtk.yaml"
        shell:
            "seqtk seq -a {input} > {output}"