rule vsearch_derep:
    input:
        "results/assembly/{sample}_{unit}/{sample}_{unit}.fasta"
    output:
        "results/assembly/{sample}_{unit}/{sample}_{unit}_derep.fasta"
    conda:
        "../envs/vsearch.yaml"
    log:
        "results/logs/{sample}_{unit}/vsearch_derep.log"
    shell:
        "vsearch --derep_fulllength {input} --output {output} --log {log}"

rule vsearch_cluster:
    input:
        "results/assembly/{sample}_{unit}/{sample}_{unit}_derep.fasta"
    output:
        "results/assembly/{sample}_{unit}/{sample}_{unit}_cluster.fasta"
    conda:
        "../envs/vsearch.yaml"
    params:
        id=config["derep"]["clustering_id"]
    threads: config["general"]["cores"]
    log:
        "results/logs/{sample}_{unit}/vsearch_cluster.log"
    shell:
        "vsearch --cluster_fast {input} --consout {output} --id {params.id} --log {log} --threads {threads}"
