if config["derep"]["centroid_selection"] == "frequency":
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
            cent = "results/assembly/{sample}_{unit}/{sample}_{unit}_cluster.fasta"
        conda:
            "../envs/vsearch.yaml"
        params:
            id=config["derep"]["clustering_id"]
        threads: config["general"]["cores"]
        log:
            "results/logs/{sample}_{unit}/vsearch_cluster.log"
        shell:
            "vsearch --cluster_size {input} --centroids {output.cent} --id {params.id} --clusters results/assembly/{wildcards.sample}_{wildcards.unit}/clust  --log {log} --threads {threads}"
            #"vsearch --cluster_fast {input} --consout {output} --id {params.id} --log {log} --threads {threads}"

elif config["derep"]["centroid_selection"] == "quality":
    rule derep:
        input:
            expand("results/assembly/{{sample}}_{{unit}}/{{sample}}_{{unit}}_assembled{filtered}.fastq",filtered=CONSTRAINT_FILTER)
        output:
            "results/assembly/{sample}_{unit}/{sample}_{unit}_derep.fastq"
        conda:
            "../envs/derep_cluster_qual.yaml"
        script:
            "../scripts/derep_quality.py"

    rule sorting:
        input:
            "results/assembly/{sample}_{unit}/{sample}_{unit}_derep.fastq"
        output:
            "results/assembly/{sample}_{unit}/{sample}_{unit}_sorted.fastq"
        conda:
            "../envs/fastq-tools.yaml"
        shell:
            "fastq-sort --mean-qual -r {input} > {output}"

    rule sorted_fasta:
        input:
            "results/assembly/{sample}_{unit}/{sample}_{unit}_sorted.fastq"
        output:
            "results/assembly/{sample}_{unit}/{sample}_{unit}_derep.fasta"
        conda:
            "../envs/seqtk.yaml"
        shell:
            "seqtk seq -a {input} > {output}"

    rule vsearch_cluster:
        input:
            "results/assembly/{sample}_{unit}/{sample}_{unit}_derep.fasta"
        output:
            cent = "results/assembly/{sample}_{unit}/{sample}_{unit}_cluster.fasta"
        conda:
            "../envs/vsearch.yaml"
        params:
            id=config["derep"]["clustering_id"]
        threads: config["general"]["cores"]
        log:
            "results/logs/{sample}_{unit}/vsearch_cluster.log"
        shell:
            "vsearch --usersort --cluster_smallmem {input} --centroid {output.cent} --id {params.id} --clusters results/assembly/{wildcards.sample}_{wildcards.unit}/clust --log {log} --threads {threads}"
            
