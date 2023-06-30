rule constraint_filtering_base:
        input: "overwrite_this"
        output: temp("overwrite_this")
        params:
            constraints = config["constraint_filtering"],
            primer_length = config["constraint_filtering"]["primer_length"],
            #sequence_length = config["constraint_filtering"]["sequence_length"],
            paired = config["general"]["paired_End"]
        log: "logs/overwrite_this.log"
        threads: config["general"]["cores"]
        conda:
            "../envs/constraint_filtering.yaml"
        script:
            "../scripts/constraint_filtering.py"

rule constraint_repair_base:
        input: "overwrite_this"
        output: "overwrite_this"
        params:
            maximum_repair_cycles=config["constraint_filtering"]["maximum_repair_cycles"],
            inplace_repair=config["constraint_filtering"]["inplace_repair"],
            repair_quality_score=config["constraint_filtering"]["repair_quality_score"],
            use_quality_mapping = config["constraint_filtering"]["use_quality_mapping"],
            sequence_length=config["constraint_filtering"]["sequence_length"],
            min_gc_content=config["constraint_filtering"]["overall_gc_content"]["gc_min"],
            max_gc_content=config["constraint_filtering"]["overall_gc_content"]["gc_max"],
            windowed_gc_content_max_gc_window_size = config["constraint_filtering"]["windowed_gc_content"]["window_size"],
            windowed_gc_content_min_gc_content = config["constraint_filtering"]["windowed_gc_content"]["gc_min"],
            windowed_gc_content_max_gc_content = config["constraint_filtering"]["windowed_gc_content"]["gc_max"],
            max_homopolymer_length=config["constraint_filtering"]["homopolymer"]["count"],
            illegal_sequences=config["constraint_filtering"]["undesired_subsequences"]["file"],
            subsequences_file= config['constraint_filtering']['undesired_subsequences']['file'] if 'undesired_subsequences' in config['constraint_filtering']['constraints'] else[],
            kmer_active=config["constraint_filtering"]["kmer_counting"]["active"],
            kmer_k=config["constraint_filtering"]["kmer_counting"]["k"],
            kmer_max_count=config["constraint_filtering"]["kmer_counting"]["upper_bound"]
        conda:
            "../envs/constraint_repair.yaml"
        #shell:
        #    "python3 setup.py install" # && python3 ./scripts/constraints/check_constraints.py"
        script:
            "../scripts/constraints/check_constraints.py"


if config['constraint_filtering']['after_assembly']:
    use rule constraint_filtering_base as constraint_filtering_3 with:
        input:
            expand("results/assembly/{{sample}}_{{unit}}/{{sample}}_{{unit}}_assembled.fastq"),
            subsequences_file=config['constraint_filtering']['undesired_subsequences'][
                'file'] if 'undesired_subsequences' in config['constraint_filtering']['constraints'] else []
        log:
            "results/logs/{sample}_{unit}/constraint_filtering_after_assembly.txt"
        output:
            expand("results/assembly/{{sample}}_{{unit}}/{{sample}}_{{unit}}_filtered.fastq"),
            expand("results/assembly/{{sample}}_{{unit}}/{{sample}}_{{unit}}_filtered_out.fastq")
        params:
            constraints=config["constraint_filtering"],
            primer_length=0,
            sequence_length=0,
            paired=False

if config['constraint_filtering']['repair_after_assembly']:
    #TODO: in this case we want to no only repair invalid centroids but ideally we would want to select:
    # the best representative of the cluster
    # with "best" being: 1) the smallest distance to the current centroid and 2) fulfilling all constraints
    use rule constraint_repair_base as constraint_repair_3 with:
        input:
            "results/assembly/{sample}_{unit}/{sample}_{unit}_assembled.fastq" if not config['constraint_filtering']['after_assembly'] else "results/assembly/{sample}_{unit}/{sample}_{unit}_filtered.fastq",
            #derep="results/assembly/{sample}_{unit}/{sample}_{unit}_derep.fasta",
            cent="results/assembly/{sample}_{unit}/{sample}_{unit}_cluster.fasta"
        output:
            expand("results/assembly/{{sample}}_{{unit}}/{{sample}}_{{unit}}_repaired.fasta"),
            expand("results/assembly/{{sample}}_{{unit}}/{{sample}}_{{unit}}_repaired_clusters.json"),
            #"results/assembly/{sample}_{unit}/{sample}_{unit}_derep_repaired.fastq",# "repaired" bases will have a lower quality score (and will be APPENDED - or replace the original ones?)
            #"results/assembly/{sample}_{unit}/{sample}_{unit}_derep_repair_mapping.json"  # mapping repaired reads to original reads

rule filtered_to_fasta:
    input:
        "results/assembly/{sample}_{unit}/{sample}_{unit}.fastq"
        #expand("results/assembly/{{sample}}_{{unit}}/{{sample}}_{{unit}}_{read}{repaired}{filtered}.fastq",read=reads ,repaired=CONSTRAINT_REPAIRED_2,filtered=CONSTRAINT_FILTER_2)
    output:
        "results/assembly/{sample}_{unit}/{sample}_{unit}.fasta"
    conda:
        "../envs/seqtk.yaml"
    shell:
        "seqtk seq -a {input} > {output}"
