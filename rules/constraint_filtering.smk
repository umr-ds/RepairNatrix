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

if config['constraint_filtering']['after_demultiplexing']:
    use rule constraint_filtering_base as constraint_filtering_1 with:
        input:
            expand("demultiplexed/{{sample}}_{{unit}}_{read}{repaired}.fastq",read=reads, repaired=CONSTRAINT_REPAIRED_1),
            subsequences_file=config['constraint_filtering']['undesired_subsequences'][
                'file'] if 'undesired_subsequences' in config['constraint_filtering']['constraints'] else []
        log:
            "results/logs/{sample}_{unit}/constraint_filtering_after_demultiplex.txt"
        output:
            temp(expand("demultiplexed/{{sample}}_{{unit}}_{read}{repaired}{filtered}.fastq",read=reads,filtered=CONSTRAINT_FILTER_1,repaired=CONSTRAINT_REPAIRED_2)),
            temp(expand("demultiplexed/{{sample}}_{{unit}}_{read}{repaired}{filtered}_out.fastq",read=reads,filtered=CONSTRAINT_FILTER_1,repaired=CONSTRAINT_REPAIRED_2))

if config['constraint_filtering']['after_quality_control']:
    use rule constraint_filtering_base as constraint_filtering_2 with:
        input:
            expand("results/assembly/{{sample}}_{{unit}}/{{sample}}_{{unit}}_{read}{repaired}.fastq",read=reads,repaired=CONSTRAINT_REPAIRED_2),
            subsequences_file=config['constraint_filtering']['undesired_subsequences'][
                'file'] if 'undesired_subsequences' in config['constraint_filtering']['constraints'] else []
        log:
            "results/logs/{sample}_{unit}/constraint_filtering_after_qc.txt"
        output:
            expand("results/assembly/{{sample}}_{{unit}}/{{sample}}_{{unit}}_{read}{repaired}{filtered}.fastq",read=reads,filtered=CONSTRAINT_FILTER_2,repaired=CONSTRAINT_REPAIRED_2),
            expand("results/assembly/{{sample}}_{{unit}}/{{sample}}_{{unit}}_{read}{repaired}{filtered}_out.fastq",read=reads,filtered=CONSTRAINT_FILTER_2,repaired=CONSTRAINT_REPAIRED_2)

if config['constraint_filtering']['after_assembly']:
    if config["derep"]["clustering"]:
        use rule constraint_filtering_base as constraint_filtering_3 with:
            input:
                expand("results/assembly/{{sample}}_{{unit}}{repaired}.fastq",repaired=CONSTRAINT_REPAIRED_3),
                subsequences_file=config['constraint_filtering']['undesired_subsequences'][
                    'file'] if 'undesired_subsequences' in config['constraint_filtering']['constraints'] else []
            log:
                "results/logs/{sample}_{unit}/constraint_filtering_after_assembly.txt"
            output:
                expand("results/assembly/{{sample}}_{{unit}}{repaired}{filtered}.fastq",repaired=CONSTRAINT_REPAIRED_3, filtered=CONSTRAINT_FILTER_3),
                expand("results/assembly/{{sample}}_{{unit}}{repaired}{filtered}_out.fastq",repaired=CONSTRAINT_REPAIRED_3, filtered=CONSTRAINT_FILTER_3)
            params:
                constraints=config["constraint_filtering"],
                primer_length=0,
                sequence_length=0,
                paired=False
    else:
        use rule constraint_filtering_base as constraint_filtering_3 with:
            input:
                expand("results/assembly/{{sample}}_{{unit}}_derep{repaired}.fast" + (
                    "a" if config["derep"]["centroid_selection"] == "frequency" else "q"),repaired=CONSTRAINT_REPAIRED_3),
                subsequences_file=config['constraint_filtering']['undesired_subsequences'][
                    'file'] if 'undesired_subsequences' in config['constraint_filtering']['constraints'] else []
            log:
                "results/logs/{sample}_{unit}/constraint_filtering_after_assembly.txt"
            output:
                expand("results/assembly/{{sample}}_{{unit}}{repaired}{filtered}.fastq",repaired=CONSTRAINT_REPAIRED_3,filtered=CONSTRAINT_FILTER_3),
                expand("results/assembly/{{sample}}_{{unit}}{repaired}{filtered}_out.fastq",repaired=CONSTRAINT_REPAIRED_3,filtered=CONSTRAINT_FILTER_3)
            params:
                constraints=config["constraint_filtering"],
                primer_length=0,
                sequence_length=0,
                paired=False

if config['constraint_filtering']['repair_after_demultiplexing']:
    use rule constraint_repair_base as constraint_repair_1 with:
        input:
            expand("demultiplexed/{{sample}}_{{unit}}_{read}.fastq",read=reads),
            #derep="results/assembly/{sample}_{unit}/{sample}_{unit}_derep.fasta",
            #cluster="results/assembly/{sample}_{unit}/{sample}_{unit}_cluster.fasta"
            cluster="results/assembly/{sample}_{unit}/{sample}_{unit}_cluster.fasta"
        output:
            expand("demultiplexed/{{sample}}_{{unit}}_{read}{repaired}.fastq",read=reads,repaired=CONSTRAINT_REPAIRED_1),
            expand("demultiplexed/{{sample}}_{{unit}}_{read}{repaired}_mapping.json", read=reads, repaired=CONSTRAINT_REPAIRED_1),
            #"results/assembly/{sample}_{unit}/{sample}_{unit}_derep_repaired.fastq",# "repaired" bases will have a lower quality score (and will be APPENDED - or replace the original ones?)
            #"results/assembly/{sample}_{unit}/{sample}_{unit}_derep_repair_mapping.json"  # mapping repaired reads to original reads


if config['constraint_filtering']['repair_after_quality_control']:
    use rule constraint_repair_base as constraint_repair_2 with:
        input:
            expand("results/assembly/{{sample}}_{{unit}}/{{sample}}_{{unit}}_{read}.fastq",read=reads),
            #derep="results/assembly/{sample}_{unit}/{sample}_{unit}_derep.fasta",
            cluster="results/assembly/{sample}_{unit}/{sample}_{unit}_cluster.fasta"
        output:
            expand("results/assembly/{{sample}}_{{unit}}/{{sample}}_{{unit}}_{read}{repaired}.fastq",read=reads,repaired=CONSTRAINT_REPAIRED_2),
            expand("results/assembly/{{sample}}_{{unit}}/{{sample}}_{{unit}}_{read}{repaired}_mapping.json", read=reads, repaired=CONSTRAINT_REPAIRED_2),
            #expand("results/assembly/{unit.sample}_{unit.unit}/{unit.sample}_{unit.unit}_{read}{repaired}.fastq",unit=units.reset_index().itertuples(),read=reads,repaired=CONSTRAINT_REPAIRED_2),
            #expand("results/assembly/{unit.sample}_{unit.unit}/{unit.sample}_{unit.unit}_{read}{repaired}_mapping.json",unit=units.reset_index().itertuples(),read=reads,repaired=CONSTRAINT_REPAIRED_2),

if config['constraint_filtering']['repair_after_assembly']:
    #TODO: in this case we want to no only repair invalid centroids but ideally we would want to select:
    # the best representative of the cluster
    # with "best" being: 1) the smallest distance to the current centroid (yet to do) and 2) fulfilling all constraints (DONE)
    if not config["derep"]["clustering"]:
        use rule constraint_repair_base as constraint_repair_3 with:
            input:
                "results/assembly/{sample}_{unit}/{sample}_{unit}_assembled.fastq",
                #derep="results/assembly/{sample}_{unit}/{sample}_{unit}_derep.fasta",
                cent="results/assembly/{sample}_{unit}/{sample}_{unit}_cluster.fasta"
            output:
                expand("results/assembly/{{sample}}_{{unit}}{repaired}.fastq",repaired=CONSTRAINT_REPAIRED_3),
                expand("results/assembly/{{sample}}_{{unit}}{repaired}_clusters.json", repaired=CONSTRAINT_REPAIRED_3),
    elif config["general"]["in_vivo"] or RES_STR != "":
        use rule constraint_repair_base as constraint_repair_3 with:
            input:
                "results/contig_assembly/assemblies/{sample}_{unit}_assembled.fasta",
                #derep="results/assembly/{sample}_{unit}/{sample}_{unit}_derep.fasta",
                #cent="results/contig_assembly/assemblies/{sample}_{unit}_cluster.fasta",
            output:
                expand("results/contig_assembly/assemblies/{{sample}}_{{unit}}_assembled{repaired}.fasta",repaired=CONSTRAINT_REPAIRED_3),
                expand("results/contig_assembly/assemblies/{{sample}}_{{unit}}_assembled{repaired}_clusters.json",repaired=CONSTRAINT_REPAIRED_3),
    else:
        use rule constraint_repair_base as constraint_repair_derep_3 with:
            input:
                "results/assembly/{sample}_{unit}_derep.fasta",
                #derep="results/assembly/{sample}_{unit}/{sample}_{unit}_derep.fasta",
                cent="results/assembly/{sample}_{unit}_cluster.fasta"
            output:
                expand("results/assembly/{{sample}}_{{unit}}_derep{repaired}.fast" + (
                "a" if config["derep"]["centroid_selection"] == "frequency" else "q"),repaired=CONSTRAINT_REPAIRED_3),
                expand("results/assembly/{{sample}}_{{unit}}_derep{repaired}_clusters.json", repaired=CONSTRAINT_REPAIRED_3),

                #expand("results/assembly/{unit.sample}_{unit.unit}/{unit.sample}_{unit.unit}{repaired}.fasta",unit=units.reset_index().itertuples(),repaired=CONSTRAINT_REPAIRED_3),
                #expand("results/assembly/{unit.sample}_{unit.unit}/{unit.sample}_{unit.unit}{repaired}_clusters.json",unit=units.reset_index().itertuples(),repaired=CONSTRAINT_REPAIRED_3),
