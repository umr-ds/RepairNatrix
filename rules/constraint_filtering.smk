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


if(config['constraint_filtering']['after_demultiplexing']):
    use rule constraint_filtering_base as constraint_filtering_1 with:
        input:
            expand("demultiplexed/{{sample}}_{{unit}}_{read}.fastq",read=reads),
            subsequences_file = config['constraint_filtering']['undesired_subsequences']['file'] if 'undesired_subsequences' in config['constraint_filtering']['constraints'] else []
        log:
            "results/logs/{sample}_{unit}/constraint_filtering_after_demultiplex.txt"
        output:
            temp(expand("demultiplexed/{{sample}}_{{unit}}_{read}{filtered}.fastq",read = reads, filtered=CONSTRAINT_FILTER_1)),
            temp(expand("demultiplexed/{{sample}}_{{unit}}_{read}{filtered}_out.fastq",read = reads, filtered=CONSTRAINT_FILTER_1))

if(config['constraint_filtering']['after_quality_control']):
    use rule constraint_filtering_base as constraint_filtering_2 with:
        input:
            expand("results/assembly/{{sample}}_{{unit}}/{{sample}}_{{unit}}_{read}.fastq",read=reads),
            subsequences_file = config['constraint_filtering']['undesired_subsequences']['file'] if 'undesired_subsequences' in config['constraint_filtering']['constraints'] else []
        log:
            "results/logs/{sample}_{unit}/constraint_filtering_after_qc.txt"
        output:
            expand("results/assembly/{{sample}}_{{unit}}/{{sample}}_{{unit}}_{read}{filtered}.fastq",read = reads, filtered=CONSTRAINT_FILTER_2),
            expand("results/assembly/{{sample}}_{{unit}}/{{sample}}_{{unit}}_{read}{filtered}_out.fastq",read = reads, filtered=CONSTRAINT_FILTER_2)

if(config['constraint_filtering']['after_assembly']):
    use rule constraint_filtering_base as constraint_filtering_3 with:
        input:
            "results/assembly/{sample}_{unit}/{sample}_{unit}_assembled.fastq",
            subsequences_file = config['constraint_filtering']['undesired_subsequences']['file'] if 'undesired_subsequences' in config['constraint_filtering']['constraints'] else []
        log:
            "results/logs/{sample}_{unit}/constraint_filtering_after_assembly.txt"
        output:
            expand("results/assembly/{{sample}}_{{unit}}/{{sample}}_{{unit}}_assembled{filtered}.fastq",filtered=CONSTRAINT_FILTER_3),
            expand("results/assembly/{{sample}}_{{unit}}/{{sample}}_{{unit}}_assembled{filtered}_out.fastq",filtered=CONSTRAINT_FILTER_3)
        params:
            constraints = config["constraint_filtering"],
            primer_length = 0,
            sequence_length = 0,
            paired = False

