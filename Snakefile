import pandas as pd
from snakemake.utils import validate

validate(config, "schema/config.schema.yaml")

units = pd.read_table(config["general"]["units"], index_col=["sample", "unit"],
    dtype=str)
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])

name_ext = config["general"]["name_ext"][:-1]

def is_single_end(sample, unit):
    return pd.isnull(units.loc[(sample,unit), "fq2"])

if config["general"]["paired_End"]:
    reads = [1,2]
else:
    reads = 1

CONSTRAINT_FILTER_1 = '_constraint_filtered' if config['constraint_filtering']['after_demultiplexing'] else ''
CONSTRAINT_FILTER_2 = '_constraint_filtered' if config['constraint_filtering']['after_quality_control'] else ''
CONSTRAINT_FILTER_3 = '_constraint_filtered' if config['constraint_filtering']['after_assembly'] else ''
CONSTRAINT_REPAIRED_1 = '_constraint_repaired' if config['constraint_filtering']['repair_after_demultiplexing'] else ''
CONSTRAINT_REPAIRED_2 = '_constraint_repaired' if config['constraint_filtering']['repair_after_quality_control'] else ''
CONSTRAINT_REPAIRED_3 = '_assembled_constraint_repaired' if config['constraint_filtering']['repair_after_assembly'] else ''
RES_STR = f"{CONSTRAINT_FILTER_1}{CONSTRAINT_FILTER_2}{CONSTRAINT_FILTER_3}{CONSTRAINT_REPAIRED_1}{CONSTRAINT_REPAIRED_2}{CONSTRAINT_REPAIRED_3}"

if config["general"]["in_vivo"]:
    rule all:
        input:
            expand("results/contig_assembly/assemblies/{unit.sample}_{unit.unit}.fasta", unit=units.reset_index().itertuples())
elif not config["derep"]["clustering"]:
    rule all:
        input:
            expand("results/assembly/{unit.sample}_{unit.unit}/{unit.sample}_{unit.unit}_derep.fasta", unit=units.reset_index().itertuples()) if config["derep"] == "frequency" else expand("results/assembly/{unit.sample}_{unit.unit}/{unit.sample}_{unit.unit}_derep.fastq", unit=units.reset_index().itertuples())
else:
    rule all:
        input:
            expand("results/assembly/{unit.sample}_{unit.unit}/{unit.sample}_{unit.unit}_fin.fasta",#{res_str}
                unit=units.reset_index().itertuples())  #res_str=RES_STR

ruleorder: assembly > prinseq

include: "rules/demultiplexing.smk"
include: "rules/read_assembly.smk"
include: "rules/dereplication.smk"
include: "rules/constraint_filtering.smk"
include: "rules/contig_assembly.smk"
