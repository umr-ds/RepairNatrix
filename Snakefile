import pandas as pd
from snakemake.utils import validate

validate(config, "schema/config.schema.yaml")

units = pd.read_table(config["general"]["units"], index_col=["sample", "unit"],
    dtype=str)
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])

name_ext = config["merge"]["name_ext"][:-1]

def is_single_end(sample, unit):
    return pd.isnull(units.loc[(sample,unit), "fq2"])

if config["merge"]["paired_End"]:
    reads = [1,2]
else:
    reads = 1

rule all:
    input:
        expand("results/assembly/{unit.sample}_{unit.unit}/{unit.sample}_{unit.unit}_assembled.fastq", unit=units.reset_index().itertuples())
       #"results/assembly/A071SE_A/A071SE_A_dada.fasta",
       # "results/assembly/A071SE_B/A071SE_B_dada.fasta",
       # "results/assembly/A172WA_A/A172WA_A_dada.fasta",
       # "results/assembly/A172WA_B/A172WA_B_dada.fasta"

ruleorder: assembly > prinseq

include: "rules/quality_control.smk"
include: "rules/read_assembly.smk"
include: "rules/dereplication.smk"
include: "rules/chim_rm.smk"
