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

rule all:
    input:
        expand("results/assembly/{unit.sample}_{unit.unit}/{unit.sample}_{unit.unit}_assembled{filtered}.fastq",
            filtered=CONSTRAINT_FILTER_3,
            unit=units.reset_index().itertuples())

ruleorder: assembly > prinseq

include: "rules/demultiplexing.smk"
include: "rules/quality_control.smk"
include: "rules/read_assembly.smk"
include: "rules/dereplication.smk"
include: "rules/chim_rm.smk"
include: "rules/constraint_filtering.smk"
