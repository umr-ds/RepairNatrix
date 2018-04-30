from glob import glob
import yaml
import pandas as pd
import numpy as np
### Create the datatable containing the samples, units and paths of all fastq files formatted correctly ###
### this is vital for the snakemake pipeline, without it, the wildcards can't be created ###

# might be better to glob abspaths
with open('../../workflow_test/amplicon_snakemake_pipeline/demultiplexer/test_data.yaml') as f_:
    config = yaml.load(f_)
file_path_list = sorted(glob('demultiplexed/' + "/*.gz"))
file_list = sorted([file_.split('/')[-1] for file_ in glob('demultiplexed/' + "/*.gz")])
df = pd.DataFrame(columns=['sample', 'unit', 'fq1', 'fq2'],
    index =range(int(len(file_list)/2)))
if config['merge']['paired_End']:
    df = pd.DataFrame(columns=['sample', 'unit', 'fq1', 'fq2'],
        index =range(int(len(file_list)/2)))
    i, j = (0, 0)
    while i < len(file_list)/2:
        df.loc[i] = file_list[j].split('_')[0], file_list[j].split('_')[1], \
                file_path_list[j][:-3], file_path_list[j+1][:-3]
        j += 2
        i += 1
else:
    df = pd.DataFrame(columns=['sample', 'unit', 'fq1', 'fq2'], 
        index = range(int(len(file_list))))
    i = 0
    while i < len(file_list):
        df.loc[i] = file_list[i].split('_')[0], file_list[i].split('_')[1], \
            file_list[i][:-3], np.nan
        i += 1
df.to_csv('units.tsv', sep='\t')
pd.DataFrame(df['sample']).to_csv('samples.tsv', sep='\t')