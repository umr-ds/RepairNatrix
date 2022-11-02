import os
#import sys
from tqdm import tqdm
#import dinopy
import subprocess

print(snakemake.params.blast_db)
blast_db_targets = snakemake.params.blast_db[str(snakemake.params.sample).split('-')[-1]]
os.makedirs(snakemake.params.spades_dir, exist_ok=True)

print(f"Processing sample")
    
sample_in_path = snakemake.input.derep

print("Filtering biological and non-biological reads")
    
sample_blast_hits_path = snakemake.output.blast_hits
sample_data_hits_path = snakemake.output.data_hits

blast_process = subprocess.run(["blastn","-db", blast_db_targets, "-query",sample_in_path, "-outfmt","7 sseqid", "-out", sample_blast_hits_path])
if blast_process.returncode != 0:
    print("Something went wrong with blast.")
    exit(1)

data_hits = set()
with tqdm(total=os.path.getsize(sample_blast_hits_path)) as progress:
    with open(sample_blast_hits_path) as file:
        current_query = ""
        for line in file:
            if line.startswith("# "):
                information = line.split(" ")
                if information[1] == "Query:":
                    current_query = information[2].rstrip()
                elif information[1] in ["BLASTN", "Database:","Fields:","BLAST"]:
                    pass
                else:
                    if int(information[1]) == 0:
                        data_hits.add(current_query)
            progress.update(len(line.encode("utf-8")))

'''
far = dinopy.FastaReader(sample_in_path)
sample_data_hits_path = snakemake.output.data_hits
faw = dinopy.FastaWriter(sample_data_hits_path)
faw.open()
    
with tqdm(total=os.path.getsize(sample_in_path)) as progress:
    for seq, name in far.reads(read_names=True):
        if name.decode() in data_hits:
            faw.write_entry((seq, name))
        progress.update(len(seq)+len(name))
   
faw.close()
'''
with tqdm(total=os.path.getsize(sample_in_path)) as progress:
    with open(sample_in_path, "r") as inp, open(sample_data_hits_path, "w") as outp:
        while True:
            name = inp.readline().rstrip()
            seq = inp.readline().rstrip()
            if not seq:
                break
            assert name[0] == ">"
            if name[1:] in data_hits:
                outp.write(name + "\n")
                outp.write(seq + "\n")
            progress.update((len(seq.encode("utf-8")) + len(name.encode("utf-8"))))


print("Assembling data reads.")
spades_out_folder = snakemake.params.spades_dir
spades_process = subprocess.run(["spades.py", "-o", spades_out_folder,"-s", sample_data_hits_path, "--only-assembler",
                                 "-k", "21,33,55,77,99"])
if spades_process.returncode != 0:
    print("Something went wrong with spades.")
    exit(1)
#
#print("Extracting largest sequence.")
#sample_contigs_path = snakemake.output.contigs
#data_assembly_path = snakemake.output.assembly
#largest_contig = []
#with open(sample_contigs_path, 'r') as contigs:
#    first = True
#    for line in contigs:
#        if line.startswith('>'):
#            if first:
#                largest_contig.append(line)
#                first = False
#            else:
#                break
#        else:
#            largest_contig.append(line)
#with open(data_assembly_path, 'w') as assembly_outfile:
#    assembly_outfile.writelines(largest_contig)
   
print("Done.")
