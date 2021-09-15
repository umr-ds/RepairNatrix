from subprocess import call

# If the sequences are single-end, they still need the read identifier
# (R1), otherwise, the string slicing should remove the last six symbols
# instead of the last 8.

if(len(snakemake.input)) == 2:
    output_edit = str(snakemake.output[0])[:-8]
    output_bad = str(snakemake.output[0])[:-8] + "_low_qual"
else:
    output_edit = str(snakemake.output[0])[:-6]
    output_bad = str(snakemake.output[0])[:-6] + "_low_qual"

call_list = [
    "prinseq-lite.pl",
    "-verbose",
    "-fastq", snakemake.input[0],
    "-ns_max_n", "0",
    "-min_qual_mean", str(snakemake.params.mq),
    "-out_good", output_edit,
    "-out_bad", output_bad,
    "-log", str(snakemake.log)
]
if(len(snakemake.input)) == 2:
    call_list.extend(["-fastq2", snakemake.input[1]])

if(snakemake.params.trim_to_length > 0):
    call_list.extend([
        "-trim_to_len", str(snakemake.params.trim_to_length)
    ])

call(call_list)