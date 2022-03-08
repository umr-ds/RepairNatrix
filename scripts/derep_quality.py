from Bio import SeqIO as sio
from Bio.SeqRecord import SeqRecord


def derep(inp):
    seqs = {}
    for record in sio.parse(inp, "fastq"):
        if record.seq not in seqs.keys():
            seqs[record.seq] = {}
            seqs[record.seq]["quality"] = record.letter_annotations["phred_quality"]
            seqs[record.seq]["size"] = 1
        else:
            new_phred = sum(record.letter_annotations["phred_quality"])/len(record.letter_annotations["phred_quality"])
            old_phred = sum(seqs[record.seq]["quality"])/len(seqs[record.seq]["quality"])
            if new_phred > old_phred:
                seqs[record.seq]["quality"] = record.letter_annotations["phred_quality"]
            seqs[record.seq]["size"] += 1
    return seqs


def get_records(seq_entries):
    for seq, vals in seq_entries.items():
        yield SeqRecord(seq, id='0', description=';size={};'.format(vals["size"]),
                        letter_annotations={"phred_quality": vals["quality"]})


derep_seqs = derep(str(snakemake.input))
with open(str(snakemake.output), "w") as outp:
    sio.write(get_records(derep_seqs), outp, "fastq")
