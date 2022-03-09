from scripts.constraints.gc_content import windowed_gc_content, overall_gc_content
from scripts.constraints.homopolymers import homopolymer
from scripts.constraints.kmer import kmer_counting
from scripts.constraints.undesired_subsequences import undesired_subsequences

constraints = [homopolymer, undesired_subsequences, overall_gc_content, windowed_gc_content, kmer_counting]
