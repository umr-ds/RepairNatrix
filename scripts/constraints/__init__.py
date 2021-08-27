from constraints.gc_content import windowed_gc_content, overall_gc_content
from constraints.homopolymers import homopolymer
from constraints.kmer import kmer_counting
from constraints.undesired_subsequences import undesired_subsequences

constraints = [homopolymer, undesired_subsequences, overall_gc_content, windowed_gc_content, kmer_counting]
