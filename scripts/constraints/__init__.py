from .gc_content import windowed_gc_content, overall_gc_content
from .homopolymers import homopolymer
from .kmer import kmer_counting
from .undesired_subsequences import UndesiredSubSequenceFinder
uFinder = UndesiredSubSequenceFinder("scripts/constraints/undesired_subsequences.txt") #scripts/constraints/

constraints = [homopolymer, uFinder.undesired_subsequences, overall_gc_content, windowed_gc_content, kmer_counting]
