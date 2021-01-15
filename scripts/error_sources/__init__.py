from scripts.error_sources.gc_content import windowed_gc_content, overall_gc_content
from scripts.error_sources.homopolymers import homopolymer
from scripts.error_sources.kmer import kmer_counting
from scripts.error_sources.undesired_subsequences import undesired_subsequences

active_rules = [overall_gc_content, windowed_gc_content, homopolymer, kmer_counting, undesired_subsequences]
