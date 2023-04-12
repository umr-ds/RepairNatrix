# RepairNatrix

![DAG of an example workflow](documentation/images/dag.svg)
*DAG of an example workflow for RepairNatrix

---
constraint_filtering can be turned on optionally  
available constraints:  
* homopolymers
* overall_gc_content
* windowed_gc_content
* undesired_subsequences
* kmer_counting
##### too harsh filtering may result in execution errors of other rules (empty files) !

---
other changes:
* added rules for vsearch derep & clustering with clustering_id option
* constraint_filtering: paired reads will be filtered together (if read 1 would be thrown away so would read 2, and vice versa)
* constraint_filtering: options for primer & sequence length, which are used in the filtering rules until assembly
* prinseq: option to trim to length (cut right side of trim)

---

#### Usage

RepairNatrix uses a yaml file to specify the input files and parameters. 
A schema to validate yaml files including a description for all field can be found [here](https://github.com/umr-ds/RepairNatrix/example_data.yaml).
An example config can be found [here](https://github.com/umr-ds/RepairNatrix/example_data.yaml).

For the full documentation see [@Natrix](https://github.com/MW55/Natrix)
