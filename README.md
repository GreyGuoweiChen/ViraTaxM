# ViraTaxM: targeting genus-level ambiguity in Viral Taxonomic classification of Metagenomic sequences
Viruses exhibit extraordinary genetic diversity. The advent of metagenomic sequencing has revolutionized viral discovery and dramatically expanded our understanding of the virosphere, revealing millions of novel viral sequences and offering unprecedented insights into viral diversity and ecology. With such a high abundance of viral particles, a central challenge in virology and bioinformatics is how to assign a given viral genomic sequence to an appropriate taxonomic label, such as family, genus, or species. Even thought International Committee on Taxonomy of Viruses (ICTV) provides an authoritative framework for virus classification, keeping over 16,000 viral species of 3,700 genera, the classification at the genus-level for metagenomic assembled viral contigs remains complicated due to three major challenges: 1) the extraordinary diversity and heterogeneity of viruses; 2) the fragmentary nature of metagenomic viral sequences; 3) the frequently updated taxonomic framework.

Here, we introduce **ViraTaxM**, an alignment-based tool designed to resolve genus-level ambiguity and enable novel genus discovery in viral classification. ViraTaxM adapts seamlessly to evolving taxonomies without requiring extensive resources. It takes viral sequences as input and predicts their taxonomic lineage down to the genus level.

ViraTaxM consists of two modules, a classification modules and a clustering modules. The classification module applies an alignment-based approach to make initial genus-level predictions, followed by two critical metrics to quantify ambiguity and determine whether the prediction is acceptable. Since metagenomic samples often contain contigs from previously uncharacterized genera, sequences that cannot be assigned to existing genera are processed by the clustering module for novel genus detection. Extensive experiments across multiple scenarios—including fragmentary data classification, out-of-distribution data classification, and open-set detection—demonstrate that ViraTaxM consistently achieves outstanding performance compared to state-of-the-art classification tools.

## Dependency:
* python 3.x
* numpy 1.23.5
* pandas 2.0.3
* scikit-learn 1.3.2
* diamond 2.0.15
* Prodigal 2.6.3
* biopython 1.83
* mcl 22.282

### Quick Installation
You could build the environment from GitHub:
```
    git clone https://github.com/GreyGuoweiChen/ViraTaxM.git
    cd ViraTaxM
    
    # Create the environment and install the dependencies using conda or mamba
    conda env create -f environment.yml
    # or
    mamba env create -f environment.yml
    
    # Activate the environment
    conda activate virataxm
    
    # Distribute ViraTaxM to your conda environment
    pip install .
```
## Usage:
### Example
```
  # you only need to distribute reference database once
  virataxm update --auto
  
  # run virataxm classification and clustering module subsequently
  virataxm predict -i virataxm/test/seq_test.fasta -o virataxm/test
```
