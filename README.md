# Measeq Genomic Summary
Quick scripts for the MeaSeq paper used to characterize and summarize the output fasta files

## Install

Create the `environment.yml` file with conda

## Scripts

### Sequence Stats
Generates simple sequence stats such as the number of Ns, genome completeness, iupacs, and indels for given samples in a directory based off of an alignment to a reference sequence. Note that it doesn't compare to the reference so it really doesn't matter which genotype samples are

#### Run
Example of how to run it with the example data

```python
cd example
python ../scripts/sequence_stats.py \
    --directory ./sequences/ \
    --reference ../reference_seqs/MK513622.1.reference.fasta 
```

#### Output

`fasta_seq_stats.tsv`: A TSV file summarizing all of the stats for the samples found in the directory

### Multi-Seq Compare
Generates a JSON file comparing an input sample and the expected sample based off of an alignment to the reference.

#### Run
Example of how to run it with the example data

```python
cd example
python ../scripts/multi_seq_compare.py \
    --sample sequences/SRR36155247.consensus.fasta \
    --external sequences/SRR36155247.ncbi.fasta \
    --reference ../reference_seqs/MH356245.1.reference.fasta  \
    --name SRR36155247 \
    >> multicomp.json
```

#### Output
Whatever JSON file the data is piped into or to stdout

## Example Data

Example data is from BioProjects: [PRJNA869081](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA869081) and [PRJNA1293457](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA1293457)
