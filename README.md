# Genotype-Phenotype Boxplots
A tool for rapidly generating swarmplots and boxplots for phenotypes across genotypes from a gzipped VCF.

```bash
usage: genotype_phenotype_boxplots.py [-h] --gzvcf GZVCF --phenotype-file PHENOTYPE_FILE --phenotype PHENOTYPE --position POSITION --output OUTPUT [--test TEST] [--groups GROUPS] [--highlights HIGHLIGHTS] [--geneid GENEID]
                                      [--sample-of-interest SAMPLE_OF_INTEREST]

This script was designed to plot the distribution of a given phenotype for all genotypes in a VCF given a specific chromosome and position.

optional arguments:
  -h, --help            show this help message and exit
  --test TEST           Test type: t-test_ind, t-test_welch, t-test_paired, Mann-Whitney, Mann-Whitney-gt, Mann-Whitney-ls, Levene, Wilcoxon, Kruskal.
  --point-size POINT_SIZE
                        Size for points in swarmplot. If size=0, then points are removed.
  --groups GROUPS       A CSV containing taxa,groupid.
  --highlights HIGHLIGHTS
                        A comma-separated list of groups to highlight. Ex. durra,milo
  --geneid GENEID       Name of gene to add to title.
  --sample-of-interest SAMPLE_OF_INTEREST
                        File containing samples of interest to print genotype information. Header (CommonName,CUSO,PI).

required arguments:
  --gzvcf GZVCF         The name of the input gzipped VCF file. Must be tbi indexed.
  --phenotype-file PHENOTYPE_FILE
                        A CSV of relevant phenotypes with taxa in the first column. Note: headers are **required**.
  --phenotype PHENOTYPE
                        The phenotype column to use in the phenotype file
  --position POSITION   Position of interest in Chr01_Position format (or what matches the VCF).
  --output OUTPUT       Name of the output figure
```
