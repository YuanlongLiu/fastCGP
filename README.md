# fastCGP: an efficient and exact method to compute gene-level *P*-values from GWAS SNP-level results

## Introduction

One essential question in gene-based analysis (e.g. pathway- and network-based analysis) is how to combine SNP-level *P*-values into a single representative gene-level *P*-value. A popular approach is to take the best SNP *P*-value, among all SNPs that are annotated to a gene, as the gene-level *P*-value. This approach is sensitive in capturing high association signals, however, is also biased to genes annotated with more SNPs. Other approach that aim to correct for bias of gene size have been proposed, such as ARTP, Fisher's combination method etc. Nonetheless, most of these approaches rely on the assumption of independence between SNP-level *P*-values, hence is inappropriate for gene-based analysis because of the well-known LD pattern between SNPs. Phenotype permutation is known as the gold standard for preservation of the LD pattern, thus is promising for correcting bias of gene size [ref]. Yet, for reason of computational intensity and the required access to raw genotype data, it was rarely applied to compute gene-level *P*-values. 

An alternative permutation strategy that accounts for the LD pattern while does not require genotype data is the Circular Genomic Permutation (CGP) strategy (https://www.ncbi.nlm.nih.gov/pubmed/22973544). Briefly, CGP considers the genome as a circle, starting from chromosome 1 to chromosome X and restarting at chromosome 1. SNP-level *P*-values of a GWAS are ordered on the circle according to the position of the SNPs. A CGP sample can be generated by rotating the *P*-values for a random position and reassigning them to each SNP. By permuting SNP-level statistics in a circular manner, CGP to keep similar patterns of correlation in the permuted data as in the original data.

Here, we took advantage of the CGP strategy and proposed a novel method named fastCGP to compute gene-level *P*-values from SNP-level *P*-values. fastCGP first annotates SNPs to genes if a SNP is located within the boundary of a gene (gene boundaries can be defined by users, for example 20Kb in both upstream and downstream of the gene coding regions). Then for each gene, the gene-level *P*-value is taken as the best SNP-level *P*-value among all SNPs mapped to the gene. This *P*-value is further corrected for bias of gene size using the concept of permutation test, in which the final gene-level *P*-value is defined as P = (k+1)/(K+1), where K is the total number of CGP samples, and k is the number of extreme samples. Specially, instead of generating some given number of CGP samples, fastCGP takes all non-repeating CGP samples into account to obtain the best obtainable *P*-value in this permutation framework. Under this specification, K becomes the total number of SNPs involved in a GWAS, while k can be obtained analytically without generating any CGP sample (details are not presented here. Reference coming soon).

The implementation of this analytical approach brings several advantages. First, the method is exact, hence does not suffer from randomness compared to many simulation-based methods where random numbers are used, such as VEGAS, ARTP, Pascal [refs]. Second, the computation is efficient. In our study of 2,370,689 SNPs and 24,120 genes, it takes around 0.5 hour on a standard PC (Intel Core i7 3.40GHz CPU, 8GB RAM), which is much faster than many simulation-based competitors [refs]. Third, the resultant *P*-values are of high precision. According to the formulation of fastCGP, the precision of each gene-level *P*-value is 1/K , where K is the total number of SNPs analyzed in a GWAS. The larger the amount of SNPs is analyzed, the higher is the precision. In our case of analyzing 2,370,689 SNPs, the precision reaches to 4.2 x 10-7.


## Installation

`git clone https://github.com/CSOgroup/fastCGP.git`

`install.packages(path_to_fastCGP, repos = NULL, type="source")` ## install from the cloned source file


## Usage

fastCGP is written in R. To perform the analysis, first load the R scripts into the R working environment. fastCGP requires three inputs:

- ```snp2gene_file```: (the path of) a space-delimited two-column text file that contains information of which SNPs are mapped to which genes. The column headers of this file should be exactly ```gene``` and ```SNP``` (order insensitive).

An example:	

	gene SNP
	KCNIP4 rs10000010
	KIAA1530 rs10000012
	BMPR1B rs10000023
	FAM114A1 rs10000037
	STK32B rs10000042
	STK32B rs10000062
	...

- ```snp_chr_pos_p_file```: (the path of) a space-delimited four-column text file that contains the information of all SNPs that are analyzed in a GWAS, including the name of a SNP, its chromosomal number, position on the chromosome, and GWAS association *P*-value. The column headers of this file should be exactly ```SNP```, ```chr```, ```pos```, and ```p``` (order insensitive). More SNPs contained in this file will lead to better precision of the result. 

An example:	

	SNP chr pos p
	rs1000000 12 126890980 0.394150
	rs10000010 4 21618674 0.821190
	rs10000012 4 1357325 0.623990
	rs10000013 4 37225069 0.081821
	rs10000017 4 84778125 0.211210
	rs1000002 3 183635768 0.305950
	...

- ```genes2compute_file``` (optional): (the path of) a one-column text file that contains the list of genes which you want to compute their *P*-values. The column header of this file should be exactly ```gene```. If this file is not provided, all genes included in the ```snp_chr_pos_p_file``` will be computed.

An example:	

	gene
	MICA
	S1PR3
	C6ORF15
	HLA-DQA1
	ZNF329
	TCF19
	...
 
Given these files, run the following code to compute gene *P*-values:
```
results = fastCGP( snp2gene_file, snp_chr_pos_p_file, genes2compute_file )
```

The results will be saved in ```results``` as a data frame, and will be saved in the ```computed_gene_p_values.tab``` file in your working directory.

## Citation

If you use fastCGP in your work, please cite: https://www.nature.com/articles/s41598-017-01058-y


## Contact information

* Author: Yuanlong LIU
* Affiliation: French National Institute of Health and Medical Research, Unit 946,  Paris, France
* Email: yliueagle@gmail.com
