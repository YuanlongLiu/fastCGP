# fastCGP: an efficient and exact method to compute gene-level P-values from GWAS SNP-level results

One essential question in gene-based analysis (e.g. pathway- and network-based analysis) is how to combine SNP-level P-values into a single representative gene-level P-value. A popular approach is to take the best SNP P-value, among all SNPs that are annotated to a gene, as the gene-level P-value. This approach is sensitive in capturing high association signals, however, is also biased to genes annotated with more SNPs. Other approach that aim to correct for bias of gene size have been proposed, such as ARTP, Fisher's combination method etc. Nonetheless, most of these approaches rely on the assumption of independence between SNP-level P-values, hence is inappropriate for gene-based analysis because of the well-known LD pattern between SNPs. Phenotype permutation is known as the gold standard for preservation of the LD pattern, thus is promising for correcting bias of gene size [ref]. Yet, for reason of computational intensity and the required access to raw genotype data, it was rarely applied to compute gene-level P-values. An alternative permutation strategy that accounts for the LD pattern while does not require genotype data is the Circular Genomic Permutation (CGP) strategy (https://www.ncbi.nlm.nih.gov/pubmed/22973544). Briefly, CGP considers the genome as a circle, starting from chromosome 1 to chromosome X and restarting at chromosome 1. SNP-level P-values of a GWAS are ordered on the circle according to the position of the SNPs. A CGP sample can be generated by rotating the P-values for a random position and reassigning them to each SNP. By permuting SNP-level statistics in a circular manner, CGP to keep similar patterns of correlation in the permuted data as in the original data.

Here, we took advantage of the CGP strategy and proposed a novel method named fastCGP to compute gene-level P-values from SNP-level P-values. fastCGP first annotates SNPs to genes if a SNP is located within the boundary of a gene (gene boundaries can be defined by users, for example 20Kb in both upstream and downstream of the gene coding regions). Then for each gene, the gene-level P-value is taken as the best SNP-level P-value among all SNPs mapped to the gene. This P-value is further corrected for bias of gene size using the concept of permutation test, in which the final gene-level P-value is defined as P = (k+1)/(K+1), where K is the total number of CGP samples, and k is the number of extreme samples. Specially, instead of generating some given number of CGP samples, fastCGP takes all non-repeating CGP samples into account to obtain the best obtainable P-value in this permutation framework. Under this specification, K becomes the total number of SNPs invovled in a GWAS, while k can be obtained analytically without generating any CGP sample (details are not presented here. Reference coming soon).

The implementation of this analytical approach brings several advantages. First, the method is exact, hence does not suffer from randomness compared to many simulation-based methods where random numbers are used, such as VEGAS, ARTP, Pascal [refs]. Second, the computation is efficient. In our study of 2,370,689 SNPs and 24,120 genes, it takes around 0.5 hour on a standard PC (Intel Core i7 3.40GHz CPU, 8GB RAM), which is much faster than many simulation-based competitors [refs]. Third, the resultant P-values are of high precision. According to the formulation of fastCGP, the precision of each gene-level P-value is 1/K , where  L is the total number of SNPs analyzed in a GWAS. The larger the amount of SNPs is analyzed, the higher is the precision. In our case of analyzing 2,370,689 SNPs, the precision reaches to 4.2 x 10-7.

## Usage

fastCGP is written in R. To perform the analysis, first load the R scripts into the R working environment. The main function ```fastCGP( snp2gene_file, snp_chr_pos_p_file, genes2compute_file )``` requires three arguments:

- ```snp2gene_file```: a tab-delimited two-column text file that contains information of which SNPs are mapped to which genes. An example:
	
For example,	

	gene	SNP
	KCNIP4	rs10000010
	KIAA1530	rs10000012	
	0.491,0.492
	0.990,0.993
	0.775,0.777
	...
	0.577,0.561

The followings are shown at standard output.


- ```snp2gene_file```: a tab-delimited four-column text le that contains the information of all SNPs that are analyzed in a GWAS, including the name of a SNP, its chromosomal number, position on the chromosome, and GWAS association p-value.
- ```snp2gene_file```: a one-column text le that contains the list of genes which you want to compute their p-values. If this le is not provided, all genes included in the "snp2gene.tab" le (SNP2gene le) will be computed.


```javascript
{
    "pref-labels": {
        "poverty": {

}

```


## Contact information
Usage and examples will be added soon

* Author: Yuanlong LIU
* Affiliation: French National Institute of Health and Medical Research, Unit 946,  Paris, France
* Email: yuanlong.liu@inserm.fr
