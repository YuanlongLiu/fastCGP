# fastCGP
An efficient and exact method to compute gene-level P-values from GWAS SNP-level results


(To be modified) One essential question in gene-based analysis (e.g. pathway- and network-based analysis) is how to combine SNP-level P-values into a single representative gene-level P-value. A popular strategy is to take the best SNP P-value as gene-level P-value, which is sensitive in capturing high association signals, however, is also biased to genes annotated with more SNPs. Many methods that aim to correct for bias of gene size have been proposed, such as ARTP, Fisher's combination method etc. [refs]. However, most of these methods rely on the assumption of independence between test statistics, which is inappropriate for gene-based analysis because of the well-known LD pattern between SNPs [ref]. Phenotype permutation is known as the gold standard for preservation of the LD pattern, thus is promising for correcting bias of gene size [ref]. Yet, for reason of computational intensity and the required access to raw genotype data, it has never been applied to compute gene-level P-values to our knowledge. An alternative permutation strategy that accounts for the LD pattern while does not require genotype data is the Circular Genomic Permutation (CGP) strategy [ref]. It permutes SNP-level statistics in a genomic manner to keep similar patterns of correlation in the permuted data as in the original data. 

We took advantage of the CGP strategy and proposed a novel method named fastCGP to compute gene-level P-values. fastCGP first annotates SNPs to genes if a SNP is located within the boundary of a gene (gene boundaries can be defined by users, for example 20Kb in both upstream and downstream of the gene coding regions). Then for each gene, the gene-level P-value is taken as the best SNP-level P-value and corrected for bias of gene size using the concept of permutation test with following steps:

permute SNP-level P-values K times
count the frequency of extreme permutations k
compute gene-level P-value: P = (k+1)/(k+1)


, defined as p = (k+1)/(K+1), where K is the total number of non-repeating permutation samples, and k is the number of extreme samples, defined as  take all non-repeating CGP samples into account to obtain the best obtainable permutation P-value. 

The using of analytical approach brings several advantages. First, the method is exact, hence does not suffer from randomness compared to many simulation-based methods where random numbers are used, such as VEGAS, ARTP, Pascal [refs]. Second, the computation is efficient. In our study of 2,370,689 SNPs and 24,120 genes, it takes around 0.5 hour on a standard PC (Intel Core i7 3.40GHz CPU, 8GB RAM), which is much faster than many simulation-based competitors [refs]. Third, the resultant P-values are of high precision. According to the formulation of fastCGP, the precision of each gene-level P-value is 1/L , where  L is the total number of SNPs analyzed in a GWAS. The larger the amount of SNPs is analyzed, the higher is the precision. In our case of analyzing 2,370,689 SNPs, the precision reaches to 4.2x10-7.

Details of how fastCGP works will be added later.



Though fastCGP is permutation-based in nature, no sample needs to be generated during computation. More particularly, instead of generating a given number of samples (100,000 for example) to get an empirical estimation,

fastCGP corrects the best SNP P-value for gene size while accounts for the correlation between SNPs. 
