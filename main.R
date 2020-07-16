	##  fastCGP: an efficient and exact method to compute gene-level P-values from GWAS SNP-level results
	##  Developed by Yuanlong LIU, contact: yuanlong.liu@inserm.fr
	##  Latest update: 26-07-2017
	##	If you this method, please cite in your publication: Liu, Y., et al. "Network-assisted analysis of GWAS data identifies a functionally-relevant gene module for childhood-onset asthma." Scientific Reports 7.1 (2017): 938.
	
	##	this sub-function constructs the 'circle' to perform fastCGP, which includes information of ordered snp p-values, and their associated genes 
	fastCGP_cir_info <- function( snp2gene, snp_chr_pos_p )
	{

		genes = unique( snp2gene$gene )
		mapped_snp  = unique( snp2gene$SNP ) ##snps in genes
		all_snps = unique( snp_chr_pos_p$SNP ) ##the complete snp list
		if( any(duplicated( all_snps )) ) stop('\nSome SNPs appear more than once in your snp_chr_pos_p file\nPlease remove duplicates first')
		
		mapped_snp_nin_all_snps = setdiff( mapped_snp, all_snps )
		diff_len = length(mapped_snp_nin_all_snps) ##SNPs mapped to gene but do not have a position in snp_chr_pos_p file
		if( diff_len > 0 ) stop('\n', len, ' SNPs are mapped to genes in your snp2gene file \nBut their chromosomal or p-value information are missed in your snp_chr_pos_p file')

		genes_count = length(genes)
		snps_count = length(all_snps)
		mapped_snp_count = length(mapped_snp)
		cat('\nNumber of genes in your snp2gene file:', genes_count, '\n')
		cat('Number of SNPs in your snp_chr_pos_p file:', snps_count, '\n')
		cat('Number of SNPs mapped to genes:', mapped_snp_count, '\n')
		
		snps_p = snp_chr_pos_p$p
		names( snps_p ) = snp_chr_pos_p$SNP ##named p values for complete SNP list
		
		snp_prk = rank( snps_p, ties.method='random' ) ##p value rank of snps

		snps_pos_table = 1:snps_count
		names(snps_pos_table) = snp_chr_pos_p[ order(snp_chr_pos_p$chr, snp_chr_pos_p$pos), ]$SNP ##place snps on the table according to their position on the chrosome
		
		table_rk = snp_prk[ names(snps_pos_table) ] ##rank of the snps on the table; rank(c(1,2)): 1 2
		
		gene_snps_poses = tapply( snps_pos_table[ snp2gene$SNP ], snp2gene$gene, function(v) {return(v)} ) ##position of snps on table of gene
		
		if(any(is.na(unlist(gene_snps_poses, recursive=TRUE)))) {cat( 'NA information from SNP table again!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' ); return(0)}
	
		lens = table( snp2gene$gene )[1:genes_count] ##number of snps in each gene
		gene_minrk = as.list(tapply( snp_prk[ snp2gene$SNP ], snp2gene$gene, min ))##best rank of each gene

		snp2gene_p = snp2gene; snp2gene_p$p = snps_p[ snp2gene$SNP ]
		gene_minP = tapply(snp2gene_p$p, snp2gene_p$gene, min) ##best rank of each gene
		
		return( list(gene_snps_poses=gene_snps_poses, gene_lens = lens, gene_minrk=gene_minrk, gene_minP=gene_minP, table_rk=table_rk) )	
	}

	fastCGP <- function( snp2gene_file, snp_chr_pos_p_file, genes2compute_file )
	{
		time_0 = proc.time()
		
		cat('\nNote: all gene names will be converted to uppercase\n')
		
		##read files:
		snp2gene = read.table( snp2gene_file, header=TRUE, stringsAsFactors=FALSE, comment.char='' )
		snp_chr_pos_p = read.table( snp_chr_pos_p_file, header=TRUE, stringsAsFactors=FALSE, comment.char='' )
		
		##lowercase to uppercase
		snp2gene$gene = toupper( snp2gene$gene )
		
		genes = unique( snp2gene$gene )
		if( missing(genes2compute_file) ) genes2compute = genes
		else { genes2compute = read.table( genes2compute_file, header=TRUE, stringsAsFactors=FALSE, comment.char='' )
		genes2compute = unique( toupper( genes2compute$gene ) ) }
		
		unmapped_genes = setdiff( genes2compute, genes )
		if( length( unmapped_genes ) > 0 ) cat( sep='', '\n', length( unmapped_genes ), ' genes in your genes2compute file do not have SNP_to_gene mapping information\nThey are removed for downstream computation\n' )
		
		cir_info = fastCGP_cir_info( snp2gene, snp_chr_pos_p )

		sorted_genes = names(sort(unlist(cir_info$gene_minrk), decreasing=TRUE))
		genes2compute = intersect( sorted_genes, genes2compute ) ##ordered according to sorted_genes
		
		x = cir_info$gene_minrk[genes2compute]
		if(!all(x==cummin(x))) stop('\nThe genes to be computed need to be ranked by minP') ##check whether genes are ranked by minP; this is required for the purpose of increasing efficiency
		table_len = length( cir_info$table_rk ) ##total number of SNPs on the circle
		pos2compare = cbind(rk=unname(cir_info$table_rk), pos=1:table_len)
		sig_pos_ini = pos2compare ## the rank at each position
		gene_lens = cir_info$gene_lens[genes2compute]
		
		genes_p =list()
		i=1
		
		cat('\n\n\n\n@@@@-------------------------------- COMPUTATION BEGINS --------------------------------@@@@\n')
		cat('\n\nNumber of genes already computed:\n')
		##compute corrected gene-level p-value gene by gene
		for(gene in genes2compute)
		{
			gene_len = gene_lens[gene]
			
			if( gene_len == 1 ) { genes_p[[ gene ]] = cir_info$gene_minP[ gene ]; next } ##if the gene only contains one SNP, the p-value will be that SNP p-value
			
			thresh = cir_info$gene_minrk[ gene ] ##Test statistics //
			sig_pos = sig_pos_ini[which(sig_pos_ini[,'rk'] <= thresh), 'pos' ] ## CONFRIMED
			if(length(sig_pos) ==1) { genes_p[[gene]]= ( gene_len + 1 )/( table_len + 1 ); next}  ##only extreme cases, that a gene contains the best snp in all_snps 
			gaps = c( diff( sig_pos ) - 1, (table_len - sig_pos[length(sig_pos)]) + (sig_pos[1] - 1) ) ##CONFRIMED
			
			big_gaps = gaps[gaps >= gene_len]
			gap_width = big_gaps - gene_len + 1 ##CONFRIMED
			genes_p[[gene]]  = 1 - sum(gap_width) / ( table_len + 1)
			
			if( i %% 100 ==0) cat(i, '|')
			sig_pos_ini = pos2compare[ sig_pos, ] ##because genes are ranked according to minP
			i=i+1
		}
		cat('\n\n\n\n@@@@--------------------------------- COMPUTATION ENDS ---------------------------------@@@@\n')

		delta_time = proc.time() - time_0
		minutes = format( round( delta_time[3] / 60, 1), nsmall=1)
		cat('\n\nTotal computational time:', minutes, 'minutes', '\n') ## output the computational time. Added this comment to see changes in git
		
		genes_p_frame = data.frame( gene=names( genes_p ), p=unname(unlist( genes_p )) )
		write.table(genes_p_frame, file='./computed_gene_p_values.tab', row.names=FALSE, sep='\t', quote=FALSE)
	
		return( genes_p_frame )
	}


