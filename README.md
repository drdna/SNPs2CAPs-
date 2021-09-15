# SNPs2CAPs-

SNPs2CAPs.pl reads a list of SNP positions between two genomes, grabs the corresponding SNP and flanking sequence from the relevant genome and identifies those where the SNP causes a restriction site polymorphism.

Script can be run using the provided resource files (gunzip the Fusarium genome file first): 

  perl SNPs2CAPs.pl RE_sites.txt Fusarium_graminearum.fasta2 Fgraminearum_snps.txt <flanking_sequence_length> (recommend ≥ 6)
 


