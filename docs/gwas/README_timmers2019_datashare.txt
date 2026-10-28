This dataset contains HRC-imputed GWAS summary statistics on parental survival, under the assumption of common effect sizes between fathers and mothers.
rsid – reference SNP cluster ID;
snpid – SNP chromosome and position ID;
chr – Chromosome;
pos – Base-pair position on chromosome (GRCh37);
a1 – the effect allele;
a0 – the reference allele;
n – Number of phenotypes analysed for the SNP. Note that combined parental residuals from UK Biobank are counted once, and LifeGen meta-analysed father and mother statistics are counted once each, and as a result, the total number of phenotypes analysed is actually greater than the number listed.
freq1 – frequency of the effect allele;
beta1 – log hazard protection ratio (negation of log hazard ratio) for carrying one copy of the a1 allele in self;
se – Standard Error;
p – the P value for the Wald test of association between imputed dosage and cox model residual; direction – direction of effect in UK Biobank and LifeGen;
info – imputation info score;
freq_se – standard error of effect allele frequency between UK Biobank and LifeGen;
min_freq1 – minimum effect allele frequency across UK Biobank and LifeGen;
max_freq1 – maximum effect allele frequency across UK Biobank and LifeGen.

Source: Edinburgh DataShare doi:10.7488/ds/2463 (file lifegen_phase2_bothpl_alldr_2017_09_18.tsv.gz).
Inspected header 2026-10-28:
rsid	snpid	chr	pos	a1	a0	n	freq1	beta1	se	p	direction	info	freq_se	min_freq1	max_freq1
