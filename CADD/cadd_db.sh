
#!/bin/sh
# Directory path containing the CADD database files
# Expected format: Path to folder containing CADD-related database files
# Input parameter: Third command-line argument ($3)

vcf=$1
cadd_db_dir=$2
outdir=$3
chr=$4

# Process VCF file by chromosome (1-22, X, Y) to separate SNVs and INDELs
# For each chromosome:
# 1. Extract SNVs (single nucleotide variants) where ref and alt are single bases
#    Format: CHROM:POS:REF:ALT
# 2. Extract INDELs (insertions/deletions) where ref or alt length > 1
#    Format: CHROM:POS:REF:ALT
# Input: ${vcf} file in VCF format
# Output: Creates chromosome-specific files for SNVs and INDELs
#         - ${vcf}.cadd.input.snv.${chr}
#         - ${vcf}.cadd.input.indel.${chr}
mkdir ${outdir}/cadd

    awk -v a="$chr" '($1==a)&&(length($4)==1)&&(length($5)==1){print $1":"$2":"$4":"$5}' ${outdir}/${vcf}.sim.multi_allele_split.annovar.input > ${outdir}/cadd/${vcf}.sim.multi_allele_split.annovar.input.cadd.input.snv.${chr}
    awk -v a="$chr" '($1==a)&&((length($4)!=1)||(length($5)!=1)){print $1":"$2":"$4":"$5}' ${outdir}/${vcf}.sim.multi_allele_split.annovar.input > ${outdir}/cadd/${vcf}.sim.multi_allele_split.annovar.input.cadd.input.indel.${chr}
echo "1" && date
#!/bin/bash
# This script processes CADD (Combined Annotation Dependent Depletion) scores for variants
# across all chromosomes (1-22, X, Y)
#
# For each chromosome, the script:
# 1. Extracts CADD scores for INDELs by matching variants with pre-computed scores
# 2. Identifies missing INDEL scores by comparing input against matched scores
# 3. Extracts CADD scores for SNVs by matching variants with pre-computed scores
# 4. Identifies missing SNV scores by comparing input against matched scores
#
# Required variables:
# - outdir: Output directory path
# - vcf: Input VCF file basename
# - cadd_db_dir: Directory containing CADD score databases
#
# Input files expected:
# - ${outdir}/cadd/${vcf}.sim.multi_allele_split.annovar.input.cadd.input.indel.${chr}
# - ${outdir}/cadd/${vcf}.sim.multi_allele_split.annovar.input.cadd.input.snv.${chr}
# - ${cadd_db_dir}/gnomad.genomes.r3.0.indel.chr${chr}.sim
# - ${cadd_db_dir}/whole_genome_SNVs.tsv.chr${chr}.sim
#
# Output files generated:
# - ${outdir}/cadd/${vcf}.sim.multi_allele_split.annovar.input.cadd.score.indel.${chr}
# - ${outdir}/cadd/${vcf}.sim.multi_allele_split.annovar.input.cadd.miss.indel.${chr}
# - ${outdir}/cadd/${vcf}.sim.multi_allele_split.annovar.input.cadd.score.snv.${chr}
# - ${outdir}/cadd/${vcf}.sim.multi_allele_split.annovar.input.cadd.miss.snv.${chr}

i=$chr

    awk 'FNR==NR{a[$1];next}($1 in a){print $0}' ${outdir}/cadd/${vcf}.sim.multi_allele_split.annovar.input.cadd.input.indel.${i} ${cadd_db_dir}/gnomad.genomes.r3.0.indel.chr${i}.sim > ${outdir}/cadd/${vcf}.sim.multi_allele_split.annovar.input.cadd.score.indel.${i}
    awk 'FNR==NR{a[$1];next}!($1 in a){print $0}' ${outdir}/cadd/${vcf}.sim.multi_allele_split.annovar.input.cadd.score.indel.${i} ${outdir}/cadd/${vcf}.sim.multi_allele_split.annovar.input.cadd.input.indel.${i} > ${outdir}/cadd/${vcf}.sim.multi_allele_split.annovar.input.cadd.miss.indel.${i}
    awk 'FNR==NR{a[$1];next}($1 in a){print $0}' ${outdir}/cadd/${vcf}.sim.multi_allele_split.annovar.input.cadd.input.snv.${i} ${cadd_db_dir}/whole_genome_SNVs.tsv.chr${i}.sim > ${outdir}/cadd/${vcf}.sim.multi_allele_split.annovar.input.cadd.score.snv.${i}
    awk 'FNR==NR{a[$1];next}!($1 in a){print $0}' ${outdir}/cadd/${vcf}.sim.multi_allele_split.annovar.input.cadd.score.snv.${i} ${outdir}/cadd/${vcf}.sim.multi_allele_split.annovar.input.cadd.input.snv.${i} > ${outdir}/cadd/${vcf}.sim.multi_allele_split.annovar.input.cadd.miss.snv.${i}    

echo "2" && date
cat ${outdir}/cadd/${vcf}.sim.multi_allele_split.annovar.input.cadd.score.indel.* >> ${outdir}/cadd/${vcf}.sim.multi_allele_split.annovar.input.cadd.score.indel_all
cat ${outdir}/cadd/${vcf}.sim.multi_allele_split.annovar.input.cadd.miss.indel.* >> ${outdir}/cadd/${vcf}.sim.multi_allele_split.annovar.input.cadd.miss.indel_all
cat ${outdir}/cadd/${vcf}.sim.multi_allele_split.annovar.input.cadd.score.snv.* >> ${outdir}/cadd/${vcf}.sim.multi_allele_split.annovar.input.cadd.score.snv_all
cat ${outdir}/cadd/${vcf}.sim.multi_allele_split.annovar.input.cadd.miss.snv.* >> ${outdir}/cadd/${vcf}.sim.multi_allele_split.annovar.input.cadd.miss.snv_all

echo "3" && date
cat ${outdir}/cadd/${vcf}.sim.multi_allele_split.annovar.input.cadd.miss.indel_all ${outdir}/cadd/${vcf}.sim.multi_allele_split.annovar.input.cadd.miss.snv_all > ${outdir}/cadd/${vcf}.sim.multi_allele_split.annovar.input.cadd.miss_all
cat ${outdir}/cadd/${vcf}.sim.multi_allele_split.annovar.input.cadd.score.indel_all ${outdir}/cadd/${vcf}.sim.multi_allele_split.annovar.input.cadd.score.snv_all > ${outdir}/cadd/${vcf}.sim.multi_allele_split.annovar.input.cadd.score_all
echo "4" && date

# python  ${annovar_dir}/annovar_db_table.py ${outdir}/${vcf}.multi_allele_split.annovar.input cadd cadd

