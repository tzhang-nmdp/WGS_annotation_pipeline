#!/bin.bash
vcf=$1
annovar_dir=$2
cadd_dir=$3
snpeff_dir=$4
ref_fasta_dir=$5
outdir=$6
con=$7
cutoff=$8

# create vcf annotation input
if [ -n "${con}" ] && ([ "${con}" == "input" ] || [ "${con}" == "all" ]); then
sh ${annovar_dir}/annovar.sh ${vcf} ${annovar_dir} ${outdir} ${ref_fasta_dir} tf tf
fi

# run annovar annotations
if [ -n "${con}" ] && ([ "${con}" == "annovar" ] || [ "${con}" == "all" ]); then
sh ${annovar_dir}/annovar.sh ${vcf} ${annovar_dir} ${outdir} ${ref_fasta_dir} all all
fi

# run annovar annotations
if [ -n "${con}" ] && ([ "${con}" == "snpeff" ] || [ "${con}" == "all" ]); then
sh ${snpeff_dir}/snpeff.sh ${vcf} ${snpeff_dir} ${outdir}
fi

# run cadd annotations
if [ -n "${con}" ] && ([ "${con}" == "cadd" ] || [ "${con}" == "all" ]); then
sh ${cadd_dir}/cadd_db.sh ${vcf} ${cadd_dir} ${outdir}
fi

# run annotation combining
if [ -n "${con}" ] && ([ "${con}" == "combine" ] || [ "${con}" == "all" ]); then
sh ${annovar_dir}/annotation_combine_filter.sh  ${vcf} ${annovar_dir} ${outdir} all
fi

# run annotation db filtering
if [ -n "${con}" ] && ([ "${con}" == "filter" ] || [ "${con}" == "all" ]); then
sh ${annovar_dir}/annotation_combine_filter.sh  ${vcf} ${annovar_dir} ${outdir} filter
fi

# run annotation vcf transform
if [ -n "${con}" ] && ([ "${con}" == "vcf_tf" ] || [ "${con}" == "all" ]); then
sh ${annovar_dir}/annotation_vcf_transform.sh  ${vcf} ${outdir} ${snpeff_dir} ${ref_fasta_dir} ${cutoff}
fi
