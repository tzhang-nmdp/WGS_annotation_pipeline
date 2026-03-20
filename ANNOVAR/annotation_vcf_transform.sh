
#!/bin/sh

vcf=$1
outdir=$2
snpeff_dir=$3
ref_fasta_dir=$4
ref_fasta=${ref_fasta_dir}/"hs38.fasta" 
annovar_dir="/home/centos/annovar"
gene_list_file="/home/centos/candidate_gene_list"
cutoff=$5
idx=7


if [ ! ${cutoff} ]; then
    echo "cutoff invalid and set to default 2"
    cutoff=2
    # exit 1
fi

# Check if required arguments are provided
if [ -z "$vcf" ] || [ -z "$outdir" ]; then
    echo "Error: Missing required arguments"
    echo "Usage: $0 <vcf> <outdir> "
    exit 1
fi

# create sim file for filtering
j=$(awk 'FNR==1{print NF}'  ${outdir}/${vcf}.sim.multi_allele_split.annovar.input.db.csv.variant_annotation_filtering.csv)
awk -v a="$j" -v b="$cutoff" '($a>=b)&&($a!=10){print $0}' ${outdir}/${vcf}.sim.multi_allele_split.annovar.input.db.csv.variant_annotation_filtering.csv > ${outdir}/${vcf}.sim.multi_allele_split.annovar.input.db.csv.variant_annotation_filtering.csv.select.${cutoff}

# create db id file for filtering
awk 'FNR>1{split($6,a,":"); print a[2]":"a[3]":"a[4]":"a[5]"\t"$42"_"$43"_"$44"_"$45}' ${outdir}/${vcf}.sim.multi_allele_split.annovar.input.db.csv.variant_annotation_filtering.csv.select.${cutoff} > ${outdir}/${vcf}.sim.multi_allele_split.annovar.input.db.csv.variant_annotation_filtering.csv.select.${cutoff}.id
awk 'FNR>1{split($6,a,":"); print a[2]":"a[3]":"a[4]":"a[5]"\t"$42"_"$43"_"$44"_"$45}' ${outdir}/${vcf}.sim.multi_allele_split.annovar.input.db.csv.variant_annotation_filtering.csv > ${outdir}/${vcf}.sim.multi_allele_split.annovar.input.db.csv.variant_annotation_filtering.csv.sim.id

# filter snpeff vcf files with lof + gnomad
awk 'FNR==NR{a["chr"$3":"$4":"$6":"$7];next} (FNR==1)||(!($10 in a)){print $0}' ${outdir}/annovarg/${vcf}.sim.multi_allele_split.annovar.input.gnomad.all_0.001ab.hg38_generic_dropped ${outdir}/${vcf}.multi_bi.vcf.sim.snpEff.id.all.vcf > ${outdir}/${vcf}.multi_bi.vcf.sim.snpEff.id.all.somatic.vcf
awk 'FNR==NR{a[$1];next} {split($8,b,"|"); if (($1~"#")||($3 in a)||(b[3]=="HIGH")) {print $0}}' ${outdir}/${vcf}.sim.multi_allele_split.annovar.input.db.csv.variant_annotation_filtering.csv.select.${cutoff}.id ${outdir}/${vcf}.multi_bi.vcf.sim.snpEff.id.all.somatic.vcf > ${outdir}/${vcf}.multi_bi.vcf.sim.snpEff.id.anno_sel_lof.somatic.${cutoff}.vcf
awk 'FNR==NR{a[$1];b[$1]=$2; next} {if ($1~"#") {print $0} else if ($3 in a) {$8=$8"|"b[$3]; print $0}}' ${outdir}/${vcf}.sim.multi_allele_split.annovar.input.db.csv.variant_annotation_filtering.csv.sim.id ${outdir}/${vcf}.multi_bi.vcf.sim.snpEff.id.anno_sel_lof.somatic.${cutoff}.vcf > ${outdir}/${vcf}.multi_bi.vcf.sim.snpEff.id.anno_sel_lof.somatic.${cutoff}.vcf.anno
mv ${outdir}/${vcf}.multi_bi.vcf.sim.snpEff.id.anno_sel_lof.somatic.${cutoff}.vcf.anno ${outdir}/${vcf}.multi_bi.vcf.sim.snpEff.id.anno_sel_lof.somatic.${cutoff}.vcf

# run maf transformation
i=$(awk '($1!~"#"){print NF}' ${outdir}/${vcf}.multi_bi.vcf.sim.snpEff.id.anno_sel_lof.somatic.${cutoff}.vcf | head -1)
sed 's% %\t%g' ${outdir}/${vcf}.multi_bi.vcf.sim.snpEff.id.anno_sel_lof.somatic.${cutoff}.vcf > ${outdir}/${vcf}.multi_bi.vcf.sim.snpEff.id.anno_sel_lof.somatic.${cutoff}.vcf.tmp
python ${snpeff_dir}/snpeff_to_maf.py ${outdir}/${vcf}.multi_bi.vcf.sim.snpEff.id.anno_sel_lof.somatic.${cutoff}.vcf.tmp ${outdir}/${vcf}.multi_bi.vcf.sim.snpEff.id.anno_sel_lof.somatic.${cutoff}.vcf.maf ${outdir}/${vcf}.multi_bi.vcf.sim.snpEff.id.anno_sel_lof.somatic.${cutoff}.vcf.error $i
rm ${outdir}/${vcf}.multi_bi.vcf.sim.snpEff.id.anno_sel_lof.somatic.${cutoff}.vcf.tmp
awk 'FNR==NR{a[$1];b[$1]=$2; next} {if (FNR==1) {print $0"\tannovar"} else if ($10 in a) {print $0"\t"b[$10]}}' ${outdir}/${vcf}.sim.multi_allele_split.annovar.input.db.csv.variant_annotation_filtering.csv.sim.id ${outdir}/${vcf}.multi_bi.vcf.sim.snpEff.id.anno_sel_lof.somatic.${cutoff}.vcf.maf > ${outdir}/${vcf}.multi_bi.vcf.sim.snpEff.id.anno_sel_lof.somatic.${cutoff}.vcf.maf.anno
mv ${outdir}/${vcf}.multi_bi.vcf.sim.snpEff.id.anno_sel_lof.somatic.${cutoff}.vcf.maf.anno ${outdir}/${vcf}.multi_bi.vcf.sim.snpEff.id.anno_sel_lof.somatic.${cutoff}.vcf.maf

# combine bi-allele vcf into multi-allele vcf
(grep "^#" ${outdir}/${vcf}.multi_bi.vcf.sim.snpEff.id.anno_sel_lof.somatic.${cutoff}.vcf; grep -v "^#" ${outdir}/${vcf}.multi_bi.vcf.sim.snpEff.id.anno_sel_lof.somatic.${cutoff}.vcf | sed 's% %\t%g' | sort -T /home/centos/tmp/ -k1,1 -k2,2n) | bgzip  -c >  ${outdir}/${vcf}.multi_bi.vcf.sim.snpEff.id.anno_sel_lof.somatic.${cutoff}.vcf.gz
tabix -p vcf ${outdir}/${vcf}.multi_bi.vcf.sim.snpEff.id.anno_sel_lof.somatic.${cutoff}.vcf.gz
bcftools norm -m +both  -a -c w --atom-overlaps . -f ${ref_fasta} ${outdir}/${vcf}.multi_bi.vcf.sim.snpEff.id.anno_sel_lof.somatic.${cutoff}.vcf.gz -o ${outdir}/${vcf}.multi_bi_multi.vcf.sim.snpEff.id.anno_sel_lof.somatic.${cutoff}.vcf

# filter multi-allele vcf with candidate gene list
grep \#\# ${outdir}/${vcf}.multi_bi.vcf.sim.snpEff.id.anno_sel_lof.somatic.${cutoff}.vcf > ${outdir}/${vcf}.multi_bi.vcf.sim.snpEff.id.anno_sel_lof.somatic.${cutoff}.vcf.tmp.header
sed 's%#CHROM%CHROM%g' ${outdir}/${vcf}.multi_bi.vcf.sim.snpEff.id.anno_sel_lof.somatic.${cutoff}.vcf | sed 's% %\t%g' > ${outdir}/${vcf}.multi_bi.vcf.sim.snpEff.id.anno_sel_lof.somatic.${cutoff}.vcf.tmp
python ${annovar_dir}/variable_mapping.py  ${outdir}/${vcf}.multi_bi.vcf.sim.snpEff.id.anno_sel_lof.somatic.${cutoff}.vcf.tmp ${gene_list_file} ${outdir}/${vcf}.multi_bi.vcf.sim.snpEff.id.anno_sel_lof.somatic.${cutoff}.candidate_gene.vcf.tmp $idx
cat ${outdir}/${vcf}.multi_bi.vcf.sim.snpEff.id.anno_sel_lof.somatic.${cutoff}.vcf.tmp.header ${outdir}/${vcf}.multi_bi.vcf.sim.snpEff.id.anno_sel_lof.somatic.${cutoff}.candidate_gene.vcf.tmp > ${outdir}/${vcf}.multi_bi.vcf.sim.snpEff.id.anno_sel_lof.somatic.${cutoff}.candidate_gene.vcf
python ${annovar_dir}/variable_mapping.py  ${outdir}/${vcf}.multi_bi.vcf.sim.snpEff.id.anno_sel_lof.somatic.${cutoff}.vcf.maf ${gene_list_file} ${outdir}/${vcf}.multi_bi.vcf.sim.snpEff.id.anno_sel_lof.somatic.${cutoff}.candidate_gene.vcf.maf $idx
awk 'FNR==NR{a[$1];b[$1]=$2; next} {if (FNR==1) {print $0"\tannovar"} else if ($10 in a) {print $0"\t"b[$10]}}' ${outdir}/${vcf}.sim.multi_allele_split.annovar.input.db.csv.variant_annotation_filtering.csv.sim.id ${outdir}/${vcf}.multi_bi.vcf.sim.snpEff.id.anno_sel_lof.somatic.${cutoff}.candidate_gene.vcf.maf > ${outdir}/${vcf}.multi_bi.vcf.sim.snpEff.id.anno_sel_lof.somatic.${cutoff}.candidate_gene.vcf.maf.anno
mv ${outdir}/${vcf}.multi_bi.vcf.sim.snpEff.id.anno_sel_lof.somatic.${cutoff}.candidate_gene.vcf.maf.anno ${outdir}/${vcf}.multi_bi.vcf.sim.snpEff.id.anno_sel_lof.somatic.${cutoff}.candidate_gene.vcf.maf
rm *.tmp*

