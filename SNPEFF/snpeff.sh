
#!/bin/bash
vcf=$1 
snpeff_dir=$2
outdir=$3

java -d64 -Xmx4g -jar ${snpeff_dir}/snpEff.jar -v GRCh38.86 ${outdir}/${vcf}.multi_bi.vcf.sim -canon > ${outdir}/${vcf}.multi_bi.vcf.sim.snpEff.vcf

awk '{if ($1~"#CHROM"){print $0} else if ($1!~"##"){$3=$1":"$2":"$4":"$5; print $0}}' ${outdir}/${vcf}.multi_bi.vcf.sim.snpEff.vcf > ${outdir}/${vcf}.multi_bi.vcf.sim.snpEff.id.vcf.rh
awk '($1~"##"){print $0}' ${outdir}/${vcf}.multi_bi.vcf.sim.snpEff.vcf > ${outdir}/${vcf}.multi_bi.vcf.sim.snpEff.id.vcf.header
j=$(awk '($1!~"#"){print NF}' ${outdir}/${vcf}.multi_bi.vcf | head -1)
awk '{if ($1~"#CHROM"){print $0} else if ($1!~"##") {$3=$1":"$2":"$4":"$5; print $0}}' ${outdir}/${vcf}.multi_bi.vcf | sed 's% %\t%g' | cut -d $'\t' -f 3,10-${j} > ${outdir}/${vcf}.multi_bi.id.vcf.rh
j1=$((j+1))
awk 'FNR==NR{a[$1]; b[$1]=$0; next}($3 in a) {$10=b[$3]; print $0}' ${outdir}/${vcf}.multi_bi.id.vcf.rh ${outdir}/${vcf}.multi_bi.vcf.sim.snpEff.id.vcf.rh | sed 's% %\t%g' | cut -d $'\t' -f 1-9,11-${j1} >  ${outdir}/${vcf}.multi_bi.vcf.sim.snpEff.id.all.vcf.rh
rm -rf  ${outdir}/${vcf}.multi_bi.id.vcf.rh
rm -rf ${outdir}/${vcf}.multi_bi.vcf.sim.snpEff.id.vcf.rh
cat ${outdir}/${vcf}.multi_bi.vcf.sim.snpEff.id.vcf.header  ${outdir}/${vcf}.multi_bi.vcf.sim.snpEff.id.all.vcf.rh > ${outdir}/${vcf}.multi_bi.vcf.sim.snpEff.id.all.vcf
rm -rf  ${outdir}/${vcf}.multi_bi.vcf.sim.snpEff.id.all.vcf.rh
i=$(awk '($1!~"#"){print NF}' ${outdir}/${vcf}.multi_bi.vcf.sim.snpEff.id.all.vcf | head -1)
echo $i
python ${snpeff_dir}/snpeff_to_maf.py ${outdir}/${vcf}.multi_bi.vcf.sim.snpEff.id.all.vcf ${outdir}/${vcf}.multi_bi.vcf.sim.snpEff.id.all.vcf.maf ${outdir}/${vcf}.multi_bi.vcf.sim.snpEff.all.vcf.error $i

