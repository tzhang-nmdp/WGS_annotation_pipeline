
#!/bin/sh
#####################################################################
# chr1	16356650	ACT	A " to " 1	16356650	16356652	ACT	A info:
#####################################################################
# 1	1207339	1207761	0	-	info:chr1:1207339:1207761:DEL:DEL0055SUR    
#####################################################################
# example  sh /home/tzhang/MDS_data/17.snv_annovar.sh NMDP188_all_somatic.snv.annovar  /home/tzhang/database/annovar /home/tzhang/MDS_data/Reanalyzed/all_variant_new/somatic rep all

#  awk '{split($8,a,";"); print a[1]"\t"$2"\t"$1}' ${vcf}.sim.multi_allele_split.annovar.input.${db} |sed 's%info:%%g'  > ${vcf}.sim.multi_allele_split.annovar.input.${db}.sim
###############################################################################
vcf=$1
annovar_dir=$2
outdir=$3
ref_fasta_dir=$4
ref_fasta=${ref_fasta_dir}/"hs38.fasta" 
con=$5
con2=$6

if [ -n "${con2}" ] && ([ "${con2}" == "tf" ] || [ "${con2}" == "all" ]); then
    # split multi-allele vcf into bi-allele vcf
    echo " bcftools norm -m -both  -a -c w --atom-overlaps . -f ${ref_fasta} ${outdir}/${vcf} -o ${outdir}/${vcf}.multi_bi.vcf"
    bcftools norm -m -both  -a -c w --atom-overlaps . -f ${ref_fasta} ${outdir}/${vcf} -o ${outdir}/${vcf}.multi_bi.vcf
    awk '{if ($1~"##") {print $0} else {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}}' ${outdir}/${vcf}.multi_bi.vcf >  ${outdir}/${vcf}.multi_bi.vcf.sim
    awk '($1~"^chr"){a=substr($1,4,length($1)); b=$2+length($4)-1; print a"\t"$2"\t"b"\t"$4"\t"$5"\tinfo:"$1":"$2":"$4":"$5}'  ${outdir}/${vcf}.multi_bi.vcf.sim | awk '!a[$6]++' > ${outdir}/${vcf}.sim.multi_allele_split.annovar.input 
   mkdir -p ${outdir}/annovargg
   mkdir -p ${outdir}/annovarr
   mkdir -p ${outdir}/annovarg
fi

if [ -n "${con2}" ] && ([ "${con2}" == "ele" ] || [ "${con2}" == "all" ]); then
    for db in gwasCatalog dbSUPER distal_hQTLs dsQTL GeneHancer local_hQTLs otherDB_eQTLs SEA  UCNEbase wgEncodeRegDnaseClustered  EnhancerAtlas JEME encRegTfbsClustered # ENCODE+RoadMap1 ENCODE+RoadMap2 EpiTensor # dgvMerged
    do
    #     if [ -n "${con}" ] && ([ "${con}" == "af" ] || [ "${con}" == "all" ]); then
    #    ${annovar_dir}/annotate_variation.pl --filter -dbtype generic -genericdbfile hg38_${db}.txt -buildver hg38   -colsWanted all -out ${outdir}/annovarg/${vcf}.sim.multi_allele_split.annovar.input.${db} ${outdir}/${vcf}.sim.multi_allele_split.annovar.input ${annovar_dir}/humandb/ 
    #     fi
        if [ -n "${con}" ] && ([ "${con}" == "ar" ] || [ "${con}" == "all" ]); then
       ${annovar_dir}/annotate_variation.pl -regionanno --buildver hg38 -dbtype ${db}_sv -colsWanted all ${outdir}/${vcf}.sim.multi_allele_split.annovar.input ${annovar_dir}/humandb/ -out ${outdir}/annovarr/${vcf}.sim.multi_allele_split.annovar.input.${db} 
        fi
    #     if [ -n "${con}" ] && ([ "${con}" == "ag" ] || [ "${con}" == "all" ]); then
    #    ${annovar_dir}/annotate_variation.pl --geneanno -dbtype ${db} -buildver hg38  -out ${outdir}/annovargg/${vcf}.sim.multi_allele_split.annovar.input.${db} ${outdir}/${vcf}.sim.multi_allele_split.annovar.input ${annovar_dir}/humandb/ 
    #     fi
    done
fi

if [ -n "${con2}" ] && ([ "${con2}" == "ele" ] || [ "${con2}" == "all" ]); then
    for db in GTEx_eQTL PancanQTL
    do
        if [ -n "${con}" ] && ([ "${con}" == "af" ] || [ "${con}" == "all" ]); then
       ${annovar_dir}/annotate_variation.pl --filter -dbtype generic -genericdbfile hg38_${db}.txt -buildver hg38   -colsWanted all -out ${outdir}/annovarg/${vcf}.sim.multi_allele_split.annovar.input.${db} ${outdir}/${vcf}.sim.multi_allele_split.annovar.input ${annovar_dir}/humandb/ 
        fi
    #     if [ -n "${con}" ] && ([ "${con}" == "ar" ] || [ "${con}" == "all" ]); then
    #    ${annovar_dir}/annotate_variation.pl -regionanno --buildver hg38 -dbtype ${db}_sv -colsWanted all ${outdir}/${vcf}.sim.multi_allele_split.annovar.input ${annovar_dir}/humandb/ -out ${outdir}/annovarr/${vcf}.sim.multi_allele_split.annovar.input.${db} 
    #     fi
    #     if [ -n "${con}" ] && ([ "${con}" == "ag" ] || [ "${con}" == "all" ]); then
    #    ${annovar_dir}/annotate_variation.pl --geneanno -dbtype ${db} -buildver hg38  -out ${outdir}/annovargg/${vcf}.sim.multi_allele_split.annovar.input.${db} ${outdir}/${vcf}.sim.multi_allele_split.annovar.input ${annovar_dir}/humandb/ 
    #     fi
    done
fi

if [ -n "${con2}" ] && ([ "${con2}" == "rep" ] || [ "${con2}" == "all" ]); then
    for db in repeatmasker genomicSuperDups simpleRepeat
    do
    #     if [ -n "${con}" ] && ([ "${con}" == "af" ] || [ "${con}" == "all" ]); then
    #    ${annovar_dir}/annotate_variation.pl --filter -dbtype generic -genericdbfile hg38_${db}.txt -buildver hg38   -colsWanted all -out ${outdir}/annovarg/${vcf}.sim.multi_allele_split.annovar.input.${db} ${outdir}/${vcf}.sim.multi_allele_split.annovar.input ${annovar_dir}/humandb/
    #     fi
        
        if [ -n "${con}" ] && ([ "${con}" == "ar" ] || [ "${con}" == "all" ]); then
       ${annovar_dir}/annotate_variation.pl -regionanno --buildver hg38 -dbtype ${db}_sv -colsWanted all ${outdir}/${vcf}.sim.multi_allele_split.annovar.input ${annovar_dir}/humandb/ -out ${outdir}/annovarr/${vcf}.sim.multi_allele_split.annovar.input.${db} 
        fi
    #     if [ -n "${con}" ] && ([ "${con}" == "ag" ] || [ "${con}" == "all" ]); then
    #    ${annovar_dir}/annotate_variation.pl --geneanno -dbtype ${db} -buildver hg38  -out ${outdir}/annovargg/${vcf}.sim.multi_allele_split.annovar.input.${db} ${outdir}/${vcf}.sim.multi_allele_split.annovar.input ${annovar_dir}/humandb/ 
    #     fi
    done
fi

if [ -n "${con}2" ] && ([ "${con2}" == "var" ] || [ "${con2}" == "all" ]); then
    for db in gnomad_all gnomad.all_0.001ab avsnp142 ICGC+TCGA+COSMIC+ClinVar1 ICGC+TCGA+COSMIC+ClinVar2 cosmic70 HGMD clinvar.benign clinvar.pathogenic # clinvar_20160302
    do
        if [ -n "${con}" ] && ([ "${con}" == "af" ] || [ "${con}" == "all" ]); then
       ${annovar_dir}/annotate_variation.pl --filter -dbtype generic -genericdbfile hg38_${db}.txt -buildver hg38   -colsWanted all -out ${outdir}/annovarg/${vcf}.sim.multi_allele_split.annovar.input.${db} ${outdir}/${vcf}.sim.multi_allele_split.annovar.input ${annovar_dir}/humandb/ 
        fi
    #     if [ -n "${con}" ] && ([ "${con}" == "ar" ] || [ "${con}" == "all" ]); then
    #    ${annovar_dir}/annotate_variation.pl -regionanno --buildver hg38 -dbtype ${db}_sv -colsWanted all ${outdir}/${vcf}.sim.multi_allele_split.annovar.input ${annovar_dir}/humandb/ -out ${outdir}/annovarr/${vcf}.sim.multi_allele_split.annovar.input.${db} 
    #     fi
    #     if [ -n "${con}" ] && ([ "${con}" == "ag" ] || [ "${con}" == "all" ]); then
    #    ${annovar_dir}/annotate_variation.pl --geneanno -dbtype ${db} -buildver hg38  -out ${outdir}/annovargg/${vcf}.sim.multi_allele_split.annovar.input.${db} ${outdir}/${vcf}.sim.multi_allele_split.annovar.input ${annovar_dir}/humandb/ 
    #     fi
    done
fi

if [ -n "${con2}" ] && ([ "${con2}" == "gene" ] || [ "${con2}" == "all" ]); then
   ${annovar_dir}/annotate_variation.pl  --geneanno --buildver hg38 ${outdir}/${vcf}.sim.multi_allele_split.annovar.input ${annovar_dir}/humandb/ -dbtype refGene -out ${outdir}/annovargg/${vcf}.sim.multi_allele_split.annovar.input.refGene
   ${annovar_dir}/annotate_variation.pl  --geneanno --buildver hg38 ${outdir}/${vcf}.sim.multi_allele_split.annovar.input ${annovar_dir}/humandb/ -dbtype knownGene -out  ${outdir}/annovargg/${vcf}.sim.multi_allele_split.annovar.input.knownGene
fi 

if [ -n "${con2}" ] && ([ "${con2}" == "sco" ] || [ "${con2}" == "all" ]); then
    for db in phyloP100way phastConsElements100way spliceai_score_wgs_exome_0.5
    do
        if [ -n "${con}" ] && ([ "${con}" == "af" ] || [ "${con}" == "all" ]); then
       ${annovar_dir}/annotate_variation.pl --filter -dbtype generic -genericdbfile hg38_${db}.txt -buildver hg38   -colsWanted all -out ${outdir}/annovarg/${vcf}.sim.multi_allele_split.annovar.input.${db} ${outdir}/${vcf}.sim.multi_allele_split.annovar.input ${annovar_dir}/humandb/ 
        fi
         if [ -n "${con}" ] && ([ "${con}" == "ar" ] || [ "${con}" == "all" ]); then
        ${annovar_dir}/annotate_variation.pl -regionanno --buildver hg38 -dbtype ${db}_sv -colsWanted all ${outdir}/${vcf}.sim.multi_allele_split.annovar.input ${annovar_dir}/humandb/ -out ${outdir}/annovarr/${vcf}.sim.multi_allele_split.annovar.input.${db} 
         fi
    #     if [ -n "${con}" ] && ([ "${con}" == "ag" ] || [ "${con}" == "all" ]); then
    #    ${annovar_dir}/annotate_variation.pl --geneanno -dbtype ${db} -buildver hg38  -out ${outdir}/annovargg/${vcf}.sim.multi_allele_split.annovar.input.${db} ${outdir}/${vcf}.sim.multi_allele_split.annovar.input ${annovar_dir}/humandb/ 
        #fi
    done
fi
 
