
#!/bin/sh
#####################################################################
# chr1	16356650	ACT	A " to " 1	16356650	16356652	ACT	A info:
#####################################################################
# 1	1207339	1207761	0	-	info:chr1:1207339:1207761:DEL:DEL0055SUR    
#####################################################################
# example  sh /home/tzhang/MDS_data/17.snv_annovar.sh NMDP188_all_somatic.snv.annovar  /home/tzhang/database/annovar /home/tzhang/MDS_data/Reanalyzed/all_variant_new/somatic rep all

#  awk '{split($8,a,";"); print a[1]"\t"$2"\t"$1}' ${vcf}.multi_allele_split.annovar.input.${db} |sed 's%info:%%g'  > ${vcf}.multi_allele_split.annovar.input.${db}.sim
###############################################################################
vcf=$1
annovar_dir=$2
outdir=$3
dbs=$4
db_type=$5  # generic, region, gene, cadd
# Check if required arguments are provided
if [ -z "$vcf" ] || [ -z "$annovar_dir" ] || [ -z "$outdir" ]; then
    echo "Error: Missing required arguments"
    echo "Usage: $0 <vcf> <annovar_dir> <outdir> [db] [db_type]"
    exit 1
fi

# awk '($1~"^chr"){if ($5~",") {split($5,a,","); for (i in a){$5=a[i];print $0}} else {print $0}}' ${outdir}/${vcf}.sim | awk '($5!="\*"){print $0}'  > ${outdir}/${vcf}.sim.multi_allele_split # sed 's%\*%-%g' > ${outdir}/${vcf}.sim.multi_allele_split
# awk '($1~"^chr"){a=substr($1,4,length($1)); b=$2+length($4)-1; print a"\t"$2"\t"b"\t"$4"\t"$5"\tinfo:"$1":"$2":"$4":"$5}'  ${outdir}/${vcf}.sim.multi_allele_split | awk '!a[$6]++' > ${outdir}/${vcf}.sim.multi_allele_split.annovar.input # if ((length($4)>1) && (length($4)>length($5))){$5="-"};

if [ -n "${dbs}" ] && [ -n "${db_type}" ] && [ "${dbs}" != "all" ] && [ "${dbs}" != "filter" ]; then
    python  ${annovar_dir}/annovar_db_table.py ${outdir}/${vcf}.sim.multi_allele_split.annovar.input ${db} ${db_type}
    cp ${outdir}/${vcf}.sim.multi_allele_split.annovar.input.db.csv ${outdir}/${vcf}.sim.multi_allele_split.annovar.input

elif [ -n "${dbs}" ] && [ "${dbs}" == "all" ]; then

    cat /home/centos/BMT/hd.tmp ${outdir}/${vcf}.sim.multi_allele_split.annovar.input > ${outdir}/${vcf}.sim.multi_allele_split.annovar.input.rh
    cp ${outdir}/${vcf}.sim.multi_allele_split.annovar.input ${outdir}/${vcf}.sim.multi_allele_split.annovar.input.tmp
    mv ${outdir}/${vcf}.sim.multi_allele_split.annovar.input.rh ${outdir}/${vcf}.sim.multi_allele_split.annovar.input

    for db in gwasCatalog dbSUPER distal_hQTLs dsQTL GeneHancer local_hQTLs otherDB_eQTLs SEA UCNEbase wgEncodeRegDnaseClustered  EnhancerAtlas JEME encRegTfbsClustered phyloP100way phastConsElements100way spliceai_score_wgs_exome_0.5
    do
    python  ${annovar_dir}/annovar_db_table.py ${outdir}/${vcf}.sim.multi_allele_split.annovar.input ${db} region
    cp ${outdir}/${vcf}.sim.multi_allele_split.annovar.input.db.csv ${outdir}/${vcf}.sim.multi_allele_split.annovar.input
    done


    for db in repeatmasker genomicSuperDups simpleRepeat
    do
    python  ${annovar_dir}/annovar_db_table.py ${outdir}/${vcf}.sim.multi_allele_split.annovar.input ${db} region
    cp ${outdir}/${vcf}.sim.multi_allele_split.annovar.input.db.csv ${outdir}/${vcf}.sim.multi_allele_split.annovar.input
    done


    for db in knownGene refGene 
    do
    python  ${annovar_dir}/annovar_db_table.py ${outdir}/${vcf}.sim.multi_allele_split.annovar.input ${db} gene
    cp ${outdir}/${vcf}.sim.multi_allele_split.annovar.input.db.csv ${outdir}/${vcf}.sim.multi_allele_split.annovar.input
    done

    for db in gnomad_all gnomad.all_0.001ab avsnp142 'ICGC+TCGA+COSMIC+ClinVar1' 'ICGC+TCGA+COSMIC+ClinVar2' cosmic70 HGMD clinvar.benign clinvar.pathogenic GTEx_eQTL PancanQTL 
    do
    python  ${annovar_dir}/annovar_db_table.py ${outdir}/${vcf}.sim.multi_allele_split.annovar.input ${db} generic
    cp ${outdir}/${vcf}.sim.multi_allele_split.annovar.input.db.csv ${outdir}/${vcf}.sim.multi_allele_split.annovar.input
    done

    for db in cadd
    do
    python  ${annovar_dir}/annovar_db_table.py ${outdir}/${vcf}.sim.multi_allele_split.annovar.input ${db} cadd
    cp ${outdir}/${vcf}.sim.multi_allele_split.annovar.input.db.csv ${outdir}/${vcf}.sim.multi_allele_split.annovar.input
    done
fi

if [ -n "${dbs}" ] && ([ "${dbs}" == "filter" ] || [ "${dbs}" == "all" ]); then
    python ${annovar_dir}/variant_annotation_filtering.py ${outdir}/${vcf}.sim.multi_allele_split.annovar.input.db.csv
fi
