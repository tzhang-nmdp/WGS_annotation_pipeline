# WGS_annotation_pipeline
## Requirement
ANNOVAR package:https://www.openbioinformatics.org/annovar/annovar_download_form.php

SNPEFF package: https://pcingola.github.io/SnpEff/download

CADD package: https://cadd.bihealth.org/download

Master script: run_variant_annotation.sh

#### ANNOVAR
sh ${annovar_dir}/annovar.sh ${vcf} ${annovar_dir} ${outdir} ${ref_fasta_dir} all all

#### SNPEFF
sh ${snpeff_dir}/snpeff.sh ${vcf} ${snpeff_dir} ${outdir}

#### CADD
sh ${cadd_dir}/cadd_db.sh ${vcf} ${cadd_dir} ${outdir}
