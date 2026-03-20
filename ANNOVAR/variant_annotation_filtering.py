
#// ...existing code...
import sys
import os
import pandas as pd

if len(sys.argv) < 1:
    raise SystemExit("Usage: python annovar_db_table.py <input_file>")

input_path = str(sys.argv[1])
input_dir = os.path.dirname(input_path)
input_base = os.path.basename(input_path)
db_file = input_path 
db_list=['Chr' , 'Pos' , 'End' , 'Ref' , 'Alt' , 'ID' , 'gwasCatalog' , 'dbSUPER' , 'distal_hQTLs' , 'dsQTL' , 'GeneHancer' , 'local_hQTLs' , 'otherDB_eQTLs' , 'SEA' , 'UCNEbase' , 'wgEncodeRegDnaseClustered' , 'EnhancerAtlas' , 'JEME' , 'encRegTfbsClustered' , 'phyloP100way' , 'phastConsElements100way' , 'repeatmasker' , 'genomicSuperDups' , 'simpleRepeat' , 'knownGene_all_func' , 'knownGene_exon_func' , 'refGene_all_func' , 'refGene_exon_func' , 'gnomad.all_0.001ab' , 'avsnp142' , 'ICGC+TCGA+COSMIC+ClinVar1' , 'ICGC+TCGA+COSMIC+ClinVar2' , 'cosmic70' , 'HGMD' , 'clinvar.benign' , 'clinvar.pathogenic' , 'GTEx_eQTL' , 'PancanQTL' , 'gnomad_all' , 'spliceai_score_wgs_exome_0.5', 'cadd']

# Read original file and database file
# var_df = pd.read_csv(input_path, sep='\t', header=0, dtype=str, na_values = 'NAN', low_memory=False)
db_df = pd.read_csv(db_file, sep='\t', header=0, dtype=str, na_values = 'NAN', low_memory=True)
print([db for db in db_list if db not in db_df.columns])

db_df['somatic_filter_key'] = 0
db_df['somatic_filter_key2'] =0
db_df['somatic_filter_key3'] = 0
db_df['somatic_status'] = 0

# filter variant subset 1 with keeping somatic/pathogenic annotation: ICGC+TCGA+COSMIC+ClinVar1 ICGC+TCGA+COSMIC+ClinVar2 cosmic70 HGMD clinvar.pathogenic 
db_df.loc[(db_df['HGMD'].notna()), 
                        'somatic_filter_key'] +=2
db_df.loc[(db_df['clinvar.pathogenic'].notna()), 
                        'somatic_filter_key'] +=2
db_df.loc[(db_df['clinvar.benign'].notna()), 
                        'somatic_filter_key'] -=2
db_df.loc[(db_df['cosmic70'].notna()), 
                        'somatic_filter_key'] +=2
# filter variants subset 2 with exonic function annotation: knownGene refGene 
db_df.loc[(db_df['knownGene_exon_func'].notna()) | 
                        (db_df['refGene_exon_func'].notna()) ,
                        'somatic_filter_key'] +=5

# filter variants subset 3 with function prediction annotation: cadd 
db_df['somatic_filter_key'] = db_df[['somatic_filter_key', 'cadd']].apply(
    lambda x: x.iloc[0] if pd.isna(x.iloc[1]) 
    else x.iloc[0] + 1 if float(x.iloc[1]) >= 15 
    else x.iloc[0] - 1 if float(x.iloc[1]) < 15 
    else x.iloc[0], axis=1)  
# db_df.loc[(db_df['cadd'].atype(int)>=15), 
#                         'somatic_filter_key'] +=1
# db_df.loc[(db_df['cadd'].atype(int)<15), 
#                         'somatic_filter_key2'] -=1  

# filter variants subset 4 with splice function prediction annotation: spliceai_score_wgs_exome_0.5         
db_df.loc[(db_df['spliceai_score_wgs_exome_0.5'].notna()), 
                        'somatic_filter_key'] +=1

# filter variant subset 1 with keeping somatic/pathogenic annotation: ICGC+TCGA+COSMIC+ClinVar1 ICGC+TCGA+COSMIC+ClinVar2
db_df.loc[(db_df['ICGC+TCGA+COSMIC+ClinVar1'].notna()) | 
                        (db_df['ICGC+TCGA+COSMIC+ClinVar2'].notna()), 
                        'somatic_filter_key2'] +=1         
# filter variants subset 4 with regulatory annotation: GTEx_eQTL PancanQTL distal_hQTLs dsQTL local_hQTLs otherDB_eQTLs gwasCatalog
db_df.loc[(db_df['GTEx_eQTL'].notna()) | 
                        (db_df['distal_hQTLs'].notna()) | 
                        (db_df['dsQTL'].notna()) | 
                        (db_df['local_hQTLs'].notna()) |                                                                         
                        (db_df['otherDB_eQTLs'].notna()) |
                        (db_df['PancanQTL'].notna()), 
                        'somatic_filter_key2'] +=1
db_df.loc[(db_df['gwasCatalog'].notna()), 
                        'somatic_filter_key2'] +=1
# filter variants subset 5 with regulatory annotation: dbSUPER GeneHancer SEA  UCNEbase wgEncodeRegDnaseClustered  EnhancerAtlas JEME encRegTfbsClustered
db_df.loc[(db_df['dbSUPER'].notna()) | 
                        (db_df['GeneHancer'].notna()) | 
                        (db_df['SEA'].notna()) |                        
                        (db_df['EnhancerAtlas'].notna()) | 
                        (db_df['JEME'].notna()),
                        'somatic_filter_key2'] +=1
db_df.loc[(db_df['UCNEbase'].notna()),
                        'somatic_filter_key2'] +=1
db_df.loc[(db_df['wgEncodeRegDnaseClustered'].notna()) |  
                        (db_df['encRegTfbsClustered'].notna()), 
                        'somatic_filter_key2'] +=1

# filter variants subset 6 with high germline frequency annotation: gnomad.all_0.001ab avsnp142
db_df.loc[(db_df['gnomad.all_0.001ab'].notna()) | 
                        (db_df['avsnp142'].notna()), 
                        'somatic_filter_key3'] -=10
# filter variants subset 7 with all germline frequency annotation: gnomad.all
db_df.loc[(db_df['gnomad_all'].notna()), 
                        'somatic_filter_key3'] -=1
# filter variants subset 6 with low complex sequence annotation: repeatmasker genomicSuperDups simpleRepeat
db_df.loc[(db_df['repeatmasker'].notna()) | 
                        (db_df['genomicSuperDups'].notna()) | 
                        (db_df['simpleRepeat'].notna()), 
                        'somatic_filter_key3'] -=5

db_df.loc[(db_df['somatic_filter_key'] == 0) & (db_df['somatic_filter_key2'] > 0) & (db_df['somatic_filter_key3'] >= -6), 'somatic_status']=1
db_df.loc[(db_df['somatic_filter_key'] == 0) & (db_df['somatic_filter_key2'] > 0) & (db_df['somatic_filter_key3'] >= -5), 'somatic_status']=2
db_df.loc[(db_df['somatic_filter_key'] == 0) & (db_df['somatic_filter_key2'] > 0) & (db_df['somatic_filter_key3'] >= -1), 'somatic_status']=2
db_df.loc[(db_df['somatic_filter_key'] == 0) & (db_df['somatic_filter_key2'] > 1) & (db_df['somatic_filter_key3'] >= -6), 'somatic_status']=2
db_df.loc[(db_df['somatic_filter_key'] == 0) & (db_df['somatic_filter_key2'] > 1) & (db_df['somatic_filter_key3'] >= -1), 'somatic_status']=3
db_df.loc[(db_df['somatic_filter_key'] == 0) & (db_df['somatic_filter_key2'] > 0) & (db_df['somatic_filter_key3'] ==0), 'somatic_status']=3
db_df.loc[(db_df['somatic_filter_key'] > 0) & (db_df['somatic_filter_key2'] > 1) & (db_df['somatic_filter_key3'] >= -6), 'somatic_status']=3
db_df.loc[(db_df['somatic_filter_key'] > 0) & (db_df['somatic_filter_key2'] > 0) & (db_df['somatic_filter_key3'] ==-1), 'somatic_status']=3
db_df.loc[(db_df['somatic_filter_key'] > 0) & (db_df['somatic_filter_key2'] > 0) & (db_df['somatic_filter_key3'] ==0), 'somatic_status']=4
db_df.loc[(db_df['somatic_filter_key'] > 0) & (db_df['somatic_filter_key3'] >= -1), 'somatic_status']=4
db_df.loc[(db_df['somatic_filter_key'] > 0) & (db_df['somatic_filter_key3'] ==0), 'somatic_status']=5
db_df.loc[(db_df['somatic_filter_key'] > 1) & (db_df['somatic_filter_key3'] >= -6), 'somatic_status']=5
db_df.loc[(db_df['somatic_filter_key'] >= 4), 'somatic_status']=6

db_df.loc[(db_df['somatic_filter_key'] > 0) & (db_df['somatic_filter_key3'] <= -11), 'somatic_status']=10
db_df.loc[(db_df['somatic_filter_key'] == 0) & (db_df['somatic_filter_key2'] > 1) & (db_df['somatic_filter_key3'] <= -11), 'somatic_status']=10

# conditions = [
#     ((db_df['somatic_filter_key']>0) & (db_df['somatic_filter_key3']<-6)),
#     ((db_df['somatic_filter_key2']>0) & (db_df['somatic_filter_key3']>=-5)),
#     ((db_df['somatic_filter_key']>1) & (db_df['somatic_filter_key3']>=-6)),
#     ((db_df['somatic_filter_key']>0) & (db_df['somatic_filter_key3']>=-1)),
#     ((db_df['somatic_filter_key']==0) & (db_df['somatic_filter_key2']>0) & (db_df['somatic_filter_key3']>=-1))
# ]
# values = [10, 10, 3, 2, 1]
# db_df['somatic_status'] = 0
# for condition, value in zip(conditions, values):
#     db_df.loc[condition & (db_df['somatic_status']==0), 'somatic_status'] = value
db_df.to_csv(db_file + '.variant_annotation_filtering.csv', sep='\t', index=False, na_rep='NAN', header=True)

