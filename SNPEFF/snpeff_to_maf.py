

############################################################

## example usage:
##
## python VCF_MAF.py -v --vcf examples/example1.vcf --maf out.maf

import os,sys,argparse,re
## output stream for logging
#_log = sys.stderr

#def main():
#  parser = argparse.ArgumentParser(description='Convert VCF files to MAF format.')
#  parser.add_argument('--vcf', dest='vcf_filenames', nargs='+', help='Filename of VCF file(s) to read')
#  parser.add_argument('--vcf', dest='vcf_filenames', help='Filename of VCF file(s) to read')
#  parser.add_argument('--maf', dest='maf_filename', help='Filename of MAF file to output') 
#  args = parser.parse_args()
#  if args.vcf_filenames:
#      in_vcf = ( open(filename, 'r') for filename in args.vcf_filenames )
#  if args.maf_filename:
#      out_maf = ( open(filename, 'w') for filename in args.maf_filename )    
#def remove_quotes(string):
#  return re.sub(r'^\s*([\'|"])(.*)\1\s*$', r'\2', string)

## snpEFF Functional class {NONE, SILENT, MISSENSE, NONSENSE}.
## snpEFF Effects:
## TODO: should be checked

#try:
#    in_vcf=open(args.vcf_filenames,'r')
#    out_maf=open(args.maf_filename,'w')
#except:
#    vcf_files = [ sys.stdin ] 
## call main method if this file is run as a script
#if __name__ == '__main__':
#    main()
in_vcf=open(sys.argv[1],'r')
out_maf=open(sys.argv[2],'w')
error_out=open(sys.argv[3],'w')
No=int(sys.argv[4])
snpeff_variant_classification = {
  "TF_binding_site_variant" :  "TF_binding_site_variant",
  "missense_variant" :  "Missense_Mutation",
  "upstream_gene_variant" :  "5'Flank",
  "intron_variant" :  "Intron",
  #"sequence_feature" :  "Missense_Mutation",
  "TFBS_ablation" :  "TFBS_ablation",
  "5_prime_UTR_premature_start_codon_gain_variant" :  "Translation_Start_Site",
  "downstream_gene_variant" :  "3'Flank",
  "splice_region_variant&intron_variant" :  "Splice_Site",
  "synonymous_variant" :  "Silent",
  "3_prime_UTR_variant" :  "3'UTR",
  "5_prime_UTR_variant" :  "3'UTR",
  "structural_interaction_variant" : "Frame_Shift_Del",
  "splice_region_variant&synonymous_variant" :  "Splice_Site",
  "splice_donor_variant&intron_variant" :  "Splice_Site",
  "intragenic_variant" :  "IGR",
  "frameshift_variant" : "Frame_Shift_Del",
  "splice_acceptor_variant&intron_variant" :  "Splice_Site",
  "disruptive_inframe_deletion" :  "In_Frame_Del",
  "splice_acceptor_variant&splice_region_variant&intron_variant" :  "Splice_Site",
  "missense_variant&splice_region_variant" :  "Splice_Site",
  "splice_region_variant&non_coding_transcript_exon_variant" :  "Splice_Site",
  "disruptive_inframe_insertion" :  "In_Frame_Ins",
  "stop_gained" :  "Nonsense_Mutation",
  "start_lost" :  "Translation_Start_Site",
  "conservative_inframe_insertion" :  "In_Frame_Ins",
  "frameshift_variant&splice_acceptor_variant&splice_region_variant&intron_variant" : "Frame_Shift_Del",
  "conservative_inframe_deletion" :  "In_Frame_Del",
  "frameshift_variant&stop_gained" : "Frame_Shift_Del",
  "frameshift_variant&splice_donor_variant&splice_region_variant&intron_variant" : "Frame_Shift_Del",
  "splice_region_variant" :  "Splice_Site",
  "protein_protein_contact" :  "Silent",
  "splice_donor_variant&splice_region_variant&intron_variant" :  "Splice_Site",
  "stop_lost&splice_region_variant" : "Nonstop_Mutation",
  "INTERGENIC" :  "IGR",
  "INTERGENIC_CONSERVED" :  "IGR",
  "INTRAGENIC" :  "Silent",
  "UTR_5_PRIME" :  "5'UTR",
  "UTR_5_DELETED" :  "5'UTR",
  "UTR_3_PRIME" :  "3'UTR",
  "UTR_3_DELETED" :  "3'UTR",
  "START_GAINED" :  "Translation_Start_Site",
  "START_LOST" :  "Translation_Start_Site",
  "SYNONYMOUS_START" :  "Translation_Start_Site",
  "NON_SYNONYMOUS_START" :  "Translation_Start_Site",
  "STOP_GAINED" :  "Nonsense_Mutation",
  "SYNONYMOUS_STOP" :  "Silent",
  "STOP_LOST" :  "Nonstop_Mutation",
  "SPLICE_SITE_ACCEPTOR" :  "Splice_Site",
  "SPLICE_SITE_DONOR" :  "Splice_Site",
  "NON_SYNONYMOUS_CODING" :  None,
  "SYNONYMOUS_CODING" :  "Silent",
  "FRAME_SHIFT" :  None,
  
  "CODON_CHANGE" :  "Missense_Mutation",
  "CODON_INSERTION" :  "In_Frame_Ins",
  "CODON_CHANGE_PLUS_CODON_INSERTION" :  "In_Frame_Ins",
  "CODON_DELETION" :  "In_Frame_Del",
  "CODON_CHANGE_PLUS_CODON_DELETION" :  "In_Frame_Del",
  "TRANSCRIPT" :  None,
  "CDS" :  None,
  "GENE" :  None,
  "EXON" :  None,
  "EXON_DELETED" :  "In_Frame_Del",
  "INTRON_CONSERVED" :  "Intron",
  "RARE_AMINO_ACID" :  "Missense_Mutation"
  }
sp_id=[]
for line in in_vcf:
    if line.startswith("##"):
        continue        
    item=line.strip().split('\t')
    if line.startswith("#CHROM") or line.startswith("CHROM"):
        for i in range(0,No):
            sp_id.append(item[i])
        out_maf.write('Chromosome'+'\t'+'Start_Position'+'\t'+'End_Position'+'\t'+'Reference_Allele'+'\t'+'Tumor_Seq_Allele2'+'\t'+'Variant_Classification'+'\t'+'Variant_Type'+'\t'+'Hugo_Symbol'+'\t'+'Tumor_Sample_Barcode'+'\t'+'var_id'+'\t'+ 'Protein_Change'+'\n')  
        continue          
    else:
        chr_No=item[0]
        posst=item[1]
        ID_list=item[0]+':'+item[1]+':'+item[3]+':'+item[4]
        if ID_list[3]=='DEL': 
            posend=ID_list[2]
            ref='TCGATCGA'
            alt='TCGA'
            vt=ID_list[3]
        elif ID_list[3]=='INV':        
            posend=ID_list[2]
            ref='TCGA'
            alt='TCGA'    
            vt=ID_list[3]
        elif ID_list[3]=='INS':
            posend=ID_list[2]
            ref='TCGA'
            alt='TCGATCGA'  
            vt=ID_list[3]      
        elif ID_list[3]=='DUP':
            posend=ID_list[2]
            ref='TCGA'
            alt='TCGATCGATCGA'
            vt=ID_list[3]      
        else:
            ref=item[3]
            alt=item[4]
            posend=int(item[1])+len(ref)-1 
            if len(ref)==1 and len(alt)==1:
                vt='SNP'
            elif len(ref)==1 and len(alt)==2:
                vt='DNP'
            elif len(ref)==3 and len(alt)==3:  
                vt='TNP'     
            elif len(ref)>len(alt):  
                vt='DEL'    
            elif len(ref)<len(alt):  
                vt='INS'                           
        info=item[7].split('|')
    #    print(len(info))
        # remainder = (len(info) % 15)
        # for va in ['MODIFIER', 'HIGH', 'MODERATE', 'LOW']:
        #     try:
        #         info_idx = info.index(va)
        #         break  # Stop after first match is found
        #     except ValueError:
        #         continue  # Move to next value if not found
        # m=int((len(info)-remainder)/15+1)
        # if remainder!=(info_idx+1):
        #     remainder=info_idx
        gene=''
        vc=''
        m=len(info)
        for i in range(1,m):
            if info[i] in ['MODIFIER', 'HIGH', 'MODERATE', 'LOW']:
                # g=int(i*15-15+remainder+1)
                # v=int(i*15-15+remainder-1)  
                # p=int(i*15-15+remainder+8)    
                gene=info[i+1]
                vc=info[i-1] 
                pc=info[i+8]
            # print(pc)
                if not pc.startswith('p.'):
                    pc=info[i+7]
            # print(pc)
                #print(chr_No+'\t'+str(posst)+'\t'+str(posend)+'\t'+str(ref)+'\t'+str(alt)+'\t'+str(vc)+'\t'+str(vt)+'\t'+str(gene)+'\t'+str(sp_id[i])+'\t'+str(item[2])+'\t'+str(pc)+'\n')
                # if gene!='':
                #     break  
                if gene=='' or vc=='':
                    error_out.write(line) 
                    continue                                  
                if vc in snpeff_variant_classification:
                    vc=snpeff_variant_classification[vc]
                for j in range(9,No):
                    geno=item[j].split(':')[0]
                    #print(item[i])
                    if not geno.startswith('0/0') and not geno.startswith('./.'):
                        out_maf.write(chr_No+'\t'+str(posst)+'\t'+str(posend)+'\t'+str(ref)+'\t'+str(alt)+'\t'+str(vc)+'\t'+str(vt)+'\t'+str(gene)+'\t'+str(sp_id[j])+'\t'+str(item[2])+'\t'+str(pc)+'\n')
            
