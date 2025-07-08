#!/usr/bin/env python3


#%% MODULES TO IMPORT 

import os
import sys
import subprocess as sbp
import gzip
import pandas as pd
import argparse

def get_genotype(genotype, ref, alt):
  out = []
  for g in genotype:
    tmp = g.split(':')[0]
    if tmp == "0/0":
      out.append(''.join((ref)))
    if tmp == "0|0":
      out.append(''.join((ref)))
    elif tmp == "0/1":
      #out.append(''.join((alt[0]))) #overestimate difference
      out.append(''.join((ref,alt[0]))) #might cause different length
    elif tmp == "0|1":
      #out.append(''.join((alt[0])))
      out.append(''.join((ref,alt[0]))) #might cause different length
    elif tmp == "1/1":
      out.append(''.join((alt[0])))
    elif tmp == "1|1":
      out.append(''.join((alt[0])))
    elif tmp == "2/2" or tmp == "2|2": #multi-allelic, won't be activated unless there's an error
      out.append(''.join((alt[1])))
      print('encounter multi, error',alt[1])
    elif tmp == "0/2" or tmp == "0|2": #multi-allelic, won't be activated
      out.append(''.join((alt[1])))
      print('encounter multi, error',alt[1])
    elif tmp == ("./."):
      out.append(''.join(('N')))
    elif tmp == (".|."):
      out.append(''.join(('N')))
  #if len(set([len(i) for i in out])) != 1: #different length
  #  add=max([len(i) for i in out])-1
  #  out=[i if len(i) == add+1 else i+'-'*add for i in out] #add to the same length as we don't do alignment after
    #print(out)
  #return ';'.join(out)
  return out


def main(vcffile):
  N=0
  NBI=0
  NMU=0
  GEN=[]
  SNP={}

  with gzip.open(vcffile, 'rb') as f:
    for line in f:
      line=bytes.decode(line)
      if line.startswith('#CHROM'):
        strains = line.split('\t')[9:]
        print('there are ',len(strains),' strains')
      elif not line.startswith('##'):
        N = N + 1
        REF = line.split('\t')[3]
        ALT = line.split('\t')[4]
        POS = int(line.split('\t')[1]) #position of the base
        if len(ALT)==1 and ALT != "*": #bi-allelic SNPs only
          NBI = NBI + 1
          l=get_genotype(line.split('\t')[9:], REF, ALT) #all strain at that position
          GEN.append(l)
          SNP[POS]=l
          #print(len(GEN),len(l))
        if len(ALT)>1:
          NMU = NMU + 1
        if N % 1000000 == 0:
          sys.stderr.write('==> Processed %i,000,000 sites\n' % (N/1000000))
        elif N % 100000 == 0:
          sys.stderr.write('==> Processed %i00,000 sites\n' % (N/100000))
        elif N % 10000 == 0:
          sys.stderr.write('==> Processed %i0,000 sites\n' % (N/10000))

  sys.stderr.write('\n==> Found %i bi-allelic SNPs.\n' % (NBI))
  sys.stderr.write('==> Found %i multi-allelic SNPs. These are ignored.\n' % (NMU))
  sys.stderr.write('==> Writing sequences to FASTA file %s\n\n' % os.path.abspath(str(vcffile) + '.fa'))
  print([len(SNP[k]) for k in SNP])
  SNP.pop(list(SNP.keys())[-1])
  df=pd.DataFrame.from_dict(SNP)
  df.index=strains
  print(df)
  def no_filter_out():
    with open(str(vcffile) + '.fa', 'a') as f2:
      for i in range(0,len(strains)):
        #print(strains[i],' there are ',len(GEN[0].split(';')),' bases')
        #try: 
          #print('>'+strains[i]+'\n'+''.join([item.split(';')[i] for item in GEN])+'\n')
       # except:
        #  print("An exception occurred")
        #f2.write('>'+strains[i]+'\n'+''.join([item.split(';')[i] for item in GEN[:-1]])+'\n') #not including the last position 
        seq=''.join([item[i] for item in GEN])
        #print(len(seq))
        f2.write('>'+strains[i]+'\n'+seq+'\n') 
    f2.close()  
  #make dataframe and filter the positions by coverage depth:  - for deletions indicated by low coverage
  def filter_coverage(SNP,bams):
    df=pd.DataFrame.from_dict(SNP)
    df.index=strains
    print(df)
    snps=df.to_dict('index')
    tbeqev=sbp.check_output(f'ls {bams}',shell=True).decode().strip('\n').split()
    for f in tbeqev[:]:
      if f.split('/')[-1] in snps.keys():
        print(f"processing {f}")
        coverage=sbp.check_output(f"samtools depth {f}|cut -f2,3",shell=True).decode().strip('\n').split('\n')
        coverage={int(i.split('\t')[0]):int(i.split('\t')[1]) for i in coverage}
        #print(coverage)
        for p in snps[f.split('/')[-1]]:
          if coverage.get(p,0)<15:
            snps[f.split('/')[-1]][p]='-'*len(snps[f.split('/')[-1]][p])
        seq=''.join(list(snps[f.split('/')[-1]].values()))
        print(seq)
            #print(p,coverage.get(p,0),snps[f.split('/')[-1]][p])
    #output file
    with open(str(vcffile).replace('.vcf.gz','') + '.dp.fa', 'a') as f2:
      for k in snps:
        seq=''.join(list(snps[k].values()))
        f2.write('>'+k+'\n'+seq+'\n') 
    f2.close() 
        
    #print(tbeqev,snps.keys())
  #no_filter_out()
  filter_coverage(SNP,'/home/zed/disk1/Tb_Fre_WGS/SNP_calling/20240122_with_Tbeveq/TbeqTbev_bam/*bam')
      
if __name__ == '__main__':
  parser = argparse.ArgumentParser(
	  description='Convert VCF to FASTA file.', \
	  usage = 'vcf2fasta.py <vcf>')
  parser.add_argument('vcf', help='VCF file', metavar='vcf')
  options = parser.parse_args()
  
  main(options.vcf)

