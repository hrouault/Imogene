import sys
import os
import shutil
import string
import glob
import datetime
import errno
import string

def mkdir(path):
   try:
      os.mkdir(path)
   except os.error, e:
      if e.errno != errno.EEXIST:
         raise

def spe2short(name):
   if name == "mus_musculus": return "MusMus"
   elif name == "pan_troglodytes":  return "PanTro"
   elif name == "pongo_abelii":  return "PonAbe"
   elif name == "macaca_mulatta":   return "MacMul"
   elif name == "homo_sapiens":  return "HomSap"
   elif name == "rattus_norvegicus":   return "RatNor"
   elif name == "equus_caballus":   return "EquCab"
   elif name == "canis_familiaris": return "CanFam"
   elif name == "bos_taurus": return "BosTau"
   elif name == "sus_scrofa": return "SusScr"
   elif name == "gorilla_gorilla":  return "GorGor"
   elif name == "callithrix_jacchus":  return "CalJac"
   else: return "NA"

def spe2num(name):
   if name == "mus_musculus": return 0
   elif name == "pan_troglodytes":  return 1
   elif name == "pongo_abelii":  return 2
   elif name == "macaca_mulatta":   return 3
   elif name == "homo_sapiens":  return 4
   elif name == "rattus_norvegicus":   return 5
   elif name == "equus_caballus":   return 6
   elif name == "canis_familiaris": return 7
   elif name == "bos_taurus": return 8
   elif name == "sus_scrofa": return 9
   elif name == "gorilla_gorilla":  return 10
   elif name == "callithrix_jacchus":  return 11
   else: return -1

def chr2str(chr):
   if chr == "1": return "chr1"
   elif chr == "2":  return "chr2"
   elif chr == "3":  return "chr3"
   elif chr == "4":  return "chr4"
   elif chr == "5":  return "chr5"
   elif chr == "6":  return "chr6"
   elif chr == "7":  return "chr7"
   elif chr == "8":  return "chr8"
   elif chr == "9":  return "chr9"
   elif chr == "10": return "chr10"
   elif chr == "11": return "chr11"
   elif chr == "12": return "chr12"
   elif chr == "13": return "chr13"
   elif chr == "14": return "chr14"
   elif chr == "15": return "chr15"
   elif chr == "16": return "chr16"
   elif chr == "17": return "chr17"
   elif chr == "18": return "chr18"
   elif chr == "19": return "chr19"
   elif chr == "X":  return "chrX"
   elif chr == "Y":  return "chrY"
   else: return "others"

def convert():
   for file in glob.glob('*.emf'):
      with open(file) as femf:
         print file
         line=femf.readline()
         while 1:
            while line[0:3] != 'SEQ':
               line=femf.readline()
            cols = [ -1 for i in range(12) ]
            spes = [ "" for i in range(12) ]
            chrs = [ "" for i in range(12) ]
            stas = [ -1 for i in range(12) ]
            stos = [ -1 for i in range(12) ]
            strs = [ 0 for i in range(12) ]
            seqs = [ "" for i in range(12) ]
            vseqs = [ [] for i in range(12) ]
            col=0
            while line[0:3] == 'SEQ':
               coords=line.rstrip().split(" ")
               if coords[1] != "ancestral_sequences":
                  spe=coords[1]
                  chr=coords[2]
                  sta=coords[3]
                  sto=coords[4]
                  str=coords[5]
                  ns=spe2num(spe)
                  cols[ns] = col
                  spes[ns] = spe2short(spe)
                  chrs[ns] = chr2str(chr)
                  stas[ns] = sta
                  stos[ns] = sto
                  strs[ns] = int(str)
               col+=1
               line=femf.readline()
            
            nref=spe2num("mus_musculus")

            if cols[nref]==-1 or chrs[nref]=="":
               continue
           
            strand=strs[nref]
            
            print chrs[nref],stas[nref],stos[nref],"length =",int(stos[nref])-int(stas[nref])+1
            
            femf.readline()
            line=femf.readline()
            
#            counter=0
            while line[0:2] != '//':
               for i in range(12):
                  if cols[i] != -1:
                     vseqs[i].append(line[cols[i]])
#               counter += 1
#               if not counter % 10000:
#                  print counter
               line=femf.readline()
            
            for i in range(12):
               if strand == -1:
                  vseqs[i].reverse()
               seqs[i]=''.join(vseqs[i])
               if strand == -1:
                  seqs[i]=seqs[i].translate(string.maketrans('AaTtCcGg','TtAaGgCc'))

            folder=chrs[nref]
            mkdir(folder)
            filename=folder+"/"+stas[nref]+"-"+stos[nref]+".fa"
            fasta=open(filename,'w')
            
            for i in range(12):
               if spes[i] != "":
                  if spes[i] == 'MusMus':
                     fasta.write('>'+spes[i]+' '+chrs[i]+' '+stas[i]+' '+stos[i]+'\n')
                  else:
                     fasta.write('>'+spes[i]+'\n')
                  
                  fasta.write(seqs[i]+'\n')
            
            line=femf.readline()
            if not line: break

