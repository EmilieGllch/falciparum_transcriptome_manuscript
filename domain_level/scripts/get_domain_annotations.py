import sys
from collections import defaultdict
import numpy as np
from mungo.fasta import FastaReader
from mungo.sequence import sixFrameTranslation


hmmfile = sys.argv[1]
prefix = sys.argv[2]
DNA_file = sys.argv[3]
PROT_file = sys.argv[4]

#Load DNA fasta
dna_seqs = {}
for h,s in FastaReader(DNA_file):
  h = h.split()[0]
  dna_seqs[h] = s

#Load Protein File
prot_seqs = {}
for h,s in FastaReader(PROT_file):
  if len(h.split()) > 2:
    h = h.split()[0] + "_" + h.split()[-1]
  else:
    h = "_".join(h.split())
  prot_seqs[h] = s

def overlap(start1, end1, start2, end2):
    """Does the range (start1, end1) overlap with (start2, end2)?"""
    return not (end1 < start2 or end2 < start1)

transcripts = defaultdict(list)
with open(hmmfile, 'r') as infile:
  for line in infile:
    if line[0]=="#":
      continue
    line=line.split()
    tran = line[0]
    if line[-1]!="-":
      tran = tran + "_" + line[-1]
    dom = line[3]
    E = float(line[12])
    from_ = int(line[17])
    to = int(line[18])

    transcripts[tran].append((dom, from_, to, E))

annotations = {}
for tran in transcripts:
  #sort by evalue
  keep = []
  transcripts[tran] = sorted(transcripts[tran], key=lambda tup: tup[3])
  for dom in transcripts[tran]:
    isOverlap=False
    for k in keep:
      if overlap(dom[1],dom[2],k[1],k[2]):
        isOverlap=True
        break
    if not isOverlap:
      keep.append(dom)
  annotations[tran] = keep

with open(prefix+"_hits.txt",'w') as outfile:
  for tran in annotations:
    annotations[tran] = sorted(annotations[tran], key=lambda tup: tup[1])
    #outfile.write(",".join([tran]+[tup[0].split("_")[-1] for tup in annotations[tran]])+"\n")
    outfile.write(tran+"\t"+",".join([str(tup) for tup in annotations[tran]])+"\n")



def findPosition(protein, dna):
  translated = sixFrameTranslation(dna)
  for frame in translated:
    st = translated[frame].find(prot)
    if st!=-1: break
  left=st
  right=left+len(protein)
  if frame<0: #we need to start from the opposite end.
    dna_start = len(dna) - (left*3 - frame-1) - len(protein)*3
  else:
    dna_start = left*3 + frame-1
  dna_end = dna_start + 3*len(protein) + 1

  translated2 = sixFrameTranslation(dna[dna_start:dna_end])
  found=False
  for frame in translated2:
    if translated2[frame]==protein:
      found=True
  if not found:
    print("Not found!")
    # print translated2
    print translated
    print translated2
    print left, right
    print protein
    sys.exit()

  return (dna_start+1, dna_end)


with open(prefix+"_domains.fasta", 'w') as outfas:
  with open(prefix+"_domains.SAF", 'w') as outSAF:
    outSAF.write("GeneID\tChr\tStart\tEnd\tStrand\n")
    for tran in annotations:
      dna_h = tran.split("_f")[0]
      dna = dna_seqs[dna_h]
      dom_num=0
      for tup in annotations[tran]:
        dom_num += 1
        domname = tran + "_" + tup[0] + "_domNum" + str(dom_num)
        outfas.write(">" + domname + "\n")
        prot = prot_seqs[tran][tup[1]-1:tup[2]]
        outfas.write(prot + "\n")
        dna_pos = findPosition(prot, dna)
        outSAF.write("\t".join([domname, dna_h, str(dna_pos[0]), str(dna_pos[1]), "-"])+"\n")



gaps = []
for tran in annotations:
  prev=("",0,0,0)
  for tup in annotations[tran]:
    g=tup[1]-prev[2]
    prev=tup
    if g>0:
      gaps.append(g)

print "Mean: ", np.mean(gaps)
print "Min: ", min(gaps)
print "Max: ", max(gaps)
print "Median: ", np.percentile(gaps,50)
