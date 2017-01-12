import sys,os
import argparse
from subprocess import check_call
# from mungo.fasta import FastaReader
import re
import numpy as np
from collections import defaultdict
from itertools import groupby
from mungo.sequence import sixFrameTranslation

ACCEPTABLE_SYMBOLS = ['*', '-', 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
ACCEPTABLE_SYMBOLS_SET = set(ACCEPTABLE_SYMBOLS)

class writeable_dir(argparse.Action):
  def __call__(self,parser, namespace, values, option_string=None):
    prospective_dir=values
    if not os.path.isdir(prospective_dir):
        raise argparse.ArgumentTypeError("writeable_dir:{0} is not a valid path".format(prospective_dir))
    if os.access(prospective_dir, os.W_OK):
        setattr(namespace,self.dest,prospective_dir)
    else:
        raise argparse.ArgumentTypeError("writeable_dir:{0} is not a writeable dir".format(prospective_dir))

class readable_file(argparse.Action):
  def __call__(self,parser, namespace, values, option_string=None):
    for v in values:
      prospective_file=v
      if not os.path.isfile(prospective_file):
          raise argparse.ArgumentTypeError("readable_dir:{0} is not a valid file path".format(prospective_file))
      if os.access(prospective_file, os.R_OK):
          setattr(namespace,self.dest,prospective_file)
      else:
          raise argparse.ArgumentTypeError("readable_dir:{0} is not a readable file".format(prospective_file))

class MyError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)


def readFasta(fastafile):

  with open(fastafile, 'rU') as fp:
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            line=line[1:]
            if name: yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))

def checkSeq(test_str):
    return set(test_str) <= ACCEPTABLE_SYMBOLS_SET

def loadMSA(msa_file
  , dna_seq_file
  , occ_theshold
  , low_occ_threshold
  , length_threshold
  , verbose):

  #load dna seq file
  dna_seqs={}
  for h,s in readFasta(dna_seq_file):
    h=h.split()[0]
    dna_seqs[h] = s

  #Load MSA and check input
  if verbose:
    print "Loading and checking alignment file..."

  msa_dict = {}
  seq_count = 0
  for h,s in readFasta(msa_file):
    if "consensus" in h: continue #get rid of consensus sequecne

    h = h.split(";sample")[0]
    h = h.split(";size")[0]
    h = re.sub(r"{.*}", "", h)
    h = "_".join(h.split())
    s = s.upper()
    s = s.replace("X","-")

    assert checkSeq(s), "MSA has an unexpected symbol in: " + s

    if seq_count == 0:
      seq_len = len(s)
    else:
      assert seq_len==len(s), "Sequences are not all of the same length"


    #Convert flanking spaces to 0's to be removed later
    if seq_count<10 and verbose:
      print s

    s = list(s)
    for i, aa in enumerate(s):
      if aa!="-":
        break
      s[i] = '0'

    for i, aa in enumerate(reversed(s)):
      if aa!="-":
        break
      s[-(i+1)] = '0'

    s="".join(s)

    if seq_count<10 and verbose:
      print s

    seq_count += 1
    msa_dict[h] = s

  if seq_count!=len(msa_dict.keys()):
    print seq_count, len(msa_dict.keys())
    print msa_dict.keys()
    raise MyError("Duplicate sequence identifiers")

  #Convert to binary matrix and calculate occupancy
  seq_names = msa_dict.keys()
  occupancy_matrix = np.zeros((seq_count, seq_len), dtype=bool)
  non_flagging = np.zeros((seq_count, seq_len), dtype=bool)
  msa_matrix = np.chararray((seq_count, seq_len))
  for i, h in enumerate(seq_names):
    occupancy_matrix[i, :] = (np.array(list(msa_dict[h]))!="-") & (np.array(list(msa_dict[h]))!="0")
    non_flagging[i, :] = np.array(list(msa_dict[h]))!="0"
    msa_matrix[i, range(seq_len)] = np.array(list(msa_dict[h]))


  col_sum = np.sum(occupancy_matrix, axis=0)
  col_sum = np.array(col_sum, dtype=float)
  col_occ = col_sum/np.sum(non_flagging, axis=0)

  #Check if gaps have sufficient insert probability
  isPrevOccCount = 0
  designation = np.zeros(len(col_occ), dtype=int)

  if verbose:
    print col_sum
    print col_occ

  for i, occ in enumerate(col_occ):
    if occ>occ_theshold:
      isPrevOccCount+=1
      if isPrevOccCount >= length_threshold:
        for k in range(length_threshold):
          designation[i-k] = 1

    else:
      isPrevOccCount = 0

  print designation


  #Write out blocks to fasta files ready for clustering
  splits = np.split(
    np.arange(len(designation)), np.where(
      designation[:-1] != designation[1:])[0]+1)

  block_dict = {}
  index=0

  prefix = os.path.splitext(os.path.basename(msa_file))[0]
  print msa_matrix[1,:]

  i=0
  with open(prefix+ "_block_safile.SAF", 'w') as saffile:
    saffile.write("GeneID\tChr\tStart\tEnd\tStrand\n")
    for s in splits:
      # print s

      if len(s) < length_threshold: continue
      i += 1

      with open(prefix + "_block"+str(i)+".fasta", 'w') as outfile:
        for row in range(seq_count):

          seq = "".join(msa_matrix[row, s])
          seq = seq.replace("-","")
          seq = seq.replace("0","") #remove flanking 0s

          if len(seq) < length_threshold: continue

          dom = "".join(msa_matrix[row, ])
          dom = dom.replace("-","")
          dom = dom.replace("0","")

          dna_header, dna_start, dna_end = findDNALocation(seq_names[row]
            , seq
            , dom
            , dna_seqs
            , verbose)

          saffile.write("\t".join([seq_names[row] + "_block" + str(i),
            dna_header, str(dna_start+1), str(dna_end), "-"]) + "\n")
          # seq = seq.replace("0","-") #remove flanking 0s



          outfile.write(">" + seq_names[row] + "_block" + str(i) + "\n")
          outfile.write(seq + "\n")

  return


def findDNALocation(header, protein, domain_protein, dna_seqs, verbose):

  dna_header = header.split()[0].split("_frm")[0].split("_frame")[0].split("_[Plasmodium")[0]

  dna_seq = dna_seqs[dna_header]

  translated = sixFrameTranslation(dna_seq)

  count = 0
  num_stops = 0
  for f in translated:
    st = translated[f].find(domain_protein)
    if st != -1:
      count += 1
      start = st
      end = start + len(domain_protein)
      frame = f

  assert count <= 1, "multiple locations found!"

  if count==0: #deal with stop codons
    for f in translated:
      st = translated[f].replace("*","").find(domain_protein)

      if st != -1:
        st = st+translated[f][:st].count("*")
        count += 1
        start = st
        end = start + len(domain_protein)
        print "translated: ", translated[f][start:end]
        num_stops = translated[f][start:end].count("*")
        tranlated_protein = translated[f][start:end]
        frame = f

  if count == 0:
    print header
    print domain_protein
    print sixFrameTranslation(dna_seq)
    raise MyError("no location found!")

  if frame<0: #we need to start from the opposite end.
    dna_start = len(dna_seq) - (start*3 - frame-1) - (len(domain_protein)+num_stops)*3
  else:
    dna_start = start*3 + frame-1

  dna_end = dna_start + 3*(len(domain_protein)+num_stops) + 1
  dna_sub = dna_seq[dna_start:dna_end]

  trans_sub = sixFrameTranslation(dna_sub)

  if domain_protein not in [f.replace("*","") for f in trans_sub.values()]:
    print header
    print num_stops
    print "protein: ", domain_protein
    print dna_start, dna_end
    print frame
    print trans_sub.values()
    print sixFrameTranslation(dna_seq)
    raise MyError("did not find domain correctly")

  st = domain_protein.find(protein)

  num_stops_block=0
  pre_stops=0
  if num_stops>0:
    num_stops_block = tranlated_protein[st:st+len(protein)+2].count("*")
    pre_stops = tranlated_protein[:st].count("*")
    st=st+pre_stops

  if frame >= 0:
    block_dna_start = dna_start+3*st
    block_dna_end = block_dna_start + 3*(len(protein)+num_stops_block) + 1
  else:
    block_dna_end = dna_end - 3*st
    block_dna_start = block_dna_end - 3*(len(protein)+num_stops_block) - 1

  dna_sub = dna_seq[block_dna_start:block_dna_end]

  trans_sub = sixFrameTranslation(dna_sub)

  # if protein not in [f.replace("*","") for f in trans_sub.values()]: #check for stop codons
  #   pre_count = block_dna_start
  #   if frame >= 0:
  #     block_dna_start = dna_start+3*st
  #     block_dna_end = block_dna_start + 3*(len(protein)) + 1
  #   else:
  #     block_dna_end = dna_end - 3*st
  #     block_dna_start = block_dna_end - 3*(len(protein)) - 1
  #   dna_sub = dna_seq[block_dna_start:block_dna_end]
  #   trans_sub = sixFrameTranslation(dna_sub)

  # print protein
  # print  [f.replace("*","") for f in trans_sub.values()]
  if protein not in [f.replace("*","") for f in trans_sub.values()]:
    print "protein: ", protein
    print block_dna_start, block_dna_end
    print dna_start, dna_end
    print frame
    print trans_sub.values()
    print sixFrameTranslation(dna_seq)
    raise MyError("did not find small prot correctly")

  return dna_header, block_dna_start, block_dna_end





def main():
  parser = argparse.ArgumentParser(description='Generate binary fasta file for phylogentic analysis.')

  parser.add_argument('-o', '--outdir'
    , dest='outdir'
    , help="location of output file."
    , required=True)

  parser.add_argument('--msa', dest='msa'
    , help="location of protein MSA generated by gismo."
    , required=True)

  parser.add_argument('--dna', dest='dnafile'
    , help="location of original dna sequence fasta file."
    , required=True)

  parser.add_argument('--length_threshold', dest='insert', type=int
    , default=3
    , help="length threshold cutoff. (default=3)")

  parser.add_argument('--occupancyThreshold', dest='occupancy', type=float
    , default=0.95
    , help="column occupancy threshold cutoff. (default=0.95)")

  parser.add_argument('--lowOcupancyThreshold', dest='low_occupancy', type=float
    , default=0.01
    , help="column occupancy threshold cutoff. (default=0.01)")

  parser.add_argument('--cpu', dest='cpu', type=int, default=1
    , help="number of cpus to use. (default=1)")

  parser.add_argument('--verbose', dest='verbose', action='store_true'
    , default=False
    , help='print verbose output (default=False)')

  args = parser.parse_args()

  args.outdir = os.path.abspath(args.outdir)
  args.msa = os.path.abspath(args.msa)


  loadMSA(args.msa
    , args.dnafile
    , args.occupancy
    , args.low_occupancy
    , args.insert
    , args.verbose)


  return


if __name__ == '__main__':
  main()