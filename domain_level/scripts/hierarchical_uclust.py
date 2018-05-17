from mungo.fasta import FastaReader
from mungo.sequence import sixFrameTranslation
import sys, os
from optparse import OptionParser
from subprocess import check_call
from collections import defaultdict

USEARCH = "/home/users/allstaff/tonkin-hill.g/var_blocks/third-party/usearch8"

def preprocess(inputfile, outdir):
  #now extract prefix from read name
  prefix = os.path.splitext(os.path.basename(inputfile))[0]

  outputfile = outdir + prefix + "_PreProcessed.fa"

  with open(outputfile, 'w') as outfile:
    for h,s in FastaReader(inputfile):
      h="_".join(h.split())
      if len(h)>250:
        h = h[0:100]+h[-150:]
      s=s.replace("*","X")
      outfile.write(">" + h + "\n")
      outfile.write(s + "\n")

  return outputfile

def run_uclust(inputfile, outdir, identities
  , minLength, verbose , cpu):
  #now extract prefix from read name
  prefix = os.path.splitext(os.path.basename(inputfile))[0]

  #sort identities from largest to smallest
  identities = sorted(identities, reverse=True)

  #iteratively run cd-hit
  prev = inputfile
  for identity in identities:

    out = outdir + prefix + "_c" + str(identity)

    cmd_usearch = (USEARCH
      + " -cluster_fast " + prev
      + " -sort length"
      + " -wordlength 2"
      #+ " -fulldp"
      #+ " -top_hits_only"
      + " -maxaccepts 1 -maxrejects 0"
      + " -alpha aa"
      + " -threads " + str(cpu)
      + " -mincols " + str(minLength)
      + " -id " + str(identity)
      + " -centroids " + out
      + " -uc " + out + ".uc")

    prev = out

    if verbose:
      print cmd_usearch

    check_call(cmd_usearch, shell=True)

  #collect results
  transcripts_clusters = defaultdict(list)
  start = True
  for identity in sorted(identities):
    clstr = outdir + prefix + "_c" + str(identity) + ".uc"
    with open(clstr, 'r') as infile:
      current_level = {}
      for line in infile:
        line = line.split()
        seq = line[8]
        cluster = line[1]
        if seq in current_level and current_level[seq]!=cluster:
          print "Problem: ", identity, seq, current_level[seq], cluster
          sys.exit(1)
        current_level[seq] = cluster
      previous_levels = {}
      if start:
        for t in current_level:
          transcripts_clusters[t] = [current_level[t]]
        start = False
      else:
        for t in transcripts_clusters:
          previous_levels[current_level[t]] = transcripts_clusters[t] + [current_level[t]]
        for t in current_level:
          transcripts_clusters[t] = previous_levels[current_level[t]]

  #write output to file
  with open(outdir + prefix + "_hierarchical_clusters.csv", 'w') as outfile:
    outfile.write(",".join(["Transcript"]
      +[str(i) for i in sorted(identities)])+ "\n")
    for t in transcripts_clusters:
      outfile.write(",".join([t]+[c for c in transcripts_clusters[t]])+ "\n")



def generate_SAF(inputfile, dnafile, outdir, verbose):
  if verbose:
    print "Generating SAF file..."

  prefix = os.path.splitext(os.path.basename(inputfile))[0]

  safFile = outdir + prefix + ".SAF"

  #load DNA sequences
  dna_seqs = {}
  for h,s in FastaReader(dnafile):
    dna_seqs[h.split()[0]] = s

  #load sequences in to dictionary
  with open(safFile , 'w') as outfile:
    outfile.write("GeneID\tChr\tStart\tEnd\tStrand\n")
    for h,s in FastaReader(inputfile):
      if "frm" in h:
        name = h.split("_frm")[0]
        oframe = int(h.split("_frm")[1].split("_")[0])
      else:
        name = h.split(" ")[0]
        oframe = int(h.split(" frame_")[1].split("_")[0])
      dna = dna_seqs[name]

      translated = sixFrameTranslation(dna)
      count = 0
      for f in [oframe]:#translated:
        st = translated[f].find(s)
        if st != -1:
          count += 1
          start = st
          end = start + len(s)
          frame = f

      if count > 1:    
        print "multiple locations found!"
        sys.exit(1)
      if count == 0:
        print "no location found!"
        sys.exit(1)

      if frame<0: #we need to start from the opposite end.
       dna_start = len(dna) - (start*3 - frame-1) - len(s)*3-1
       dna_start = max(dna_start, 0)
      else:
        dna_start = start*3 + frame-1
      dna_end = dna_start + 3*len(s) + 1

      dna_sub = dna[dna_start:dna_end]


      trans_sub = sixFrameTranslation(dna_sub)
      #print s, trans_sub.values()

      if s not in trans_sub.values():
        print "frame, ", frame
        print "did not find correctly"
        sys.exit(1)

      seg_name = "_".join(h.split())
      if len(seg_name) > 250:
        seg_name = seg_name[0:100]+seg_name[-150:]
      outfile.write("\t".join([seg_name
        ,name,str(dna_start),str(dna_end),"-"])+"\n")

  return safFile



def main():
    parser = OptionParser()

    parser.add_option("-o", "--outputdir", dest="outputdir",
        help=("the output directory"))

    parser.add_option("-i", "--inputSequences", dest="inputfile",
        help="A fasta file with sequences to be clustered")

    parser.add_option("", "--minLength", type=int, dest="minLength"
      , default=50
      , help="minimum number of columns in alignment (default=50)")

    parser.add_option("-c", "--identityThresholds", type=str, dest="identities",
        default="0.97,0.95,0.9,0.85,0.8,0.75,0.7,0.65,0.6,0.5",
        help=("A list of identities to perform each level of the hierarchical "
          + "clustering sperated by commas. "
          + "(default=0.97,0.95,0.9,0.85,0.8,0.75,0.7,0.65,0.6,0.5)"))

    parser.add_option("","--verbose", action="store_true", dest="verbose"
        , default=False, help="turns on more detailed output")

    parser.add_option("","--SAF", dest="saf", default=None
        , help="generates SAF file for featurecounts. Requires original dna seq")

    parser.add_option("", "--cpu", type=int, dest="num_cpu", default=1,
        help="number of cores to use")


    (options, args) = parser.parse_args()

    options.identities = options.identities.strip(",").split(",")

    preProcessed = preprocess(options.inputfile, options.outputdir)

    run_uclust(preProcessed, options.outputdir, options.identities
      , options.minLength
      , options.verbose , options.num_cpu)

    if options.saf != None:
      generate_SAF(options.inputfile, options.saf, options.outputdir
        , options.verbose)

if __name__ == '__main__':
    main()

