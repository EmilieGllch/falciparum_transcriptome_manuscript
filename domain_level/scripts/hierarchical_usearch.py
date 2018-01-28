from mungo.fasta import FastaReader
from mungo.sequence import sixFrameTranslation
import sys, os
import shutil
from optparse import OptionParser
from subprocess import check_call
from collections import defaultdict
import subprocess


USEARCH = "/home/users/allstaff/tonkin-hill.g/var_blocks/third-party/usearch8"


def run_usearch(inputfile, outdir, minI, cpu, verbose):

  prefix = os.path.splitext(os.path.basename(inputfile))[0]

  #copy file to output directory and rename transcripts with temprary
  # place folders
  names = {}
  segCount = 0
  with open(outdir + prefix + "_renamed.fa", 'w') as outfile:
    for h,s in FastaReader(inputfile):
      names["block_" + str(segCount)] = "_".join(h.split())
      outfile.write(">" + "block_" + str(segCount) + "\n")
      s = s.replace("*", "X")
      outfile.write(s + "\n")
      segCount += 1

  inputfile = outdir + prefix + "_renamed.fa"

  #read in sequences
  sequences = []
  for h,s in FastaReader(inputfile):
    sequences.append((h,s))

  num_sequences = len(sequences)
  sub_size = num_sequences/cpu + 1

  #write out into seperate files
  files = []
  for i in range(0,cpu):
    filename = outdir + prefix + "_query" + str(i) + ".fa"
    with open(filename, 'w') as outfile:
      for seq in sequences[i*sub_size:(i+1)*sub_size]:
        outfile.write(">"+seq[0]+"\n")
        outfile.write(seq[1]+"\n")
    files.append(filename)

  #now run usearch with each query file
  processes = set()
  #-usearch_global ./test.fa -db ./param_p1e-5_k5_radiusThresh1.fa -id 0.2 -maxaccepts 0 -maxrejects 0
  #-fastapairs alignFasta.out -uc uc.out  -qmask none -dbmask none -wordlength 2 -alpha aa -fulldp -threads 20
  for f in files:
    out = f.strip(".fa") + "_usearch.tab"
    search_cmd = (USEARCH
      + " -usearch_global " + f
      + " -db " + inputfile
      + " -id " + str(minI)
      + " -maxaccepts 0 -maxrejects 0"
      + " -qmask none -dbmask none"
      + " -wordlength 2"
      + " -mincols 5"
      + " -fulldp"
      + " -alpha aa"
      + " -threads 5"
      + " -blast6out " + out)

    if verbose:
      print search_cmd
    #append system call
    processes.add(subprocess.Popen(search_cmd,  shell=True))
  #Check if all the child processes were closed
  for p in processes:
    p.wait()

  #now combine output into a results file and rename to original names
  with open(outdir + prefix + "combinedSearch.tab", 'w') as outfile:
    known = {}
    for f in files:
      f = f.strip(".fa") + "_usearch.tab"
      with open(f, 'r') as infile:
        for line in infile:
          line=line.split()
          match = tuple(sorted([line[0],line[1]]))
          if match in known:
            known[match] = min(1-float(line[2])/100.0, known[match])
          else:
            known[match]=1-float(line[2])/100.0

    for edge in known:
      nodeA = names[edge[0]]
      nodeB = names[edge[1]]
      if len(nodeA)>250:
        nodeA = nodeA[0:100]+nodeA[-150:]
      if len(nodeB)>250:
        nodeB = nodeB[0:100]+nodeB[-150:]
      outfile.write("\t".join(
        [nodeA,nodeB, str(known[edge])])+"\n")

  return

def main():
    parser = OptionParser()

    parser.add_option("-o", "--outputdir", dest="outputdir",
        help=("the output directory"))

    parser.add_option("-i", "--inputSequences", dest="inputfile",
        help="A fasta file with sequences to be clustered")

    parser.add_option("", "--minI", type=float, dest="minI", default=0.5,
        help="minimum identity threshold (default=0.5)")

    parser.add_option("","--verbose", action="store_true", dest="verbose"
        , default=False, help="turns on more detailed output")

    parser.add_option("", "--cpu", type=int, dest="num_cpu", default=1,
        help="number of cores to use")


    (options, args) = parser.parse_args()

    run_usearch(options.inputfile, options.outputdir
      , options.minI
      , options.num_cpu, options.verbose)



if __name__ == '__main__':
    main()


