import glob

mult=1
with open("combined_blocks_clustered.csv", 'w') as outfile:
  for f in glob.glob("*hierarchical_clusters.csv"):
    print f
    with open(f, 'rU') as infile:
      if mult==1:
        outfile.write(infile.next())
      else:
        infile.next()
      for line in infile:
        line = line.strip().split(",")
        outfile.write(",".join([line[0]] + [str(int(t)+10000*mult) for t in line[1:]])+"\n")
    mult+=1

mult=0
with open("combined_blocks.SAF", 'w') as outfile:
  for f in glob.glob("*safile.SAF"):
    mult+=1
    print f
    with open(f, 'rU') as infile:
      if mult==1:
        outfile.write(infile.next())
      else:
        infile.next()
      for line in infile:
        outfile.write(line)
