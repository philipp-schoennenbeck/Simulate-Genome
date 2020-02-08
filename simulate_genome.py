import csv
import random
import datetime
import argparse
import subprocess

parser = argparse.ArgumentParser()

parser.add_argument("closestMaskedCoordinates", help="A sorted tab seperated table with the start position in the second, the end position in the third and the closest next coordinates in the fourth column (in the first column should be the name/id)")
parser.add_argument("-v", "--verbose", help="Displays the progress of the program", action="store_true")
parser.add_argument("-a","--adjust", help="Adjust the repeat sequence sizes to the size of the simulated genome", action="store_false")
parser.add_argument("genomeSize", help="The size of the simulated genome in Mb")
parser.add_argument("-o", "--output", help="Path (and name) to the output file")
args = parser.parse_args()

res = subprocess.check_output("awk '{if($1 != atm) {sum += last_end;} atm = $1; last_end = $3 } END {sum += last_end; print sum}' " + args.closestMaskedCoordinates, shell=True)
res = int(res.decode('ascii'))
print(res)
gc_content = 0.5
repeatlist = []
genomesize_mb_sim = int(10**6*float(args.genomeSize))
genomesize_reference = 440159624
factor = 1
if (args.adjust):
  factor = genomesize_mb_sim/genomesize_reference
factor = min(factor, 1)

def createSequence(length):
    #print(length)
    sequence = []
    nucleotides = ["A", "T", "G", "C"]
    for i in range(length):
        if random.random() <= gc_content:
            sequence.append(nucleotides[random.randint(0,1)+2])
        else:
            sequence.append(nucleotides[random.randint(0,1)])
    return sequence

last = 0
if args.verbose : print("reading file...")
with open(args.closestMaskedCoordinates) as tsvfile:
  reader = csv.reader(tsvfile, delimiter='\t')
  for row in reader:
    last = row[2]
    repeatlist.append([int(abs((int(row[1])-int(row[2]))*factor)),int(abs(int(row[3])*factor)), createSequence(random.randint(1,10))])
if args.verbose : print("Done!")



random.shuffle(repeatlist)
if (args.adjust):
  repeatlist = repeatlist[0:int(factor*len(repeatlist))]


genome = []
percent = 0
if args.verbose : print("Simulating genome...")
while (len(genome) < genomesize_mb_sim):
  if (str(int(len(genome)/genomesize_mb_sim*100)) != percent):
    percent = str(int(len(genome)/genomesize_mb_sim*100))
    if args.verbose : print(percent, "%")
  dummy = repeatlist[random.randint(0,len(repeatlist)-1)]
  genome.extend(dummy[2]*int(dummy[0]/len(dummy[2])))
  genome.extend(createSequence(dummy[1]))

date = str(datetime.datetime.now()).replace(" ","_")
if args.output:
  outfile = open(args.output,"w")
else:
  outfile = open(str(genomesize_mb_sim) + "Mbp_"+ date,"w")

outfile.write(">SimulationGenome_"+str(genomesize_mb_sim)+"Mbp_" + str(date)+"\n")
counter = 0

if args.verbose : print("writing to file...")
for i in range(len(genome)):
    outfile.write(genome[i])
    counter +=1
    if (counter==60):
        outfile.write("\n")
        counter = 0
if args.verbose : print("Done!")