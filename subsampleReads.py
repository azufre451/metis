from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import argparse,sys,random

parser = argparse.ArgumentParser(description="Create multiple random samples of reads from FASTX files.")
parser.add_argument("pick", help="how many reads to pick? (multiple values will result in multiple subsamples)", type=int,nargs='+')
parser.add_argument("--labels", help="list of labels for the output files (length must be paired with the number of subsamples)",nargs='+')
parser.add_argument("--input", help="file containing the sequences")
parser.add_argument("--output_prefix", help="prefix for outputing the sequences")
parser.add_argument("--type", default="fasta", help="input file type (fasta|fastq)")
parser.add_argument("--verbose", action="store_true")
args=parser.parse_args()
 
if args.input: stream = args.input
else: stream = sys.stdin

seqList= []
i =0
for seq_record in SeqIO.parse(stream, args.type):
	seqList.append(seq_record) 
	i+=1
	if i % 10000 == 0: print str(i/10000)+"k seqs (",args.input,")"
 
print "sequences read.",len(seqList)

for reduction,label in zip(args.pick,args.labels):
	if reduction <= len(seqList):
		fname = (args.output_prefix if args.output_prefix else '')+'_'+label+'.'+args.type
		print 'outputing '+fname+' ( '+str(reduction)+' )'
		SeqIO.write(random.sample(seqList,reduction),fname,args.type)
		print '\tdone '+fname
print "complete!"
