import logging
import fnmatch
import os
import numpy as np
import glob
import argparse
import pysam
import pandas
import concurrent.futures

def get_seq_len_from_bam(samfile):
	temp = []
	for i, dic in enumerate(samfile.header['SQ']):
		try:
			temp.append({'Seq':dic['SN'], 'Length':dic['LN']})
		except Exception:
			print "i: {0} d: {1}".format(i,dic)
			raise
	return temp


def coverage_vectors(contigs_size):
	coverage = {}
	for i in contigs_size:
		temp = {}
		temp["positions"] = np.zeros(i["Length"])
		temp["nb_reads"] = 0
		temp["nb_bp"] = 0
		coverage[i["Seq"]] = temp
	return coverage

def parse(samfile,coverage):
	for l in samfile.fetch():
		try:
			coverage[samfile.getrname(l.tid)]["nb_reads"]+=1
			if args.pv=="0.7.5":
				# pysam 0.7.5
				begin = l.pos+1
				end = l.pos+l.alen
				coverage[samfile.getrname(l.tid)]["positions"][begin:end] = 1
				coverage[samfile.getrname(l.tid)]["nb_bp"]+=l.alen
			else:
				# pysam 0.8.4
				begin = l.reference_start+1
				end = l.reference_start+l.reference_length
				coverage[samfile.getrname(l.reference_id)]["positions"][begin:end] = 1
				coverage[samfile.getrname(l.tid)]["nb_bp"]+=l.reference_length
				
		except Exception:
			print "line: {0}".format(l)
			raise
	return coverage

def process(inputfile,outputfile,pv):
	if pv=="0.7.5":
		samfile = pysam.Samfile(inputfile, "rb") # for older pysam version on OSC
	else:
		samfile = pysam.AlignmentFile(inputfile, "rb")
	print("Reading the length of sequences from bam file header")
	contigs_size=get_seq_len_from_bam(samfile)
	print("Getting the coverage tables ready")
	coverage = coverage_vectors(contigs_size)
	print("Getting coverage values")
	coverage = parse(samfile,coverage)
	print("Calculating sequence coverage ratios and writing the output")
	coverage_prop = {}
	for contig,vector in coverage.items():
		temp = {}
		for i in contigs_size:
			if contig == i["Seq"]:
				temp["length"] = i["Length"]
		temp["length_covered"] = np.sum(vector["positions"])/float(len(vector["positions"]))*100
		temp["number_reads"] = vector["nb_reads"]
		temp["number_bp"] = vector["nb_bp"]
		temp["coverage"] = vector["nb_bp"]/temp["length"]
		if vector["nb_reads"]>1 :
			coverage_prop[contig] = temp
	print("Now writing the output")
	output=pandas.DataFrame(coverage_prop).transpose()
	output = output[['length', 'length_covered', 'coverage', 'number_bp', 'number_reads']]
	output=output.sort(['number_bp','length_covered'],ascending=[0,0])
	output.to_csv(outputfile)
	samfile.close()

def main():
	parser = argparse.ArgumentParser(description='Parse a directory of bam files all mapped to the same reference database.')
	parser.add_argument('--bam_dir','-b', dest='bam_dir', required=True, 
				   help='directory of input bam files (sorted and indexed)')
	parser.add_argument('--out','-o', dest='outdir', required=True, 
				   help='output directory')
	parser.add_argument("--verbose",'-v', dest= 'v', help="verbose mode",
					action="store_true")
	parser.add_argument("--pysam_version",'-pv', dest= 'pv', help="used to specify the version of pysam available (specific variable names for version 0.7.5)",default="0.8.4")
	global args
	args = parser.parse_args()
	print args
	list = glob.glob(args.bam_dir+'*_sorted.bam')
	print list
	for file in list:
		outfi = args.outdir+"/"+os.path.splitext(os.path.basename(file))[0]+"_coverage.csv";
		if os.path.isfile(outfi):
			print "{0} is already there, we don't recompute".format(outfi)
		else:
			print "we will process file {0} -> will go to {1}".format(file,outfi)
			arguments=(file,outfi,args.pv)
			with concurrent.futures.ProcessPoolExecutor(max_workers=1) as executor:
				result = executor.submit(process, *arguments).result()
				print result


if __name__ == "__main__":
	output = main()
