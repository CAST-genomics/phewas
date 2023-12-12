#!/usr/bin/env python3
"""
Perform phewas on all Pan-UKB phenotypes

Example:
./panukb_phewas.py --loc 9:27573483-27573544 --name C9ORF72 --numpts 1 --outdir panukb-test
"""

import argparse
import pandas as pd
import os
import sys

# TODO fail if the command returns non-zero exit status
# TODO don't spit all the output to the terminal
def RunCmd(cmd):
	os.system(cmd)

def MSG(msg):
	sys.stderr.write("[panukb_phewas.py]:"+msg.strip()+"\n")

def ParseCoords(chrom_loc, window):
	chrom_list = chrom_loc.replace('-',':').split(':')
	chrom_num, start, end = chrom_list
	start = int(start) - window
	end = int(end) + window
	chrom_range = f'{chrom_num}:{start}-{end}'
	return chrom_range

def AddTraitCode(outfile, trait, newoutfile):
	newf = open(newoutfile, "w")
	header = True
	with open(outfile, "r") as f:
		for line in f:
			if header:
				newf.write("phenocode\t"+line.strip()+"\n")
				header = False
			else:
				newf.write(trait+"\t"+line.strip()+"\n")
	newf.close()

def main():
	parser = argparse.ArgumentParser(__doc__)
	parser.add_argument("--loc", help="chr:start-end of target variant", type=str, required=True)
	parser.add_argument("--name", help="Name of locus", type=str, required=True)
	parser.add_argument("--window", help="Window size (bp)", type=int, default=50000)
	parser.add_argument("--manifest", help="Path to manifest file", type=str, default="panukb-manifest.csv")
	parser.add_argument("--numpts", help="Number of phenotypes to process (for debugging)", type=int, default=-1)
	parser.add_argument("--outdir", help="Output directory to store results", type=str, required=True)
	args = parser.parse_args()

	chrom_loc = args.loc
	chrom_name = args.name
	chrom_range = ParseCoords(chrom_loc, args.window)
	MSG("Processing locus {} ({})".format(chrom_range, chrom_name))

	manifest = pd.read_csv(args.manifest)

	# set up output directories
	out_dir = os.path.join(args.outdir, chrom_name)
	if not os.path.isdir(args.outdir):
		os.mkdir(args.outdir)
	if not os.path.isdir(out_dir):
		os.mkdir(out_dir)

	###### Process each phenotype one at a time #####
	numpts = 0
	traitinfo_path_list = []
	for index, row in manifest.iterrows():
		trait = row['phenocode']
		MSG("Processing {}".format(trait))
		outfile = os.path.join(out_dir, "{}_{}.tab".format(trait, chrom_name))
		newoutfile = os.path.join(out_dir, "{}_{}_withinfo.tab".format(trait, chrom_name))
		# Download header line and tabix index
		RunCmd("wget -O - {} | zcat | head -n 1 > {}".format(row['aws_link'], outfile))
		# TODO - download these to somewhere else. also, check if already exists
		# Extract tabix info to a file
		RunCmd("tabix {} {} >> {}".format(row['aws_link'], chrom_range, outfile))
		# Add column with trait code
		AddTraitCode(outfile, trait, newoutfile)
		traitinfo_path_list.append(newoutfile)
		# Done processing phenotype. Increment counter
		numpts += 1
		if numpts >= args.numpts and args.numpts > 0: break

	###### Combine data from all phenotypes #######
	final_out_file = os.path.join(args.outdir, "{}_phewas.tab".format(chrom_name))
	MSG("Combining all phenotype data to {}".format(final_out_file))
	combined_df = pd.concat([pd.read_csv(f, sep='\t') for f in traitinfo_path_list])
	combined_df["phenocode"] = combined_df["phenocode"].apply(str)
	combined_df = pd.merge(manifest[["phenocode", "trait_type", "description"]], combined_df, on="phenocode")
	combined_df.to_csv(final_out_file, sep="\t", index=None)

if __name__ == "__main__":
	main()