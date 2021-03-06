#!/usr/bin/env python
import sys,time,argparse

def main(args):
	sys.stdout.write("Start analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S") + "\n")
	sys.stdout.flush()
	make_fusion_gpd(args.input,args.output)
	sys.stdout.write("Finish analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S") + "\n")
	sys.stdout.flush()

def make_fusion_gpd(input_gpd,output_gpd):
	dic_read_first_part = {}
	for line in input_gpd:
		read_id = line.strip().split("\t")[0]
		if read_id not in dic_read_first_part.keys():
			dic_read_first_part[read_id] = line.strip().split("\t")[2:11]
		else:
			first_part = dic_read_first_part[read_id]
			first_part_shc = [int(i) for i in dic_read_first_part[read_id][5].split("_")]
			second_part = line.strip().split("\t")[2:11]
			second_part_shc = [int(i) for i in second_part[5].split("_")]
			if first_part[1] == "+":
				if first_part_shc[2] == 0 and first_part_shc[3] != 0:
					print >>output_gpd, "\t".join([read_id,"\t".join(first_part),"\t".join(second_part)])
				elif first_part_shc[2] != 0 and first_part_shc[3] == 0:
					print >>output_gpd, "\t".join([read_id,"\t".join(second_part),"\t".join(first_part)])
				else:
					sys.stderr.write("Abnormal hard clip for the read below:\n"+read_id+"\n")
					
			else:
				if first_part_shc[2] != 0 and first_part_shc[3] == 0:
					print >>output_gpd, "\t".join([read_id,"\t".join(first_part),"\t".join(second_part)])
				elif first_part_shc[2] == 0 and first_part_shc[3] != 0:
					print >>output_gpd, "\t".join([read_id,"\t".join(second_part),"\t".join(first_part)])
				else:
					sys.stderr.write("Abnormal hard clip for the read below:\n"+read_id+"\n")
			del dic_read_first_part[read_id]

	input_gpd.close()
	output_gpd.close()
	
def do_inputs():
	output_gpd_format = '''
1. read id

2-10: gpd information for first part of chimera alignment

2. chromosome id
3. strand
4. start site of alignment
5. end site of alignment
6. MAPQ 
7. Number of nucleotides that are soft- and hard- clipped by GMAP aligner; SOFTleft_SOFTright_HARDleft_HARDright
8. exon count
9. exon start set
10. exon end set

11-21: gpd information for second part of chimera alignment (structure is same with first part)'''

	parser = argparse.ArgumentParser(description="Function: make fusion-specific gpd file",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-i','--input',type=argparse.FileType('r'),required=True,help="Input: gpd file (generated by 'py_isoseqfusion_polish.py')")
	parser.add_argument('-o','--output',type=argparse.FileType('w'),required=True,help="Output: fusion-specific gpd file")
	args = parser.parse_args()
	return args

if __name__=="__main__":
	args = do_inputs()
	main(args)
