#!/usr/bin/env python
import sys,time,argparse

def main(args):
	sys.stdout.write("Start analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S") + "\n")
	sys.stdout.flush()
	convert_gpd2gtf(args.input,args.output)
	sys.stdout.write("Finish analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S") + "\n")
	sys.stdout.flush()

def convert_gpd2gtf(input_gpd,output_gtf):
	for line in input_gpd:
		fusion_id,full_length_LR_count,LR_count,same_chromosome,overlap_bewteen_two_parts,fst_derived_gene,fst_derived_isoform,fst_fusion_site_info,sec_derived_gene,sec_derived_isoform,sec_fusion_site_info,chr1,strand1,tss1,tts1,sequence_length1,exon_number1,start_exon_set1,end_exon_set1,chr2,strand2,tss2,tts2,sequence_length2,exon_number2,start_exon_set2,end_exon_set2 = line.strip().split("\t")
		read_count_info = 'full_length_LR_count "' + full_length_LR_count + '"; LR_count "' + LR_count + '"; same_chr "' + '"; overlap_two_parts "' + '";'
		fst_id = 'gene_id "' + fusion_id + '"; transcript_id "' + fusion_id + '_First_Part' + '"; '
		sec_id = 'gene_id "' + fusion_id + '"; transcript_id "' + fusion_id + '_Second_Part' + '"; '
		fst_info = 'fst_derived_gene "' + fst_derived_gene + '"; fst_derived_isoform "' + fst_derived_isoform + '"; fst_fusion_site_info "' + fst_fusion_site_info + '"; fst_sequence_length "' + sequence_length1 + '"; fst_exon_number "' + exon_number1 + '"; '
		sec_info = 'sec_derived_gene "' + sec_derived_gene + '"; sec_derived_isoform "' + sec_derived_isoform + '"; sec_fusion_site_info "' + sec_fusion_site_info + '"; sec_sequence_length "' + sequence_length2 + '"; sec_exon_number "' + exon_number2 + '"; '
		# print first part
		fst_attribute_info = "".join([fst_id,fst_info,sec_info,read_count_info])
		print >>output_gtf, "\t".join([chr1,".","transcript",str(int(tss1)+1),str(tts1),".",strand1,".",fst_attribute_info])
		if strand1 == "+":
			for i in range(0,int(exon_number1)):
				exon_id = ' exon_number "' + str(i+1) + '";'
				fst_attribute_info = "".join([fst_id,fst_info,sec_info,read_count_info,exon_id])

				print >>output_gtf, "\t".join([chr1,".","exon",str(int(start_exon_set1.split(",")[i])+1),end_exon_set1.split(",")[i],".",strand1,".",fst_attribute_info])


		else:
			for i in range(0,int(exon_number1)):
				exon_id = ' exon_number "' + str(int(exon_number1)-i) + '";'
				fst_attribute_info = "".join([fst_id,fst_info,sec_info,read_count_info,exon_id])
				print >>output_gtf, "\t".join([chr1,".","exon",str(int(start_exon_set1.split(",")[i])+1),end_exon_set1.split(",")[i],".",strand1,".",fst_attribute_info])
		# print second part
		sec_attribute_info = "".join([sec_id,sec_info,fst_info,read_count_info])
		print >>output_gtf, "\t".join([chr1,".","transcript",str(int(tss1)+1),str(tts1),".",strand1,".",sec_attribute_info])
		if strand2 == "+":
			for i in range(0,int(exon_number2)):
				exon_id = ' exon_number "' + str(i+1) + '";'
				sec_attribute_info = "".join([sec_id,sec_info,fst_info,read_count_info,exon_id])
				print >>output_gtf, "\t".join([chr2,".","exon",str(int(start_exon_set2.split(",")[i])+1),end_exon_set2.split(",")[i],".",strand2,".",sec_attribute_info])
		else:
			for i in range(0,int(exon_number2)):
				exon_id = ' exon_number "' + str(int(exon_number2)-i) + '";'
				sec_attribute_info = "".join([sec_id,sec_info,fst_info,read_count_info,exon_id])
				print >>output_gtf, "\t".join([chr2,".","exon",str(int(start_exon_set2.split(",")[i])+1),end_exon_set2.split(",")[i],".",strand2,".",sec_attribute_info])

	input_gpd.close()
	output_gtf.close()

def do_inputs():
	parser = argparse.ArgumentParser(description="Convert gpd to gtf",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-i','--input',type=argparse.FileType('r'),required=True,help="Input: gpd file")
	parser.add_argument('-o','--output',type=argparse.FileType('w'),required=True,help="Output: gtf file")
	args = parser.parse_args()
	return args

if __name__=="__main__":
	args = do_inputs()
	main(args)
