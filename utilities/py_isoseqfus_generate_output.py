#!/usr/bin/env python
import sys,time,argparse
import numpy as np

def main(args):
	sys.stdout.write("Start analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S") + "\n")
	sys.stdout.flush()
	make_matrix(args.input,args.output,args.fl_read,args.tt_read,args.chr,args.overlap)
	sys.stdout.write("Finish analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S") + "\n")
	sys.stdout.flush()

def decide_flag(gpd_info_list):
	chr1,strand1,tss_set1,tts_set1,score1,shc1,exon_number1,exon_start1,exon_end1,chr2,strand2,tss_set2,tts_set2,score2,shc2,exon_number2,exon_start2,exon_end2 = gpd_info_list
	chr_flag,overlap_flag = ["No","No"]
	if chr1 == chr2:
		chr_flag = "Yes"
		tss1 = min([int(i) for i in tss_set1.split("|")])
		tts1 = max([int(i) for i in tts_set1.split("|")])
		tss2 = min([int(i) for i in tss_set2.split("|")])
		tts2 = max([int(i) for i in tts_set2.split("|")])
		if (tts1 > tss2 and tss1 < tts2) or (tts2 > tss1 and tss2 < tts1):
			overlap_flag = "Yes"
	else:
		chr_flag = "No"
	return chr_flag,overlap_flag

def make_matrix(input_gpd,output_gpd,min_fl_read_rq,min_tt_read,chr_flag_rq,overlap_flag_rq):
	fusion_id_idx = 0
	for line in input_gpd:
		read_set,chr1,strand1,tss_set1,tts_set1,score1,shc1,exon_number1,exon_start1,exon_end1,chr2,strand2,tss_set2,tts_set2,score2,shc2,exon_number2,exon_start2,exon_end2,dev_gene1,dev_iso1,dev_gene2,dev_iso2 = line.strip().split("\t")[:23]
		full_length_read = 0
		total_read= 0
		for read in read_set.split(","):
			total_read += 1
			if read.endswith("F1P1T1"):
				full_length_read += 1
		if full_length_read < min_fl_read_rq or total_read < min_tt_read: continue # check read count

		gpd_info_list = line.strip().split("\t")[1:-4]
		chr_flag,overlap_flag = decide_flag(gpd_info_list)
			
		if chr_flag_rq == "yes" and chr_flag == "Yes": continue #  exclude fusion transcripts if two parts are same chromosome
		if overlap_flag_rq == "yes" and overlap_flag == "Yes": continue # exclude fusion transcripts with overlap between two parts

		if int(shc1.split("_")[-1]) == 0 and int(shc1.split("_")[-2]) != 0: # fusion site for first part
			fusion_site_list1 = [int(i) for i in tss_set1.split("|")]
		elif int(shc1.split("_")[-2]) == 0 and int(shc1.split("_")[-1]) != 0:
			fusion_site_list1 = [int(i) for i in tts_set1.split("|")]
		else:
			sys.stderr.write("Error: The both ends of first part are hard-clipped. Please check the line below.\n"+line)

		if int(shc2.split("_")[-1]) == 0 and int(shc2.split("_")[-2]) != 0: # fusion site for second part
			fusion_site_list2 = [int(i) for i in tss_set2.split("|")]
		elif int(shc2.split("_")[-2]) == 0 and int(shc2.split("_")[-1]) != 0:
			fusion_site_list2 = [int(i) for i in tts_set2.split("|")]
		else:
			sys.stderr.write("Error: The both ends of first part are hard-clipped. Please check the line below.\n"+line)

		fusion_site_min1 = min(fusion_site_list1)
		fusion_site_max1 = max(fusion_site_list1)
		fusion_site_mode1 = max(set(fusion_site_list1),key=fusion_site_list1.count)
		fusion_site_mode_frq1 = fusion_site_list1.count(fusion_site_mode1)
		fusion_site_median1 = int(round(np.median(fusion_site_list1)))
		fusion_site_mean1 = int(round(np.mean(fusion_site_list1)))
		fusion_site_std1 = int(round(np.std(fusion_site_list1)))
		fusion_site_set1 = "_".join(str(i) for i in [fusion_site_min1,fusion_site_max1,fusion_site_mode1,fusion_site_mode_frq1,fusion_site_median1,fusion_site_mean1,fusion_site_std1])
		fusion_site_min2 = min(fusion_site_list2)
		fusion_site_max2 = max(fusion_site_list2)
		fusion_site_mode2 = max(set(fusion_site_list2),key=fusion_site_list2.count)
		fusion_site_mode_frq2 = fusion_site_list2.count(fusion_site_mode2)
		fusion_site_median2 = int(round(np.median(fusion_site_list2)))
		fusion_site_mean2 = int(round(np.mean(fusion_site_list2)))
		fusion_site_std2 = int(round(np.std(fusion_site_list2)))
		fusion_site_set2 = "_".join(str(i) for i in [fusion_site_min2,fusion_site_max2,fusion_site_mode2,fusion_site_mode_frq2,fusion_site_median2,fusion_site_mean2,fusion_site_std2])

		fusion_id_idx += 1
		fusion_id = "Fusion_transcript_" + str(fusion_id_idx)

		fst_part_length = 0
		for i in range(int(exon_number1)):
			fst_part_length += (int(exon_end1.split(",")[i]) - int(exon_start1.split(",")[i]))
		sec_part_length = 0
		for i in range(int(exon_number2)):
			sec_part_length += (int(exon_end2.split(",")[i]) - int(exon_start2.split(",")[i]))

		print >>output_gpd, "\t".join([fusion_id,str(full_length_read),str(total_read),chr_flag,overlap_flag,dev_gene1,dev_iso1,fusion_site_set1,dev_gene2,dev_iso2,fusion_site_set2,chr1,strand1,exon_start1.split(",")[0],exon_end1.split(",")[0],str(fst_part_length),exon_number1,exon_start1,exon_end1,chr2,strand2,exon_start2.split(",")[0],exon_end2.split(",")[0],str(sec_part_length),exon_number2,exon_start2,exon_end2])

	input_gpd.close()
	output_gpd.close()

def do_inputs():
	output_gpd_format = '''
1. Fusion transcript id
2. Number of support full-length long reads
3. Number of support long reads (including both full-length and non-full-length)
4. If two parts of fusion transcript are from same chromosome
5. If two parts of fusion transcript have overlap
6. Derived fusion gene of first part
7. Derived fusion isoform of first part
8. Fusion site information of first part. [min, max, mode, mode_frequency, median, mean, standard deviation]
9. Derived fusion gene of second part
10. Derived fusion isoform of second part
11. Fusion site information of second part. [min, max, mode, mode_frequency, median, mean, standard deviation]

12-19: structue for first part of fusion transcript
12. chromosome
13. strand
14. tss (+)
15. tts (+)
16. sequence length
17. exon number
18. start exon set
19. end exon set

20-27: structue for second part of fusion transcript
20. chromosome
21. strand
22. tss (+)
23. tts (+)
24. sequence length
25. exon number
26. start exon set
27. end exon set'''

	parser = argparse.ArgumentParser(description="Function: make matrix for final fusion transcript",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-i','--input',type=argparse.FileType('r'),required=True,help="Input: fusion-specific gpd file generated by 'py_isoseqfusion_add_anno.py'")
	parser.add_argument('-o','--output',type=argparse.FileType('w'),required=True,help="Output: fusion matrix")
	parser.add_argument('--fl_read',type=int,default=0,help="Minimal number of support full-length long reads")
	parser.add_argument('--tt_read',type=int,default=1,help="Minimal number of support total long reads")
	parser.add_argument('--chr',type=str,choices=['yes','no'],default='no',help="Exclude it if two parts of fusion transcript are from same chromosome")
	parser.add_argument('--overlap',type=str,choices=['yes','no'],default='yes',help="Exclude it if two parts of fusion transcript have any overlap")
	args = parser.parse_args()
	return args

if __name__=="__main__":
	args = do_inputs()
	main(args)
