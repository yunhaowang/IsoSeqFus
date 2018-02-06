#!/usr/bin/env python
import sys,time,argparse

def main(args):
	sys.stdout.write("Start analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S") + "\n")
	sys.stdout.flush()
	construction(args.input,args.input_prefix,args.output)
	sys.stdout.write("Finish analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S") + "\n")
	sys.stdout.flush()


def construction(input_gpd,input_prefix,output_gpd):
	input_fs = open(input_prefix+".FS.sort.concat.gpd","r")
	input_fm = open(input_prefix+".FM.concat.gpd","r")
	input_ss = open(input_prefix+".SS.sort.concat.gpd","r")
	input_sm = open(input_prefix+".SM.concat.gpd","r")
	dic_qname_gpd = {}
	for line in input_gpd:
		qname = line.strip().split("\t")[0]
		dic_qname_gpd[qname] = line.strip().split("\t")
	dic_sec_gpd = {}
	for line in input_ss:
		read_set = line.strip().split("\t")[0]
		dic_sec_gpd[read_set] = line.strip()
	for line in input_sm:
		read_set = line.strip().split("\t")[0]
		dic_sec_gpd[read_set] = line.strip()

	for line in input_fs:
		read_set = line.strip().split("\t")[0]
		read_count = int(line.strip().split("\t")[6])
		if read_count == 1:
			print >>output_gpd, "\t".join(dic_qname_gpd[read_set])
		else:
			fst_read_list = read_set.split(",")
			for i in range(read_count):
				if len(fst_read_list) == 0:
					break
				for sec_read_set in dic_sec_gpd.keys():
					sec_read_list = sec_read_set.split(",")
					intersect_list = list(set(fst_read_list).intersection(sec_read_list))
					if len(intersect_list) != 0:
						new_read_set,new_fst_tss,new_fst_tts,new_fst_mapq,new_fst_shc,new_sec_tss,new_sec_tts,new_sec_mapq,new_sec_shc = [],[],[],[],[],[],[],[],[]
						for ist_read in intersect_list:
							new_read_set.append(ist_read)
							new_fst_tss.append(dic_qname_gpd[ist_read][3])
							new_fst_tts.append(dic_qname_gpd[ist_read][4])
							new_fst_mapq.append(dic_qname_gpd[ist_read][5])
							new_fst_shc.append(dic_qname_gpd[ist_read][6])
							new_sec_tss.append(dic_qname_gpd[ist_read][12])
							new_sec_tts.append(dic_qname_gpd[ist_read][13])
							new_sec_mapq.append(dic_qname_gpd[ist_read][14])
							new_sec_shc.append(dic_qname_gpd[ist_read][15])

							fst_read_list.remove(ist_read)
						gpd_info_list = dic_qname_gpd[intersect_list[0]]
						gpd_info_list[0] = ",".join(new_read_set)
						gpd_info_list[3] = "|".join(new_fst_tss)
						gpd_info_list[4] = "|".join(new_fst_tts)
						gpd_info_list[5] = "|".join(new_fst_mapq)
						gpd_info_list[6] = "|".join(new_fst_shc)
						gpd_info_list[8] = str(min([int(i) for i in new_fst_tss])) + ","
						gpd_info_list[9] = str(max([int(i) for i in new_fst_tts])) + ","
						gpd_info_list[12] = "|".join(new_sec_tss)
						gpd_info_list[13] = "|".join(new_sec_tts)
						gpd_info_list[14] = "|".join(new_sec_mapq)
						gpd_info_list[15] = "|".join(new_sec_shc)
						if int(gpd_info_list[16]) == 1:
							gpd_info_list[17] = str(min([int(i) for i in new_sec_tss])) + ","
							gpd_info_list[18] = str(max([int(i) for i in new_sec_tts])) + ","
						else:
							gpd_info_list[17] = str(min([int(i) for i in new_sec_tss])) + "," + ",".join(gpd_info_list[17].split(",")[1:])
							gpd_info_list[18] = ",".join(gpd_info_list[18].split(",")[:-2]) + "," + str(max([int(i) for i in new_sec_tts])) + ","
						print >>output_gpd, "\t".join(gpd_info_list)
	for line in input_fm:
		read_set = line.strip().split("\t")[0]
		read_count = int(line.strip().split("\t")[6])
		if read_count == 1:
			print >>output_gpd, "\t".join(dic_qname_gpd[read_set])
		else:
			fst_read_list = read_set.split(",")
			for i in range(read_count):
				if len(fst_read_list) == 0:
					break
				for sec_read_set in dic_sec_gpd.keys():
					sec_read_list = sec_read_set.split(",")
					intersect_list = list(set(fst_read_list).intersection(sec_read_list))
					if len(intersect_list) != 0:
						new_read_set,new_fst_tss,new_fst_tts,new_fst_mapq,new_fst_shc,new_sec_tss,new_sec_tts,new_sec_mapq,new_sec_shc = [],[],[],[],[],[],[],[],[]
						for ist_read in intersect_list:
							new_read_set.append(ist_read)
							new_fst_tss.append(dic_qname_gpd[ist_read][3])
							new_fst_tts.append(dic_qname_gpd[ist_read][4])
							new_fst_mapq.append(dic_qname_gpd[ist_read][5])
							new_fst_shc.append(dic_qname_gpd[ist_read][6])
							new_sec_tss.append(dic_qname_gpd[ist_read][12])
							new_sec_tts.append(dic_qname_gpd[ist_read][13])
							new_sec_mapq.append(dic_qname_gpd[ist_read][14])
							new_sec_shc.append(dic_qname_gpd[ist_read][15])

							fst_read_list.remove(ist_read)
						gpd_info_list = dic_qname_gpd[intersect_list[0]]
						gpd_info_list[0] = ",".join(new_read_set)
						gpd_info_list[3] = "|".join(new_fst_tss)
						gpd_info_list[4] = "|".join(new_fst_tts)
						gpd_info_list[5] = "|".join(new_fst_mapq)
						gpd_info_list[6] = "|".join(new_fst_shc)
						gpd_info_list[8] = str(min([int(i) for i in new_fst_tss])) + "," + ",".join(gpd_info_list[8].split(",")[1:])
						gpd_info_list[9] = ",".join(gpd_info_list[9].split(",")[:-2]) + "," + str(max([int(i) for i in new_fst_tts])) + ","
						gpd_info_list[12] = "|".join(new_sec_tss)
						gpd_info_list[13] = "|".join(new_sec_tts)
						gpd_info_list[14] = "|".join(new_sec_mapq)
						gpd_info_list[15] = "|".join(new_sec_shc)
						if int(gpd_info_list[16]) == 1:
							gpd_info_list[17] = str(min([int(i) for i in new_sec_tss])) + ","
							gpd_info_list[18] = str(max([int(i) for i in new_sec_tts])) + ","
						else:
							gpd_info_list[17] = str(min([int(i) for i in new_sec_tss])) + "," + ",".join(gpd_info_list[17].split(",")[1:])
							gpd_info_list[18] = ",".join(gpd_info_list[18].split(",")[:-2]) + "," + str(max([int(i) for i in new_sec_tts])) + ","
						print >>output_gpd, "\t".join(gpd_info_list)

	input_gpd.close()
	output_gpd.close()
	input_fs.close()
	input_fm.close()
	input_ss.close()
	input_sm.close()
	
def do_inputs():
	output_gpd_format = '''
1. read id list, split by ","

2-10: gpd information for first part of chimera alignment

2. chromosome id
3. strand
4. start site of alignment list, split by "|"
5. end site of alignment list, split by "|"
6. MAPQ list, split by "|"
7. number of nucleotides that are soft- and hard- clipped by GMAP aligner; SOFTleft_SOFTright_HARDleft_HARDright, list, split by "|"
8. exon count
9. exon start set
10. exon end set

11-19: gpd information for second part of chimera alignment (structure is same with first part)'''

	parser = argparse.ArgumentParser(description="Function: split fusion gpd file into 4 files including 1) first part, singleton alignment; 2) first part, multi-exon alignment; 3) second part, singleton alignment; and 4) second part, multi-exon alignment.",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-i','--input',type=argparse.FileType('r'),required=True,help="Input: fusion-specific gpd file generated by 'py_isoseqfusion_make_fusion_gpd.py'")
	parser.add_argument('-p','--input_prefix',type=str,default="isoseqfus",help="Input prefix including 4 files (generated by 'py_isoseqfusion_concat_sgt.py': input_prefix.FS.sort.concat.gpd and input_prefix.SS.sort.concat.gpd; and 'py_isoseqfusion_concat_mlt.py': input_prefix.FM.concat.gpd and input_prefix.SM.concat.gpd)")
	parser.add_argument('-o','--output',type=argparse.FileType('w'),required=True,help="Output: constructed fusion isoforms")
	args = parser.parse_args()
	return args

if __name__=="__main__":
	args = do_inputs()
	main(args)
