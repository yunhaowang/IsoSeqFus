#!/usr/bin/env python
import sys,time,argparse

def main(args):
	sys.stdout.write("Start analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S") + "\n")
	sys.stdout.flush()
	dic_chr_iso_spliceSite,dic_chr_iso_info = extract_info_from_annotation(args.anno)
	add_anno(args.input,args.output,dic_chr_iso_spliceSite,dic_chr_iso_info)
	sys.stdout.write("Finish analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S") + "\n")
	sys.stdout.flush()


#=== extract splice site, junction set, and gene region informations from known gene annotation library ===
def extract_info_from_annotation(anno_gpd):
	dic_chr_iso_spliceSite = {}
	dic_chr_iso_info = {}
	for line in anno_gpd:
		gene,iso,chr,strand,tss,tts,cds_start,cds_end,exon_number,exon_start,exon_end = line.rstrip("\n").split("\t")[:11]
		if chr not in dic_chr_iso_info.keys():
			dic_chr_iso_info[chr] = {}
			dic_chr_iso_info[chr][iso] = gene + "&" + tss + "&" + tts # gene + tss + tts combination
			dic_chr_iso_info[chr]["iso_list"] = []
			dic_chr_iso_info[chr]["iso_list"].append(iso) # add isoform by the genomic coordinate
			dic_chr_iso_spliceSite[chr] = {}
			if int(exon_number) > 1: # spliced alignment
				dic_chr_iso_spliceSite[chr][iso] = exon_start.split(",")[1:-1] + exon_end.split(",")[:-2]
			else: # singleton alignment
				dic_chr_iso_spliceSite[chr][iso] = []
		else:
			dic_chr_iso_info[chr][iso] = gene + "&" + tss + "&" + tts
			dic_chr_iso_info[chr]["iso_list"].append(iso)
			if int(exon_number) > 1:
				dic_chr_iso_spliceSite[chr][iso] = exon_start.split(",")[1:-1] + exon_end.split(",")[:-2]
			else:
				dic_chr_iso_spliceSite[chr][iso] = []
	
	anno_gpd.close()
	return dic_chr_iso_spliceSite,dic_chr_iso_info
		
#=== determine derived gene and isoform ===

def determine_derived_loci(gpd_info_list,dic_chr_iso_spliceSite,dic_chr_iso_info):
	chr,strand,tss_set,tts_set,score,shc,exon_number,exon_start,exon_end = gpd_info_list
	tss = min([int(i) for i in tss_set.split("|")])
	tts = max([int(i) for i in tts_set.split("|")])
	if int(exon_number) == 1: # singleton alignemt
		dic_over_pct = {}
		for iso in dic_chr_iso_info[chr]["iso_list"]:
			if int(tts) <= int(dic_chr_iso_info[chr][iso].split("&")[1]): # break the loop when the TTS of reads is < the TSS of annotation
				break
			else:
				if int(tss) < int(dic_chr_iso_info[chr][iso].split("&")[2]) and int(tts) > int(dic_chr_iso_info[chr][iso].split("&")[1]): # with overlap
					over_pct = float(min(int(tts),int(dic_chr_iso_info[chr][iso].split("&")[2]))-max(int(tss),int(dic_chr_iso_info[chr][iso].split("&")[1])))/(int(tts)-int(tss))
					if over_pct not in dic_over_pct.keys(): # collect isoforms with overlap
						dic_over_pct[over_pct] = []
						dic_over_pct[over_pct].append(iso)
					else:
						dic_over_pct[over_pct].append(iso)
		if dic_over_pct != {}: # choose the maximal overlap
			gene_over_set = set()
			max_pct = max(dic_over_pct.keys())
			derived_iso = ",".join(dic_over_pct[max_pct])
			for iso in dic_over_pct[max_pct]:
				gene_over_set.add(dic_chr_iso_info[chr][iso].split("&")[0])
			derived_gene = ",".join(list(gene_over_set))
		else: # no overlap
			derived_iso = "Novel_loci_" + "_".join([chr,str(tss),str(tts)])
			derived_gene = derived_iso
	else: # spliced alignment
		splice_site_list = exon_start.split(",")[1:-1]+exon_end.split(",")[:-2]
		if chr in dic_chr_iso_spliceSite.keys():
			dic_over_pct = {}
			dic_intersect_num = {}
			for iso in dic_chr_iso_info[chr]["iso_list"]:
				if int(tts) <= int(dic_chr_iso_info[chr][iso].split("&")[1]): # break the loop when the TTS of reads is < the TSS of annotation
					break
				else:
	
					if int(tss) < int(dic_chr_iso_info[chr][iso].split("&")[2]) and int(tts) > int(dic_chr_iso_info[chr][iso].split("&")[1]): # with overlap
						intersect_num = len(set(splice_site_list).intersection(dic_chr_iso_spliceSite[chr][iso])) # splice site
						if intersect_num not in dic_intersect_num.keys(): # test if splice sites are included by annotated isoforms
							dic_intersect_num[intersect_num] = []
							dic_intersect_num[intersect_num].append(iso)
						else:
							dic_intersect_num[intersect_num].append(iso)
						over_pct = float(min(int(tts),int(dic_chr_iso_info[chr][iso].split("&")[2]))-max(int(tss),int(dic_chr_iso_info[chr][iso].split("&")[1])))/(int(tts)-int(tss))
						if over_pct not in dic_over_pct.keys(): # overlap percentage
							dic_over_pct[over_pct] = []
							dic_over_pct[over_pct].append(iso)
						else:
							dic_over_pct[over_pct].append(iso)
			if dic_intersect_num != {} and max(dic_intersect_num.keys()) != 0: # splice sites are included by known isoforms
				derived_iso = ",".join(dic_intersect_num[max(dic_intersect_num.keys())])
				derived_gene_set = set()
				for iso in dic_intersect_num[max(dic_intersect_num.keys())]:
					derived_gene_set.add(dic_chr_iso_info[chr][iso].split("&")[0])	
				derived_gene = ",".join(list(derived_gene_set))
			else: # splice sites are not included by any known isoform
				if dic_over_pct != {}: # choose the maximal overlap
					gene_over_set = set()
					max_pct = max(dic_over_pct.keys())
					derived_iso = ",".join(dic_over_pct[max_pct])
					for iso in dic_over_pct[max_pct]:
						gene_over_set.add(dic_chr_iso_info[chr][iso].split("&")[0])
					derived_gene = ",".join(list(gene_over_set))
				else: # no overlap
					derived_iso = "Novel_loci_" + "_".join([chr,str(tss),str(tts)])
					derived_gene = derived_iso
		else:
			derived_iso = "Novel_loci_" + "_".join([chr,str(tss),str(tts)])
			derived_gene = derived_iso

	return derived_gene,derived_iso

def add_anno(fusion_gpd,output_gpd,dic_chr_iso_spliceSite,dic_chr_iso_info):
	for line in fusion_gpd:
		gpd_info_list1 = line.strip().split("\t")[1:10] # first part
		gpd_info_list2 = line.strip().split("\t")[10:19] # second part
		derived_gene1,derived_iso1 = determine_derived_loci(gpd_info_list1,dic_chr_iso_spliceSite,dic_chr_iso_info)
		derived_gene2,derived_iso2 = determine_derived_loci(gpd_info_list2,dic_chr_iso_spliceSite,dic_chr_iso_info)

		print >>output_gpd, line.strip() + "\t" + "\t".join([derived_gene1,derived_iso1,derived_gene2,derived_iso2])

	fusion_gpd.close()
	output_gpd.close()

def do_inputs():
	output_gpd_format = '''
1. read id list, split by ","

2-10: gpd information for first part of chimera alignment

2. chromosome id
3. strand
4. start site of alignment list, split by "|"
5. end site of alignment, split by "|"
6. MAPQ list, split by "|"
7. number of nucleotides that are soft- and hard- clipped by GMAP aligner; SOFTleft_SOFTright_HARDleft_HARDright, split by "|"
8. exon count
9. exon start set
10. exon end set

11-19: gpd information for second part of chimera alignment (structure is same with first part)

20. Derived fusion gene of first part
21. Derived fusion isoform of first part
22. Derived fusion gene of second part
23. Derived fusion isoform of second part'''

	parser = argparse.ArgumentParser(description="Function: add annotation informations (including derived gene and isoform) for fusion isoforms ",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-i','--input',type=argparse.FileType('r'),required=True,help="Input: fusion-specific gpd file generated by 'py_isoseqfusion_construct.py'")
	parser.add_argument('-a','--anno',type=argparse.FileType('r'),required=True,help="Annotation libray (gpd format), should be sorted by 'sort -k3,3 -k5,5n -k6,6n'")
	parser.add_argument('-o','--output',type=argparse.FileType('w'),required=True,help="Output: fusion isoform with annotation information")
	args = parser.parse_args()
	return args

if __name__=="__main__":
	args = do_inputs()
	main(args)
