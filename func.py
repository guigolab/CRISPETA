#sys and regular functions
import argparse, tempfile, os, re, math, sys, resource, time, datetime, _mysql
from subprocess import call
from os import listdir
from config import *

#Biopithon
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from operator import itemgetter

#graphics
import numpy as np
try:
	import plotly
	_plotly=1
except ImportError:
	_plotly=0
if _plotly==1:
	import plotly.plotly as py
	import plotly.graph_objs as go
	from plotly.graph_objs import Scatter, Layout

#import matplotlib.pyplot as plt
#from matplotlib.backends.backend_pdf import PdfPages


gibson_5 = 'atcttGTGGAAAGGACGAAACACCg'
constant_h1 = 'GTTTTAGAGCTAGAAGAGACGGAATTCCTAGGATCCCGTCTCTCTGTATGAGACCACTCTTTCCC'
gibson_3 = 'gttttagagctaGAAAtagcaagttaaaataaggc'

def check_user_options(options):
	'''
	Check user defined options
	'''
	#-i input
	if not os.path.isfile(options.infile):
		sys.exit("  Input file (-i [Path/File]) does not exist... Exiting the program.")

	#-g Genome
	if not os.path.isfile(options.genome):
		sys.exit('  Genome file (-g [Path/File]) does not exist... Exiting the program.')

	#-pmask Positive mask
	if options.pmask != None and not os.path.isfile(options.pmask):
		sys.exit('  Positive mask file (-pmask [Path/File]) does not exist... Exiting the program.')

	#-nmask Negative mask
	if options.nmask != None and not os.path.isfile(options.nmask):
		sys.exit('  Negative mask file (-m [Path/File]) does not exist... Exiting the program.')

	#-off Off targets
	#if not os.path.isfile(options.off_targets):
	#	sys.exit('  Off-target file (-off [Path/File]) does not exist... Exiting the program.')


def try_except(success_list,success_position, failure, *exceptions):
    try:
        return success_list[success_position]
    except exceptions:
        return failure if callable(failure) else failure


def try_number(value):
	try:
		return int(value)
	except ValueError:
		return 0


def check_input(input_file,tmp_dir):
	infile = open(input_file, "r")
	outfile = open(tmp_dir+'/checked_tmp_bed_file','w')
	ucsc=re.compile("^[a-zA-Z0-9]+:[0-9,\.]+-[0-9,\.]+:.")
	n = 1
	ids = []

	for line in infile:
		line = line.strip()
		n += 1
		if ucsc.match(line):
			strand = line[-1]
			line = line.split(":")
			line = [line[0]]+line[1].split("-")
			chrom = line[0]
			start = line[1].replace(",","").replace(".","")
			end = line[2].replace(",","").replace(".","")
			ID = 'unknown'+str(n)
			strand = '+'

		else:
			line = line.split()
			if len(line < 3): continue 
			chrom = line[0]
			start = line[1]
			end = line[2]
			ID = try_except(line,3,'unknown'+str(n),IndexError)
			strand = try_except(line,5,'+',IndexError)





def target_bed(input_file, design_up, exclude_up, exclude_down, design_down, tmp_dir):
	'''
	Reads input BED file and creates a new bed file with regions
	modified with input options (up and down ranges)
	'''
	infile = open(input_file, "r")
	outfile = open(tmp_dir+'/tmp_bed_file','w')
	
	m=0 		#Total number of regions in input BED file
	n=0 		#Number of modified regions
	ids = {}	#Id's of each modified region

	for line in infile:
		m += 1
		arguments = line.strip().split()
				
		#check distande and ranges
		if len(arguments) < 3:
			print('WARNING: arguments missing in region '+line)
			continue
		
		else:
			n+=1
			chrom = try_except(arguments,0,'chr?',IndexError)
			start = try_except(arguments,1,'',IndexError)
			end = try_except(arguments,2,'',IndexError)
			ID = try_except(arguments,3,'',IndexError) if try_except(arguments,3,'',IndexError) != '.' else 'unknown'+str(n)
			ID = 'unknown'+str(n) if ID == '' else ID
			if try_except(arguments,5,'+',IndexError) == '.':
				strand = '+'
			else:
				strand = try_except(arguments,5,'+',IndexError)
			ids['>'+chrom+':'+str(int(start)-design_up-exclude_up)+'-'+str(int(end)+exclude_down+design_down)]=(ID)

			if strand == '-':
				design_up, design_down = design_down, design_up
				exclude_up, exclude_down = exclude_down, exclude_up
				ids['>'+chrom+':'+str(int(start)-design_up-exclude_up)+'-'+str(int(end)+exclude_down+design_down)]=(ID)

			outfile.write(chrom+'\t'+
				str(int(start)-design_up-exclude_up)+'\t'+
				str(int(end)+design_down+exclude_down)+'\t'+
				str(ID)+'\t'+'.'+'\t'+
				strand+'\n')

	infile.close()
	outfile.close()
	if n == 0 or os.path.getsize(tmp_dir+'/tmp_bed_file') == 0:
		sys.exit('No regions to analyze... Exiting the program.')	#If modified regions = 0
	else:
		return([m,n,ids])


def get_sequences(genome_file_name, input_file, output_file, tmp_dir):
	'''
	Create a new fasta file by obtaining the sequences from the specified genome fasta file using bedtools
	'''

	genome = genome_file_name
	input_bed_file = tmp_dir+'/'+input_file
	output_fasta_file = tmp_dir+'/'+output_file
	command_line='bedtools getfasta -fi '+genome+' -bed '+input_bed_file+' -s -fo '+output_fasta_file
	os.system(command_line)
	with open(output_fasta_file) as f:
	    for i, l in enumerate(f):
	        pass
	os.remove(tmp_dir+'/'+input_file)
	if os.path.getsize(output_fasta_file) == 0:
		sys.exit('No regions to analyze... Exiting the program.')	#If 0 sequences are obtained
	else:
		return i + 1


def get_gRNAs(arguments, seq, reverse_seq, design_up, design_down, exclude_up, exclude_down, strand, iscore, tmp_dir, genome, off_target_file, conn):
	'''
	Analyze sequence to find all posibles sgRNAs and filter the undesired ones
	'''
	gRNAs_upstream_downstream_filtered = [0,0,[0,0,0,0,0]]		#sgRNAs upstream and downstream in forward chain and sgRNAs filtered
	direction = "forward" if strand == '+' else "reverse"		#Direction of the positive chain
	reverse_strand = '-' if strand == '+' else '+'				#Sign of reverse strand of the region in the input BED file

	new_files=open(tmp_dir+'/upstream_bed','w')
	new_files=open(tmp_dir+'/downstream_bed','w')
	new_files.close()

	for l in [seq,strand],[reverse_seq,reverse_strand]:
		n = search_gRNA_in_sequence(l[0], arguments, design_up, design_down, exclude_up, exclude_down, direction, iscore, off_target_file, l[1], tmp_dir, conn)
		gRNAs_upstream_downstream_filtered[0] += n[0]			#Add number of positive sgRNAs upstream  in chain 
		gRNAs_upstream_downstream_filtered[1] += n[1]			#Add number of positive sgRNAs downstream in chain
		gRNAs_upstream_downstream_filtered[2][0] += n[2][0]		#Add number of 
		gRNAs_upstream_downstream_filtered[2][1] += n[2][1]		#Add number of 
		gRNAs_upstream_downstream_filtered[2][2] += n[2][2]		#Add number of 
		gRNAs_upstream_downstream_filtered[2][3] += n[2][3]		#Add number of 
		gRNAs_upstream_downstream_filtered[2][4] += n[2][4]		#Add number of 

	file_list = ['upstream_bed', 'downstream_bed']				#Sorting BED files with sgRNAs info
	sort_bed_files(tmp_dir, file_list,1)						#Sorting BED files with sgRNAs info

	return(gRNAs_upstream_downstream_filtered)


def search_gRNA_in_sequence(sequence, arguments, design_up, design_down, exclude_up, exclude_down, direction, iscore, off_target_file, strand_sign, tmp_dir, conn):

	grnas = [0,0,[0,0,0,0,0]]
	for m in re.finditer(r'(?=(.{24}.GG...))',str(sequence)):															#All sgRNAs that fits the pattern
		
		filters = gRNA_filter(arguments, m.start(), design_up, design_down, direction, iscore, off_target_file, m.group(1), conn)	#Filter out sgRNAs
		grnas[2][0] += filters[1][0]
		grnas[2][1] += filters[1][1]
		grnas[2][2]	+= filters[1][2]
		grnas[2][3] *= filters[1][3]
		grnas[2][4] += filters[1][4]
		if filters[0]:
			continue

		n = match_to_bed(arguments,							#arguments = [chr,start,end]
			m.start(),										#m_start
			design_up, design_down,							#ranges
			exclude_up, exclude_down,						#distance from TSS
	 		m.group(1),										#sgRNA
			strand_sign, tmp_dir)							#sgRNA strand, tmp_dir
	
		grnas[0] += n[0]
		grnas[1] += n[1]

	return grnas


def gRNA_filter(arguments, start, design_up, design_down, direction, iscore, off_target_file, group, conn):
	'''
	Filter out sgRNAs
	'''
	if 'TTTTT' in group[4:-3]:															#remove sgRNAs with TTTTT stop signal for polimerase III
		return [1,[1,0,0,0,0]]
	if 'N' in group[4:-3]:
		return [1,[0,1,0,0,0]]
	if check_grna_position(arguments, start, design_up, design_down, direction):		#Avoid sgRNAs in targeting region
		return [1,[0,0,1,0,0]]
	if calc_score(group) < iscore:														#remove sgRNAs with score < min individual score
		return [1,[0,0,0,1,0]]
	if check_grna_off_targets(group[4:-6], off_target_file, conn) == 0:					#Filter sgRNAs by off-targets
		return [1,[0,0,0,0,1]]
	return [0,[0,0,0,0,0]]

def check_grna_position(arguments, match_position, design_up, design_down, direction):
	'''
	Checks if the sgRNA falls into a non interesting region
	'''
	if direction == "forward":
		cut_position = int(arguments[1])+match_position+21
	else:
		cut_position = int(arguments[2])-match_position-21

	undesired_region = range(int(arguments[1])+1+design_up,int(arguments[2])-design_down)
	if cut_position in undesired_region:
		return 1


def calc_score(s):
    '''
    Doench score
    '''
    s_list = list(s)
    if len(s_list) != 30:
    	return 0
    s_20mer = s[4:24]
    nuc_hash = {'A':0, 'T':1, 'C':2, 'G':3}
    score = 0.597636154
    gc = s_20mer.count('G')+s_20mer.count('C')
    gc_low = -0.202625894
    gc_high = -0.166587752
    if gc < 10:
        gc_val = abs(gc-10)
        score = score+(gc_val*gc_low)
    elif gc > 10:
        gc_val = gc-10
        score = score+(gc_val*gc_high)
    #rows[1-30]cols['ATCG']
    sing_nuc_hash = {'G2':-0.275377128,'A3':-0.323887456,'C3':0.172128871,'C4':-0.100666209,'C5':-0.20180294, \
                    'G5':0.245956633,'A6':0.036440041,'C6':0.098376835,'C7':-0.741181291,\
                    'G7':-0.393264397,'A12':-0.466099015,'A15':0.085376945,'C15':-0.013813972,\
                    'A16':0.272620512,'C16':-0.119022648,'T16':-0.285944222,'A17':0.097454592,\
                    'G17':-0.17554617,'C18':-0.345795451,'G18':-0.678096426,'A19':0.22508903,\
                    'C19':-0.507794051,'G20':-0.417373597,'T20':-0.054306959,'G21':0.379899366,\
                    'T21':-0.090712644,'C22':0.057823319,'T22':-0.530567296,'T23':-0.877007428,\
                    'C24':-0.876235846,'G24':0.278916259,'T24':-0.403102218,'A25':-0.077300704,\
                    'C25':0.287935617,'T25':-0.221637217,'G28':-0.689016682,'T28':0.117877577,\
                    'C29':-0.160445304,'G30':0.386342585}
    #score_mat = np.matrix('0 0 0 0;0 0 0 -0.275377128;-0.323887456 0 0.172128871 0;0 0 -0.100666209 0;0 0 -0.20180294 0.245956633;0.036440041 0 0.098376835 0;0 0 -0.741181291 -0.393264397;0 0 0 0;0 0 0 0;0 0 0 0;0 0 0 0;-0.466099015 0 0 0;0 0 0 0;0 0 0 0;0.085376945 0 -0.013813972 0;0.272620512 -0.285944222 -0.119022648 0;0.097454592 0 0 -0.17554617;0 0 -0.345795451 -0.678096426;0.22508903 0 -0.507794051 0;0 -0.054306959 0 -0.417373597;0 -0.090712644 0 0.379899366;0 -0.530567296 0.057823319 0;0 -0.877007428 0 0;0 -0.403102218 -0.876235846 0.278916259;-0.077300704 -0.221637217 0.287935617 0;0 0 0 0;0 0 0 0;0 0.117877577 0 -0.689016682;0 0 -0.160445304 0;0 0 0 0.386342585')
    dinuc_hash = {'GT2':-0.625778696,'GC5':0.300043317,'AA6':-0.834836245,'TA6':0.760627772,'GG7':-0.490816749,'GG12':-1.516907439,'TA12':0.7092612,'TC12':0.496298609,'TT12':-0.586873894,'GG13':-0.334563735,'GA14':0.76384993,'GC14':-0.53702517,'TG17':-0.798146133,'GG19':-0.66680873,'TC19':0.353183252,'CC20':0.748072092,'TG20':-0.367266772,'AC21':0.568209132,'CG21':0.329072074,'GA21':-0.836456755,'GG21':-0.782207584,'TC22':-1.029692957,'CG23':0.856197823,'CT23':-0.463207679,'AA24':-0.579492389,'AG24':0.649075537,'AG25':-0.077300704,'CG25':0.287935617,'TG25':-0.221637217,'GT27':0.117877577,'GG29':-0.697740024}
    for i,nuc in enumerate(s_list):
        key = nuc+str(i+1)
        if sing_nuc_hash.has_key(key):
            nuc_score = sing_nuc_hash[key]
        else:
            nuc_score = 0
        #nuc_score = score_mat[i,nuc_hash[nuc]]
        score = score+nuc_score
        if i<29:
            dinuc = nuc+s[i+1]+str(i+1)
            if dinuc in dinuc_hash.keys():
                score = score+dinuc_hash[dinuc]
    partial_score = math.e**-score
    final_score = 1/(1+partial_score)
    return final_score


def use_database(host='localhost', dbuser='crispeta', database='crispeta', dbpass='pwd'):
	#DB_HOST = host
	#DB_USER = dbuser
	#DB_NAME = database 
	#DB_PASS = dbpass
	from config import config
	dic = config['mysql_db']
	#conn = _mysql.connect(user=DB_USER,passwd=DB_PASS,host=DB_HOST,db=DB_NAME)
	conn = _mysql.connect(user=dic['user'],passwd=dic['passwd'],host=dic['host'],db=dic['db'])
	return conn

def check_grna_off_targets(gRNA, off_target_file, conn):
	'''
	Checks the number of offtargets for a sgRNA
	'''
	from config import config
	t=off_target_file.split(',')
	for i,j in enumerate(t):
		if j=='x':
			t[i]=float("inf")
		else:
			t[i]=int(j)
	table = config['table']
	grna = gRNA
	conn.query("SELECT * FROM "+table+" where grna = '"+grna+"';")
	result = conn.store_result()
	off = result.fetch_row()
	if len(off) == 0:
		return 0
	if int(off[0][1]) <= t[0] and int(off[0][2]) <= t[1] and int(off[0][3]) <= t[2] and int(off[0][4]) <= t[3] and int(off[0][5]) <= t[4]:
		return 1
	else:
		return 0


def match_to_bed(arguments, m_start, design_up, design_down, exclude_up, exclude_down, gRNA, strand, tmp_dir):
	'''
	Add sgRNAs to a BED file
	'''	
	gRNAs_upstream = 0
	gRNAs_downstream = 0

	if arguments[3] == '+':
		cut_position_up =  int(arguments[1]) + design_up - 21
		cut_position_down =  int(arguments[2]) - design_down - 21
	else:
		cut_position_up = int(arguments[2]) - design_up - 9
		cut_position_down = int(arguments[1]) + design_down - 9
	
	if strand == '+':
		position = int(arguments[1]) + m_start
	else:
		position = int(arguments[2]) - m_start - 31
		
	if position < cut_position_up and arguments[3] == '+':
		gRNAs_upstream+=1
		name = 'upstream_bed'
	elif position > cut_position_down and arguments[3] == '+':
		gRNAs_downstream+=1		
		name = 'downstream_bed'
	elif position < cut_position_down and arguments[3] == '-':
		gRNAs_downstream+=1	
		name = 'downstream_bed'
	elif position > cut_position_up and arguments[3] == '-':
		gRNAs_upstream+=1
		name = 'upstream_bed'
	else:
		return([0,0])

	write_in_bed(arguments, name, position, gRNA, strand, tmp_dir)
	return([gRNAs_upstream,gRNAs_downstream])


def write_in_bed(arguments, file_name, position, gRNA, strand, tmp_dir,):
	'''
	Writes the positions of sgRNA in BED format
	'''
	score=round(calc_score(gRNA),3)
	with open(tmp_dir+'/'+file_name, 'a') as file_bed:
		file_bed.write(arguments[0]+'\t'+
			str(position+4)+'\t'+
			str(position+27)+'\t'+
			gRNA[4:27]+'\t'+
			str(score)+'\t'+
			strand+'\t'+'0'+'\t'+'0'+'\n')
	file_bed.close()


def sort_bed_files(tmp_dir, file_list,key_item):
	'''
	sort all bed files in the file_list
	'''
	# for element in file_list:
	# 	if os.path.exists(element):
	# 		call('sort -k1,1 -k2,2n '+element+' -o '+element, shell=True)
	# 	else:
	# 		pass
	for element in file_list:
		if os.path.exists(tmp_dir+'/'+element):								#If file exist
			with open(tmp_dir+'/'+element) as fin:							#Create a list with the lines in the file
				lines = [line.split() for line in fin]
    		lines.sort(key=itemgetter(key_item))							#Sort list by start position
    		with open(tmp_dir+'/output.txt', 'w') as fout:
    			for el in lines:
        			fout.write('{0}\n'.format('\t'.join(el)))				#Write sorted list into new file
        	os.rename(tmp_dir+'/'+'output.txt', tmp_dir+'/'+element)		#Overwrite the file with old name


def new_design_regions(n, design_up, design_down):
	if n[0] == 0 and n[1] == 0:
		return [design_up*2, design_down*2]
	elif n[0] == 0 and n[1] != 0:
		return [design_up*2, 0]
	elif n[0] != 0 and n[1] == 0:
		return [0, design_down*2]
	else:
		return [0,0]


def strict_bed_file(arguments, up_range, new_up_range, down_range, new_down_range, tmp_dir):
	
	of = open(tmp_dir+'/strict_bed', 'w')
	
	start = int(arguments[1])+up_range-new_up_range
	end = int(arguments[2])-down_range+new_down_range
	of.write(arguments[0]+'\t'+str(start)+'\t'+str(end)+'\t'+'.'+'\t'+'.'+'\t'+arguments[3]+'\n')

	of.close()


def strict_sequence(input_file, tmp_dir):
	'''
	Returns first sequence in a fasta file
	'''
	fi = open(tmp_dir+'/'+input_file, 'r')
	for line in fi:
		if line.startswith('>'):
			strand = line[-3]
			args = line[1:-4].replace(":",",").replace("-",",",1).strip().split(",") + [strand]
		else:
			seq = line.strip()
	fi.close()
	os.remove(tmp_dir+'/'+input_file)
	return([seq]+args)


def write_browser_target_regions(tmp_dir, output_name, arguments, new_design_up, exclude_up, new_design_down, exclude_down):

	fb = open(tmp_dir+'/'+output_name, 'a')
	chrom = arguments[0]
	start1 = str(int(arguments[1])+new_design_up+exclude_up)
	start2 = str(int(arguments[1])+new_design_up)
	start3 = str(int(arguments[1]))
	end1 = str(int(arguments[2])-new_design_down-exclude_down)
	end2 = str(int(arguments[2])-new_design_down)
	end3 = str(int(arguments[2]))
	ID = arguments[4]
	tab = '\t'
	fb.write(chrom+tab+start1+tab+end1+tab+ID+tab+'0	+	0	0	0,0,0'+'\n')
	fb.write(chrom+tab+start2+tab+end2+tab+ID+tab+'0	+	0	0	180,180,130'+'\n')
	fb.write(chrom+tab+start3+tab+end3+tab+ID+tab+'0	+	0	0	40,175,40'+'\n')
	fb.close()


def get_masks(arguments, pmask, nmask, iscore, tmp_dir, element_name):
	'''
	Compare sgRNA location with DNAse hypersensitivity regions
	and re-sctructure the files in order to clean all unnecesary data
	'''
	if os.path.isfile(tmp_dir+'/'+element_name):
		if pmask == None:							#Positive mask
			pass
		else:
			get_pmask(element_name, pmask, tmp_dir)

		if nmask == None:							#Negative mask
		 	pass
		else:
			get_nmask(element_name, nmask, tmp_dir)

		rank_mask(element_name, tmp_dir)



def get_pmask(element, pmask, tmp_dir):
	'''
	Does intersection between two bed files and creates a new bed file with the info obtained
	'''
	awk = "awk '{OFS=\"\t\";print $1,$2,$3,$4,$5,$6,$(NF),$8}' > "
	call('intersectBed -a '+tmp_dir+'/'+element+' -b '+pmask+' -sorted -f 1 -wao | '+awk+tmp_dir+'/x.bed', shell=True)	
	os.rename(tmp_dir+'/x.bed', tmp_dir+'/'+element)


def get_nmask(element, nmask, tmp_dir):
	'''
	Intersect sgRNA bed file with another bed file with undesired regions
	'''
	awk = "awk '{OFS=\"\t\";print $1,$2,$3,$4,$5,$6,$7,$(NF)}' > "
	call('intersectBed -a '+tmp_dir+'/'+element+' -b '+nmask+' -sorted -f 1 -wao | '+awk+tmp_dir+'/x.bed', shell=True)	
	os.rename(tmp_dir+'/x.bed', tmp_dir+'/'+element)


def rank_mask(element, tmp_dir):
	'''
	Make the combination of masks scores
	'''
	fi = open(tmp_dir+'/'+element, 'r')
	fo = open(tmp_dir+'/'+element+'_x', 'a')

	for line in fi:
		arg = line.strip().split('\t')
		pmask = int(arg[-2])
		nmask = 1 if arg[-1] == "0" else 0
		if pmask == 1 and nmask == 1:
			combined_mask = 2
		elif pmask == 1 or nmask == 1:
			combined_mask = 1
		else:
			combined_mask = 0
		for i in arg[:-3]+[str(combined_mask)]:
			fo.write(i+'\t')
		fo.write('\n')
	os.rename(tmp_dir+'/'+element+'_x', tmp_dir+'/'+element)


def make_pairs(arguments, file_list, number_results, tmp_dir, outfile, variety_results, pscore, design_up, design_down, combined_method, construct_method, rank_method, seq_id, strict):
   	'''
	Reads files with sgRNAs up/down stream and make sgRNA pairs
	'''
	n_pairs = 0
	n_score = 0
	n_no_G = 0
	n_ok = 0
	arg = arguments

	f_up = open(tmp_dir+'/'+file_list[0], 'r')

	with open(tmp_dir+'/gRNA_pairs', 'a') as of:

		for line_up in f_up:

			up = line_up.strip().split('\t')
			f_down = open(tmp_dir+'/'+file_list[1], 'r')

			for line_down in f_down:
				
				down = line_down.strip().split('\t')

				#Make and filter pairs
				n_list = write_pairs(combined_method, up, down, pscore, arg, of, construct_method, seq_id)

				n_pairs += n_list[0]
				n_score += n_list[1]
				n_no_G += n_list[2]
				n_ok += n_list[3]
					
			f_down.close()
		
		if 	n_pairs == 0:
			print('0 pairs found for '+seq_id)
		elif n_score == n_pairs and n_pairs != 0:
			print('All pair scores are lower than the minimum pair score allowed ('+str(pscore)+') for '+seq_id)
		elif n_no_G == n_pairs and n_pairs != 0:
			print('All pair discarded by DECKO filter for U6 promotor for '+seq_id)
		elif n_no_G + n_score == n_pairs and n_pairs != 0:
			print('All pair discarded by score or DECKO filter for '+seq_id)
	
	f_up.close()
	of.close()
	
	#with open(outfile,'a') as output:
	#	output.write(arg[0]+':\t'+arg[1]+'-'+arg[2]+'\t'+'('+arg[3]+') '+arg[4]+'\n')
	#output.close()
	#call('echo ">'+arg[0]+':\t'+arg[1]+'-'+arg[2]+'\t'+'('+arg[3]+') '+arg[4]+'" >> '+outfile, shell=True)
	rank_pairs(tmp_dir, rank_method)

	returned_pairs = pairs_filter(arguments, tmp_dir+'/gRNA_pairs', outfile, number_results, variety_results, design_up, design_down, construct_method, tmp_dir, strict)
	#call('cat '+tmp_dir+'/gRNA_pairs >> check.txt', shell=True)
	if os.path.isfile(tmp_dir+'/gRNA_pairs'): os.remove(tmp_dir+'/gRNA_pairs')
	if os.path.isfile(tmp_dir+'/'+file_list[0]): os.remove(tmp_dir+'/'+file_list[0])	
	if os.path.isfile(tmp_dir+'/'+file_list[1]): os.remove(tmp_dir+'/'+file_list[1])

	return([n_pairs, n_score, n_no_G, returned_pairs[0],returned_pairs[1]])


def write_pairs(combined_method, up, down, pscore, arg, of, method, seq_id):
	'''
	Analyze 2 gRNAs and write them into a file
	'''
	#Filter by DECKO construction U6 promotor
	if method == "DECKO":
		if up[3][0]!='G' and down[3][0]!='G':
			return[1,0,1,0]

	#Paired score method
	if combined_method == "+":
		score = round(float(up[4])+float(down[4]),3)
	elif combined_method == "x":
		score = round(float(up[4])*float(down[4]),3)
			
	#Pairs discarded by paired score			
	if score < pscore:
		return [1,1,0,0]

	#Distance between gRNAs
	if arg[3] == '+':
		dist = str(int(down[1])-int(up[2]))
	else:
		dist = str(int(up[1])-int(down[2]))

	new_line = [seq_id]+up[:5]+down[:5]+[dist]+[str(score)]+[str(int(up[-1])+int(down[-1]))]
	for element in new_line:
		of.write(element+'\t')
	of.write('\n')
	
	return [1,0,0,1]


def rank_pairs(tmp_dir, rank_method):
	'''
	Sort pairs file depending on the rank method
	'''
	if rank_method == 'score':
		with open(tmp_dir+'/gRNA_pairs') as fin:							#Create a list with the lines in the files
			lines = [line.split() for line in fin]
		lines.sort(key=itemgetter(13,12),reverse=True)									#Sort list by start position
		with open(tmp_dir+'/'+'output.txt', 'w') as fout:
			for el in lines:
				fout.write('{0}\n'.format('\t'.join(el)))					#Write sorted list into new file
		os.rename(tmp_dir+'/output.txt', tmp_dir+'/gRNA_pairs')				#Overwrite the file with old name

	elif rank_method == 'dist':
		with open(tmp_dir+'/gRNA_pairs') as fin:							#Create a list with the lines in the files
			lines = [line.split() for line in fin]
		lines.sort(key=itemgetter(13,11),reverse=True)									#Sort list by start position
		with open(tmp_dir+'/'+'output.txt', 'w') as fout:
			for el in lines:
				fout.write('{0}\n'.format('\t'.join(el)))					#Write sorted list into new file
		os.rename(tmp_dir+'/output.txt', tmp_dir+'/gRNA_pairs')				#Overwrite the file with old name

	#call('sort -k13,13 -k11,11n '+tmp_dir+'/gRNA_pairs -o '+tmp_dir+'/gRNA_pairs', shell=True)
	#call('sort -k13,13 -k12,12nr '+tmp_dir+'/gRNA_pairs -o '+tmp_dir+'/gRNA_pairs', shell=True)


def pairs_filter(arguments, input_file, output_file, number_results, variety, up_range, down_range, method, tmp_dir, strict):
	'''
	Remove those gRNA pairs that appears more times than minimum variety
	'''

	par1_dict = dict()
	par2_dict = dict()
	i = 0
	v = 0
	result_lines = []
	browser_lines = []

	fi = open(input_file, 'r')
	if os.path.exists(output_file):
		fo = open(output_file, 'a')
	else:
		fo = open(output_file, 'a')
		tab="\t"
		fo.write('Sequence_ID(#pair)'+tab+'chromosome'+tab+'start'+tab+'end'+tab+'sgRNA_1+PAM'+tab+'score_1'+tab+'chromosome'+tab+'start'+tab+'end'+tab+'sgRNA_2+PAM'+tab+'score_2'+tab+'distance_to_exclude_up_region'+tab+'distance_to_exclude_down_region'+tab+'distance_between_sgRNAs'+tab+'paired_score'+tab+'mask_score'+tab+'oligo\n')

	fbpairs = open(tmp_dir+'/browser_pairs', 'a')

	for line in fi:
		arg = line.strip().split('\t')
		seq_id = arg.pop(0)
		if arg[3] in par1_dict:
			par1_dict[arg[3]] += 1.0
		else:
			par1_dict[arg[3]] = 1.0

		if arg[8] in par2_dict:
			par2_dict[arg[8]] += 1.0
		else:
			par2_dict[arg[8]] = 1.0

		if (par1_dict[arg[3]] / number_results) > variety or (par2_dict[arg[8]] / number_results) > variety:
			v += 1
			pass
		else:

			if i >= number_results:
				fo.close()
				break
			else:

				i += 1
				if method == "DECKO":
					if arg[3][0] == 'G':
						oligo = gibson_5+arg[3][:-3]+constant_h1+arg[8][:-3]+gibson_3
					else: 
						oligo = gibson_5+arg[8][:-3]+constant_h1+arg[3][:-3]+gibson_3
				else:
					oligo = '.'

				#Distance till excluded region
				if arguments[3] == '+':
					dist_up_TSS = int(arg[2]) - (int(arguments[1]) + up_range)
					dist_down_TSS = int(arg[6]) - (int(arguments[2]) - down_range)
				else:
					dist_down_TSS = int(arg[7]) - (int(arguments[1]) + down_range)
					dist_up_TSS = int(arg[1]) - (int(arguments[2]) - up_range)

				if arguments[3] == '+':
					result_lines.append(seq_id+'('+str(i)+')'+'\t'+arg[0]+'\t'+arg[1]+'\t'+arg[2]+'\t'+arg[3]+'\t'+arg[4]+'\t'
								+arg[5]+'\t'+arg[6]+'\t'+arg[7]+'\t'+arg[8]+'\t'+arg[9]+'\t'
				 				+str(dist_up_TSS)+'\t'+str(dist_down_TSS)+'\t'
				 				+arg[10]+'\t'+arg[11]+'\t'+arg[12]+'\t'+oligo+'\n')
				else:
					result_lines.append(seq_id+'('+str(i)+')'+'\t'+arg[5]+'\t'+arg[6]+'\t'+arg[7]+'\t'+arg[8]+'\t'+arg[9]+'\t'
								+arg[0]+'\t'+arg[1]+'\t'+arg[2]+'\t'+arg[3]+'\t'+arg[4]+'\t'
								+str(dist_down_TSS)+'\t'+str(dist_up_TSS)+'\t'
								+arg[10]+'\t'+arg[11]+'\t'+arg[12]+'\t'+oligo+'\n')

				positions = [arg[1],arg[7]] if int(arg[1]) < int(arg[7]) else [arg[6],arg[2]]
				start = str(positions[0])
				end = str(positions[1])
				browser_lines.append(arg[0]+'\t'+start+'\t'+end+'\t'+arguments[4]+'('+str(i)+')'+'\t'+arg[11]+'\t'+'.'+
						'\t'+start+'\t'+end+'\t'+str(round((float(arg[11])*255.0)/2.0,0))+',0,0'+'\t'+'2'+'\t'+
						'23,23'+'\t'+'0,'+str(int(end)-int(start)-23)+'\n')				

	if strict and len(result_lines) < number_results:
		return([i,v])

	for result in result_lines:
		fo = open(output_file, 'a')
		fo.write("%s" % result)
		fo.close()

	for browser in browser_lines:
		fbpairs = open(tmp_dir+'/browser_pairs', 'a')
		fbpairs.write("%s" % browser)
		fbpairs.close()

	if i < number_results and i != 0:
		print('Only '+str(i)+' sgRNA pair(s) found for '+seq_id)
	fi.close()
	fo.close()
	fbpairs.close()
	return([i,v])


def write_browser_file(tmp_dir, of):

	
	browser_file = open(of+'_browser','a')

	browser_file.write('track name="'+of+'_Target_Regions" description=" " visibility=1 itemRgb="On"\n')
	try:
		with open(tmp_dir+'/browser_target_regions') as infile:
			browser_file.write(infile.read())
		infile.close()
	except IOError:
		return

	browser_file.write('track name="'+of+'_sgRNA_pairs" description=" " visibility=3 itemRgb="On"\n')
	try:
		with open(tmp_dir+'/browser_pairs') as infile:
			browser_file.write(infile.read())
		infile.close()
	except IOError:
		return

	browser_file.close()


def comput_means(outfile):

	fi = open(outfile, 'r')
	oldID = ''

	pairs, pairmeans = 0, []
	distances, dmeans = [], []
	iscores, imeans = [], []
	pscores, pmeans = [], []

	for line in fi:
		if line.startswith('Sequence_ID(#pair)'):
			continue
		arg  = line.strip().split('\t')
		newID = arg[0].split('(')[0]

		if  newID != oldID:
			if oldID == '':
				pairs += 1
				distances.append(int(arg[13]))
				iscores.append(float(arg[5]))
				iscores.append(float(arg[10]))
				pscores.append(float(arg[14]))
				oldID = newID	
				continue	

			pairmeans.append(pairs)
			dmeans.append(np.mean(distances))
			imeans.append(np.mean(iscores, dtype=np.float32))
			pmeans.append(np.mean(pscores, dtype=np.float32))	

			distances, iscores, pscores, pairs = [], [], [], 0
			pairs += 1
			distances.append(int(arg[13]))
			iscores.append(float(arg[5]))
			iscores.append(float(arg[10]))
			pscores.append(float(arg[14]))
		
		else:
			pairs += 1
			distances.append(int(arg[13]))
			iscores.append(float(arg[5]))
			iscores.append(float(arg[10]))
			pscores.append(float(arg[14]))		

		oldID = newID

	pairmeans.append(pairs)
	dmeans.append(np.mean(distances))
	imeans.append(np.mean(iscores, dtype=np.float32))
	pmeans.append(np.mean(pscores, dtype=np.float32))

	values = [np.mean(pairmeans), np.mean(imeans, dtype=np.float32), np.mean(pmeans, dtype=np.float32), np.mean(dmeans), pairmeans]

	return values


def write_log_results(log, end_time, dt, pair_list, grna_list, outfile, n):

	if not os.path.isfile(outfile):
		return
	
	days, seconds = dt.days, dt.seconds
	hours = days * 24 + seconds // 3600
	minutes = (seconds % 3600) // 60
	seconds = seconds % 60
	n = int(n)

	means = comput_means(outfile)

	g_up_min,g_up_mean,g_up_max=min(grna_list[0]),sum(grna_list[0])/len(grna_list[0]),max(grna_list[0])
	g_dwn_min,g_dwn_mean,g_dwn_max=min(grna_list[1]),sum(grna_list[1])/len(grna_list[1]),max(grna_list[1])
	t_min,t_mean,t_max=min(grna_list[2]),sum(grna_list[2])/len(grna_list[2]),max(grna_list[2])
	n_min,n_mean,n_max=min(grna_list[3]),sum(grna_list[3])/len(grna_list[3]),max(grna_list[3])
	pos_min,pos_mean,pos_max=min(grna_list[4]),sum(grna_list[4])/len(grna_list[4]),max(grna_list[4])
	s_min,s_mean,s_max=min(grna_list[5]),sum(grna_list[5])/len(grna_list[5]),max(grna_list[5])
	ot_min,ot_mean,ot_max=min(grna_list[6]),sum(grna_list[6])/len(grna_list[6]),max(grna_list[6])

	tp_min,tp_mean,tp_max=min(pair_list[0]),sum(pair_list[0])/len(pair_list[0]),max(pair_list[0])
	rp_min,rp_mean,rp_max=min(pair_list[3]),sum(pair_list[3])/len(pair_list[3]),max(pair_list[3])
	ps_min,ps_mean,ps_max=min(pair_list[1]),sum(pair_list[1])/len(pair_list[1]),max(pair_list[1])
	no_g_min,no_g_mean,no_g_max=min(pair_list[2]),sum(pair_list[2])/len(pair_list[2]),max(pair_list[2])
	v_min,v_mean,v_max=min(pair_list[11]),sum(pair_list[11])/len(pair_list[11]),max(pair_list[11])

	n_n = means[4].count(n)
	n_i = len(means[4])-means[4].count(n)
	n_0 = int(pair_list[4]) - n_i - n_n

	log.write('\n\nTarget regions summary:\n'
				'\tNumber of regions in input_file: {0}\n'.format(pair_list[4])+
				'\tModified regions to analyze: {0}\n'.format(pair_list[5])+
				'\tNumber of regions failed extracting sequence to analyze: {0}\n'.format(pair_list[7])+
				'\tNumber of sequences analyzed: {0}\n'.format(pair_list[6])+
				'\tNumber of complete regions (#pairs = n): {0}\n'.format(n_n)+
				'\tNumber of incomplete regions (#pairs < n): {0}\n'.format(n_i)+
				'\tNumber of regions with 0 results (#pairs = 0): {0}\n\n'.format(n_0)+
				
				'Individual sgRNAs summary:\n'
				'\t#Positive sgRNAs upstream:\t\t\t\tmin:{0} mean:{1} max:{2}\n'.format(g_up_min,g_up_mean,g_up_max)+
				'\t#Positive sgRNAs downstream:\t\t\t\tmin:{0} mean:{1} max:{2}\n'.format(g_dwn_min,g_dwn_mean,g_dwn_max)+
				'\t#sgRNAs excluded by sgRNApol III stop signal (TTTTT):\tmin:{0} mean:{1} max:{2}\n'.format(t_min,t_mean,t_max)+
				'\t#sgRNAs excluded by "N" nucleotide in sgRNA sequence:\tmin:{0} mean:{1} max:{2}\n'.format(n_min,n_mean,n_max)+
				'\t#sgRNAs excluded by score:\t\t\t\tmin:{0} mean:{1} max:{2}\n'.format(s_min,s_mean,s_max)+
				'\t#sgRNAs excluded by off-target analysis:\t\t\tmin:{0} mean:{1} max:{2}\n'.format(ot_min,ot_mean,ot_max)+
				'\t#mean of individual sgRNAs scores:\t\t\t{0}\n\n'.format(str(means[1]))+

				'Paired sgRNAs summary:\n'
				'\t#Total pairs per region:\t\t\t\tmin:{0} mean:{1} max:{2}\n'.format(tp_min,tp_mean,tp_max)+
				'\t#Returned pairs per region:\t\t\t\tmin:{0} mean:{1} max:{2}\n'.format(rp_min,rp_mean,rp_max)+
				'\t#Pairs excluded by paired score:\t\t\tmin:{0} mean:{1} max:{2}\n'.format(ps_min,ps_mean,ps_max)+
				'\t#Pairs excluded by variety:\t\t\t\tmin:{0} mean:{1} max:{2}\n'.format(v_min,v_mean,v_max)+
				'\t#Pairs excluded by DECKO construction method:\t\tmin:{0} mean:{1} max:{2}\n'.format(no_g_min,no_g_mean,no_g_max)+
				'\t#mean of distances(nt):\t\t\t\t\t'+str(means[3])+'\n'+
				'\t#mean of paired sgRNAs scores:\t\t\t\t'+str(means[2])+'\n\n\n'+

				'CRISPETA execution time: {0}h{1}m{2}s\n'.format(hours,minutes,seconds)+
				'CRISPETA finished on: '+end_time.strftime("%A, %d. %B %Y %I:%M%p"))
	

def memory_usage_resource():
    
    rusage_denom = 1024.
    if sys.platform == 'darwin':
        # ... it seems that in OSX the output is different units ...
        rusage_denom = rusage_denom * rusage_denom
    mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / rusage_denom
    return mem


def drange(start, stop, step):
	r = start    
	result = []
	while r < stop:
		result.append(r)
		r += step
	return result

def plots(file_name, results_regions):

	#if not os.path.exists(file_name+'_graphics'):
	#	os.makedirs(file_name+'_graphics')
	if not os.path.isfile(file_name):
		return 
	f = open(file_name,'r')

	pairs = 0
	regions = 0
	switch = 0
	n_pairs = []
	score1 = []
	score2 = []
	pscore = []
	dist = []
	previous_ID = ''

	for line in f:
		if line.startswith('Sequence_ID(#pair)'):
			continue
		row = line.strip().split('\t')
		current_ID = row[0].split('(')[0]

		if current_ID == previous_ID:
			pairs += 1
		else:
			if switch == 1:
				n_pairs.append(int(pairs))
			regions += 1
			pairs = 1

		switch = 1
		score1.append(float(row[5]))
		score2.append(float(row[10]))
		pscore.append(float(row[14]))
		dist.append(int(row[13]))
		previous_ID = current_ID

	n_pairs.append(int(pairs))
	iscore=score1+score2

	#create html codes for each hist
	fig = def_fig(n_pairs+[0]*results_regions[0], 'Pairs per regions', '# pairs', 'Frequency')
	html1 = plotly.offline.plot(fig, auto_open=False, show_link=False, link_text='',output_type='div')

	fig = def_fig(iscore, 'Individual sgRNA scores', 'scores [0-1]', 'Frequency')
	html2 = plotly.offline.plot(fig, auto_open=False, show_link=False, link_text='',output_type='div')

	fig = def_fig(pscore, 'Paired sgRNA scores', 'scores [0-2]', 'Frequency')
	html3 = plotly.offline.plot(fig, auto_open=False, show_link=False, link_text='',output_type='div')

	fig = def_fig(dist, 'Pair distances', 'distance (nt)', 'Frequency')
	html4 = plotly.offline.plot(fig, auto_open=False, show_link=False, link_text='',output_type='div')

	chartpie = def_chartpie(n_pairs+[0]*results_regions[0])
	html5 = plotly.offline.plot(chartpie, auto_open=False, show_link=False, link_text='',output_type='div')

	#Create html file with all images using 4 htmml codes (1 code from each hist)
	fo = open(file_name+'_images.html','w')
	fo.write('\n'.join([html1,html5,html2,html3,html4,html5]))

	#Converts to pdf
	try:
		import pdfkit
		_pdfkit=1
	except ImportError:
		_pdfkit=0
		print "\n\tpdfkit module not installed. Install it to obtain pdf with graphics."
	
	if _pdfkit==1:
		pdfkit.from_file(file_name+'_images.html',file_name+'images.pdf',options={'quiet':''})


def def_chartpie(values):

	zeros = values.count(0)
	tens = values.count(10)
	rest = len(values)-zeros-tens
	fig = {
	'data': [{	'labels': ['n=0', '0<n<10', 'n>=10'],
				'values': [zeros,rest,tens],
				'type': 'pie'}],
	'layout': {	'title': 'Number of pairs along the analyzed regions'}
	}

	return fig

def def_fig(data, title, xname, yname):

	data = [go.Histogram(x=data)]
	layout = def_layout(title, xname, yname)
	fig = go.Figure(data=data, layout=layout)

	return fig


def def_layout(title, xname, yname, xfont='Courier New, monospace', yfont='Courier New, monospace', xsize=18, ysize=18, xcolor='#7f7f7f', ycolor='#7f7f7f'):

	layout = go.Layout(	title=title,
					xaxis=dict(	title=xname,
								titlefont=dict(
									family=xfont,
									size=xsize,
									color=xcolor)
								),
					yaxis=dict(	title=yname,
								titlefont=dict(
									family=yfont,
									size=ysize,
									color=ycolor)
								))
	return layout


def get_conn(dictionary):
	conn = _mysql.connect(**dictionary)													#connects to mysql
	return conn


def create_database(conn, dbname):	
	conn.query("CREATE DATABASE IF NOT EXISTS "+dbname)									#creates database


def create_table(conn,tablename,cnames):

	table = "CREATE TABLE IF NOT EXISTS "+tablename+" ("+\
			cnames['col1']+" VARCHAR(20) NOT NULL,"+\
			cnames['col2']+" INT NOT NULL,"+\
			cnames['col3']+" INT NOT NULL,"+\
			cnames['col4']+" INT NOT NULL,"+\
			cnames['col5']+" INT NOT NULL,"+\
			cnames['col6']+" INT NOT NULL,"\
			"PRIMARY KEY (grna));"
		
	conn.query(table)


def insert_into_table(conn,infile,table):

	query = "LOAD DATA LOCAL INFILE '"+infile+\
			"' INTO TABLE "+table+\
			" FIELDS TERMINATED BY ','"\
			" LINES TERMINATED BY '\n'"

	conn.query(query)