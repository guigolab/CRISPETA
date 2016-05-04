"""
DEKO_tool
"""
from func import *

parser = argparse.ArgumentParser(description = "DECKOtool is a really flexible tool to obtain CRISPR paired gRNAs on desired regions of a genome.\
                                                Using as input a BED format file DECKOtool is able to find, analyze, and score all posible pairs.\
                                                As a result the program returns a file with {n} paired gRNAs/region,\
                                                a BED file with genome coordinates of paired gRNAs ready to be uploaded to Genome Browser\
                                                and a summary of the execution.")

parser.add_argument('-i', '--input',
                    required = True,
                    dest = 'infile',
                    action = 'store',
                    help = 'Path to input bed file with target features')

parser.add_argument('-g', '--genome',
                    dest = 'genome',
                    action = 'store',
                    required = True,
                    default = '../genomes/hg19/Genome_v19.fasta',
                    help = 'Genome fasta or multifasta file path')

parser.add_argument('-t', '--off-targets',
                    dest = 'off_targets',
                    action = 'store',
                    default = '1,0,0,x,x',
                    help = 'Maximum number of off-targets allowed per mismatch')

parser.add_argument('-o', '--output',
                    dest = 'outfile',
                    action = 'store',
                    default = 'gRNA_pairs.txt',
                    help = 'Output file path (default="gRNA_pairs.txt")')

parser.add_argument('-n', '--number',
                    dest = 'number_results',
                    action = 'store',
                    default = 10,
                    type= int,
                    help = 'Number of paired gRNAs per feature retrieved by the program')

parser.add_argument('-v', '--diversity',
                    dest = 'diversity_results',
                    action = 'store',
                    default = 0.5,
                    type = float,
                    help = 'Percentage of max variety in gRNA pairs')

parser.add_argument('-eu', '--exclude_up',
                    dest = 'exclude_up',
                    action = 'store',
                    default = 100,
                    type= int,
                    help = 'Distance upstream from target region for the program to start to look for gRNAs')

parser.add_argument('-ed', '--exclude_down',
                    dest = 'exclude_down',
                    action = 'store',
                    default = 100,
                    type= int,
                    help = 'Distance downstream from target region for the program to start to look for gRNAs')

parser.add_argument('-du', '--design_up',
                    dest = 'design_up',
                    action = 'store',
                    default = 500,
                    type= int,
                    help = 'Range upstream from the excluded region to look for gRNAs (value must be positive)')

parser.add_argument('-dd', '--design_down',
                    dest = 'design_down',
                    action = 'store',
                    default = 500,
                    type= int,
                    help = 'Range downstream from the excluded region look for gRNAs (value must be positive)')

parser.add_argument('-mp', '--positive_mask',
                    dest = 'pmask',
                    action = 'store',
                    default = None, 
                    help = 'Path to BED file with information about positive regions in genome.')

parser.add_argument('-mn', '--negative_mask',
                    dest = 'nmask',
                    action = 'store',
                    default = None, 
                    help = 'Path to BED file with information about negative regions in genome.')

parser.add_argument('-si', '--min_iscore',
                    dest = 'individual_score',
                    action = 'store',
                    default = 0.2,
                    type= float,
                    help = 'Minimum individual gRNA score allowed')

parser.add_argument('-sp', '--min_pscore',
                    dest = 'paired_score',
                    action = 'store',
                    default = 0.4,
                    type= float,
                    help = 'Minimum gRNA paired score allowed')

parser.add_argument('-sc', '--score_combination',
                    dest = 'score_combination',
                    action = 'store',
                    default = "+",
                    type = str,
                    choices=['+', 'x'],
                    help = 'Method applied to paired score: "+" (sum) or "x" (product)')

parser.add_argument('-c', '--construct_method',
                    dest = 'construct_method',
                    action = 'store',
                    default = None,
                    type = str,
                    choices=[None,'DECKO'],
                    help = 'Method applied when making gRNA pairs and oligo construction')

parser.add_argument('-r', '--rank_method',
                    dest = 'rank',
                    action = 'store',
                    default = 'score',
                    type = str,
                    choices=['score','dist'],
                    help = 'Method applied when ranking gRNA pairs [score/dist]')

parser.add_argument('-str', '--strict',
                    dest = 'strict',
                    default = False,
                    help = 'Searching windows will not be increased if 0 gRNAs are found in a regions')


#Storing user options
options = parser.parse_args()

#Checking user defined parameters
check_user_options(options)

#starting variables
complete_regions,incomplete_regions, regions_0 = 0,0,0
n_pairs, n_paired_score, n_no_G, returned_pairs, n_variety_pairs = [], [], [], [], []
p_grna_up, p_grna_dwn, polIII, n_in_seq, bad_position, bad_score, ot = [],[],[],[],[],[],[]
starting_time=datetime.datetime.now()

#Creating a temporary directory for temporary files and opening log file
tmp_dir = tempfile.mkdtemp(prefix='gRNA')
log = open(options.outfile+'.log', 'a')
log.write("CRISPETA started on: "+starting_time.strftime("%A, %d. %B %Y %I:%M%p\n")+'Options:\n')
for k in vars(options).keys():
     log.write('\t-'+k+': '+str(vars(options)[k])+'\n')

#Transforming input bed file with new positions exclude/design up/down
print('Getting new BED regions...')
n_regions_ids = target_bed(options.infile,
                            options.design_up, options.exclude_up,
                            options.exclude_down, options.design_down,
                            tmp_dir)

total_regions = n_regions_ids[0]    #Total number of regions in input BED file
n_regions = n_regions_ids[1]        #Number of modified regions
ids = n_regions_ids[2]              #Id's of modified regions

#Extracting sequences using modified bed file regions
print('Getting sequences from new regions... '+'('+str(n_regions)+' regions)')
n_sequences = get_sequences(options.genome, 'tmp_bed_file', 'tmp_fasta_file', tmp_dir)/2

conflictive_regions = total_regions - n_sequences   #Those regions that bedtools could not retrieve a sequence

#Exploring the sequences
up_n_gRNA = []              #Number of gRNAs upstream per regions
down_n_gRNA = []            #Number of gRNAs downstream per regions
error_regions = 0           #Number of regions without gRNAs
i=-1
switch = 1                  #Switcher for mask analysis
print('Searching sgRNAs in sequences.. ('+str(n_sequences)+' sequences)')

sequences_file = open(tmp_dir+"/tmp_fasta_file")

conn = use_database()
for line in sequences_file:
     if line.startswith(">"):    #>chr11:65266175-65266683(+)
          i+=1.0
          strand = line[-3]
          seq_id = ids[line[:-4]]
          arguments = line[1:-4].replace(":",",").replace("-",",",1).strip().split(",") + [strand] + [seq_id] #Region features [chr, start, end, strand, ID]

     else:
          #Sequence and reverse sequence to analyze
          seq = Seq(line.strip(), generic_dna)
          reverse_seq = seq.reverse_complement()
        
          if arguments[3] == '-':     #Reverse up/down stream values if strand = "-"
               design_up, design_down = options.design_down, options.design_up
               exclude_up, exclude_down = options.exclude_down, options.exclude_up
          else:
               design_up, design_down = options.design_up, options.design_down
               exclude_up, exclude_down = options.exclude_up, options.exclude_down

          #Search for all posible gRNA patterns and creating a bed file with info for each match
          n = get_gRNAs(arguments, seq, reverse_seq,
                         design_up, design_down,
                         exclude_up, exclude_down,                         
                         strand, options.individual_score,
                         tmp_dir, options.genome,
                         options.off_targets, conn)

          #increase windows for gRNA search if [-strict] == True and pairs < n; max 3 times
          counter=0
          new_design_up, new_design_down = design_up, design_down

          if n[0]*n[1] < 1:                                      #No gRNAs to make pairs
               start = str(int(arguments[2])-new_design_up)
               end = str(int(arguments[2])+new_design_down)
               regions_0 += 1
               print('0 sgRNA pair(s) found for '+seq_id)
               continue

          write_browser_target_regions(tmp_dir, 'browser_target_regions', arguments, new_design_up, exclude_up, new_design_down, exclude_down)


          #Analyze the gRNAs and compare their positions to masks positions
          for element_name in ['upstream_bed', 'downstream_bed']:
               get_masks(arguments, options.pmask, options.nmask, options.individual_score, tmp_dir, element_name)

          statistics = make_pairs(arguments, ['upstream_bed', 'downstream_bed'], options.number_results,
                              tmp_dir, options.outfile, options.diversity_results,
                              options.paired_score, design_up, design_down,
                              options.score_combination, options.construct_method, options.rank, seq_id, options.strict)

          if statistics[0] == 0:                                 #0 pairs returned
               regions_0 += 1
               print('0 sgRNA pair(s) found for '+seq_id)
               continue
          elif statistics[3] < options.number_results:           #pairs returned < n
               if options.strict:
                    regions_0 += 1
                    continue
               else:
                    incomplete_regions += 1
          else:
               complete_regions += 1


          n_pairs.append(statistics[0])
          n_paired_score.append(statistics[1])
          n_no_G.append(statistics[2])
          returned_pairs.append(statistics[3])
          n_variety_pairs.append(statistics[4])
          p_grna_up.append(n[0])
          p_grna_dwn.append(n[1])
          polIII.append(n[2][0])
          n_in_seq.append(n[2][1])
          bad_position.append(n[2][2])
          bad_score.append(n[2][3])
          ot.append(n[2][4])


#Plot histogram with summary results

results_regions = [regions_0,incomplete_regions, complete_regions]
write_browser_file(tmp_dir, options.outfile)
sequences_file.close()

for file_name in listdir(tmp_dir):
     os.remove(tmp_dir+'/'+file_name)
os.rmdir(tmp_dir)

contig_pairs_list = [n_pairs, n_paired_score, n_no_G, returned_pairs,
                     total_regions, n_regions, n_sequences, conflictive_regions,
                     incomplete_regions, regions_0, complete_regions, n_variety_pairs]
contig_grnas_list = [p_grna_up, p_grna_dwn, polIII, n_in_seq, bad_position, bad_score, ot]
end_time=datetime.datetime.now()
dt = end_time-starting_time

print('Writing log of CRISPETA designs')
write_log_results(log, end_time, dt, contig_pairs_list, contig_grnas_list, options.outfile, options.number_results)

print('Making plots')
try:
     import plotly
     _plotly=1
except ImportError:
     _plotly=0
     print "\n\tplotly module not installed. Install it to obtain graphics."

if _plotly==1:
     plots(options.outfile, results_regions)

     
print('CRISPETA finished without errors =D!')