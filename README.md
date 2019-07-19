#**About CRISPETa**
		                            
CRISPETa is a flexible tool to design optimal pairs of sgRNAs for deletion of desired genomic regions. Using as input a BED format file CRISPETa is able to find, analyze, and score all possible sgRNAs. As a result the program returns:

1. A file with information of n ranked pairs of sgRNAs for every target region (n, the desired number of targets, can be selected by the user): position in the genome, sequence of the sgRNA+PAM, individual and paired scores and distance between sgRNAs.
2. A log file with summary of run settings: number of analyzed regions, mean of sgRNA scores (individual and paired), sgRNAs filtered out for each of the filters, etc.
3. A BED file of designed sgRNAs ready to be uploaded to UCSC Genome Browser for visualization of target regions and sgRNA pairs related to them.
4. PDF with graphics based on results: histogram with pairs per target region, histograms of individual and paired scores distribution and pie chart with distributions of complete and incomplete designs.
5. HTML with same graphics as the PDF above.

The code can be found on github: https://github.com/guigolab/CRISPETA or on our web server:
http://crispeta.crg.eu


##**Requirements**

* python 2.7
    * Numpy
    * BioPython
    * python-mysqldb
    * Plotly (optional)
    * pdfkit (optional)
* BEDtools
* MySQL (tested on v5.1 and v5.5)
	

##**Before Starting**

**WARNING: CRISPETa uses MySQL. Options to connect to MySQL database can be found in config.py file. Change parameters if necessary in order to connect to your MySQL season. If you change options names while creating the database remember to change theese values in config.py file**

Before running CRISPETa the user must create a database in MySQL to store off-target information for sgRNAs. This step can take a while depending on the size of the database and computer resources (more than 1 hour for human). Files with precomputed off-target information for some organisms can be directly download from our web server (http://crispeta.crg.eu/download).

####**Create database**

The off-target database can be created:

1. Using module crispeta_mysql.py
2. Manually using MySQL.

####1. Using crispeta_mysql module:

crispeta_mysql can use comma separated files downloaded from the web site to create the data base. By default crispeta_mysql will use the following values to create and access to database:	

* user name: "crispeta"
* password: "crispeta"
* host: "localhost"
* database name: "crispeta"
* table name: "crispeta"
* column names: "gnra", "off0", "off1", "off2", "off3", "off4"

The user can modify theese parameters, except column names, directly from the command line. If you change theese values remember to change them also in the config.py file.

Eexample:

	$ python crispeta_mysql -i [coma_separated_file.txt] -u [user_name] -p [pwd]

Make sure that all required files (crispeta_mysql.py, config.py and func.py) are in the same directory

####2. Manually using MySQL:

The comma separated file can be loaded directly to MySQL from the terminal using the following commands:
		
	mysql> CREATE DATABASE crispeta;
	mysql> USE crispeta;
	mysql> CREATE TABLE [genome_name] (
	->	grna VARCHAR(20) NOT NULL,
	->	off0 INT NOT NULL,
	->	off1 INT NOT NULL,
	->	off2 INT NOT NULL,
	->	off3 INT NOT NULL,
	->	off4 INT NOT NULL,
	->	PRIMARY KEY (grna));
	mysql>	LOAD DATA LOCAL INFILE '[coma_separated_file.txt]' INTO TABLE  [genome_name]
	-> FIELDS TERMINATED BY ',' LINES TERMINATED BY '\n';
	

##**CRISPETA**
**Running example**


	$ python CRISPETA.py -i [file.bed] -g [genome.fasta] -o results.txt 

Make sure that all program files (CRISPETA.py, config.py and func.py) are in the same directory
	
**Parameters** (Also see Table 1 of the CRISPETa manuscript for details)

"dir" = directory e.g: path/to/file; "int" = integer e.g: 10;  "float" = decimal number e.g: 0.5;  "string" = text e.g: this is a string;  "bool" = boolean e.g: True or False
  
- **-i dir:** Path to input BED file.
- **-g dir:** Path to genome fasta file
- **-t string (default: 1,0,0,x,x):** Off-targets: String with maximum number of off-targets allowed with 0,1,2,3 and 4 mismatches (x: no limit). Text must have five integers (or x), comma-separated, with no spaces.
- **-o dir (default: ./sgRNA_pairs):** Path/prefix of output files
- **-n int (default: 10):** Maximum number of pairs to be returned for each target.
- **-du int (default: 500):** Upstream design region: Length of upstream region for sgRNAs search
- **-dd int (default: 500):** Downstream design region: Length of downstream region for sgRNAs search
- **-eu int (default: 100):** Exclude upstream: Length of upstream region adjacent to target excluded from sgRNAs search
- **-ed int (default: 100):** Exclude downstream: Length of downstream region adjacent to target excluded from sgRNAs search
- **-v float (default: 0.5):** Diversity measure: the maximum fraction of returned pairs for each target that contain the same sgRNA
- **-si float (default: 0.2):** Individual score cutoff: The minimum score individual sgRNAs must have to be considered
- **-sp float (default: 0.4):** Paired score cutoff: The minimum combined score that a sgRNAs pair must have to be considered
- **-sc string (default: +):** Score combination: Method by which individual scores are combined to yield pair score: addition ("sum") or multiplied ("product")
- **-r string (default: score):** Ranking method: Criteria for ranking protospacer pairs ("score"or "dist")
- **-c string (default: None):** Construction method: Method applied when making sgRNAs pairs and oligo construction: “none” or “DECKO” (only returns pairs where first protospacer starts with G)
- **-mp dir:** Positive mask: File with favoured regions from genome, in BED format
- **-mn dir:** Negative mask: File with disfavoured regions from genome, in BED format

**Input** 

The input data should be specified using a tab separated file (BED format) and passing it to the pipeline command with the option -i. Here is an example of the file format:

	chr1	29557453	29557454	region1	0	-
	chr2	32716839	32716840	region2	0	+
	chr7	151129139	151129140	region3	0	+
	chr12	151138423	151138424	region4	0	-

The fields in the file correspond to:

1. Chromosome name
2. Start position
3. End position
4. Unique ID
5. Score (irrelevant here)
6. Strand
    
**Outputs**

1. DESIGN file: sgRNA pairs found by CRISPETa for each target regions. Each line  correspond to one pair

		Sequence_ID(#pair)	chromosome	start	end	sgRNA_1+PAM	score_1	chromosome	start	end	sgRNA_2+PAM	score_2	distance_to_exclude_up_region	distance_to_exclude_down_region	distance_between_sgRNAs	paired_score	mask_score	oligo
		region1(1)	chr1	29557261	29557284	GCTTGTCTATGGGCACCACGGGG	0.944	chr1	29557844	29557867	CGTGTACTCTCCTCAGTGTAGGG	0.572	-69	290	560	1.516	2	.
		region1(2)	chr1	29557261	29557284	GCTTGTCTATGGGCACCACGGGG	0.944	chr1	29557805	29557828	CCTATGCCGTTACATGGTAGTGG	0.542	-69	251	521	1.486	2	.
		region1(3)	chr1	29557261	29557284	GCTTGTCTATGGGCACCACGGGG	0.944	chr1	29557576	29557599	GACTGCGTGTGGGCCCCGGAGGG	0.501	-69	22	292	1.445	2	.

2. DESIGN BED file: Pairs and target regions in BED format ready to upload to UCSC GenomeBrowser as a custom trak. The new tracks will show target regions and pairs for each region passed as input for CRISPETa

		track name="Target_Regions" description="Regions_for_sgRNA_searching" visibility=1 itemRgb="On"
		chr1	29557453	29557454	region1	0	+	0	0	0,0,0
		chr1	29557353	29557554	region1	0	+	0	0	180,180,130
		chr1	29556853	29558054	region1	0	+	0	0	40,175,40
		...
		track name="sgRNA_pairs" description="sgRNA_pairs" visibility=3 itemRgb="On"
		chr1	29557261	29557867	region1(1)	1.516	.	29557261	29557867	193.0,0,0	2	23,23	0,583
		chr1	29557261	29557828	region1(2)	1.486	.	29557261	29557828	189.0,0,0	2	23,23	0,544
		chr1	29557261	29557599	region1(3)	1.445	.	29557261	29557599	184.0,0,0	2	23,23	0,315
		...

3. DESIGN Settings & Statistics: A summary of the design performance. It contains the number of regions analyzed, mean of individual and paired scores, mean of pair distances, number of sgRNAs excluded by filters, etc.

4. DESIGN Plots: Some graphics are plotted in html and pdf format using designs information.


##CRISPETA Data

    CRISPETA.py -> Main script.
    crispeta_mysql.py -> Module to load off-target information to MySQL database.
    config.py -> Options and values for MySQL configuration.
    func.py -> necessary functions for CRISPETA.py and crispeta_mysql.py to work.
    README.md -> Markdown file with information about CRISPETA.

