#**About CRISPETA**
		                            
CRISPETA is a really flexible tool to obtain CRISPR paired sgRNAs on desired regions of a genome. Using as input a BED format file CRISPETA is able to find, analyze, and score all posible sgRNAs. As a result the program returns:

1. A file with information of {n} paired sgRNAs per region: position in the genome, sequence of the sgRNA+PAM, individual and paired scores, distance between sgRNAs and oligo construction if DECKO method is selected
2. A log file with summary of the execution. Contains statistics about the execution: number of analyzed regions, mean of sgRNA scores (individual and paired), sgRNAs filtered out for each of the filters, etc.
3. A BED file with genome coordinates of sgRNAs ready to be uploaded to Genome Browser for visualization of the target regions and sgRNA pairs related to them.
4. PDF with graphics based on summary numbers: histogram with pairs per target region, pie chart with distributions of regions with n=0,n=10 or 0<n<10 , etc.
5. HTML with same graphics than the PDF above to open in browser and interact with them.

The code for the tool can be found on github: https://github.com/guigolab/CRISPETA
Or on our web server: http://crispeta.crg.eu


##**Requirements**

* python 2.7
    * Numpy
    * BioPython
    * python-mysqldb
    * Plotly (optional)
    * pdfkit (optional)
* BEDtools
* MySQL (tested on v5.1 and v5.5)
	

##**Previous Steps**

**WARNING: CRISPETA uses MySQL. Options to connect to MySQL database can be found in config.py file. Change parameters if necessary in order to connect to your MySQL season. If you change options names while creating the database remember to change this values on config.py file in order to allow CRISPETA.py to run without errors**

Before running CRISPETA the user must create a database in MySQL to store off-target information for sgRNAs in target genome (this step can take a while depending on the size of the database and computer resources. More than 1 hour for human off-targets). Files with precomputed off-target information for some organisms can be directly download from our web server (http://crispeta.crg.eu).

####**Create database**

This database can be created:

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

The user can modify this parameters, except column names, directly from the command line. If you change this values remember to change them also in the config.py file.

Running example:

	$ python crispeta_mysql -i [coma_separated_file.txt] -u [user_name] -p [pwd]

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
	mysql>	LOAD DATA LOCAL INFILE '<dir>/offtarget_analysis.txt' INTO TABLE  [genome_name]
	-> FEILDS TERMINATED BY ',' LINES TERMINATED BY '\n';
	
##**CRISPETA **
**Running example**

	$ python CRISPETA.py -i [file.bed] -g [genome.fasta] -o results.txt 
	
**Parameters**

"dir" = directory e.g: path/to/file; "int" = integer e.g: 10;  "float" = decimal number e.g: 0.5;  "string" = text e.g: this is a string;  "bool" = boolean e.g: True or False
  
- **-i dir:** Path to input BED file.
- **-g dir:** Path to genome fasta file
- **-o dir (default: ./sgRNA_pairs):** Path where outputs will be placed
- **-t string (default: 1,0,0,x,x):** String without spaces containing the maximum number of off-targets allowed with 0, 1, 2, 3 and 4 mismatches separated by coma. "x" means infinite number of off-targets
- **-n int (default: 10):** Number of results per regions
- **-v float (default: 0.5):** Percentage of maximum variety in sgRNA pairs 
- **-eu int (default: 100):** Distance upstream from target region for the program to start to look for sgRNAs
- **-ed int (default: 100):** Distance downstream from target region for the program to start to look for sgRNAs
- **-du int (default: 500):** Range upstream from the excluded region to look for sgRNAs (value must be positive)
- **-dd int (default: 500):** Range downstream from the excluded region to look for sgRNAs (value must be positive)
- **-mp dir:** Path to BED file with positive regions. sgRNAs overlapping positive regions will be favored when ranking result pairs
- **-mn dir:** Path to BED file with negative regions. sgRNAs overlapping negative regions will be disfavored when ranking result pairs
- **-si float (default: 0.2):** Minimum individual sgRNA score allowed
- **-sp float (default: 0.4):** Minimum sgRNA paired score allowed
- **-sc string (default: +):** Method applied to paired score: "+" (sum) or "x" (product) of individual sgRNA scores. Multiplication method will favor those sgRNA pairs where both sgRNAs have high scores e.g: sum 0.8+0.2=1, product 0.8*0.2=0.16.
- **-c string (default: None):** Method applied when making sgRNA pairs: "None" or "DECKO". DECKO method will create pairs where at least one of the sgRNAs starts with G.
- **-r string (default: score):** Key value for sorting results [score/dist]

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
5. Score
6. Strand
    
**Outputs**

1. Main output: Pairs found by CRISPETA for target regions in Input file

		Sequence_ID(#pair)	chromosome	start	end	sgRNA_1+PAM	score_1	chromosome	start	end	sgRNA_2+PAM	score_2	distance_to_exclude_up_region	distance_to_exclude_down_region	distance_between_sgRNAs	paired_score	mask_score	oligo
		region1(1)	chr1	29557261	29557284	GCTTGTCTATGGGCACCACGGGG	0.944	chr1	29557844	29557867	CGTGTACTCTCCTCAGTGTAGGG	0.572	-69	290	560	1.516	2	.
		region1(2)	chr1	29557261	29557284	GCTTGTCTATGGGCACCACGGGG	0.944	chr1	29557805	29557828	CCTATGCCGTTACATGGTAGTGG	0.542	-69	251	521	1.486	2	.
		region1(3)	chr1	29557261	29557284	GCTTGTCTATGGGCACCACGGGG	0.944	chr1	29557576	29557599	GACTGCGTGTGGGCCCCGGAGGG	0.501	-69	22	292	1.445	2	.

2. GenomeBrowser file: Pairs and target regions in BED format ready to upload in GenomeBrowser as a custom traks

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

3. Execution summary: A brief explanation of execution is returned in txt format. Also some graphics are plotted in html and pdf format using summary results.


##CRISPETA Data

    CRISPETA.py -> Main script with pipeline of crispeta.
    crispeta_mysql.py -> module to load off-target information to MySQL in a easy way.
    config.py -> Options and values for MySQL configuration.
    func.py -> necessary functions for CRISPETA.py and crispeta_mysql.py to work.
    README.md -> Marckdown file with information about CRISPETA.
