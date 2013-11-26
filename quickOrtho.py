#! /usr/bin/python
# Python 2.7.5, requires Biopython.

# 1. Opens .xml BLAST result file
# 2. filters out records based on e-value
# 3. order records based on e-value
# 4. Take top 50 unique gene id numbers.
# 5. Query GenBank for records that match gene ids.
# 6. Output records, including sequences, as a multifasta file.

# Version 1. Darcy Jones, November 2013.


### Import modules
import argparse; #for command line flags and input.
import sys;
import xml.etree.ElementTree as etree; #for generic xml handling
from Bio import SeqIO; #For .fasta output
from Bio import Entrez; #For .fasta output
import re; #used for extracting gi numbers from complex strings

### Functions

def gi_extract(hit_object): #takes <hit> object, looks for gid in <Hit_id> then <Hit_def>, returns gid
	regex = re.compile("^gi\|(\w+)");
	hit_id = regex.findall(hit_object.find("Hit_id").text.strip());
	if hit_id == []:
		hit_id = regex.findall(hit_object.find("Hit_def").text.strip());
	return hit_id[0];

def iteration_gi_extract(iteration_object): #takes iteration number, finds gid for query id in either <Iteration_query-def> or <Iteration_query-ID>
	regex = re.compile("^gi\|(\w+)");
	query_id = regex.findall(iteration_object.find("Iteration_query-ID").text.strip());
	if query_id == []:
		query_id = regex.findall(iteration_object.find("Iteration_query-def").text.strip());
		if query_id == []: #handles the event that the query has no gid, ie it may be an unpublished sequence.
			query_id = [iteration_object.find("Iteration_query-def").text.replace(" ", "_")]; #note that this is as a list because the regex.findall() terms give lists as output, it took fewer steps change this.
	return query_id[0];

def query_hit_table(dictionary): #
	query_heading = [x for x in dictionary];
	table='';
	for query in query_heading:
		if len(dictionary[query])!=0:
			if len(query)<=23:
				row = '{:>26}'.format(str(query+" ||"))
			else:
				row = '{:>26}'.format(str(query[:20]+"... ||"));
			i=0;
			for hit in dictionary[query]:
				row+=''.join(' gi|{:<10} e|{:<10e};'.format(hit[0], hit[1]));
				i+=1;
				if i>= number_unique_gids:
					break;
			table+= row+"\n";
	return table;
def write_sequences_fast(): #faster but more memory intensive output system, less appropriate for retrieving large sequences. Stores all records retrieved from entrez in a single object then writes either to stdout or file.
	try:
		Entrez_handle = Entrez.efetch(db=database, rettype="fasta", retmode="text", id = gene_list_master); 
		if file_dir_output == sys.stdout:
			Entrez_record = "".join(Entrez_handle);
			sys.stdout.write(Entrez_record);
		else:
			Entrez_record = SeqIO.parse(Entrez_handle, "fasta");
			SeqIO.write(Entrez_record, file_output, "fasta");
		Entrez_handle.close();
	except: #Catches all exceptions
		print("Error retrieving or writing from Entrez, please check fasta file/xml and try again.");
		print("Hint: Check that gi|'number' is present in xml file");
		print(repr(sys.exc_info()));
		raise; #If you want the loop to continue running after errors, comment out the raise (this line).
	if not quiet:
		print("Successfully wrote sequences to file: {}.".format(file_dir_output));


def write_sequences_slow(): #Slower output system, using multiple entrez queries and progressive writing to file or stdout. Should be less memory instensive than write_sequences_fast()
	records_written=0;
	for gid in gene_list_master:
		try:
			Entrez_handle = Entrez.efetch(db=database, rettype="fasta", retmode="text", id = gid); 
			if file_dir_output == sys.stdout:
				Entrez_record = "".join(Entrez_handle);
				sys.stdout.write(Entrez_record);
			else:
				Entrez_record = SeqIO.read(Entrez_handle, "fasta");
				SeqIO.write(Entrez_record, file_output, "fasta");
			Entrez_handle.close();
			records_written+=1;
			if not quiet: 
				print("Retrieving record {} of {}. gi|{}.".format(records_written, len(gene_list_master), gid));
		except: #Catches all exceptions
			print("Error retrieving or writing {}, please check fasta file/xml and try again.".format(repr(gid)));
			print("Hint: Check that gi|'number' is present in xml file");
			print(repr(sys.exc_info()));
			raise; #If you want the loop to continue running after errors, comment out the raise (this line).	
	if not quiet:
		print("Successfully wrote {} sequences to file: {}.".format(records_written, file_dir_output));

### Argument handling
arg_parser = argparse.ArgumentParser(description='Takes the top x number of unique gids for each query from xml, outputs non-redundant .fas file containing sequences. Give command as $ python quickOrtho.py openFile.xml database yourEmail@address.com -n integer -e number -o outputFile.fas');
arg_parser.add_argument("file_dir", help="Directory to NCBI .xml file. To use stdin (for piped input) enter '-'");
arg_parser.add_argument("database", help="Which database to direct entrez query to.", choices=['protein','nucleotide']);
arg_parser.add_argument("email", help="Email for entrez record retrieval, tells NCBI who you are.");
arg_parser.add_argument("-n", "--number_unique_gids", type=int, default=50, help="Number of unique gids to extract for each query");
arg_parser.add_argument("-e", "--e_value_threshold", type=float, default=1e-20, help="Maximum e-value allowed in screening, enter as decimal or in scientific notation (eg. 1e-20)");
arg_parser.add_argument("-o", "--file_dir_output", default=None, help="Directory/name of output file. To use stdout (for piped output) enter '-'. NB, quiet operation will automatically be enabled for stdout output");
arg_parser.add_argument("-t", "--table", action="store_true", help="Toggles option to create another .txt file showing a table with the top n hits and their evalues for each query. ie. which hits came from which queries. Nb best viewed without text wrapping.")
arg_parser.add_argument("-q", "--quiet", action="store_true", help="Toggles option, to run without printing running feedback");
args = arg_parser.parse_args();

### Variable definitions/declarations
file_dir = args.file_dir; #directory of .xml file to be parsed.
if file_dir=='-':
	file_dir = sys.stdin;

file_dir_output = args.file_dir_output;
if file_dir_output is None:
	if args.file_dir =='-':
		file_dir_output = "stdin_quickOrthoResults.fas";
	else:
		file_dir_output = file_dir.split('.xml')[0]+"_quickOrthoResults.fas";
elif file_dir_output=='-':
	file_dir_output = sys.stdout;

number_unique_gids = args.number_unique_gids; #number of unique gids to be fetched for each query id.
e_value_threshold = args.e_value_threshold; #E-value threshold.
gene_dict = {}; #dictionary keyed by query id, with value = list of tuples (gid, evalue) that match evalue requirement
gene_list_master = []; #list of x ordered unique gids for entrez record retrieval
Entrez.email = args.email; #Required for request 
temp_hit_set = set(); #temporary set for handling duplicated results between query ids/dictionary keys
database = args.database;
quiet = args.quiet; #Boolean toggle for realtime feedback. Default is not quiet, ie verbose, ie print realtime feedback.
if file_dir_output==sys.stdout:
	quiet = 'true'
table_output=args.table; #Boolean tobble for outputting summary table for hits to a different .txt file.
table_dir_output = str(file_dir_output).split('.fas')[0]+'_summaryTable.txt';

### Code
if not quiet:
	print("#################### Begin quickOrtho ####################");
	print("");
	print("Extracting records from {}.".format(file_dir));

xml_root = etree.parse(file_dir).getroot();
xml_iterations = xml_root.find("BlastOutput_iterations").findall("Iteration"); #list of xml elements <iteration>
for iteration in xml_iterations: #Loops through each alignment query.
	xml_iterations_queryID = iteration_gi_extract(iteration); 
	gene_list = [];
	xml_hits = iteration.find("Iteration_hits").findall("Hit"); #List of xml hits for each iteration.
	if xml_iterations_queryID not in gene_dict: #Prevents looping through alignments already done.
		for hit in xml_hits: #Loops through alignment hits in xml.
			hit_id = gi_extract(hit);
			hit_evalue = float(hit.find("Hit_hsps").find("Hsp").findall("Hsp_evalue")[0].text);
			if hit_evalue <= e_value_threshold and hit_id != xml_iterations_queryID: #adds hits e_value less than required to list, excludes hits already in list.
				gene_list.append((hit_id, hit_evalue));
		gene_dict[xml_iterations_queryID] = gene_list;

if not quiet:
	print("Sorting extracted results and retrieving top {} for each query.".format(str(number_unique_gids)));

for key in gene_dict:
	seen = set(); #temporary set for handling duplications within query id lists
	gene_dict[key] = [x for x in gene_dict[key] if x not in seen and not seen.add(x)]; #Removes duplicates from the list of tuples for each query id
	gene_dict[key].sort(key=lambda tup: tup[1]); # Sort in place by evalue (tup[1]), small to big
	if len(gene_dict[key]) > number_unique_gids: #Selects the requested number of lowest e-values from gene_list then outputs to gene_list_master
		i=0;
		while (i<number_unique_gids):
			if gene_dict[key][i][0] not in temp_hit_set: #prevents duplicate gi instances from being added to gene_list_master 
				gene_list_master.append(gene_dict[key][i][0]);
				temp_hit_set.add(gene_dict[key][i][0]);
			i+=1;
	else: # Handler for event that there are fewer list items than requested
		for tu in gene_dict[key]: # where tu = the tuple (gid, evalue)
			if tu[0] not in temp_hit_set: #prevents duplicate gi instances from being added to gene_list_master 
				gene_list_master.append(tu[0]);
				temp_hit_set.add(tu[0]);
	#if key not in temp_hit_set: #adds query id to gene_list_master
	#	gene_list_master.append(key);
	#	temp_hit_set.add(key);
temp_hit_set.clear();

if not quiet:
	print("Found {} queries in xml, containing total {} unique (to query) gid's with e-value < {}.".format(len(gene_dict), sum([len(gene_dict[key]) for key in gene_dict]), str(e_value_threshold)));
	print("Fetching {} unique records from NCBI.".format(len(gene_list_master)));

if file_dir_output==sys.stdout:
	write_sequences_fast();
else:
	with open(file_dir_output, 'w') as file_output: #opens output fasta file, can also be append mode if required?
		write_sequences_fast();

if table_output: #if table is flagged to be output. 
	with open(table_dir_output, 'w') as file_output:
		table=query_hit_table(gene_dict);
		file_output.write(table);
	if not quiet:
		print("Successfully wrote hit summary table to file: {}.".format(table_dir_output));

if not quiet:
	print("");
	print("#################### End quickOrtho ####################");

