#! /usr/local/bin/python
# Python 2.7.5, requires Biopython.

# 1. Opens .xml BLAST result file
# 2. filters out records based on e-value
# 3. order records based on e-value
# 4. Take top 50 unique gene id numbers.
# 5. Query GenBank for records that match gene ids.
# 6. Output records, including sequences, as a multifasta file.

# Version 1. Darcy Jones, November 2013.


### Import modules
import sys; #for command line argument handling
import argparse; 
import xml.etree.ElementTree as etree; #for generic xml handling
from Bio import SeqIO; #For .fasta output
from Bio import Entrez; #For .fasta output
import time; #For output file name
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
	return query_id[0];

### Argument handling
arg_parser = argparse.ArgumentParser(description='Takes the top x number of unique gids for each query from xml, outputs non-redundant .fas file containing sequences. Give command as $ python quickOrtho.py openFile.xml database yourEmail@address.com -n integer -e number -o outputFile.fas');
arg_parser.add_argument("file_dir", help="Directory to NCBI .xml file.");
arg_parser.add_argument("database", help="Which database to direct entrez query to.", choices=['protein','nucleotide']);
arg_parser.add_argument("email", help="Email for entrez record retrieval, tells NCBI who you are.");
arg_parser.add_argument("-n", "--number_unique_gids", type=int, default=50, help="Number of unique gids to extract for each query");
arg_parser.add_argument("-e", "--e_value_threshold", type=float, default=1e-20, help="Maximum e-value allowed in screening, enter as decimal or in scientific notation (eg. 1e-20)");
arg_parser.add_argument("-o", "--file_dir_output", default=arg_parser.parse_args().file_dir[:-4]+"_"+time.strftime("%d%m%y")+"_records.fas", help="Directory/name of output file");

args = arg_parser.parse_args();

### Variable definitions/declarations
file_dir = args.file_dir; #directory of .xml file to be parsed.
number_unique_gids = args.number_unique_gids; #number of unique gids to be fetched for each query id.
e_value_threshold = args.e_value_threshold; #E-value threshold.
gene_dict = {}; #dictionary keyed by query id, with value = tuple (gid, evalue) that match evalue requirement
gene_list_master = []; #list of x ordered unique gids for entrez record retrieval
Entrez.email = args.email; #Required for request 
temp_hit_set = set(); #temporary set for handling duplicated results between query ids/dictionary keys
file_dir_output = args.file_dir_output;
database = args.database;

### Code
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

for key in gene_dict:
	seen = set(); #temporary set for handling duplications within query id lists
	gene_dict[key] = [x for x in gene_dict[key] if x not in seen and not seen.add(x)] #Removes duplicates from the list of tuples for each query id
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

print "Fetching "+`len(gene_list_master)`+" records from NCBI."

records_written=0; #counts how many files retrieved successfully 
with open(file_dir_output, 'w') as file_output: #opens output fasta file, can also be append mode if required?
	for gid in gene_list_master:
		try:
			Entrez_handle = Entrez.efetch(db=database, rettype="fasta", retmode="text", id = gid); 
			Entrez_record = SeqIO.read(Entrez_handle, "fasta");
			Entrez_handle.close();
			SeqIO.write(Entrez_record, file_output, "fasta");
			records_written+=1
		except: #Catches all exceptions
			print "Error retrieving or writing "+repr(gid)+", please check fasta file or try again."
			print `sys.exc_info()`+"\n";

print "Successfully wrote "+repr(records_written)+" sequences to file: "+ file_dir_output +"."
