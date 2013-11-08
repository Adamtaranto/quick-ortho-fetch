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
import xml.etree.ElementTree as etree; #for generic xml handling
from Bio import SeqIO; #For .fasta output
from Bio import Entrez; #For .fasta output
import time; #For output file name
import re;




### Variable definitions/declarations
e_value_threshold = 1e-20;
number_unique_gids = 100;
queryID = ""; #gid of query being evaluated
file_dir_output = ""; #name of multifasta file to be created.
gene_list = []; #list of tuples (gid, evalue) that match evalue requirement
temp_iteration_set = set();
temp_hit_set = set(); #temporary set for handling duplicated results, can also add query gid to this so that we avoid self hits.
gene_list_master = []; #list of x ordered unique gids
Entrez.email = 'da8jones@students.latrobe.edu.au'; #Required for request 
database = "protein"; # nb conditional handler later using blast type in xml to decide whether nucleotide db is better.

### Code
file_dir=sys.argv[1];
xml_root = etree.parse(file_dir).getroot();
queryID = sys.argv[1][:-4];

def gi_extract(gi_string):
	regex = re.compile("^(gi\|\w+)");
	return regex.findall(gi_string.strip());

xml_iterations = xml_root.find("BlastOutput_iterations").findall("Iteration"); #list of xml elements <iteration>
for iteration in xml_iterations: #Loops through each alignment query.
	xml_iterations_queryID = iteration.find("Iteration_query-ID").text; 
	temp_hit_set.add(xml_iterations_queryID); #prevents self hits
	xml_hits = iteration.find("Iteration_hits").findall("Hit"); #List of xml hits for each iteration.
	if xml_iterations_queryID not in temp_iteration_set: #prevents looping through alignments already done.
		for hit in xml_hits: #Loops through alignment hits in xml.
			hit_id = hit.find("Hit_id").text;
			hit_evalue = float(hit.find("Hit_hsps").find("Hsp").findall("Hsp_evalue")[0].text);
			if hit_evalue <= e_value_threshold and hit_id not in temp_hit_set: #adds hits e_value less than required to list, excludes hits already in list.
				gene_list.append((gi_extract(hit_id), hit_evalue));
				temp_hit_set.add(hit_id);
		temp_iteration_set.add(xml_iterations_queryID); #to ignore processing an iteration twice

gene_list.sort(key=lambda tup: tup[1]); # Sort in place by evalue (tup[1]), small to big

if len(gene_list) > number_unique_gids: #Selects the requested number of lowest e-values from gene_list then outputs to gene_list_master
	i=0;
	while (i<number_unique_gids):
		gene_list_master.append(gene_list[i][0]);
		i+=1;
else: # Handler for event that there are fewer list items than requested
	for tu in gene_list: # where tu = the tuple (gid, evalue)
		gene_list_master.append(tu[0]);

print `len(gene_list)`+" unique GenBank ids found with e-value < "+`e_value_threshold`+". Fetching "+`len(gene_list_master)`+" records from NCBI."

file_dir_output = queryID+"_"+time.strftime("%d%m%y")+"_records.fas" #name of multifasta file to be created. Probably want to change this

if xml_root.find("BlastOutput_program").text == "blastn" or xml_root.find("BlastOutput_program").text == "tblastn" or xml_root.find("BlastOutput_program").text == "tblastx": # database is protein by default
	database = "nucleotide";

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

