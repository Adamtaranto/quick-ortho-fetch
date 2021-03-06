Quick-Ortho-Fetch
================================

What: Provides raw material for building phyogenies from paraphyletic gene families.  

How: Extracts a non-redundant list of best hits from a multi-query blast and fetches matching fasta seqs from NCBI.


Example
-------------------------

    python quickOrtho.py myMultiBlast.xml protein myemail@gmail.com -n 40 -e 1e-5 -o outputFile.fas -t -q

Takes the top 40 unique gids for each query from xml, outputs non-redundant .fas file containing sequences.

Options
-------------------------

usage: quickOrtho.py file_dir {protein,nucleotide} email [-h] [-n NUMBER_UNIQUE_GIDS] [-e E_VALUE_THRESHOLD] [-o OUTPUT_DIR] [-t table][-q quiet] 

<table>

  <tr>
	<th>Arg</th><th>Name</th><th>Description</th>
  </tr>
  
  <tr>
	<td>[1]</td><td>file_dir</td><td>Directory to NCBI .xml file.</td>
  </tr>
  
  <tr>
	<td>[2]</td><td>{protein,nucleotide}</td><td>Which database to direct entrez query to.</td>
  </tr>

  <tr>
	<td>[3]</td><td>email</td><td>Email for entrez record retrieval, tells NCBI who you are.</td>
  </tr>

  <tr>
	<td>-e</td><td>E_VALUE_THRESHOLD</td><td>Maximum e-value allowed in screening, enter as decimal or in scientific notation (eg. 1e-20). Default = 1e-20</td>
  </tr>

  <tr>
    <td>-n</td><td>NUMBER_UNIQUE_GIDS</td><td>Number of unique gids to extract for each query. Default = 50</td>
  </tr>

  <tr>
    <td>-o</td><td>OUTPUT_DIR</td><td>Set name of output fasta file. Default = "'input_dir'_quickOrthoResults.fas"</td>
  </tr>
  <tr>
    <td>-t</td><td>table</td><td>Creates a new .txt table summarising top hits for each query. Writes file to "'output_dir'_summaryTable.txt". Recommend viewing this file in a text editor without text wrapping.</td>
  </tr>
  <tr>
    <td>-q</td><td>quiet</td><td>Runs the program in quiet mode, with no running feedback</td>
  </tr>
  <tr>
    <td>-h</td><td>help</td><td>Print help message and exit</td>
  </tr>

</table>
