# Сommand-line parameters

	
**POST-INSTALL DATA AND CONFIGURATION**

- `--data`  
Creates the *uorf4u_data* folder in the current working directory.
The folder will contain the adjustable configuration file templates, palettes, tables as well as the necessary sample.

- `--linux`  
**All Linux user should run it only once after installation.**  
Replaces the tools paths in the premade config files from the MacOS' version [default] to the Linux'.
 
- `--blastp_path`  
Update the blastp path in the pre-made config files.  
**Required for using local blastp databases with** `-lbdb` **parameter.**  


**MANDATORY ARGUMENTS**

- `-an` *accession_number*  
Protein's RefSeq accession number.

	OR

- `-hl` *accession_number1 [accession_number2, ...]*  
Space separated list of proteins accession numbers which will be used as list of homologous.

	OR

- `-hlf` *file.txt*  
Path to a file with list of accession numbers. File format: one accession number per line, no header.

	OR

- `-fa` *file.fa*  
    Path to a fasta file with upstream sequences.


- `-c` *bacteria|eukaryotes|file.cfg*  
Path to a configuration file [default: internal].


**OPTIONAL ARGUMENTS**

- `-bdb` *<efseq_select|refseq_protein*  
Online blastp database to perform blastp searching for homologues.  
[default: from config; *refseq_select* for bacteria, *refseq_protein* for eukaryotes]

- `-lbdb` *path to a database*  
Local blastp database to perform blastp searching for homologues.  
Note: You have to specify path to a blastp with --blastp_path command before using this argument.

- `-bh` *number_of_hits*  
Max number of blastp hits in homologous searching.

- `-bid` *identity_cutoff [0-1]*  
BlastP searching cutoff for hit's identity to your query protein.

- `-mna` *number_f_assemblies*  
Max number of assemblies to take into analysis for each protein. If there are more sequences in the identical protein database then random sampling will be used.

- `-al` *path_to/assemblies_list.tsv*  
Path to an assemblies list file. During each run of uorf4u, a tsv table with information about assemblies (from identical protein database, ncbi) for each protein is saved to your output folder (output_dir_name/assemblies_list.tsv). There are cases with multiple assemblies for one protein accession numbers (up to thousands). In case to control assemblies included in the analysis this table can be filtered (simply by removing rows) and then used with this parameter as part of input to the next run.  
In addition, config file (see config parameters section) has max_number_of_assemblies parameter. It can be used to limit max number of assemblies included in the analysis. In case number of assemblies is more than the cutoff, random sampling will be used to take only subset of them.

- `-annot`  
Retrieve sequences annotation (to be sure that annotated uORFs is not overlapped with a known CDS. 

- `-ul` *length*     
Length of upstream sequences to retrieve.

- `-dl` *length*    
Length of downstream sequences to retrieve.

- `-asc`   
Include alternative start codons in uORF annotation step. List of alternative start codons are taken from the ncbi genetic code.

- `-nsd`  
Deactivate filtering ORFs by SD sequence presence. [default: True for 'prokaryotes' config and False for 'eukaryotes' config].

- `-at` *aa|nt*  
Alignment type used by uorf4u for conserved ORFs searching [default: aa]. 

- `-pc` *cutoff [0-1]*  
A cutoff of presence (number of ORFs in a list/number of sequences) for an ORFs set to be called conserved and returned [default: 0.4, set in the config].

- `-fast`  
Fast searching mode with less accuracy (>~300 sequences or >~2000 ORFs).

- `-o` *dirname*  
Output dirname. It will be created if it's not exist. All output dirs will be then created in this folder [default: uorf4u_{current_date}; e.g. uorf4u_2022_07_25-20_41].


**MISCELLANEOUS ARGUMENTS**

- `-h`, `--help`  
Show help message and exit.

- `-v`, `--version`  
Show program version.

- `--debug`  
Provide detailed stack trace for debugging purposes.

- `--quiet`  
Don't show progress messages.
