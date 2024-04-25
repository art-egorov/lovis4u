# Version log
* **Ver 0.9.5** - 4 May 2023
	* Local blastp searching function now also returns the query accession number.
* Ver 0.9.4 - 19 April 2023
	* A minor bug with `-al` parameter is fixed.
* Ver 0.9.3 - 18 April 2023
	* New cmd parameter `-pdb` for choosing an online blastp database.
	* New cmd parameter `-lpdb` for using local blastp database.
* Ver 0.9.2 - 19 February 2023
	* Blastp searching is now available for protein sequences. 
* Ver 0.9.1 - 19 February 2023
	* A major bug of 0.9.0 version in conservation analysis was fixed.
* Ver 0.9.0 - 17 February 2023
	* Conservation analysis algorithm was updated.
	* -c parameter is now always required.
	* Logs messages were updated.
* Ver 0.8.7 - 20 January 2023
	* Logs messages were updated.
* Ver 0.8.6 - 18 January 2023
	* Sequences in MSA plot now are ordered according to their similarity. (--reorder MAFFT parameter)
* Ver 0.8.5 - 17 January 2023
	* Update for web version
* Ver 0.8.4 - 30 November 2022
	* Report files were updated.
* Ver 0.8.3 - 30 November 2022
	* Minor bugs were fixed.
	* Warning messages were updated.
	* NCBI database parsing was optimised.
* Ver 0.8.2 - 7 November 2022
	* xml files and assemblies annotation bugs were fixed.
	* Annotation parsing was optimised.
* Ver 0.8.1 - 4 November 2022
	* Large assemblies annotation bug was fixed.
* Ver 0.8.0 - 2 November 2022
	* New exceptions control.
	* New cmd parameter (*-pc*).
* Ver 0.7.0 - 1 November 2022
	* NCBI database parsing was optimised and became ~10 times faster.
* Ver 0.6.4 - 31 October 2022
	* MAFFT version 7.505 was replaced with v. 7.490 since it's more stable. 
* Ver 0.6.3 - 29 October 2022
	* Entrez.email for the NCBI requests was set. 
* Ver 0.6.2 - 26 October 2022  
	* A problem with xml file writing was fixed.
* Ver 0.6.1 - 25 October 2022
	* **!** After the NCBI API update all previous version have a bug with identical protein database parsing. The bug was fixed in this version. 
* Ver 0.6.0 - 23 October 2022
	* New implementation of MSA visualisation
* Ver 0.5.4 - 13 October 2022
	* Annotation visualisation' bug was fixed 
* Ver 0.5.3 - 12 October 2022
	* Minor bugs with pathes were fixed  
* Ver 0.5.2 - 7 October 2022
	* MSA tool's path bug was fixed.
	* Fast searching now is set as 'auto'.
* Ver 0.5.1 - 6 October 2022
	* Now compatible with python3.7 (previous versions were compatible only with python3.8+).
* Ver 0.5.0 - 4 October 2022
	* 'Eukaryotes' and 'Prokaryotes' mode were introduced.
	* MSA now perfomed with MAFT.
	* New cmd and configs parameters were added.
	* Minor bugs were fixed.
* Ver 0.4.0 - 30 August 2022
	* Visualisation of loci annotation was added.
	* Minor bugs were fixed.
* Ver 0.3.1 - 17 August 2022
	* New cmd and configs parameters were added.
	* Annotation of uORFs overlapped with the main CDSs was added.
* Ver 0.3.0 - 7 August 2022
	* Algorithm of conserved ORFs searching was updated.
	* New configs parameteres were added.
* Ver 0.2.1 - 5 August 2022
	* Annotation parsing bug was fixed.
* Ver 0.2.0 - 5 August 2022  
	* New cmd and configs parameters were added.
	* New classes and methods were developed.
* Ver 0.1.5 - 31 July 2022
	* MSA visualisation functions updated.
	* Bugs were fixed. 
	* New cmd and configs parameters were added.
* Ver 0.1 - 27 July 2022 - Initial release. 

