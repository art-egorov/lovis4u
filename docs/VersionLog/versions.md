# Version log

* **Ver 0.1.5** - 13 May 2025
	- New command line parameters to control mmseqs run are added.

* Ver 0.1.4.1 - 14 March 2025
	- An issue with reading files without any features is fixed.

* Ver 0.1.3 - 5 March 2025
	- A bug with CDS group type and colour assignment was fixed. 

* Ver 0.1.2 - 17 January 2025
	- A bug with reading non-integer count bedGraph files is fixed.
	- New parameter `-smp, --set-mmseqs-path <path>` to specify user's mmseqs path is introduced.

* Ver 0.1.1 - 14 January 2025
	- Support of bigWig files for coverage profiles input with `-bw` parameter.

* Ver 0.1.0 - 4 January 2025
	- Visualisation of genomic signal tracks from sequencing expeirments is now supported!
	- Multiple features and styles for single locus visualisation are introduced (GC content, GC skew, novel styles for label)
	- Multiple parameters are introduced, including `-w, --windows` that allows to specify coordinates for visualisation more easily.

* Ver 0.0.14 - 29 December 2024
	- An issue with visualisation of feature lables for single locus without mmseqs clustering is was fixed. 

* Ver 0.0.13 - 19 December 2024
	- Bug with parsing folder containing genbank files is fixed.

* Ver 0.0.12 - 4 December 2024
	- Minor bug related to non-coding features without specified orientation is fixed.

* Ver 0.0.11 - 24 October 2024
	- Functional annotation of proteins with pyhmmer hmmscan is now supported!
	- A set of parameters to control hmmscan annotation is addded. 
	- Now single file path input instead of folder name with '-gff' and '-gb' parameters is supported.

* Ver 0.0.10.1 - 22 October 2024
	- `-snl, --show-noncoding-labels` and `-sfnl, --show-first-noncoding-label` paraemetrs are added for controlling non-coding feature labels.

* Ver 0.0.10 - 22 October 2024
	- Visualisation and parsing of non-coding features (tRNAs, tmRNAs, pseudogene..) is added!
	- New parameter `-alip, --add-locus-id-prefix` is added in order to process the duplication feature problem.

* Ver 0.0.9.3 - 1 October 2024 
	- Linux path adjustment bug is fixed.
	- Genbank files parsing adjusted for incomplete features.

* **Ver 0.0.9** - 27 September 2024 
	- Set of new configuration files and styles for layouts of A4 page (one-,two-columns, and landscape)
	- New design and locus label positions.
	- New config and cmd parameters.
	- Update of the home page with new uiser guide and the gallery page.

* Ver 0.0.8 - 18 September 2024 
	- Suffix requirments for input files are removed.
	- Subset of files with incorrect format will not lead to the error, only to the warning message.

* Ver 0.0.7 - 4 September 2024 
	-  Minor API update. 
	- `--one-cluster` parameter is added.

* Ver 0.0.6 - 2 July 2024 
	-  Update of parameter names.

* Ver 0.0.5 - 27 June 2024 
	-  Minor fixes.

* Ver 0.0.4 - 13 June 2024 
	-  Minor fixes.

* Ver 0.0.3 - 3 May 2024 - first public release. 
