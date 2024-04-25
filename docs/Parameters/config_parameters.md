# Configuration file parameters


uorf4u configuration file allows detailed customization of the tool's parameters.

***Note:***   
uorf4u has two pre-made configuration files: *uorf4u_eukaryotes.cfg* and *uorf4u_bacteria* located in: ./uorftu/uorf4u_data/ folder (internal). By default, 'bacteria' config file is used if no path or name of premade file is specified by a cmd parameter: `-c  bacteria|eukaryotes|<file.cfg>`.  

You can copy the *uorf4u_data* folder that contains the config files to your wiking directory with `uorf4u --data` command and safely edit and use them without affecting 'internal' set of configs. If you want to use a copied config file, use `-c path_to_config`.

---

***;[General]***  
**ncbi_genetic_code_name = Bacterial**; *the ncbi genetic code name ('Standard' for eukaryotes' config)*   
**ncbi_entrez_email = uorf4u@gmail.com**; *e-mail for the NCBI API.*   
**upstream_region_length = 500**; *[int or 'all'] Length of upstream region to retrieve. 'all' value is set for eukaryotes config file since by default it uses only mRNAs sequences. (can be overriden by `-ul` cmd parameter).*  
**minimal_upstream_region_length = 300**; *[int] minimal upstream region length for sequences to retrieve. If available sequence length to retrieve is shorter then this record won't be taken in the analysis*.    
**downstream_region_length = 100**; *[int] downstream region (overlapped with CDS) length to retrieve. (can be overriden by `-dl` cmd parameter).*  
**filter_refseq_sequences_by_regex = True**; *[bool] use or not regex parameter (below) for filtering the NCBI RefSeq sequences to retrieve.*    
**refseq_sequences_regex = ^[N][ZCW]_.***; *[regex] that will be used to filter the NCBI RefSeq sequnces. For eukaryotes set as '^[NX]M_.*' that means that only sequences that start with NM_ or XM_ (mRNAs) will be taken in the analysis.*  
**max_number_of_assemblies = 1**; *[int] max number of assemblies to take into analysis for each protein. If there are more sequences in the identical protein database then random sampling will be used. (can be overriden by `-mna` cmd parameter).*


***;[blastp homologous searching]***  
**blastp_database = refseq_select**; *[refseq_select or refseq_protein] Online blastp database to perform blastp searching for homologues. (can be overriden by `-bdb` cmd parameter).*  
**blastp_evalue_cutoff = 1e-5**; *[float] blastp e-value cutoff during the searching for homologs against the RefSeq database.*  
**blastp_hit_list_size = 200**; *[int] max number of blastp hits to take in the analysis.*  
**blastp_max_number_of_alignments = 1000**; *[int] max number of alignments during the searching (there could be several alignments for 1 hit, see blastp documentation)*    
**blastp_pident_to_query_length_cutoff = 0.5**; *[float: 0-1] cutoff for hit's identity to your query protein.*  

***;[ORF annotation]***  
**alternative_start_codons = False**; *[bool] use or not set of alternative start codons.*   
**main_start_codon = ATG**;  *[str]*        
**min_orf_length = 9**; *[int] cutoff for ORFs length during annotation*    
**filter_by_sd = True**; *[bool] filter annotated ORFs by Shine-Dalgarno sequence prersence. Filtering based on calculation of binding energy between aSD sequence (UCCUCC) and putative SD sequence in an upstream to uAUG window. Energy calculation performed as described here: [Yang et.al, 2016](10.1534/g3.116.032227)*  
**sd_energy_cutoff = -3**; *[float] cutoff for aSD-SD binding energy.*    
**sd_window_length = 20**; *[int] length of a region for SD sequence search.*    
**check_assembly_annotation = False**; [bool] retrieve or not the NCBI sequences annotation to be sure that annotated uORFs are not overlapped with known CDSs (can be overriden by `-annot` cmd parameter).

***;[conserved ORFs searching]***  
**fast_searching = auto**; *[bool(true or false) or 'auto'] use or not fast searching mode with less accuracy (needed for >~200 sequences or >~2000 ORFs). Can be also set as auto [default]. (can be overriden by `-fast` cmd parameter).*  
**fast_searching_fraction_of_initial_genomes = 0.3**; *[bool] fraction of input sequences that will be used as initial step in algorithm searching. Applied only if the fast_searching parameter is True.*    
**orf_length_group_range = 0.25**; *[float or int] orf's lengths window within conserved uORFs set can be annotated. If it's a float value [0-1] then the radius of window is a set percentage of the claster's length, while if it's int then the window radius is fixed.*    
**orfs_presence_cutoff = 0.5**; *[float] a set of ORFs will be returned only if they were found in a fraction of input sequences larger than this cutoff.*    
**paths_identity_cutoff = 0.5**; *[float] if two sets of found ORFs are ovelapped more than this cutoff, then only a set with a higher. score will be returned. (Helps to remove duplicates).*    
**max_num_of_initial_genome_iteration = 100**; *[int] similar to the fast_searching_fraction_of_initial_genomes parametr, but used with a normal mode for optimisation.*    

***;[Pairwise alignment]***  
**alignment_type = aa**; *[nt or aa] alignment type of uORFs during conservation analysis.*    
; *Below listed global alignments parametersduring conservation analysis. uorf4u uses Bio.Align. package to perfome pairwise alignment of uORFs.*    
**global_match_score = 1**; *[float]*  
**global_mismatch_score = -1**; *[float]*    
**global_open_gap_score = -1**; *[float]*  
**global_extend_gap_score = -1**; *[float]*    
**global_target_end_gap_score = -1**; *[float]*    
**global_query_end_gap_score = -1**; *[float]*    
**alignment_score_cutoff = 0**; *[float] if a pairwise alignment score is larger than this cutoff then two uORFs are considered as aligned.*    

***;[Multiple Sequence Alignment]***  
**consensus_threshold = 0.7**;  *[float] treshold for MSA position to consider a nucleotide/amino acid as conserved in consensus sequence building.*  

***;[Paths]***   
;*Pathes to scripts and files used by the tool. {internal} means a folder uorf4u/uorf4u_data in the tool location.*  
**ref_energies = {config_path}/energyRef-CCTCCT.json**; *aSD-SD energy table downloaded from: [Yang et.al, 2016](10.1534/g3.116.032227)*   
**maft_binary = {config_path}/bin/mafft-mac/mafft.bat**    
**palette = {config_path}/palette.txt**    
**colors_nt = {config_path}/colors_nt.txt**    
**colors_aa = {config_path}/colors_aa.txt**    

***;[Output]***  
**sequences_to_write = nt, aa, sd**; *[list] type of sequences results for that (logos, MSAs, fasta files) will writetn. nt - nucleotide seqs of uORFs, aa - amino acid seqs, sd - SD seqs (sd is not available for 'eukaryotes' mode)*    
**logo_type = both**; *[probability, information or both] type of logo, see logomaker package documentation.*       
**output_dir = uorf4u_{current_date}**; *[str] default name of the output dir. default: uorf4u_{current_date}; e.g. uorf4u_2022_07_25-20_41. (can be overriden by `-o` cmd parameter).*  

;------------------------  
***;Annotation visualisation***  
;------------------------  
***;[General figure parameters]***  
margin = 0.1  
gap = 0.03  
label_gap = 0.07  
orf_height = 0.15  
annotation_width = auto  
mm_per_nt = 0.04  
font_regular = {config_path}/fonts/Lato-Regular.ttf  
font_bold = {config_path}/fonts/Lato-Bold.ttf  
font_mono = {config_path}/fonts/RobotoMono-Regular.ttf  

***;[Sequence labels]***  
label_color = #3D3D3D  
label_color_alpha = 1  
label_height_to_orf_height = 0.65  

***;[Axis tics]***  
axis_tics_font_size = auto  
axis_tics_line_width = 0.3  

***;[Loci annotations]***  
upstream_seq_line_color = #CECECE  
upstream_seq_line_color_alpha = 1  
upstream_seq_line_width = 0.5  
cds_seq_stroke_color = #489143  
cds_seq_stroke_color_alpha = 0.8  
cds_seq_fill_color = #9ee19b  
cds_seq_fill_color_alpha = 0.03  
orf_line_width = 0.5  
conserved_uorfs_stroke_color = #4e4e4e  
conserved_uorfs_stroke_color_alpha = 1  
conserved_uorfs_fill_color = #ee8fb1  
conserved_uorfs_fill_color_alpha = 0.6  
other_uorfs_stroke_color = #CECECE  
other_uorfs_stroke_color_alpha = 1  
annotated_orf_stroke_color = #3d6f8e  
annotated_orf_stroke_color_alpha = 1  


;------------------------  
**;MSA plot**  
;------------------------  
tile_size = 0.1  
tile_stroke = 0.05  
char_size = 0.7  
label_size = 0.6  
tile_alpha = 0.8  
logo_alpha = 0.8  
