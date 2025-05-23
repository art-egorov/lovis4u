# Сommand-line parameters



### Post-install steps

`--data`  
:    Creates the 'lovis4u_data' folder in the current working directory.
     The folder contains adjustable configuration files used by lovis4u
     (e.g. config, palettes...)

`--linux`
:    Replaces the mmseqs path in the pre-made config file from the MacOS
     version [default] to the Linux version.

`--mac`
:    Replaces the mmseqs path in the pre-made config file from the Linux

`-smp, --set-mmseqs-path <path>`
:	Specify mmseqs path that will be used by LoVis4u. Can be either full path
    to the binary mmseqs, "default_mac", or "default_linux". In case if mmseqs
    is installed in the system and available without specifying its path, you can
    run "-smp mmseqs" and this should work.

`--get-hmms`
:    Download HMMs (hmmscan format) of defence, anti-defence, virulence,
     and AMR proteins from our server [data-sharing.atkinson-lab.com]

### Mandatory arguments
`-gff <folder>`
:    Path to a folder containing extended gff files.
     Each gff file should contain corresponding nucleotide sequence.
     (designed to handle pharokka produced annotation files).

 OR

`-gb <folder>`
:    Path to a folder containing genbank files.

### Optional arguments | Data processing
`-w, --windows <locus_id1:start1:end1:strand [locus_id1:start1:end1:strand ...]>`
:    Specify window of visualisation (coordinates) for a locus or multiple loci

`-bg, --bedgraphs <bedGraph_file1 [bedgGaph_file2 ...]>`
:    Space separated list of paths to bedGraph files to plot coverage profiles.
     (>=1 file) ! Can be applied only for single locus.

`-bw, --bigwigs <bigWig_file1 [bigWig_file2 ...]>`
:    Space separated list of paths to bigWig files to plot coverage profiles.
     (>=1 file) ! Can be applied only for single locus.

`-bgl, --bedgraph-labels <bedgraph_label1 [bedgraph_label2 ...]>`
:    List of labels for bedgraph/bigwig tracks (order the same as order of bedgraph/bigwig files)
     By default basename of files will be used.

`-bgc, --bedgraph-colours <bedgraph_colour1 [bedgraph_colour2 ...]>`
:    List of colours for bedgraph tracks (order the same as order of bedgraph files)
    Each value can be either HEX code of colour name (e.g. pink, blue, etc (from the palette file))

`-gc, --gc-track`
:    Show GC content track. ! Can be applied only for single locus.

`-gc_skew, --gc_skew-track`
:    Show GC skew track. ! Can be applied only for single locus.

`-ufid, --use-filename-as-id`
:    Use filename (wo extension) as track (contig) id instead  
     of the contig id written in the gff/gb file.

`-alip, --add-locus-id-prefix`
:    Add locus id prefix to each feature id.

`-laf, --locus-annotation-file <file path>`
:    Path to the locus annotation table.
     (See documentation for details)

`-faf, --feature-annotation-file <file path>`
:    Path to the feature annotation table.
     (See documentation for details)

`-mmseqs-off, --mmseqs-off`
:   Deactivate mmseqs clustering of proteomes of loci.

`-mmsi, --mmseqs-min-seq-id <float>`
:    MMSeqs2 parameter for minimal sequence identity during clustering.

`-mc, --mmseqs-coverage <float>`
:    MMSeqs2 parameter for minimal coverage during clustering.

`-hmmscan, --run-hmmscan`
:    Run hmmscan search for additional functional annotation.

`-dm, --defence-models <DefenseFinder|PADLOC|both>`
:    Choose which defence system database to use for hmmscan search
     [default: both (DefenseFinder and PADLOC)]

`-hmm, --add-hmm-models <folder_path [name]>`
:    Add your own hmm models database for hmmscan search. Folder should
     contain files in HMMER format (one file per model). Usage: -hmm path [name].
     Specifying name is optional, by default it will be taken from them folder name.
     If you want to add multiple hmm databases you can use this argument several
     times: -hmm path1 -hmm path2.

`-omh, --only-mine-hmms`
:    Force to use only models defined by user with -hmm, --add-hmm-models parameter.

`-kdn, --keep-default-name`
:    Keep default names and labels for proteins that have hits with
     hmmscan search. [default: name is replaced with target hmm model name]

`-kdc, --keep-default-category`
:    Keep default category for proteins that have hits with hmmscan
     search. [default: category is replaced with database name]

`-salq, --show-all-labels-for-query`
:    Force to show all labels for proteins that have hits to any database with hmmscan search.
    [default: False]

`-cl-owp, --cluster-only-window-proteins`
:    Cluster only proteins that are overlapped with  
     the visualisation windows, not all.

`-fv-off, --find-variable-off`
:    Deactivate annotation of variable or conserved protein clusters.

`-cl-off, --clust_loci-off`
:    Deactivate defining locus order and using similarity based hierarchical
    clustering of proteomes.

`-oc, --one-cluster`
:    Consider all sequences to be members of one cluster but use clustering
    dendrogram to define the optimal order.

`-reorient_loci, --reorient_loci`
:    Auto re-orient loci (set new strands) if they are not matched.
     (Function tries to maximise co-orientation of homologous features.)

### Optional arguments | Locus visualisation
`-sgc-off, --set-group-colour-off`
:    Deactivate auto-setting of feature fill and stroke colours.
     (Pre-set colours specified in feature annotation table will be kept.)

`-sgcf, --set-group-colour-for <feature_group1 [feature group2 ...]>`
:    Space-separated list of feature groups for which colours should be set.
     [default: variable, labeled]

`-scc, --set-category-colour`
:    Set category colour for features and plot category colour legend.

`-cct, --category-colour-table <file path>`
:    Path to the table with colour code for categories.
     Default table can be found in lovis4u_data folder.

`-lls, --locus-label-style <id|description|full>`
:    Locus label style based on input sequence annotation.

`-llp, --locus-label-position <left|bottom>`
:    Locus label position on figure.

`-safl, --show-all-feature-labels`
:    Display all feature labels.

`-sflf, --show-feature-label-for  <feature_group1 [feature group2 ...]>`
:    Space-separated list of feature groups for which label should be shown.
     [default: variable, labeled]

`-sfflf, --show-first-feature-label-for <feature_group1 [feature group2 ...]>`
:    Space-separated list of feature group types for which label will be displayed
      only for the first occurrence of feature homologues group.
     [default: shell/core]

`-snl, --show-noncoding-labels`
:    Show all labels for non-coding features. [default: False]

`-sfnl, --show-first-noncoding-label`
:    Show labels only for the first occurrence for non-coding features.
     [default: False]

`-ifl, --ignored-feature-labels <feature_label1 [feature_label2 ...]>`
:    Space-separated list of feature names for which label won't be shown.
     [default: hypothetical protein, unknown protein]

`-sxa, --show-x-axis`
:    Plot individual x-axis for each locus track.

`-hix, --hide-x-axis`
:    Do not plot individual x-axis for each locus track.

`-dml, --draw-middle-line`
:    Draw middle line for each locus.

`-mm-per-nt, --mm-per-nt <float value>`
:   Scale which defines given space for each nt cell on canvas.
     [default: 0.05]

`-fw, --figure-width <float value>`
:    Output figure width in mm.

### Optional arguments | Additional tracks
`-hl, --homology-links`
:    Draw homology link track.

`-slt, --scale-line-track`
:    Draw scale line track.


### Optional arguments | others
`-o <name>`
:    Output dir name. It will be created if it does not exist.
  	 [default: lovis4u_{current_date}; e.g. uorf4u_2022_07_25-20_41]

`--pdf-name <name>`
:    Name of the output pdf file (will be saved in the output folder).
     [default: lovis4u.pdf]

`-c <standard|<file.cfg>`
:    Path to a configuration file or name of a premade config file
     [default: standard]

### Miscellaneous arguments
`-h, --help`
:    Show this help message and exit.

`-v, --version`
:    Show program version.

`--debug`
:    Provide detailed stack trace for debugging purposes.

`--parsing-debug`
:    Provide detailed stack trace for debugging purposes   
     for failed reading of gff/gb files.

`-q, --quiet`	
:    Don't show progress messages.

	
