# Example-driven guide

Here we show usage examples of lovis4u command-line interface. Through this guide we will show step-by-step how you can optimise your visualisation starting from default parameters.
 
**Before start:** The necessary sample data as well as adjustable tool configuration files are provided by lovis4u at the post-install step:    
`lovis4u --data` which copies *lovis4u_data* folder to your working directory. 
  
**If you work on a Linux machine** after installation you should run: `uorf4u --linux`  
This command replaces the tools paths (MMseqs2) in the pre-made config files from the MacOS version (default) to the Linux.  
If you run it for fun and want to change it back you can use `uorfu --mac`.

^^For demonstration we will use pharokka generated gff files with sequences of 5 Enterobacteria phages.  
Gff files are stored at: *lovis4u_data/guide/gff_files*.^^      
The main difference of pharokka generated gff files from regular gff3 (for ex. which you can download from the NCBI) is that in addition to annotation rows they contain corresponding to the annotation nucleotide sequence. 

---

## Example run with default parameters

Let's start with running lovis4u without using any optional argument. The only mandatory argument is a folder path containing pharokka generated gff (`-gff`) files or genbank files (`-gb`).   

```
lovis4u -gff lovis4u_data/guide/gff_files
``` 

As results of running this command an output folder named lovis4u_{current_date} (e.g. lovis4u_2024_04_28-16_36) will be created.  
Name of the output folder can be changed with `-o <output_folder_name>`.

**Output folder structure:**

- *lovis4u.pdf* - vector graphic output (file name can be changed with `--pdf-name <filename>` parameter)  
- *loci_annotation_table.tsv* - table containing annotation (sequence_id, length, coordinates, etc..) for each locus.  
- *features_annotation_table.tsv* - table containing annotation (feature_id, locus_id, coordinates, etc..) for each feature (e.g. CDS)  
- *mmseqs* (folder)  
	- *DB* - folder with mmseqs' databases.  
	- *mmmseqs_clustering.tsv* - table with proteomes clustering results.  
	- *mmmseqs_(stdout/stderr).txt* - mmseqs logs.  
	- *input_proteins.fa* - fasta file with all annotated protein sequences (input to mmseqs).  
- *proteome_similarity_matrix.tsv* - pairwise proteome similarity scores indicating fraction of shared proteins homologues.  

*Visualisation results:*
<img  src="img/lovis4u_default_1.png" width="100%"/>

**By default, lovis4u utilises the following data preprocessing steps:**

1. Full length of each locus are taken for analysis *(this can be adjusted using loci annotation table, see below)*. 
2. All proteins annotated on input sequences are used as input for MMseqs2 clustering *(can be deactivated with `--mmseqs-off` parameter)*. MMseqs2 arguments can be adjusted using config file. Proteins clustered together are considered as a set of homologues. Based on that the "group" attribute of each CDS is set. 
3. Taking into account information about set of homologues from the previous step, lovis4u applies similarity based hierarchical clustering of corresponding proteomes with which it finds optimal order for visualisation and set "group" attribute for each locus. The purpose is to group together only related proteomes (keeping average proteome set similarity > ~80%). This step can be skipped with (`-cl-off` or `--clust_loci-off`) parameter.
4. Defining feature attribute "group_type" which allows to apply visualisation parameter targeting particular set of feature groups (e.g. set color or show labels only for "group_type" = "variable"). By default it sets group_type "variable" for CDS features that found in less than 0.25 of loci within the loci group and "shell/core" for others. 
5. Setting feature color based on feature "group" attribute *(can be deactivated with `-sgc-off` or `--set-group-color-off`)*. Be default it sets distinct colours only for features with group_type "variable". But, you can change it with `-sgcf` or `--set-group-color-for`. For example, if you want to set color only for features with group_type "shel/core", run `--set-group-color-for shell/core`. 
6. Defining labels to be shown. By default lovis4u shows all labels for "variable" features and only first occurrence for "shell/core" features. You can show all labels with `--show-all-feature-labels` or specify group types for which all labels will be shown with `-sflf` or `--show-feature-label-for`. Additionally, by default lovis4u ignores labels:  hypothetical protein, unknown protein. The list of ignored labels can be set with `-ifl, --ignored-feature-labels <feature_label1 [feature_label2 ...]>`. The list with this argument can be left empty to not filter out labels by their name.

**We also consider corresponding parameters as highly useful for basic runs:**

- `--reorient_loci` - Auto re-orient loci (set new strands) if they are not matched.  Function tries to maximise co-orientation of homologous features.
- `-hl`, `--homology-links` - Draw homology link track.
- `-o <name>` - Output dir name.  

While loci in our test set are already correctly orientated, let's add -hl parameter to update the output.

```
lovis4u -gff lovis4u_data/guide/gff_files  -hl -o lovis4u_output
```
<img  src="img/lovis4u_default_hl.png" width="100%"/>

## Using loci annotation table 

As it was already mentioned, full length of each locus is taken for visualisation by default. However, you can specify coordinates of multiple regions for each locus to be shown. This coordinates together with other information about each locus can be specified in loci annotation table and used as input with `-laf` or `--loci-annotation-file` parameter. 

Additionally, after each run lovis4u saves the *loci_annotation_table.tsv* with annotation parameters used in this particular run. If no table was specified by input then all annotation columns are set with default values.  

*Default table generated from previous runs:*

{{ read_table('loci_annotation_table.tsv', sep = '\t') }} 

After default run we can take the output loci_annotation_table and edit information we want to change.  
**Important to note** that we can use as input a table only with subset of columns (only *sequence_id* column is essential), for other columns or empty cells, lovis4u will set default values.

For example, let's use this table as input:

{{ read_table('loci_annotation_table_alt.tsv', sep = '\t') }} 

Here we specified only the coordinates, order and group for each locus. Order and group are also specified and kept since it's logical to use clustering results defined on full locus length run and to turn off new attempt to cluster sequences with `-cl-off, --clust_loci-off` parameter.   

**Format for coordinates specification:** comma-separated list of start\:end:strand. Start and end are in 1-based format, strand: 1 for plus strand and -1 for minus.

The table can be found in the guide folder: *lovis4u_data/guide/loci_annotation_table_alt.tsv*

Now we can run: 
```
lovis4u -gff lovis4u_data/guide/gff_files  -hl -o lovis4u_output \
	--loci-annotation-file lovis4u_data/guide/loci_annotation_table_alt.tsv -cl-off
```
<img  src="img/lovis4u_window_1.png" width="100%"/>


## Using features annotation table 

Similar way of adjusting features visualisation parameters is implemented. After each run a *features_annotation_table.tsv* file is saved in your output folder. Below you can see header of default table created with default parameter run.

{{ read_table('features_annotation_table_alt_header.tsv', sep = '\t') }} 

And here we also can use in input with parameter `-faf` or `--features-annotation-file` a table that contains only subset of annotation columns. For example, you can specify only new label and color for a particular CDS.

{{ read_table('features_annotation_table_alt.tsv', sep = '\t') }} 

This table can be found in the guide folder: *lovis4u_data/guide/features_annotation_table_alt.tsv*

```
lovis4u -gff lovis4u_data/guide/gff_files  -hl -o lovis4u_output \
	--loci-annotation-file lovis4u_data/guide/loci_annotation_table_alt.tsv -cl-off \
	--features-annotation-file lovis4u_data/guide/features_annotation_table_alt.tsv
```
<img  src="img/lovis4u_updated_name.png" width="100%"/>

Here also **important to note** that if you don't use feature attribute *"group_type"* defined by full loci run, then it will be re-defined by mmseqs run on subset of proteins found only in specified regions (you can deactivate it with `--find-variable-off` and `--mmseqs-off`) . It can result in situation that ORF called in that case "variable" in reality are encoded by each proteome but outside of their shown coordinates. 

## Other features

### Category color and annotation 

Using parameter `--set-category-color` you can use category annotation column for features. By default it was designed to parse PHROGs category annotation for proteins and retrieve information about category in "function" qualifiers in Genbank or GFF files *(used qualifiers can be changed in config file)*. However, you can set up category for each CDS using described above features annotation table with "category" column. Additionally, you can set up color code for your categories using `--category-color-table`. For categories not found in a table a random color will be set. By default, lovis4u uses pre-made color table which location is *lovis4u_data/category_colors.tsv*. 

```
lovis4u -gff lovis4u_data/guide/gff_files -hl --set-category-color
```
<img  src="img/lovis4u_category_colour.png" width="100%"/>


### Scale line instead of individual x-axis

Using parameter `--hide-x-axis` you can deactivate visualisation of individual x-axis for each locus track and instead of them with `-slt` or `--scale-line-track` you can draw a scale line track below.

```
lovis4u -gff lovis4u_data/guide/gff_files -hl --hide-x-axis --scale-line-track
``` 
<img  src="img/lovis4u_scale.png" width="100%"/>  

### Highlighting conserved genes instead of variable

For many analysis purposes (e.g. conserved neighbourhood visualisation) we need to colorise conserved gene clusters instead of variable. It can be easily switched in lovis4u using `--set-group-color-for` parameter. By default it's set as "variable" but using `--set-group-color-for shell/core` will change it to the opposite mode.  
**Note** that if you have other feature group set in your features annotation table and want to set auto-colorising for them as well you can specify them in space separated list with this argument (e.g. `--set-group-color-for shell/core your_group_1 your_group_2`.

```
lovis4u -gff lovis4u_data/guide/gff_files -hl --set-group-color-for shell/core 
```

<img  src="img/lovis4u_conserved_colorised.png" width="100%"/>  

**Note:** By default colours for groups are randomly set for each group using [seaborn husl palette](https://seaborn.pydata.org/tutorial/color_palettes.html). In config file you can change to more intense hsl palette or change the desaturation parameter.

### Specifying figure width

Lovis4u tries to set an optimal figure width taking into account nucleotide size of visualisation window. You can change it in two ways:  
1) Using `--mm-per-nt <float value>` argument changing scale which defines given space for each nt cell on canvas. Default: 0.0022 - 0.02, depending on window size.  
2) With `-fw, --figure-width <float value [cm]>` parameter which defines total output figure width.  
We demonstrate usage by plotting compact visualisation of full loci together with `--show-first-feature-label-for` argument with empty list deactivating showing first occurrence label for shell/core genes.

```
lovis4u -gff lovis4u_data/guide/gff_files -hl -o width_test --show-first-feature-label-for --figure-width 7
```
<img  src="img/lovis4u_set_width.png" width="100%"/>  

