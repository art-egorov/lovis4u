# Gallery

Here we show usage examples of the LoVis4u command-line interface. 
For description of parameters and modes shown on this page visit the [Example-driven guide page](https://art-egorov.github.io/lovis4u/ExampleDrivenGuide/cmd_guide/). 

## Sample data & Configuration files 

The necessary sample data as well as adjustable tool configuration files are provided by lovis4u at the post-install step with the command:    
`lovis4u --data` This copies the *lovis4u_data* folder to your working directory.


Also important to describe the set of editable configuration files with parameters already adjusted for different tasks and page layout:

- **A4p1** - adjusted for one-column portrait A4 page layout. Figure width - 90mm.
- **A4p2** - adjusted for two-column portrait A4 page layout. Figure width - 190mm.
- **A4L** - adjusted for landscape A4 page layout. Figure width - 240mm.
- **standard** - flexible width length, scale - 0.05 mm per nucleotide with minimal width - 80mm. 

???+ card_hearts "**A4 page layout examples**" 
	[***Corresponding PDF file***](gallery/img/A4_page.pdf){target="_blank"}
	![a4page](gallery/img/A4_page.png){loading=lazy width="100%" } 
	

## The full-length visualisation of loci

### Miltiple sequences

**Visualisation with homology lines and category colours. Figure size: two-column of A4 page publication layout (190mm width).**

```sh
lovis4u -gff lovis4u_data/guide/gff_files -hl --set-category-colour -c A4p2
```

![f2](gallery/img/lovis4u_default_hl_A4p2.png){loading=lazy width="100%" }  

---

**Applying additional functional annotation of proteins using pyhmmer hmmscan. Figure size: two-column of A4 page publication layout (190mm width).**

```sh
lovis4u -gff lovis4u_data/guide/gff_files -hl --set-category-colour --run-hmmscan -c A4p2
```

![f2](gallery/img/lovis4u_hmmscan.png){loading=lazy width="100%" }  


--- 

**Visualisation with homology lines and category colours. Compact one-column A4 page figure without labelling of conserved proteins (90mm width).**


```sh
lovis4u -gff lovis4u_data/guide/gff_files -hl --set-category-colour \
	--show-first-feature-label-for  -c A4p1
```

![f2](gallery/img/lovis4u_full_length_a4p1.png){loading=lazy width="500px" }  

---

### One sequence

**Visualisation of a single sequence with highlighted protein groups of homologues by different colours and additional functional annotation. Two-column A4 page publication layout (190mm width).**

```sh
lovis4u -gff lovis4u_data/guide/gff_files/NC_001895.1.gff -hl --set-category-colour -c A4p2 \
	--set-group-colour-for conserved  --run-hmmscan --defence-models PADLOC
```

![f2](gallery/img/lovis4u_single_suquence.png){loading=lazy width="100%" }  

**Visualisation of a single sequence with GC content and GC skews tracks. Two-column A4 page publication layout (190mm width).**

```sh
 lovis4u -gff lovis4u_data/guide/gff_files/NC_001895.1.gff -hl --set-category-colour -c A4p2 \
    -gc -gc_skew
```
![f212](gallery/img/lovis4u_single_suquence_with_GC.png){loading=lazy width="100%" } 

## Visualisation of user-defined regions

**Visualisation of user-defined regions, one region per locus. One-column A4 page layout (90mm width).**

```sh
lovis4u -gff lovis4u_data/guide/gff_files -c A4p1 --locus-annotation-file lovis4u_data/guide/locus_annotation_table_demo_2.tsv \ 
	-cl-off -hl --set-category-colour
```

![f2](gallery/img/lovis4u_region_a4p1.png){loading=lazy width="500px" }  

**Visualisation of user-defined regions, multiple regions per locus.**


```sh
lovis4u -gff lovis4u_data/guide/gff_files  -hl -o lovis4u_output_m_loci  \
    --locus-annotation-file lovis4u_data/guide/locus_annotation_table_demo_1.tsv -cl-off -llp left  -slt
```

![f2](gallery/img/lovis4u_regions_multiple.png){loading=lazy width="800px" }  

### Non-coding features

**Visualisation of loci with non-coding features (tRNAs) and displaying their labels.**

```sh
lovis4u -gff lovis4u_data/guide/BaselSubset -laf lovis4u_data/guide/locus_annotation_table_trnas.tsv \
     -hl --set-category-colour -c A4p1 --show-first-noncoding-label
``` 
 
![fn](gallery/img/lovis4u_trnas_first_label.png){loading=lazy width="530px" }  

---

## Visualisation of genomic signal tracks from sequencing experiments

Starting with version 0.1.0, LoVis4u supports visualisation of genomic signal tracks from sequencing experiments using bedGraph file as input for a single track. We will demonstrate this functionality using three DNA sequencing coverage profiles for BASEL collection phage Bas01.

**Visualisation of a single genome with GC content, GC skews tracks and three bedgraph track profiles. Two-column A4 page publication layout (190mm width).**

```sh
 lovis4u -gff lovis4u_data/guide/BaselCollection/Bas01.gff -c A4p2 -scc -gc -gc_skew \ 
	-bg lovis4u_data/guide/bedgraphs/Bas01_DNA-Seq_WT.bedGraph lovis4u_data/guide/bedgraphs/Bas01_DNA-Seq_1.bedGraph lovis4u_data/guide/bedgraphs/Bas01_DNA-Seq_2.bedGraph \
	-bgl "Bas01 WT | DNA-Seq" "Bas01 E1 | DNA-Seq" "Bas01 E2 | DNA-Seq"
```

![f_bg2](gallery/img/lovis4u_bg_files_advanced.png){loading=lazy width="100%"}

**Visualisation of a single genome with GC content, GC skews tracks and three bedgraph track profiles. Multiple windows per image. Two-column A4 page publication layout (190mm width).**

```sh
lovis4u -gff lovis4u_data/guide/BaselCollection/Bas01.gff -c A4p2 -scc -gc -gc_skew \
 -bg lovis4u_data/guide/bedgraphs/Bas01_DNA-Seq_WT.bedGraph lovis4u_data/guide/bedgraphs/Bas01_DNA-Seq_1.bedGraph lovis4u_data/guide/bedgraphs/Bas01_DNA-Seq_2.bedGraph \ 
-bgl "Bas01 WT | DNA-Seq" "Bas01 E1 | DNA-Seq" "Bas01 E2 | DNA-Seq" \ 
 -w Bas01:1:15000:1,Bas01:16000:31000:1 
```
![f_bg2](gallery/img/lovis4u_bg_files_advanced_windows.png){loading=lazy width="100%"}


---

## Visualisation of the Basel collection of phages

For demonstration of LoVis4u's capabilities (and because it is fun), I applied it to the entire Basel collection of 78 *Escherichia coli* phages (see [Maffei et.al. PLOS Biology](https://doi.org/10.1371/journal.pbio.3001424)). In addition, we activate pyhmmer functional annotation of proteins.

```sh
lovis4u -gff lovis4u_data/guide/BaselCollection  -hl --set-category-colour -c A4p2 -fw 500  --run-hmmscan
```

[***LoVis4u PDF file of the visualisation***](gallery/img/lovis4u_Basel.pdf){target="_blank"}

??? card_hearts "**Raster image**" 
	![a4page1](gallery/img/lovis4u_Basel.png){loading=lazy width="100%" } 

This run also highlights one of LoVis4u key features—its ability to detect clusters of similar proteomes within the input data. Beyond visualisation, LoVis4u generates a matrix of pairwise proteome similarity (see the pipeline description page). The visualisation of this matrix for the Basel collection appears as follows:

![f2](gallery/img/Heatmap.png){loading=lazy width="650px" }  





