# Configuration file parameters


lovis4u configuration file allows detailed customization of the tool's parameters.

***Note:***  
You can copy the *lovis4u_data* folder that contains the config files to your wiking directory with `lovis4u --data` command and safely edit and use them without affecting 'internal' config. If you want to use a copied config file, use `-c path_to_config`.


##### Description is under construction...

---

**;[Data Processing]**  
default_transl_table = 11  
gff_description_source = description  
gff_CDS_name_source = product 
gff_CDS_category_source = function  
genbank_id_source = protein_id  
genbank_id_alternative_source = ID,id,PROTEIN_ID  
genbank_description_source = annotations:organism  
genbank_CDS_name_source = product  
genbank_CDS_category_source = function  
CDS_is_variable_cutoff = 0.25  
keep_predefined_groups = True  
keep_predefined_colors = True  

**;[mmseqs parameters]**  
mmseqs_cluster_mode = 0  
mmseqs_cov_mode = 0  
mmseqs_min_seq_id = 0.35  
mmseqs_c = 0.8    
mmseqs_s = 6.5  

**;[Paths]**  
palette = {internal}/palette.txt  
mmseqs_binary = {internal}/bin/mmseqs_mac/bin/mmseqs    
category_colors = {internal}/category_colors.tsv  
font_italic = {internal}/fonts/Lato-Italic.ttf  
font_light = {internal}/fonts/Lato-Light.ttf  
font_light_italic = {internal}/fonts/Lato-LightItalic.ttf  
font_regular = {internal}/fonts/Lato-Regular.ttf  
font_bold = {internal}/fonts/Lato-Bold.ttf  
font_mono = {internal}/fonts/RobotoMono-Regular.ttf  

**;[Output]**  
verbose = True  
debug = False  
output_dir = lovis4u_{current_date}  

**;[General figure parameters]**  
margin = 0.1  
gap = 0.1  
mm_per_nt = auto  
figure_width = None  

**;[Locus track primary settings]**  
locus_label_style = full  
locus_label_description_font_face = italic  
locus_label_id_font_face = regular  
draw_individual_x_axis = True  
draw_middle_line = False  
feature_default_fill_color = fancy_grey  
feature_default_stroke_color = grey  
set_feature_stroke_color_based_on_fill_color = True  
feature_stroke_color_relative_lightness = 0.4  
feature_group_types_to_set_color = variable,labeled  
groups_fill_color_palette_lib = seaborn  
groups_fill_color_seaborn_palette = husl  
groups_fill_color_seaborn_desat = 0.95  
feature_labels_to_ignore = hypothetical protein, unknown protein  
feature_group_types_to_show_label = variable,labeled  
feature_group_types_to_show_label_on_first_occurrence = shell/core  
show_all_feature_labels = False  
category_color_seaborn_palette = hls  
category_color_seaborn_desat = 0.9  

**;[Locus track advanced settings]**  
feature_height = 0.3  
feature_stroke_width = 0.3  
feature_arrow_length = 0.6  
feature_bottom_gap = 0.04  
locus_label_size = 0.7  
locus_label_color = black  
locus_label_color_alpha = 1  
gap_after_locus_label = 0.1  
x_axis_ticks_height = 0.07  
x_axis_ticks_labels_height = 0.1  
x_axis_ticks_labels_font_face = light  
x_axis_line_width = 0.3  
x_axis_line_color = black  
x_axis_line_color_alpha = 1  
category_annotation_alpha = 1  
category_annotation_line_width = 0.07  
groups_fill_colors_pastel_factor = 0.4  
feature_stroke_color_alpha = 1  
feature_fill_color_alpha = 0.7  
groups_stroke_colors_alpha = 1  
feature_label_size = 0.5  
feature_label_gap = 0.5  
feature_label_font_face = regular  
gap_between_regions = 0.5  

**;[Homology_track]**  
homology_fill_color = lightgrey  
homology_fill_color_alpha = 1  
homology_line_width = 0.1  
homology_stroke_color = grey  
homology_stroke_color_alpha = 0  

**;[Scale line track]**  
scale_line_relative_size = 0.3  
scale_line_width = 0.3  
scale_line_tics_height = 0.05  
scale_line_color = black  
scale_line_color_alpha = 1  
scale_line_label_height = 0.13  
scale_line_label_font_face = regular  

**;[Category color legend]**  
color_legend_label_size = 0.15  
color_legend_line_height = 0.07  
color_legend_font_face = regular  
color_legend_label_color = black  
color_legend_label_color_alpha = 1  


