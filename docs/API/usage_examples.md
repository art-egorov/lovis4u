# Short example-drived guide to LoVis4u API  

LoVis4u has a simple API allowing it programmatic usage from within a python program. Below we describe Python snippets that mimic results of command-line calls.

**See detailed description of each class and method in the "Library" section.**





```python
import lovis4u

# Creating a parameters object and loading config
parameters = lovis4u.Manager.Parameters()
parameters.load_config("standard")

# Example of changing a particular parameter
parameters.args["output_dir"] = "API_output_example"

# To turn off progress messages
parameters.args["verbose"] = False

# Creating a loci object and loading gff files  
loci = lovis4u.DataProcessing.Loci(parameters=parameters)

# Loading locus and feature annotation tables (optional)
loci.load_locus_annotation_file("file_path")
loci.load_feature_annotation_file("file_path")

# Loading folder with gff files (for ex. the example-driven guide)
gff_folder = "lovis4u_data/guide/gff_files" 
loci.load_loci_from_extended_gff(gff_folder)

# Running mmseqs on all encoded proteins and processing results
mmseqs_clustering_results = loci.mmseqs_cluster()
loci.define_features_groups(mmseqs_clustering_results)

# Cluster loci (optional)
loci.cluster_sequences(mmseqs_clustering_results)

# Find variable protein groups (optional)
loci.find_variable_feature_groups(mmseqs_clustering_results)

# Reoirent loci (optional)
loci.reorient_loci()

# Set colours  (optional)
loci.set_feature_colours_based_on_groups()
loci.set_category_colours()

# Defining labels to be shown (optional)
loci.define_labels_to_be_shown()

# Saving annotation tables (optional)
loci.save_feature_annotation_table()
loci.save_locus_annotation_table()

# Visualisation steps
# Creating a canvas manager object
canvas_manager = lovis4u.Manager.CanvasManager(parameters)
canvas_manager.define_layout(loci)

# Adding tracks. The only mandatory: loci
canvas_manager.add_loci_tracks(loci)

# We can add scale line on the bottom (optional)
canvas_manager.add_scale_line_track()

# Category colours (optional)
canvas_manager.add_categories_colour_legend_track(loci)

# And homology line track (optional)
canvas_manager.add_homology_track()

# Finally, plotting results and saving the pdf file
canvas_manager.plot(filename="example.pdf")
 
```
