# Short example-drived guide to uorf4u API.  

uorf4u has a simple API allowing it programmatic usage from within a python program. Below we descrive several Python snippets that mimic results of command-line calls.

**Using a single RefSeq protein accession number as input in bacteria mode**
```python
import uorf4u

parameters = uorf4u.manager.Parameters()
parameters.load_config("bacteria")  # Load config (bacteria, eukaryotes or path)
parameters.arguments["output_dir"] = "ErmCL"  # to change parameters

protein = uorf4u.data_processing.RefSeqProtein(accession_number="WP_001003263.1",
                                               parameters=parameters)
homologues_list = protein.blastp_searching_for_homologues()
homologues = uorf4u.data_processing.Homologues(homologues_list, parameters)
upstream_sequences_records = homologues.get_upstream_sequences()

upstream_seqs = uorf4u.data_processing.UpstreamSequences(upstream_sequences_records, parameters)
upstream_seqs.annotate_orfs()
upstream_seqs.filter_orfs_by_sd_annotation()
upstream_seqs.save_annotated_orfs()
upstream_seqs.conserved_orf_searching()
upstream_seqs.filter_out_similar_paths()
upstream_seqs.run_msa()
upstream_seqs.save_orfs_sequences()
upstream_seqs.save_msa()
upstream_seqs.save_results_summary_table()
upstream_seqs.plot_annotation()
upstream_seqs.plot_logo_figs()
upstream_seqs.plot_msa_figs()

```  

**Will be updated with more examples...**