# FAQ


 __*Q:* Something weird is happening, e.g. program is crashing with my data. What should I do?__  
 *A:* Please add `--debug` to the command line (or `--parsing-debug` if some files were not read correctly) and post the details in [Issues](https://github.com/art-egorov/lovis4u/issues?q=).



 __*Q:* GenBank files "was not read properly and skipped".__  
 *A:* Firstly, use the `--parsing-debug` parameter. Old prokka version generated GenBank files with incorrect locus header which then cannot be read by biopython library properly. If it is the case, you will find *"ValueError: Did not recognise the LOCUS line layout:"* line in the traceback. In that case you can use GFF files that are produced by prokka in parallel, or update prokka which should solve the problem.




