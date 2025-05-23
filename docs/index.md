
# ![logo](img/lovis4u_logo.png){loading=lazy width="265" }  

## Description

**LoVis4u** is a bioinformatics tool for **Lo**cus **Vis**ualisation and coverage profiles from sequencing experiments.

**Supported genome annotation input**: Genbank, gff3 with sequence    
**Supported coverage profile input**: bedGraph, bigWig  
**Supported output**: Static vector graphics (pdf)   
**Programming language:** Python3   
**OS:** MacOS, Linux  
**Python dependencies:** biopython, bcbio-gff, scipy, configs, pandas, distinctipy, matplotlib, seaborn, reportlab, pyhmmer, progress, requests  
**Python version:** >= 3.8  
**OS-level dependencies:** MMseqs2, bigWigToBedGraph (included in the package)  
**License:** WTFPL  
**Version:** 0.1.5 (May 2025)



**Pipeline**

![pipeline](img/lovis4u_pipeline.png){loading=lazy width="100%" }  


**Visualisation example (comparative gemomics)**

![visex1](img/lovis4u_hmmscan.png){loading=lazy width="100%" }  


**Visualisation example (genomic signal tracks from sequencing experiments)**

![visex2](img/lovis4u_bg_files_advanced.png){loading=lazy width="100%" }  


See the [gallery page](https://art-egorov.github.io/lovis4u/Gallery/gallery/) for more examples

## Installation

- LoVis4u can be installed directly from pypi:

```
python3 -m pip install lovis4u
```

- The development version is available at github :

```
git clone https://github.com/art-egorov/lovis4u.git
cd lovis4u
python3 -m pip install --upgrade pip
python3 -m pip install setuptools wheel
python3 setup.py sdist
python3 -m pip install .
```

**!** If you're a linux user, run `lovis4u --linux` post-install command once to update paths in the premade config files that set by default for MacOS users.


## Reference 

If you find LoVis4u useful, please cite:  
Artyom. A. Egorov, Gemma C. Atkinson, **LoVis4u: a locus visualization tool for comparative genomics and coverage profiles**, *NAR Genomics and Bioinformatics, Volume 7, Issue 1, March 2025, lqaf009, doi: [10.1093/nargab/lqaf009](https://doi.org/10.1093/nargab/lqaf009)*


## Contact 

Please contact us by e-mail _artem**dot**egorov**AT**med**dot**lu**dot**se_ or use [Issues](https://github.com/art-egorov/lovis4u/issues?q=) to report any technical problems.  
You can also use [Discussions section](https://github.com/art-egorov/lovis4u/discussions) for sharing your ideas or feature requests! 

## Authors 

LoVis4u is developed by Artyom Egorov at [the Atkinson Lab](https://atkinson-lab.com), Department of Experimental Medical Science, Lund University, Sweden. We are open for suggestions to extend and improve LoVis4u functionality. Please don't hesitate to share your ideas or feature requests.

