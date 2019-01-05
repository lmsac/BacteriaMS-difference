# BacteriaMS-difference
Differential protein analysis of MALDI-TOF mass spectra of bacteria.


## Requirement
* [R] (https://www.r-project.org/). As an alternative, the latest version of [Microsoft R Open] (https://mran.microsoft.com/open/) should be fine.

* [RStudio] (https://www.rstudio.com/) is recommended but optional.


## Usage
If you have everything installed, you can run identification for a sample spectrum as follows:

1. Run [main.R] (main.R) to load the functions.

2. Peak alignment.
[alignment.R](scripts/alignment.R) outputs peak intensity table compatible with statistical analysis tools, e.g. [MetaboAnalyst](https://www.metaboanalyst.ca/).

3. Differential protein identification.
[matching.R](scripts/matching.R) matches m/z of peaks to protein/peptide sequences from fasta files.


## Publications
Dongxue Zhang, Yi Yang, Qin Qin, Juan Xu, Bing Wang, Jianwei Chen, Baohong Liu, Weijia Zhang and Liang Qiao, "MALDI-TOF characterization of protein expression mutation during morphological changes of bacteria under the impact of antibiotics", Submitted.


## License
BacteriaMS-difference is distributed under a BSD license. See the LICENSE file for details.


## Contacts
Please report any problems directly to the github issue tracker. Also, you can send feedback to liang_qiao@fudan.edu.cn.