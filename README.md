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
1. Yang, Y., Lin, Y., Chen, Z., Gong, T., Yang, P., Girault, H., Liu, B., Qiao, L. Bacterial whole cell typing by mass spectra pattern matching with bootstrapping assessment. *Anal Chem* **89**, 12556–12561 (2017). https://doi.org/10.1021/acs.analchem.7b03820.
2. Zhang, D., Yang, Y., Qin, Q., Xu, J., Wang, B., Chen, J., Liu, B., Zhang, W., Qiao, L. MALDI-TOF characterization of protein expression mutation during morphological changes of bacteria under the impact of antibiotics. *Anal Chem* **91**, 2352–2359 (2019). https://doi.org/10.1021/acs.analchem.8b05080.

## License
BacteriaMS-difference is distributed under a BSD license. See the LICENSE file for details.


## Contacts
Please report any problems directly to the github issue tracker. Also, you can send feedback to liang_qiao@fudan.edu.cn.
