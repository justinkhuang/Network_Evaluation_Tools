# Network Evaluation Tools

Network Evaluation Tools is a Python package with corresponding examples for evaluating a network's ability to group a given node set in network proximity. This package was originally developed as a part of the work done in [Huang and Carlin et al. 2018](http://www.cell.com/cell-systems/fulltext/S2405-4712(18)30095-4). This version of the package is being updated as needed to be compatible with Python 3.6+. For the Python 2.7 compatible version of this package that was used in the aforementioned paper, please visit the [Ideker Lab GitHub](https://github.com/idekerlab/Network_Evaluation_Tools).

## Updates to this package
  - ```print``` statements fixed to be compatible with Python 3
  - Most usages of ```.ix``` in the packages updated to be compatible with pandas 0.20.0+
  - ```run_network_evaluation.py``` script in ```Network Evaluation Examples``` has been slightly updated
  - Jupyter Notebooks have *NOT* been updated to reflect changes to the code base at this time

## Modules in this package
  - _data_import_tools_ - This module contains functions for helping import network files and gene set files for analysis.
  - _gene_conversion_tools_ - This module contains functions for helping convert, filter, and save networks from their raw database form. Used in the Network Processing Jupyter Notebooks.
  - _miscellaneous_functions_ - This module contains various functions developed to help with analysis along the way. These functions are not well tested and may contain bugs. These functions were generally used to determine other network performance metrics on network recovery of gene sets.
  - _network_evaluation_functions_ - This module contains many of the core functions of the set-based network evaluation algorithm.
  - _network_propagation_ - This module contains functions to help with network propagation steps used in the set-based network evaluation algorithm.

## Version and Dendencies
Currently, the network_evaluation_tools package has been deployed in a Python 3.6 environment. However, assuming the dependencies are updated accordingly, there is an expectation that this package will be backwards compatible with Python 2.7. Additional testing is still required for other functions not involved in the command line script deployment of the package.
network_evaluation_tools requires: 
  - Argparse >= 1.1
  - NetworkX >= 2.3
  - Numpy >= 1.16.3
  - Matplotlib >= 3.0.3
  - Pandas >= 0.24.2
  - Requests >= 2.21.0
  - Scipy >= 1.2.1
  - Scikit-learn >= 0.20.3

## Installation
1. Clone the repository 
2. cd to new respository
3. Execute following command:  
```python setup.py install```

## Network analysis
1. If the network needs to be normalized to a particular naming scheme:<br>
A Jupyter Notebook describing how each network was processed from the raw download file in the original [paper](Link) can be found in the ```Network Processing Notebooks``` folder.<br>
2. There are two ways to perform the network evaluation on a gene set:<br>
The following network analyses can be performed either from a Jupyter Notebook or from the command line (see ```Network Evaluation Examples``` folder). Jupyter notebooks are documented within the notebook and the documentation for the python scripts can be seen using the command ```python [script_name].py -h```. <br>

## Data provided in this repository (see ```Data``` Folder)
 - Database Citations - An Excel file containing details about all of the networks used in the original paper's analysis and affiliated citations for all of the databases used.
 - _DisGeNET / Oncogenic Component Gene Sets_ - Two tab separated files, each line containing a gene set from either DisGeNET or the Oncogenic Component collection. The first column of each file is the name of the gene set followed by the list of genes associated with that given gene set on the same line.
 - _Network performance (AUPRCs) on DisGeNET / Oncogenic Component Gene Sets_ - Two csv files containing the raw Z-normalized AUPRC scores (network performance scores) of each network analyzed on each gene set analyzed from DisGeNET or the Oncogenic Component gene set collection.
 - _Network performance effect sizes on DisGeNET / Oncogenic Component Gene Sets_ - Two csv files containing the relative performance gain of each network's AUPRC score over the median null AUPRC score for each gene set analyzed from DisGeNET or the Oncogenic Component gene set collection.

## Issues
Please feel free to post issues/bug reports. Questions can be sent to jkh013@ucsd.edu

- In Pandas v0.20.0+, the ```.ix``` indexer has been deprecated. Usages of ```.ix``` in the package are being updated as needed. The package should still function normally despite the warnings currently.

## License
See the [LICENSE](https://github.com/huangger/Network_Evaluation_Tools/blob/master/LICENSE.txt) file for license rights and limitations (MIT).


