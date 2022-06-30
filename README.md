# TaxonTableTools

[![Downloads](https://pepy.tech/badge/taxontabletools)](https://pepy.tech/project/taxontabletools)

## Background

TaxonTableTools (TTT) aims to provide easy to use tools for biologists and non-bioinformaticians to analyse and visualize their DNA metabarcoding data quickly and reproducible via a graphical user interface.

TaxonTableTools is an evolving software and there will be bugs and issues at few points. If so, please leave the report in the git repository or drop me an email. Furthermore, new content and functions will be gradually added. Suggestions and recommendations for new features are always welcome!

## Requirements

* Python version 3.6 or 3.7 (should also work with 3.8 and 3.9)
* Pip python package installer
* OS independent (tested on Windows 10, Ubuntu 16 & 18 and MacOS Catalina)
  * Some issues were reported during the [installation on Windows](https://github.com/TillMacher/TaxonTableTools/issues/1)
  * Some issues were reported during the [installation on MacOS](https://github.com/TillMacher/TaxonTableTools/issues/4)

## Installation

TaxonTableTools requires Python version 3.6 or 3.7 and can be easily installed using pip.

First, make sure you run the correct pip version via:

`python3.9 -m pip --version`

Which should return "pip 21.3.1 from /Library/Frameworks/Python.framework/Versions/3.9/lib/python3.9/site-packages/pip (python 3.9)" or similar.

Then install numpy and TaxonTableTools via pip:

`python3.9 -m pip install numpy`

`python3.9 -m pip install taxontabletools`

TaxonTableTools can then be started via:  

`python3.9 -m taxontabletools`

Updates can be installed via:

`python3.9 -m pip install --upgrade taxontabletools`

The installation of [Krona tools](https://github.com/marbl/Krona/wiki) is optional and needs to be done separately. Note that Krona tools is currently not supported on Windows.

## Tutorial

The easiest way to learn to how to use TaxonTableTools is to check out the [tutorial](_data/TaxonTableTools_tutorial.pdf)!

## Quick start

* When first launched TTT will ask to define an output directory. This is where all your projects and respective output files will be stored. A new folder "Projects" will be created in this directory.

* The next window will always appear when launching TTT. Here you can create new projects or load an already existing one. The output directory can also be changed here.

* Create a new project by typing the desired name or leave blank to create a "Default_project". Click on "Create new" to proceed.

* A blank template TaXon table is available [here](_data/template_taxon_table.xlsx).

## Graphical user interface
![](/_data/ttt_gui.png "figure*1")

## Examples
![](/_data/ttt_figure_1.png "figure*1") ![](/_data/ttt_figure_2.png "figure*2")

## Ecological traits

OTU-specific traits can now be added to the TaXon table. Traits can e.g. be:
- habitat
- diet type
- distribution
- ecological traits

*Requirements for trait source tables (.xlsx):*
- Column A: One of the TTT-specific taxonomic level (Phylum, Class, Order, Family, Genus, Species)
- Column B: Name of the trait

<details><summary> Example for aquatic invertebrates </summary>

![image](https://user-images.githubusercontent.com/48299746/174554773-2aa8e5d5-3c39-4285-a0a0-7b1997bc472d.png)

</details>

*Requirements for manual import:*
- The traits columns must be located between the 'Status' and 'seq' columns!
- Do not alter the location of the other TaXon table columns.

There is no limit of OTU traits columns that can be added. A TaXon table containing OTU traits is compatible with
all modules (OTU traits will simply be ignored if not required).

<details><summary> Example for aquatic invertebrates </summary>

![image](https://user-images.githubusercontent.com/48299746/174555954-a495755c-723f-40c4-9051-fe879242639a.png)

</details>

## Change log

### v 1.4.5
  TTT change log v 1.4.5

  * Added trait import to the data conversion modules. 


### v 1.4.4
  TTT change log v 1.4.4

  Bug fixes
  * Fixed incorrect text boxes.

  NMDS
  * Changed back NMDS to raw stress until sclearn updates to normalized stress.
  * TTT will be updated as soon as normalized stress is available!

### v 1.3.1
  * Metadata entries for specific samples can now be left blank to exclude the sample from the analyses. This avoids the need to create new metadata tables for each subset of a TaXon table
  * Taxon lists have been overhauled and now include presence/absence indications for samples or metadata groups
  * Fixed incorrect y-axis title in the per taxon plot
  * Enhanced y-axis title of shared OTUs per bin plot
  * Fixed Qiime2 conversion
  * Fixed OTU and read correlation layout bugs


### v 1.3.0
  * Enhanced TaXon table check. Missing taxonomy is now marked with a flag
  * User preferences can now be saved
  * User will now be informed about available updates
  * User can choose what the clustering units (e.g. OTUs, ASVs) shall be called
  * Genus and species names are displayed in italic
  * Added additional replicate analysis plot
  * TaXon table filter tools print the number of removed OTUs
  * Metadata tables can be opened in TTT
  * Added new tool to subtract reads present in negative controls from all samples
  * Added option to distinguish different taxa in the rarefaction plot
  * Added site occupancy heatmap option
  * Added parallel category chart
  * Added option to display single samples in one Krona chart
  * Added non.metric multidimensional scaling (NMDS) analysis
  * Added bray-curtis option to diversity and ordination analyses
  * Added vernacular names download in the taxon list
  * Log files can be opened in TTT
  * Several layout changes
  * Various minor bug fixes

 ### v. 1.2.4
 * Fixed "TaXon table consuming" bug on some Windows 10 machines

 ### v. 1.2.3
 * Fixed alpha diversity bug
 * Missing file in installation path fixed.
 * Corrected minor spelling and layout mistakes.

 ### v. 1.2.0
 Major update:
 * The identifier "IDs" was changed to "ID" in the TaXon table header!
  * Please adjust old TaXon tables accordingly!
  * Otherwise TTT will crash!

 New features:
 * Metadatables can now directly be adjusted from within TTT.
 * Browsing TaXon tables will directly show the current project folder.
 * The OTU prompt is more relaxed. Hash OTUs are accepted now.
 * Qiime2 feature tables can now be converted to the TTT format.
 * Additional options for customizing colors have been added.
 * Font sizes can now be adjusted.
 * Diversity and ordiniation analyses can now be performed on different taxonomic levels.
 * Diversity and ordiniation analyses automatically transform data to boolean values. Thus, also non presence/absence converted tables can be used.
 * Moved various options to drop down menus.

 Bug fixes:
 * Fixed "shared OTUs" tool. More than 3 replicates are now supported.
 * Changed x*axis value to read proportions in the per taxon analysis.
 * Added missing package 'statsmodels' to installation.

## How to cite

If you use TTT:
* Macher, T.*H., Beermann, A. J., & Leese, F. (2021). TaxonTableTools—A comprehensive, platform*independent graphical user interface software to explore and visualise DNA metabarcoding data. Molecular Ecology Resources. doi: https://doi.org/10.1111/1755*0998.13358

If you create Krona charts, please also cite:
* Ondov, B. D., Bergman, N. H., & Phillippy, A. M. (2011). Interactive metagenomic visualization in a Web browser. BMC Bioinformatics, 12(1), 385. doi: 10.1186/1471*2105*12*385
