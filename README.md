# TaxonTableTools

## Background

TaxonTableTools (TTT) aims to provide easy-to-use tools for biologists and non-bioinformaticians to analyse and visualize their metabarcoding data quickly and reproducible via a graphical user interface.

TaxonTableTools is an evolving software and there will be bugs and issues at few points. If so, please leave the report in the git repository or drop me an email. Furthermore, new content and functions will be gradualy added. Suggestions and recommodations for new features are always welcome!

## Requirements

* Python version 3.6 or 3.7
* Pip
* OS independent (tested on Windows 10, Ubuntu 16 & 18 and MacOS Catalina)
  * Installation failed on some Windows machines, but could be fixed by a second installation attempt
  * In other cases some Windows machines required the installation of C++ Build Tools

## Installation

TaxonTableTools requires Python version 3.6 and 3.7 and can be easily installed using pip.

First, make sure you run the correct pip version via:

`pip3 --version`

Which should return python3.6 or 3.7. Otherwise specify your pip using `pip3.6` or `pip3.7`

Then install TaxonTableTools via pip:

`pip3 install taxontabletools`

TaxonTableTools can then be started via:  

`python3 -m taxontabletools`

Updates can be installed via:

`pip3 install --upgrade taxontabletools`

The installation of [Krona tools](https://github.com/marbl/Krona/wiki) is optional and needs to be done separately. Note that Krona tools is currently not supported on Windows.

## Graphical user interface
![](/_data/ttt_gui.png "figure-1")

## Examples
![](/_data/ttt_figure_1.png "figure-1") ![](/_data/ttt_figure_2.png "figure-2")

## Quick start

* When first launched TTT will ask to define an output directory. This is where all your projects and respective output files will be stored. A new folder "Projects" will be created in this directory.

* The next window will always appear when launching TTT. Here you can create new projects or load an already existing one. The output directory can also be changed here.

* Create a new project by typing the desired name or leave blank to create a "Default_project". Click on "Create new" to proceed.

* Initially, the taxonomy table and read table must be converted to a TaXon table, the standard input format for TTT (see [manual](_data/TaxonTableTools_manual.pdf) for definitions and requirements). A template TaXon table is available [here](_data/template_taxon_table.xlsx).

* Check out the [tutorial](_data/TaxonTableTools_tutorial.pdf) to get a more detailed impression of the TTT workflow.
