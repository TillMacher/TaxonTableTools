# TaXonTableTools

## Installation

TaXonTableTools requires Python version 3.6 and 3.7 and can be easily installed using pip.

First, make sure you run the correct pip version via:

`pip3 --version`

Which should return python3.6 oder 3.7. Otherwise specify your pip using `pip3.6` or `pip3.7`

Then install TaXonTableTools via pip:

`pip3 install taxontabletools`

TaXonTableTools can then be started via:  

`python3 -m taxontabletools`

Updates can be installed via:

`pip3 install --upgrade taxontabletools`


## Introduction

The sequencing of DNA metabarcoding data has drastically increased over the past decade and many datasets are produced quickly nowadays. The analysis of this massive amounts of data and their translation into biological meaningful facts is often limited especially for non-experts and bioinformatics beginners. However, it is the biologists that need to work with the data and interpret these. This was the rationale for developing TaXonTableTools (TTT is the following) as part of the GeDNA project (eDNA metabarcoding in regulatory biomonitoring in Germany). The program aims to provide easy-to-use tools for biologists and non-bioinformaticians to analyse and visualize their data quickly and reproducible via a graphical user interface. Thus, the dependency on self-written R or python scripts, which can cause confusion and errors particularly when working with different datasets, is reduced. TaXonTableTools is not aiming to replace those scripts (since a specific script always outperforms a general script), but rather provide tools to quickly assess data and generate information as basis for further, more dataset-specific analyses. Furthermore, the quick data visualization integrated in TTT always comes in handy for presenting first preliminary results of a dataset.

## Manual

Find more information of how to use TaXonTableTools in the [manual](_data/TaXonTableTools_manual.docx) and the [tutorial](_data/TaXonTableTools_tutorial.pptx).
