# Gexplora

## How to start

* Download Gexplora_V1.rar, unzip it.
* Navigate to Gexplora_V1/Gexplora
* Execute it with double-click: Gexplora.exe

## Description

plotting the density of genetic elements along a chromosome using GTF files (ideally from Ensembl) as input
and using REST API calls to Ensembl, OMA and STRING database (to be done).

It allows to use a GTF file as input and then to display mRNA, exon (in the second
column). Maximum threshold per BIN from 200 non-overlapping bins respective
to largest chromosomes defines highest intensity.

Visualizing all genes from an input GTF file and Ensembl Webservice calls to retrieve sequence for a
gene as well as gene tree for a particular human gene.
![Gexplora](https://github.com/nthomasCUBE/Gexplora/blob/master/pix/fig1A.png)

Or searching for candidate genes by uploading the candidate genes

![Gexplora](https://github.com/nthomasCUBE/Gexplora/blob/master/pix/fig1B.png)

The GTF can be taken from Ensembl and then any of the element types can be selected and visualized.
![Gexplora](https://github.com/nthomasCUBE/Gexplora/blob/master/pix/fig1C.png)

## Ensembl Plant

Ensembl provide a bunch of different endpoints to obtain data. We use some of these endpoints 
in Gexplora (https://rest.ensembl.org/).

## Use case - barley genome

The annotation was downloaded from Ensembl ![Ensembl-Barley](ftp://ftp.ensemblgenomes.org/pub/plants/release-46/gff3/hordeum_vulgare).
and gene distribution was shown while also all other types can be shown in the density visualisation.
![Gexplora](https://github.com/nthomasCUBE/Gexplora/blob/master/pix/Gexplora_Barley.png)

## Additional required Python Packages

In Gexplora, we use tkinter for the graphical user interface, while
we use numpy to be faster in calculating overlaps between elements of two arrays while
requests is used to obtain data from endpoints of Ensembl.

- numpy
- requests
- tkinter
- json

For installation under Windows, you can use the "pip.exe" that can be normally
found in the subdirectory "scripts" within the python installation.
You can install Python packages using ``pip.exe install requests''-


##

To generate the executable to run:
pip install https://github.com/pyinstaller/pyinstaller/archive/develop.zip





