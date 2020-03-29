# Gexplorr (under development)

## How to start

download Gexplorr_V1.rar, unzip.
Navigate to Gexplorr_V1/Gexplorr
then execute it with double-click: Gexplorr.exe

## Description

plotting the density of genetic elements along a chromsome using GTF files (ideally from Ensembl) as input
and using REST API calls to Ensembl, OMA and StringDB (to be done).

It allows to use a GTF file as input and then to display mRNA, exon (in the second
column). Maximum threshold per BIN from 200 non-overlapping bins respective
to largest chromosomes defines highest intensity.

Visualizing all genes from an input GTF file and Ensembl Webservice calls to retrieve sequence for a
gene as well as gene tree for a particular human gene.
![Gexplorr](https://github.com/nthomasCUBE/Gexplorr/blob/master/pix/fig1A.png)

Or searching for candidate genes by uploading the candidate genes

![Gexplorr](https://github.com/nthomasCUBE/Gexplorr/blob/master/pix/fig1B.png)

The GTF can be taken from Ensembl and then any of the element types can be selected and visualized.
![Gexplorr](https://github.com/nthomasCUBE/Gexplorr/blob/master/pix/fig1C.png)

Gexplorr as a graphical user interface

## Ensembl Plant

Ensembl provide a bunch of different endpoints to obtain data. We use some of these endpoints 
in Gexplorr (https://rest.ensembl.org/).

## Additional required Python Packages

In Gexplorr, we use tkinter for the graphical user interface, while
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






