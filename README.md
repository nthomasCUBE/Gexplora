# chromoWIZv2 (under development)

plotting the gene distribution along the Brachypodium distachyon genome.

It allows to use a GTF file as input and then to display mRNA, exon (in the second
column) and to adapt the maximum threshold per BIN (using 200 non-overlapping bins respective
to largest chromosome).

Visualizing all genes from an input GTF file and Ensembl Webservice calls to retrieve sequence for a
gene as well as gene tree for a particular human gene.
![chromoWIZv2](https://github.com/nthomasCUBE/chromoWIZv2/blob/master/chromoWIZv2i_1.png)

Or searching for candidate genes by uploading the candidate genes

![chromoWIZv2](https://github.com/nthomasCUBE/chromoWIZv2/blob/master/chromoWIZv2i_2.png)

The GTF can be taken from Ensembl and then any of the element types can be selected and visualized.
![chromoWIZv2](https://github.com/nthomasCUBE/chromoWIZv2/blob/master/chromoWIZv2i_3.png)

chromoWIZ as a graphical user interface

## Ensembl Plant

Ensembl provide a bunch of different endpoints to obtain data. We use some of these endpoints 
in chromoWIZpy (https://rest.ensembl.org/).

## Additional required Python Packages

In chromoWIZpy, we use tkinter for the graphical user interface, while
we use numpy to be faster in calculating overlaps between elements of two arrays while
requests is used to obtain data from endpoints of Ensembl.

- numpy
- requests
- tkinter

For installation under Windows, you can use the "pip.exe" that can be normally
found in the subdirectory "scripts" within the python installation.
You can install Python packages using ``pip.exe install requests''-




