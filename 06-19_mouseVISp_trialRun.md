Single Nucleus Sequencing Trial Quality Control and Analysis
================

–

The starting point for the analysis is the filtered output form the
cellranger software. The software assigns all counts with the same
barcode to one cell. In addition, it discards reads from barcodes that
occur only a few number of times, as these likely correspond to read
errors and not cells.

## Quality Control

I did not remove any more cells based on low numbers of detected genes.
I also did not remove cells as potential doublets based on a high
numbers of detected genes alone.