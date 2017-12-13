RANKS: Robust Analytics and Normalization for Knockout Screens

Calculates gene scores, p-values and FDR values for CRISPR knockout screens.

Input files should be tab-delimited text files.

Read count file format: 
First column: read count; Second column: sgRNA ID

sgRNA to gene mapping file format:
First column: sgRNA ID; Subsequent column(s): gene/feature ID(s)

Control sgRNA file format
First column: sgRNA ID


Before running RANKS, you must first run once in the same folder the program which generates the control gene scores by typing:
	perl control-distribution.pl


Example:
perl control-distribution.pl        #Note: you only need to run this command once.

perl ranks.pl day0reads day15reads -lib mapping -ctrl control-sgrnas -d > output.tsv




This new version of RANKS (V2.0) produces a single score per gene which captures both sgRNA depletion (negative scores) and sgRNA enrichment (positive scores) instead of producing a separate depletion and enrichment score. It is important to run the new version of the control-distribution.pl script before using the new version of RANKS because the score distribution is different.
