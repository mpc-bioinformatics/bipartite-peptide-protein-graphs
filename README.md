# Bipartite peptide-protein-graphs

In bottom-up proteomics, intact proteins are digested to peptides via trypsin
before measurement by mass spectrometry based techniques. The relationship 
between peptides and proteins can naturally be presented as bipartite graphs, 
where there is one node type for proteins and one for peptides. An edge is drawn
if the peptide belongs to the peptide. Analysing the structure of these bipartite
graphs give a hint about the difficulty of protein quantification, especially
in presence of shared peptides. 

This R code belongs to the following manuscript:

Schork K, Turewicz M, Uszkoreit J, Rahnenführer J, Eisenacher M (2022) Characterization of peptide-protein relationships in protein ambiguity groups via bipartite graphs. PLOS ONE 17(10): e0276401. https://doi.org/10.1371/journal.pone.0276401 

The different R scripts have to be executed in the numbered order (00 to 05).

The "Figures_and_Tables_for_Paper" folder contains Code for alle tables and figures 
of the above mentioned manuscript.

The data underlying this analysis are available in the ProteomXChange Consortium 
at  http://proteomecentral.proteomexchange.org via the PRIDE partner repository, 
and can be accessed with the dataset identifiers PXD024684 and PXD030577 and PXD030603.

Currently, development is continued in the bppg package. The current development version can be found here:
https://github.com/mpc-bioinformatics/bppg
