# RBDmap Lite

## An R Shiny app for detecting RNA-binding sites within proteins

![UI Screenshot](screenshots/rbdmap-lite-screenshot.png)

RBDmap Lite is an R Shiny app used to identify RNA-binding sites within proteins from LC-MS/MS data after performing cross-linking with mass spectrometry (XL-MS) using the RBDmap method. For more information about the RBDmap workflow and how it was previously used to directly identify an RNA-binding binding site in an unbiased manner within the PRC2 protein complex, please see the references at the bottom.

# Necessary input files:
* A .fasta file of your protein or protein complex of interest (included example: "prc2_5m.fasta")
* .txt files contained tab-separated data generated from analysing raw LC-MS/MS output files using MaxQuant (included example files in "LS-MSMS-data" directory). Files must have the following keywords in the file name: "UV", "ArgC" or "LysC", and a digit based on the repeat number (e.g. "Repeat01_UV_ArgC_Eluate.txt")

# Output generated:
* A graph (or graphs if a protein complex) based on the number of hits per repeat of RBDmap. Each individual subunit within a protein complex can be selected via a drop-down menu.
* A colour-coded sequence of the protein of interest, where the number of hits of residues identified to bind RNA are colour-coded.

Functionally, RBDmap Lite is much more simple and humble version of the crisscrosslinkeR R package by [egmg726](https://github.com/egmg726/crisscrosslinker) [2]. If you're interested in analysing RNA-protein and protein-protein interactions, it would be good to look into that package. RBDmap Lite will look at RNA-protein binding interactions only and also doesn't automatically generate PyMOL files with residues coloured in as with crisscrosslinkeR. However, RBDmap Lite does provide a very user-friendly UI in the form of an online app, as well as a coloured version of the protein sequence with residue letters.

# References:
1. Zhang, Q., McKenzie, N.J., Warneford-Thomson, R., Gail, E.H., Flanigan, S.F., Owen, B.M., Lauman, R., Levina, V., Garcia, B.A., Schittenhelm, R.B. and Bonasio, R., 2019. RNA exploits an exposed regulatory site to inhibit the enzymatic activity of PRC2. Nature structural & molecular biology, 26(3), pp.237-247.
1. Gail, E.H., Shah, A.D., Schittenhelm, R.B. and Davidovich, C., 2020. crisscrosslinkeR: identification and visualization of protein–RNA and protein–protein interactions from crosslinking mass spectrometry. Bioinformatics.