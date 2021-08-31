# sc_immunorepertoire_annotator

This script in python3 cleans the BCR and TCR data from all_contig_annotation.csv provided by 10XGenomics.

It will create a more compact .csv file :

    - Only full length (full_length == True)
    - Only cells (is_cell == True)
    - One row = one cell.
   
