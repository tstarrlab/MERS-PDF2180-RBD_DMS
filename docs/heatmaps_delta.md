---
layout: heatmaps_delta
permalink: /RBD-heatmaps_delta/
---

---

*NOTE data are preliminary*

### Overview

You can use this tool to explore the experimentally determined impacts of amino acid mutations on hDPP4 binding affinity (delta-log10Ka) and expression (delta-logMFI) in merbecovirus receptor-binding domains (RBD). 

#### Instruction

To use this tool, select the merbecovirus variants and the metrics that you wish to display in the heatmap (change in hDPP4 binding affinity (delta-log10Ka), or change in mean fluorescence intensity (delta-MFI)) for RBD expression, by selecting that metric in the corresponding drop down menu. Hover over individual mutations to see exact numerical details. Click and drag on the site zoom bar to make it easier to scroll along the linear sequence.

*Note that sites are plotted according to aligned MERS-CoV spike indexing for both backgrounds, instead of their own indices (self numbering in tooltips). Columns also do not align perfectly (though tiptools/hover cursor does align to corresponding site)*

#### Technical Details

The impact on hDPP4 receptor-binding and expression level of every single amino-acid mutation in each RBD, as determined by high-throughput FACS-seq assays. Wildtype amino acids are indicated by an 'x', and gray squares indicate missing mutations from each library. The number of internally replicated barcodes with which a mutation was measured is visible as `barcode count` in the tooltips, where higher numbers indicate higher-confidence measurements.


### Data

Raw data  can be found [here](https://github.com/jbloomlab/MERS-PDF2180-RBD_DMS_/blob/main/results/final_variant_scores/final_variant_scores.csv). The code used to make these plots can be found [here](https://github.com/jbloomlab/MERS-PDF2180-RBD_DMS_/blob/main/RBD-Heatmaps-Interactive-Visualization_delta.ipynb).
