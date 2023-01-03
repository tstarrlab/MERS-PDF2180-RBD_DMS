---
layout: heatmaps
permalink: /RBD-heatmaps/
---

---

*NOTE data are preliminary*

### Overview

You can use this tool to explore the experimentally determined impacts of amino acid mutations on DPP4-binding affinity and expression in MERS and PDF-2180 receptor-binding domain (RBD) variants. 

#### Instruction

To use this tool, select the RBD variants that you wish to display in each heatmap by selecting a variant from the drop down menu corresponding to each plot. Then, you can select the metric that you wish to display in the heatmap (either change in DPP4 binding affinity (-log10 $$K_D$$), or change in RBD expression (log(MFI)) by selecting that metric in the corresponding drop down menu. Hover over individual mutations to see exact numerical details.

#### Technical Details

The impact on DPP4 receptor-binding affinity ($$\Delta$$ log10 $$K_D$$) or RBD expression ($$\Delta$$ log(MFI)) of every single amino-acid mutation in MERS and PDF-2180 CoV RBDs, as determined by high-throughput titration assays. Wildtype amino acids are indicated by an 'x', and gray squares indicate missing mutations from each library. The number of internally replicated barcodes with which a mutation was measured is visible as `Barcode Count` in the tooltips, where higher numbers indicate higher-confidence measurements.

### Data

Raw data  can be found [here](https://github.com/jbloomlab/MERS-PDF2180-RBD_DMS/blob/main/results/final_variant_scores/final_variant_scores.csv). The code used to make these plots can be found [here](https://github.com/jbloomlab/MERS-PDF2180-RBD_DMS/blob/main/RBD-Heatmaps-Interactive-Visualization.ipynb).
