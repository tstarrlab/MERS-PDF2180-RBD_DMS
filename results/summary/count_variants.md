# Count variants in each sample
This Python Jupyter notebook counts occurrences of each barcode in each sample from Illumina barcode sequencing, and adds these counts to the codon variant table.

## Set up analysis
### Import Python modules.
Use [plotnine](https://plotnine.readthedocs.io/en/stable/) for ggplot2-like plotting.

The analysis relies heavily on the Bloom lab's [dms_variants](https://jbloomlab.github.io/dms_variants) package:


```python
import itertools
import multiprocessing
import multiprocessing.pool
import os
import warnings

import alignparse
import alignparse.targets

import dms_variants.codonvarianttable
from dms_variants.constants import CBPALETTE
import dms_variants.illuminabarcodeparser
import dms_variants.utils
import dms_variants.plotnine_themes

from IPython.display import display, HTML

import pandas as pd

from plotnine import *

import yaml
```

Set [plotnine](https://plotnine.readthedocs.io/en/stable/) theme to the gray-grid one defined in `dms_variants`:


```python
theme_set(dms_variants.plotnine_themes.theme_graygrid())
```

Versions of key software:


```python
print(f"Using alignparse version {alignparse.__version__}")
print(f"Using dms_variants version {dms_variants.__version__}")
```

    Using alignparse version 0.6.0
    Using dms_variants version 1.4.2


Ignore warnings that clutter output:


```python
warnings.simplefilter('ignore')
```

### Parameters for notebook
Read the configuration file:


```python
with open('config.yaml') as f:
    config = yaml.safe_load(f)
```

Make output directory if needed:


```python
os.makedirs(config['counts_dir'], exist_ok=True)
os.makedirs(config['figs_dir'], exist_ok=True)
```

## Input variant tables
Initialize the table of barcode-variant pairs from the respective `process_ccs` notebooks for each background.


```python
variants = pd.read_csv(config['codon_variant_table_file_MERS'], na_filter=None)
variants = variants.append(pd.read_csv(config['codon_variant_table_file_PDF2180'], na_filter=None))

variants = variants.reset_index(drop=True)

display(HTML(variants.head().to_html(index=False)))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>target</th>
      <th>library</th>
      <th>barcode</th>
      <th>variant_call_support</th>
      <th>codon_substitutions</th>
      <th>aa_substitutions</th>
      <th>n_codon_substitutions</th>
      <th>n_aa_substitutions</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>MERS</td>
      <td>lib51_53</td>
      <td>AAAAAAAAACATTCGT</td>
      <td>1</td>
      <td>AAC145TTG</td>
      <td>N145L</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <td>MERS</td>
      <td>lib51_53</td>
      <td>AAAAAAAAAGACTTTC</td>
      <td>1</td>
      <td>GTA158AAA</td>
      <td>V158K</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <td>MERS</td>
      <td>lib51_53</td>
      <td>AAAAAAAAAGCGATAG</td>
      <td>4</td>
      <td>TCC128---</td>
      <td>S128-</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <td>MERS</td>
      <td>lib51_53</td>
      <td>AAAAAAAAATTGAGGT</td>
      <td>2</td>
      <td>CCT139CAA</td>
      <td>P139Q</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <td>MERS</td>
      <td>lib51_53</td>
      <td>AAAAAAAACAATCCCG</td>
      <td>4</td>
      <td>ACA188TAT</td>
      <td>T188Y</td>
      <td>1</td>
      <td>1</td>
    </tr>
  </tbody>
</table>


Are there any barcodes in the same library that are shared across targets?
If so, we need to get rid of those as they will be confounded in barcode parsing:


```python
dup_barcodes = (
    variants
    .groupby(['library', 'barcode'])
    .size()
    .rename('duplicate_count')
    .reset_index()
    .query('duplicate_count > 1')
    )

print('Here are duplicated barcodes:')
display(HTML(dup_barcodes.head().to_html(index=False)))

print(f"\nRemoving the {len(dup_barcodes)} duplicated barcodes."
      f"Started with {len(variants)} barcodes:")
variants = (
    variants
    .merge(dup_barcodes, on=['library', 'barcode'], how='outer')
    .query('duplicate_count.isnull()', engine='python')
    )
print(f"After removing duplicates, there are {len(variants)} barcodes.")
```

    Here are duplicated barcodes:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>library</th>
      <th>barcode</th>
      <th>duplicate_count</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>lib51_53</td>
      <td>AAACTTCTGAGAAAAA</td>
      <td>2</td>
    </tr>
    <tr>
      <td>lib51_53</td>
      <td>GCTCAATATATATATC</td>
      <td>2</td>
    </tr>
    <tr>
      <td>lib51_53</td>
      <td>TCTGTATAGACTTAAC</td>
      <td>2</td>
    </tr>
    <tr>
      <td>lib52_54</td>
      <td>ACAAAGAATTCAAATT</td>
      <td>2</td>
    </tr>
    <tr>
      <td>lib52_54</td>
      <td>GGACTAATGAACCATG</td>
      <td>2</td>
    </tr>
  </tbody>
</table>


    
    Removing the 5 duplicated barcodes.Started with 245760 barcodes:
    After removing duplicates, there are 245750 barcodes.


Pull out a target sequence for matching to the barcode and flanking sequence regions. Note, in this pipeline this is ok because our different backgrounds don't have differing flanks or other features within the actual N16 region covered in Illumina sequencing. If ever placing in-line barcodes here in the future, we would need to modify this.


```python
# get wildtype gene sequence for primary target
targets = alignparse.targets.Targets(seqsfile=config['amplicons_MERS'],
                                     feature_parse_specs=config['feature_parse_specs_MERS'])
```

## Setup to parse barcodes
Read data frame with list of all barcode runs.


```python
# barcode runs with R1 files by semicolon string split
barcode_runs = (pd.read_csv(config['barcode_runs'])
                .assign(R1=lambda x: x['R1'].str.split('; '))
                )

display(HTML(barcode_runs.to_html(index=False)))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>library</th>
      <th>replicate</th>
      <th>sample</th>
      <th>sample_type</th>
      <th>sort_bin</th>
      <th>concentration</th>
      <th>date</th>
      <th>number_cells</th>
      <th>R1</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>lib51_53</td>
      <td>A</td>
      <td>TiteSeq_hDPP4_01_bin1</td>
      <td>TiteSeq_hDPP4</td>
      <td>1</td>
      <td>1.0</td>
      <td>221213</td>
      <td>2979526</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221213_MERS-PDF2180_bc-seq/221213_s1_b1_S1_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib51_53</td>
      <td>A</td>
      <td>TiteSeq_hDPP4_01_bin2</td>
      <td>TiteSeq_hDPP4</td>
      <td>2</td>
      <td>1.0</td>
      <td>221213</td>
      <td>1403409</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221213_MERS-PDF2180_bc-seq/221213_s1_b2_S2_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib51_53</td>
      <td>A</td>
      <td>TiteSeq_hDPP4_01_bin3</td>
      <td>TiteSeq_hDPP4</td>
      <td>3</td>
      <td>1.0</td>
      <td>221213</td>
      <td>1288376</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221213_MERS-PDF2180_bc-seq/221213_s1_b3_S3_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib51_53</td>
      <td>A</td>
      <td>TiteSeq_hDPP4_01_bin4</td>
      <td>TiteSeq_hDPP4</td>
      <td>4</td>
      <td>1.0</td>
      <td>221213</td>
      <td>1262113</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221213_MERS-PDF2180_bc-seq/221213_s1_b4_S4_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib51_53</td>
      <td>A</td>
      <td>TiteSeq_hDPP4_02_bin1</td>
      <td>TiteSeq_hDPP4</td>
      <td>1</td>
      <td>2.0</td>
      <td>221213</td>
      <td>3407489</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221213_MERS-PDF2180_bc-seq/221213_s2_b1_S5_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib51_53</td>
      <td>A</td>
      <td>TiteSeq_hDPP4_02_bin2</td>
      <td>TiteSeq_hDPP4</td>
      <td>2</td>
      <td>2.0</td>
      <td>221213</td>
      <td>1526544</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221213_MERS-PDF2180_bc-seq/221213_s2_b2_S6_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib51_53</td>
      <td>A</td>
      <td>TiteSeq_hDPP4_02_bin3</td>
      <td>TiteSeq_hDPP4</td>
      <td>3</td>
      <td>2.0</td>
      <td>221213</td>
      <td>1153617</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221213_MERS-PDF2180_bc-seq/221213_s2_b3_S7_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib51_53</td>
      <td>A</td>
      <td>TiteSeq_hDPP4_02_bin4</td>
      <td>TiteSeq_hDPP4</td>
      <td>4</td>
      <td>2.0</td>
      <td>221213</td>
      <td>1004697</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221213_MERS-PDF2180_bc-seq/221213_s2_b4_S8_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib51_53</td>
      <td>A</td>
      <td>TiteSeq_hDPP4_03_bin1</td>
      <td>TiteSeq_hDPP4</td>
      <td>1</td>
      <td>3.0</td>
      <td>221213</td>
      <td>4374876</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221213_MERS-PDF2180_bc-seq/221213_s3_b1_S9_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib51_53</td>
      <td>A</td>
      <td>TiteSeq_hDPP4_03_bin2</td>
      <td>TiteSeq_hDPP4</td>
      <td>2</td>
      <td>3.0</td>
      <td>221213</td>
      <td>1178077</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221213_MERS-PDF2180_bc-seq/221213_s3_b2_S10_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib51_53</td>
      <td>A</td>
      <td>TiteSeq_hDPP4_03_bin3</td>
      <td>TiteSeq_hDPP4</td>
      <td>3</td>
      <td>3.0</td>
      <td>221213</td>
      <td>1107868</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221213_MERS-PDF2180_bc-seq/221213_s3_b3_S11_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib51_53</td>
      <td>A</td>
      <td>TiteSeq_hDPP4_03_bin4</td>
      <td>TiteSeq_hDPP4</td>
      <td>4</td>
      <td>3.0</td>
      <td>221213</td>
      <td>719585</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221213_MERS-PDF2180_bc-seq/221213_s3_b4_S12_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib51_53</td>
      <td>A</td>
      <td>TiteSeq_hDPP4_04_bin1</td>
      <td>TiteSeq_hDPP4</td>
      <td>1</td>
      <td>4.0</td>
      <td>221213</td>
      <td>4615914</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221213_MERS-PDF2180_bc-seq/221213_s4_b1_S13_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib51_53</td>
      <td>A</td>
      <td>TiteSeq_hDPP4_04_bin2</td>
      <td>TiteSeq_hDPP4</td>
      <td>2</td>
      <td>4.0</td>
      <td>221213</td>
      <td>1149682</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221213_MERS-PDF2180_bc-seq/221213_s4_b2_S14_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib51_53</td>
      <td>A</td>
      <td>TiteSeq_hDPP4_04_bin3</td>
      <td>TiteSeq_hDPP4</td>
      <td>3</td>
      <td>4.0</td>
      <td>221213</td>
      <td>710331</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221213_MERS-PDF2180_bc-seq/221213_s4_b3_S15_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib51_53</td>
      <td>A</td>
      <td>TiteSeq_hDPP4_04_bin4</td>
      <td>TiteSeq_hDPP4</td>
      <td>4</td>
      <td>4.0</td>
      <td>221213</td>
      <td>123303</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221213_MERS-PDF2180_bc-seq/221213_s4_b4_S16_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib51_53</td>
      <td>A</td>
      <td>TiteSeq_hDPP4_05_bin1</td>
      <td>TiteSeq_hDPP4</td>
      <td>1</td>
      <td>5.0</td>
      <td>221213</td>
      <td>5269498</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221213_MERS-PDF2180_bc-seq/221213_s5_b1_S17_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib51_53</td>
      <td>A</td>
      <td>TiteSeq_hDPP4_05_bin2</td>
      <td>TiteSeq_hDPP4</td>
      <td>2</td>
      <td>5.0</td>
      <td>221213</td>
      <td>875320</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221213_MERS-PDF2180_bc-seq/221213_s5_b2_S18_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib51_53</td>
      <td>A</td>
      <td>TiteSeq_hDPP4_05_bin3</td>
      <td>TiteSeq_hDPP4</td>
      <td>3</td>
      <td>5.0</td>
      <td>221213</td>
      <td>133991</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221213_MERS-PDF2180_bc-seq/221213_s5_b3_S19_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib51_53</td>
      <td>A</td>
      <td>TiteSeq_hDPP4_05_bin4</td>
      <td>TiteSeq_hDPP4</td>
      <td>4</td>
      <td>5.0</td>
      <td>221213</td>
      <td>4635</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221213_MERS-PDF2180_bc-seq/221213_s5_b4_S20_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib51_53</td>
      <td>A</td>
      <td>TiteSeq_hDPP4_06_bin1</td>
      <td>TiteSeq_hDPP4</td>
      <td>1</td>
      <td>6.0</td>
      <td>221213</td>
      <td>5868162</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221213_MERS-PDF2180_bc-seq/221213_s6_b1_S21_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib51_53</td>
      <td>A</td>
      <td>TiteSeq_hDPP4_06_bin2</td>
      <td>TiteSeq_hDPP4</td>
      <td>2</td>
      <td>6.0</td>
      <td>221213</td>
      <td>311628</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221213_MERS-PDF2180_bc-seq/221213_s6_b2_S22_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib51_53</td>
      <td>A</td>
      <td>TiteSeq_hDPP4_06_bin3</td>
      <td>TiteSeq_hDPP4</td>
      <td>3</td>
      <td>6.0</td>
      <td>221213</td>
      <td>5760</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221213_MERS-PDF2180_bc-seq/221213_s6_b3_S23_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib51_53</td>
      <td>A</td>
      <td>TiteSeq_hDPP4_06_bin4</td>
      <td>TiteSeq_hDPP4</td>
      <td>4</td>
      <td>6.0</td>
      <td>221213</td>
      <td>586</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221213_MERS-PDF2180_bc-seq/221213_s6_b4_S24_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib51_53</td>
      <td>A</td>
      <td>TiteSeq_hDPP4_07_bin1</td>
      <td>TiteSeq_hDPP4</td>
      <td>1</td>
      <td>7.0</td>
      <td>221213</td>
      <td>6185457</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221213_MERS-PDF2180_bc-seq/221213_s7_b1_S25_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib51_53</td>
      <td>A</td>
      <td>TiteSeq_hDPP4_07_bin2</td>
      <td>TiteSeq_hDPP4</td>
      <td>2</td>
      <td>7.0</td>
      <td>221213</td>
      <td>182421</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221213_MERS-PDF2180_bc-seq/221213_s7_b2_S26_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib51_53</td>
      <td>A</td>
      <td>TiteSeq_hDPP4_07_bin3</td>
      <td>TiteSeq_hDPP4</td>
      <td>3</td>
      <td>7.0</td>
      <td>221213</td>
      <td>5522</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221213_MERS-PDF2180_bc-seq/221213_s7_b3_S27_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib51_53</td>
      <td>A</td>
      <td>TiteSeq_hDPP4_07_bin4</td>
      <td>TiteSeq_hDPP4</td>
      <td>4</td>
      <td>7.0</td>
      <td>221213</td>
      <td>303</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221213_MERS-PDF2180_bc-seq/221213_s7_b4_S28_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib51_53</td>
      <td>A</td>
      <td>TiteSeq_hDPP4_08_bin1</td>
      <td>TiteSeq_hDPP4</td>
      <td>1</td>
      <td>8.0</td>
      <td>221213</td>
      <td>6432873</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221213_MERS-PDF2180_bc-seq/221213_s8_b1_S29_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib51_53</td>
      <td>A</td>
      <td>TiteSeq_hDPP4_08_bin2</td>
      <td>TiteSeq_hDPP4</td>
      <td>2</td>
      <td>8.0</td>
      <td>221213</td>
      <td>175979</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221213_MERS-PDF2180_bc-seq/221213_s8_b2_S30_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib51_53</td>
      <td>A</td>
      <td>TiteSeq_hDPP4_08_bin3</td>
      <td>TiteSeq_hDPP4</td>
      <td>3</td>
      <td>8.0</td>
      <td>221213</td>
      <td>2807</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221213_MERS-PDF2180_bc-seq/221213_s8_b3_S31_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib51_53</td>
      <td>A</td>
      <td>TiteSeq_hDPP4_08_bin4</td>
      <td>TiteSeq_hDPP4</td>
      <td>4</td>
      <td>8.0</td>
      <td>221213</td>
      <td>234</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221213_MERS-PDF2180_bc-seq/221213_s8_b4_S32_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib51_53</td>
      <td>A</td>
      <td>TiteSeq_hDPP4_09_bin1</td>
      <td>TiteSeq_hDPP4</td>
      <td>1</td>
      <td>9.0</td>
      <td>221213</td>
      <td>6131811</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221213_MERS-PDF2180_bc-seq/221213_s9_b1_S33_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib51_53</td>
      <td>A</td>
      <td>TiteSeq_hDPP4_09_bin2</td>
      <td>TiteSeq_hDPP4</td>
      <td>2</td>
      <td>9.0</td>
      <td>221213</td>
      <td>573092</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221213_MERS-PDF2180_bc-seq/221213_s9_b2_S34_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib51_53</td>
      <td>A</td>
      <td>TiteSeq_hDPP4_09_bin3</td>
      <td>TiteSeq_hDPP4</td>
      <td>3</td>
      <td>9.0</td>
      <td>221213</td>
      <td>814</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221213_MERS-PDF2180_bc-seq/221213_s9_b3_S35_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib51_53</td>
      <td>A</td>
      <td>TiteSeq_hDPP4_09_bin4</td>
      <td>TiteSeq_hDPP4</td>
      <td>4</td>
      <td>9.0</td>
      <td>221213</td>
      <td>119</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221213_MERS-PDF2180_bc-seq/221213_s9_b4_S36_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib51_53</td>
      <td>A</td>
      <td>SortSeq_bin1</td>
      <td>SortSeq</td>
      <td>1</td>
      <td>NaN</td>
      <td>221213</td>
      <td>11620651</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221213_MERS-PDF2180_bc-seq/221213_lib51-53_bin1_1_S53_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221213_MERS-PDF2180_bc-seq/221213_lib51-53_bin1_2_S54_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221213_MERS-PDF2180_bc-seq/221213_lib51-53_bin1_3_S55_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib51_53</td>
      <td>A</td>
      <td>SortSeq_bin2</td>
      <td>SortSeq</td>
      <td>2</td>
      <td>NaN</td>
      <td>221213</td>
      <td>4448236</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221213_MERS-PDF2180_bc-seq/221213_lib51-53_bin2_1_S56_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221213_MERS-PDF2180_bc-seq/221213_lib51-53_bin2_2_S57_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib51_53</td>
      <td>A</td>
      <td>SortSeq_bin3</td>
      <td>SortSeq</td>
      <td>3</td>
      <td>NaN</td>
      <td>221213</td>
      <td>3402367</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221213_MERS-PDF2180_bc-seq/221213_lib51-53_bin3_1_S58_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221213_MERS-PDF2180_bc-seq/221213_lib51-53_bin3_2_S59_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib51_53</td>
      <td>A</td>
      <td>SortSeq_bin4</td>
      <td>SortSeq</td>
      <td>4</td>
      <td>NaN</td>
      <td>221213</td>
      <td>2678000</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221213_MERS-PDF2180_bc-seq/221213_lib51-53_bin4_1_S60_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221213_MERS-PDF2180_bc-seq/221213_lib51-53_bin4_2_S61_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib52_54</td>
      <td>A</td>
      <td>SortSeq_bin1</td>
      <td>SortSeq</td>
      <td>1</td>
      <td>NaN</td>
      <td>221213</td>
      <td>7023038</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230303_MERS-2180-expr-addl-reads/lib52_54_b1_1_S41_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230303_MERS-2180-expr-addl-reads/lib52_54_b1_2_S42_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230303_MERS-2180-expr-addl-reads/lib52_54_b1_3_S43_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib52_54</td>
      <td>A</td>
      <td>SortSeq_bin2</td>
      <td>SortSeq</td>
      <td>2</td>
      <td>NaN</td>
      <td>221213</td>
      <td>3977028</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221213_MERS-PDF2180_bc-seq/221213_lib52-54_bin2_1_S65_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221213_MERS-PDF2180_bc-seq/221213_lib52-54_bin2_2_S66_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib52_54</td>
      <td>A</td>
      <td>SortSeq_bin3</td>
      <td>SortSeq</td>
      <td>3</td>
      <td>NaN</td>
      <td>221213</td>
      <td>3640404</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230303_MERS-2180-expr-addl-reads/lib52_54_b3_1_S44_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2023/230303_MERS-2180-expr-addl-reads/lib52_54_b3_2_S45_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib52_54</td>
      <td>A</td>
      <td>SortSeq_bin4</td>
      <td>SortSeq</td>
      <td>4</td>
      <td>NaN</td>
      <td>221213</td>
      <td>4523054</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221213_MERS-PDF2180_bc-seq/221213_lib52-54_bin4_1_S69_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221213_MERS-PDF2180_bc-seq/221213_lib52-54_bin4_2_S70_R1_001.fastq.gz]</td>
    </tr>
  </tbody>
</table>


Make sure library / sample combinations are unique:


```python
assert len(barcode_runs) == len(barcode_runs.groupby(['library', 'replicate', 'sample']))
```

Make sure the the libraries for which we have barcode runs are all in our variant table:


```python
unknown_libs = set(barcode_runs['library']) - set(variants['library'])
if unknown_libs:
    raise ValueError(f"Libraries with barcode runs not in variant table: {unknown_libs}")
```

Now we initialize an [IlluminaBarcodeParser](https://jbloomlab.github.io/dms_variants/dms_variants.illuminabarcodeparser.html#dms_variants.illuminabarcodeparser.IlluminaBarcodeParser) for each library.

First, get the length of the barcode from the alignment target after making sure the same length for all targets:


```python
bclen = len(targets.targets[0].get_feature('barcode').seq)

assert (bclen == len(target.get_feature('barcode').seq) for target in targets.targets)

print(f"Barcodes of length {bclen}")
```

    Barcodes of length 16


The other barcode parsing params come from the config file:


```python
parser_params = config['illumina_barcode_parser_params']

display(HTML(
    pd.Series(parser_params, name='value')
    .rename_axis(index='parameter')
    .reset_index()
    .to_html(index=False)
    ))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>parameter</th>
      <th>value</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>upstream</td>
      <td>GGCCGC</td>
    </tr>
    <tr>
      <td>downstream</td>
      <td></td>
    </tr>
    <tr>
      <td>minq</td>
      <td>20</td>
    </tr>
    <tr>
      <td>upstream_mismatch</td>
      <td>1</td>
    </tr>
    <tr>
      <td>downstream_mismatch</td>
      <td>0</td>
    </tr>
  </tbody>
</table>


The parser needs to know the set of valid barcodes, which are stored in the variant table and are different for each library.
So we create a different parser for each library using these valid barcode sets:


```python
# create dict keyed by library, value is parser for library
parsers = {lib: dms_variants.illuminabarcodeparser.IlluminaBarcodeParser(
                    bclen=bclen,
                    valid_barcodes=variants.loc[variants['library']==lib]['barcode'],
                    **parser_params)
           for lib in set(variants['library'])}

print('Number of valid barcodes searched for by each parser:')
display(HTML(
    pd.DataFrame([(lib, len(p.valid_barcodes)) for lib, p in parsers.items()],
                 columns=['library', 'number of valid barcodes'])
    .to_html(index=False)
    ))
```

    Number of valid barcodes searched for by each parser:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>library</th>
      <th>number of valid barcodes</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>lib51_53</td>
      <td>111372</td>
    </tr>
    <tr>
      <td>lib52_54</td>
      <td>134378</td>
    </tr>
  </tbody>
</table>


## Parse barcodes
We now parse the barcodes.
Since this will take a while, we utilize multiple CPUs via the Python [multiprocessing](https://docs.python.org/3.6/library/multiprocessing.html) module.
First, determine how many CPUs to use.
We use the minimum of the user-specified number hardcoded below and the number actually available.
(If you are running *interactively* on the Hutch cluster, you may need to reduce the number below in order to avoid an error as there is an enforced CPU limit on the home `rhino` nodes):


```python
ncpus = min(config['max_cpus'], multiprocessing.cpu_count())
print(f"Using {ncpus} CPUs")
```

    Using 8 CPUs


Parse the barcodes in parallel via a [multiprocessing.Pool](https://docs.python.org/3.6/library/multiprocessing.html#multiprocessing.pool.Pool) using all the available CPUs to get a list of the data frames with barcode counts / fates for each sample:


```python
def process_func(parser, r1files, library, replicate, sample):
    """Convenience function to be starmapped to multiprocessing pool."""
    return parser.parse(r1files, add_cols={'library': library, 'replicate': replicate, 'sample': sample})

# parallel computation of list of data frames
with multiprocessing.pool.Pool(processes=ncpus) as pool:
    bclist = pool.starmap(
                process_func,
                [(parsers[run.library], run.R1, run.library, run.replicate, run.sample)
                  for run in barcode_runs.itertuples()],
                )
```

Now concatenate the list into data frames of barcode counts and barcode fates:


```python
counts = pd.concat([samplecounts for samplecounts, _ in bclist],
                   sort=False,
                   ignore_index=True)

print('First few lines of counts data frame:')
display(HTML(counts.head().to_html(index=False)))

fates = pd.concat([samplefates for _, samplefates in bclist],
                  sort=False,
                  ignore_index=True)

print('First few lines of fates data frame:')
display(HTML(fates.head().to_html(index=False)))
```

    First few lines of counts data frame:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>barcode</th>
      <th>count</th>
      <th>library</th>
      <th>replicate</th>
      <th>sample</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>GAATAAACCGAAAACA</td>
      <td>2844</td>
      <td>lib51_53</td>
      <td>A</td>
      <td>TiteSeq_hDPP4_01_bin1</td>
    </tr>
    <tr>
      <td>AGGATGCCACAATCCC</td>
      <td>2782</td>
      <td>lib51_53</td>
      <td>A</td>
      <td>TiteSeq_hDPP4_01_bin1</td>
    </tr>
    <tr>
      <td>AAGACCGAAAAAAATC</td>
      <td>2304</td>
      <td>lib51_53</td>
      <td>A</td>
      <td>TiteSeq_hDPP4_01_bin1</td>
    </tr>
    <tr>
      <td>TTGCGCAATTACAAGA</td>
      <td>2001</td>
      <td>lib51_53</td>
      <td>A</td>
      <td>TiteSeq_hDPP4_01_bin1</td>
    </tr>
    <tr>
      <td>CTATCTACCAAGTAAG</td>
      <td>1997</td>
      <td>lib51_53</td>
      <td>A</td>
      <td>TiteSeq_hDPP4_01_bin1</td>
    </tr>
  </tbody>
</table>


    First few lines of fates data frame:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>fate</th>
      <th>count</th>
      <th>library</th>
      <th>replicate</th>
      <th>sample</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>valid barcode</td>
      <td>12232264</td>
      <td>lib51_53</td>
      <td>A</td>
      <td>TiteSeq_hDPP4_01_bin1</td>
    </tr>
    <tr>
      <td>invalid barcode</td>
      <td>2532373</td>
      <td>lib51_53</td>
      <td>A</td>
      <td>TiteSeq_hDPP4_01_bin1</td>
    </tr>
    <tr>
      <td>low quality barcode</td>
      <td>2431862</td>
      <td>lib51_53</td>
      <td>A</td>
      <td>TiteSeq_hDPP4_01_bin1</td>
    </tr>
    <tr>
      <td>unparseable barcode</td>
      <td>200265</td>
      <td>lib51_53</td>
      <td>A</td>
      <td>TiteSeq_hDPP4_01_bin1</td>
    </tr>
    <tr>
      <td>failed chastity filter</td>
      <td>0</td>
      <td>lib51_53</td>
      <td>A</td>
      <td>TiteSeq_hDPP4_01_bin1</td>
    </tr>
  </tbody>
</table>


## Examine fates of parsed barcodes
First, we'll analyze the "fates" of the parsed barcodes.
These fates represent what happened to each Illumina read we parsed:
 - Did the barcode read fail the Illumina chastity filter?
 - Was the barcode *unparseable* (i.e., the read didn't appear to be a valid barcode based on flanking regions)?
 - Was the barcode sequence too *low quality* based on the Illumina quality scores?
 - Was the barcode parseable but *invalid* (i.e., not in our list of variant-associated barcodes in the codon variant table)?
 - Was the barcode *valid*, and so will be added to variant counts.
 
First, we just write a CSV file with all the barcode fates:


```python
fatesfile = os.path.join(config['counts_dir'], 'barcode_fates.csv')
print(f"Writing barcode fates to {fatesfile}")
fates.to_csv(fatesfile, index=False)
```

    Writing barcode fates to results/counts/barcode_fates.csv


Next, we tabulate the barcode fates in wide format:


```python
display(HTML(fates
             .pivot_table(columns='fate',
                          values='count',
                          index=['library', 'replicate', 'sample'])
             .to_html()
             ))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th></th>
      <th>fate</th>
      <th>failed chastity filter</th>
      <th>invalid barcode</th>
      <th>low quality barcode</th>
      <th>unparseable barcode</th>
      <th>valid barcode</th>
    </tr>
    <tr>
      <th>library</th>
      <th>replicate</th>
      <th>sample</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th rowspan="40" valign="top">lib51_53</th>
      <th rowspan="40" valign="top">A</th>
      <th>SortSeq_bin1</th>
      <td>0</td>
      <td>17281763</td>
      <td>16397713</td>
      <td>1316213</td>
      <td>78651457</td>
    </tr>
    <tr>
      <th>SortSeq_bin2</th>
      <td>0</td>
      <td>1826820</td>
      <td>1785020</td>
      <td>148965</td>
      <td>8815303</td>
    </tr>
    <tr>
      <th>SortSeq_bin3</th>
      <td>0</td>
      <td>1137765</td>
      <td>1102475</td>
      <td>89522</td>
      <td>5491238</td>
    </tr>
    <tr>
      <th>SortSeq_bin4</th>
      <td>0</td>
      <td>763445</td>
      <td>768410</td>
      <td>57834</td>
      <td>3726930</td>
    </tr>
    <tr>
      <th>TiteSeq_hDPP4_01_bin1</th>
      <td>0</td>
      <td>2532373</td>
      <td>2431862</td>
      <td>200265</td>
      <td>12232264</td>
    </tr>
    <tr>
      <th>TiteSeq_hDPP4_01_bin2</th>
      <td>0</td>
      <td>695298</td>
      <td>675814</td>
      <td>55475</td>
      <td>3408326</td>
    </tr>
    <tr>
      <th>TiteSeq_hDPP4_01_bin3</th>
      <td>0</td>
      <td>686533</td>
      <td>647398</td>
      <td>50435</td>
      <td>3134845</td>
    </tr>
    <tr>
      <th>TiteSeq_hDPP4_01_bin4</th>
      <td>0</td>
      <td>1459466</td>
      <td>1373487</td>
      <td>100941</td>
      <td>6467063</td>
    </tr>
    <tr>
      <th>TiteSeq_hDPP4_02_bin1</th>
      <td>0</td>
      <td>2463091</td>
      <td>2394143</td>
      <td>195357</td>
      <td>12089172</td>
    </tr>
    <tr>
      <th>TiteSeq_hDPP4_02_bin2</th>
      <td>0</td>
      <td>1042558</td>
      <td>1002824</td>
      <td>83453</td>
      <td>5029427</td>
    </tr>
    <tr>
      <th>TiteSeq_hDPP4_02_bin3</th>
      <td>0</td>
      <td>757552</td>
      <td>679441</td>
      <td>51852</td>
      <td>3441227</td>
    </tr>
    <tr>
      <th>TiteSeq_hDPP4_02_bin4</th>
      <td>0</td>
      <td>650211</td>
      <td>575371</td>
      <td>43778</td>
      <td>2882348</td>
    </tr>
    <tr>
      <th>TiteSeq_hDPP4_03_bin1</th>
      <td>0</td>
      <td>3387560</td>
      <td>3298836</td>
      <td>266659</td>
      <td>16561480</td>
    </tr>
    <tr>
      <th>TiteSeq_hDPP4_03_bin2</th>
      <td>0</td>
      <td>615381</td>
      <td>590142</td>
      <td>46604</td>
      <td>2886426</td>
    </tr>
    <tr>
      <th>TiteSeq_hDPP4_03_bin3</th>
      <td>0</td>
      <td>738437</td>
      <td>678073</td>
      <td>52926</td>
      <td>3345677</td>
    </tr>
    <tr>
      <th>TiteSeq_hDPP4_03_bin4</th>
      <td>0</td>
      <td>304706</td>
      <td>269195</td>
      <td>20380</td>
      <td>1340113</td>
    </tr>
    <tr>
      <th>TiteSeq_hDPP4_04_bin1</th>
      <td>0</td>
      <td>2977675</td>
      <td>3024685</td>
      <td>239356</td>
      <td>14509064</td>
    </tr>
    <tr>
      <th>TiteSeq_hDPP4_04_bin2</th>
      <td>0</td>
      <td>447903</td>
      <td>427782</td>
      <td>32812</td>
      <td>2105546</td>
    </tr>
    <tr>
      <th>TiteSeq_hDPP4_04_bin3</th>
      <td>0</td>
      <td>439386</td>
      <td>400236</td>
      <td>30570</td>
      <td>1947145</td>
    </tr>
    <tr>
      <th>TiteSeq_hDPP4_04_bin4</th>
      <td>0</td>
      <td>296624</td>
      <td>265645</td>
      <td>16812</td>
      <td>1296289</td>
    </tr>
    <tr>
      <th>TiteSeq_hDPP4_05_bin1</th>
      <td>0</td>
      <td>3778630</td>
      <td>3737872</td>
      <td>294364</td>
      <td>18295799</td>
    </tr>
    <tr>
      <th>TiteSeq_hDPP4_05_bin2</th>
      <td>0</td>
      <td>2742</td>
      <td>8199</td>
      <td>253</td>
      <td>11641</td>
    </tr>
    <tr>
      <th>TiteSeq_hDPP4_05_bin3</th>
      <td>0</td>
      <td>320680</td>
      <td>290979</td>
      <td>19248</td>
      <td>1405930</td>
    </tr>
    <tr>
      <th>TiteSeq_hDPP4_05_bin4</th>
      <td>0</td>
      <td>5043</td>
      <td>6258</td>
      <td>241</td>
      <td>23866</td>
    </tr>
    <tr>
      <th>TiteSeq_hDPP4_06_bin1</th>
      <td>0</td>
      <td>3033132</td>
      <td>3007919</td>
      <td>231431</td>
      <td>14529960</td>
    </tr>
    <tr>
      <th>TiteSeq_hDPP4_06_bin2</th>
      <td>0</td>
      <td>305956</td>
      <td>289093</td>
      <td>21381</td>
      <td>1418299</td>
    </tr>
    <tr>
      <th>TiteSeq_hDPP4_06_bin3</th>
      <td>0</td>
      <td>14637</td>
      <td>16298</td>
      <td>1239</td>
      <td>64510</td>
    </tr>
    <tr>
      <th>TiteSeq_hDPP4_06_bin4</th>
      <td>0</td>
      <td>1148</td>
      <td>2305</td>
      <td>86</td>
      <td>5133</td>
    </tr>
    <tr>
      <th>TiteSeq_hDPP4_07_bin1</th>
      <td>0</td>
      <td>1328832</td>
      <td>1288094</td>
      <td>108424</td>
      <td>6354521</td>
    </tr>
    <tr>
      <th>TiteSeq_hDPP4_07_bin2</th>
      <td>0</td>
      <td>201000</td>
      <td>194369</td>
      <td>15323</td>
      <td>955709</td>
    </tr>
    <tr>
      <th>TiteSeq_hDPP4_07_bin3</th>
      <td>0</td>
      <td>6364</td>
      <td>6534</td>
      <td>479</td>
      <td>26936</td>
    </tr>
    <tr>
      <th>TiteSeq_hDPP4_07_bin4</th>
      <td>0</td>
      <td>978</td>
      <td>1240</td>
      <td>111</td>
      <td>3533</td>
    </tr>
    <tr>
      <th>TiteSeq_hDPP4_08_bin1</th>
      <td>0</td>
      <td>2870263</td>
      <td>2764395</td>
      <td>220157</td>
      <td>13667600</td>
    </tr>
    <tr>
      <th>TiteSeq_hDPP4_08_bin2</th>
      <td>0</td>
      <td>247519</td>
      <td>246283</td>
      <td>19127</td>
      <td>1184999</td>
    </tr>
    <tr>
      <th>TiteSeq_hDPP4_08_bin3</th>
      <td>0</td>
      <td>8658</td>
      <td>8372</td>
      <td>347</td>
      <td>39572</td>
    </tr>
    <tr>
      <th>TiteSeq_hDPP4_08_bin4</th>
      <td>0</td>
      <td>1010</td>
      <td>836</td>
      <td>13</td>
      <td>906</td>
    </tr>
    <tr>
      <th>TiteSeq_hDPP4_09_bin1</th>
      <td>0</td>
      <td>2328526</td>
      <td>2260485</td>
      <td>183967</td>
      <td>11166836</td>
    </tr>
    <tr>
      <th>TiteSeq_hDPP4_09_bin2</th>
      <td>0</td>
      <td>278963</td>
      <td>270566</td>
      <td>21532</td>
      <td>1347292</td>
    </tr>
    <tr>
      <th>TiteSeq_hDPP4_09_bin3</th>
      <td>0</td>
      <td>1860</td>
      <td>5387</td>
      <td>144</td>
      <td>9514</td>
    </tr>
    <tr>
      <th>TiteSeq_hDPP4_09_bin4</th>
      <td>0</td>
      <td>145</td>
      <td>1337</td>
      <td>10</td>
      <td>760</td>
    </tr>
    <tr>
      <th rowspan="4" valign="top">lib52_54</th>
      <th rowspan="4" valign="top">A</th>
      <th>SortSeq_bin1</th>
      <td>0</td>
      <td>10985679</td>
      <td>7046891</td>
      <td>527493</td>
      <td>28563890</td>
    </tr>
    <tr>
      <th>SortSeq_bin2</th>
      <td>0</td>
      <td>1942843</td>
      <td>1187040</td>
      <td>95742</td>
      <td>5203681</td>
    </tr>
    <tr>
      <th>SortSeq_bin3</th>
      <td>0</td>
      <td>2153397</td>
      <td>1415562</td>
      <td>101423</td>
      <td>5784647</td>
    </tr>
    <tr>
      <th>SortSeq_bin4</th>
      <td>0</td>
      <td>2460868</td>
      <td>1591440</td>
      <td>114237</td>
      <td>6490595</td>
    </tr>
  </tbody>
</table>


Now we plot the barcode-read fates for each library / sample, showing the bars for valid barcodes in orange and the others in gray.
We see that the largest fraction of barcode reads correspond to valid barcodes, and most of the others are invalid barcodes (probably because the map to variants that aren't present in our variant table since we didn't associate all variants with barcodes). The exception to this is lib2 Titeseq_03_bin3; the PCR for this sample in the original sequencing run failed, so we followed it up with a single MiSeq lane. We did not filter out the PhiX reads from this data before parsing, so these PhiX reads will deflate the fraction of valid barcode reads as expected, but does not indicate any problems.


```python
barcode_fate_plot = (
    ggplot(
        fates
        .assign(sample=lambda x: pd.Categorical(x['sample'],
                                                x['sample'].unique(),
                                                ordered=True),
                fate=lambda x: pd.Categorical(x['fate'],
                                              x['fate'].unique(),
                                              ordered=True),
                is_valid=lambda x: x['fate'] == 'valid barcode'
                ), 
        aes('fate', 'count', fill='is_valid')) +
    geom_bar(stat='identity') +
    facet_grid('sample ~ library + replicate') +
    facet_grid('sample ~ library + replicate') +
    scale_fill_manual(CBPALETTE, guide=False) +
    theme(figure_size=(1.4 * (1 + fates['library'].nunique() + fates['replicate'].nunique()),
                       1.7 * (1.2 + fates['sample'].nunique())),
          axis_text_x=element_text(angle=90),
          panel_grid_major_x=element_blank()
          ) +
    scale_y_continuous(labels=dms_variants.utils.latex_sci_not,
                       name='number of reads')
    )

_ = barcode_fate_plot.draw()
```


    
![png](count_variants_files/count_variants_42_0.png)
    


## Output csv of barcode counts in variant-barcode lookup table


```python
print(f"Writing variant counts to {config['variant_counts_file']}")
counts.to_csv(config['variant_counts_file'], index=False)
```

    Writing variant counts to results/counts/variant_counts.csv


The [CodonVariantTable](https://jbloomlab.github.io/dms_variants/dms_variants.codonvarianttable.html#dms_variants.codonvarianttable.CodonVariantTable) has lots of nice functions that can be used to analyze the counts it contains.
However, we do that in the next notebook so we don't have to re-run this entire (rather computationally intensive) notebook every time we want to analyze a new aspect of the counts.


```python

```
