<h1>Table of Contents<span class="tocSkip"></span></h1>
<div class="toc"><ul class="toc-item"><li><span><a href="#Process-CCSs" data-toc-modified-id="Process-CCSs-1">Process CCSs</a></span><ul class="toc-item"><li><span><a href="#Setup" data-toc-modified-id="Setup-1.1">Setup</a></span></li><li><span><a href="#PacBio-amplicons" data-toc-modified-id="PacBio-amplicons-1.2">PacBio amplicons</a></span></li><li><span><a href="#CCS-stats-for-PacBio-runs" data-toc-modified-id="CCS-stats-for-PacBio-runs-1.3">CCS stats for PacBio runs</a></span></li><li><span><a href="#Align-CCSs-to-amplicons" data-toc-modified-id="Align-CCSs-to-amplicons-1.4">Align CCSs to amplicons</a></span></li><li><span><a href="#Write-valid-CCSs" data-toc-modified-id="Write-valid-CCSs-1.5">Write valid CCSs</a></span></li></ul></li></ul></div>

This Python Jupyter notebook processes the PacBio circular consensus sequences (CCSs) to extract barcodes and call mutations in the gene. It then builds consensus sequences for barcoded variants from the mutations, to build a barcode-variant lookup table.

# Process CCSs
First, process the PacBio CCSs to extract barcodes and call mutations in the gene.

Define the background for CCS mapping


```python
background = "panmerbeco"
```

## Setup

Import Python modules

Plotting is done with [plotnine](https://plotnine.readthedocs.io/en/stable/), which uses ggplot2-like syntax.

The analysis uses the Bloom lab's [alignparse](https://jbloomlab.github.io/alignparse) and [dms_variants](https://jbloomlab.github.io/dms_variants) packages.


```python
import collections
import math
import os
import re
import time
import warnings

import alignparse
import alignparse.ccs
from alignparse.constants import CBPALETTE
import alignparse.minimap2
import alignparse.targets
import alignparse.consensus
import alignparse.utils

import dms_variants
import dms_variants.codonvarianttable
import dms_variants.plotnine_themes
import dms_variants.utils

from IPython.display import display, HTML

import numpy

import pandas as pd

from plotnine import *

import yaml

%matplotlib inline
```

Set [plotnine](https://plotnine.readthedocs.io/en/stable/) theme to the one defined in [dms_variants](https://jbloomlab.github.io/dms_variants):


```python
theme_set(dms_variants.plotnine_themes.theme_graygrid())
```

Versions of key software:


```python
print(f"Using alignparse version {alignparse.__version__}")
print(f"Using dms_variants version {dms_variants.__version__}")
print(f"Using pandas version {pd.__version__}")
```

    Using alignparse version 0.6.0
    Using dms_variants version 1.4.3
    Using pandas version 1.3.5


Ignore warnings that clutter output:


```python
warnings.simplefilter('ignore')
```

Read the configuration file:


```python
with open('config.yaml') as f:
    config = yaml.safe_load(f)
```

Make output directory for figures:


```python
os.makedirs(config['process_ccs_dir'], exist_ok=True)
```

## PacBio reads
Get the amplicons sequenced by PacBio as the alignment target along with the specs on how to parse the features:


```python
pacbio_runs = (
    pd.read_csv(config['pacbio_runs'], dtype=str)
    .drop(columns=['ccs'])
    .assign(name=lambda x: x['library'] + '_' + x['bg'] + '_' + x['run'],
            fastq=lambda x: config['ccs_dir'] + '/' + x['name'] + '_ccs.fastq.gz'
            )
    )
 
pacbio_runs.drop(pacbio_runs[pacbio_runs['bg'] != background].index,inplace=True)

display(HTML(pacbio_runs.to_html(index=False)))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>library</th>
      <th>bg</th>
      <th>run</th>
      <th>name</th>
      <th>fastq</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>lib82</td>
      <td>panmerbeco</td>
      <td>240206</td>
      <td>lib82_panmerbeco_240206</td>
      <td>results/ccs/lib82_panmerbeco_240206_ccs.fastq.gz</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>panmerbeco</td>
      <td>240206</td>
      <td>lib83_panmerbeco_240206</td>
      <td>results/ccs/lib83_panmerbeco_240206_ccs.fastq.gz</td>
    </tr>
  </tbody>
</table>


## PacBio amplicons
Get the amplicons sequenced by PacBio as the alignment target along with the specs on how to parse the features:


```python
print(f"Reading amplicons from {config['amplicons_' + background]}")
print(f"Reading feature parse specs from {config['feature_parse_specs_' + background]}")

targets = alignparse.targets.Targets(
                seqsfile=config['amplicons_' + background],
                feature_parse_specs=config['feature_parse_specs_' + background])
```

    Reading amplicons from data/PacBio_amplicon_panmerbeco.gb
    Reading feature parse specs from data/feature_parse_specs_panmerbeco.yaml


Draw the target amplicons:


```python
fig = targets.plot(ax_width=7,
                   plots_indexing='biopython',  # numbering starts at 0
                   ax_height=2,  # height of each plot
                   hspace=1.2,  # vertical space between plots
                   )

plotfile = os.path.join(config['process_ccs_dir'], 'amplicons_'+background+'.pdf')
print(f"Saving plot to {plotfile}")
fig.savefig(plotfile, bbox_inches='tight')
```

    Saving plot to results/process_ccs/amplicons_panmerbeco.pdf



    
![png](process_ccs_panmerbeco_files/process_ccs_panmerbeco_22_1.png)
    


Write out the specs used to parse the features (these are the same specs provided as `feature_parse_specs` when initializing `targets`, but with defaults filled in):


```python
print(targets.feature_parse_specs('yaml'))
```

    PaGB01: &id001
      query_clip5: 5
      query_clip3: 5
      termini5:
        filter:
          clip5: 4
          mutation_nt_count: 1
          mutation_op_count: null
          clip3: 0
        return: []
      gene:
        filter:
          mutation_nt_count: 90
          mutation_op_count: null
          clip5: 0
          clip3: 0
        return:
        - mutations
        - accuracy
      spacer:
        filter:
          mutation_nt_count: 1
          mutation_op_count: null
          clip5: 0
          clip3: 0
        return: []
      barcode:
        filter:
          mutation_nt_count: 0
          mutation_op_count: null
          clip5: 0
          clip3: 0
        return:
        - sequence
        - accuracy
      termini3:
        filter:
          clip3: 4
          mutation_nt_count: 1
          mutation_op_count: null
          clip5: 0
        return: []
    BatCoV_SWI2: *id001
    AncHKU25d_alt: *id001
    AncHKU25f_alt: *id001
    AncMERSHKU4_alt: *id001
    hMERS_AKM76239: *id001
    AncHKU25b_alt: *id001
    AncHKU4_alt: *id001
    AncHKU4c_alt: *id001
    NL140455: *id001
    AncHKU25h_MAP_rpk: *id001
    AncHKU25b_alt_rpk: *id001
    AncHKU4b_alt_rpk: *id001
    AncHKU5a_MAP: *id001
    AncMerbeco1_MAP: *id001
    AncHKU4c_MAP: *id001
    hMERS_EMC2012: *id001
    AncNeoGhEri_MAP: *id001
    AncNeo_MAP: *id001
    AncEri_MAP: *id001
    AncMerbeco1_alt: *id001
    HKU5r_19S: *id001
    AncHKU4b_MAP: *id001
    HKU5r_2: *id001
    AncHKU25d_MAP: *id001
    AncMERSHKU4_MAP: *id001
    HKU4r_CZ01: *id001
    AncHKU5b_alt_rpk: *id001
    AncNeoEri_MAP: *id001
    AncHKU5_MAP: *id001
    AncNeo_alt: *id001
    AncHKU4d_MAP: *id001
    HKU4r_BtCoV_1332005: *id001
    BtCoV_B12080: *id001
    HKU5r_31S: *id001
    VsCoVkj15: *id001
    AncHKU25a_MAP: *id001
    AncNeoEri_alt: *id001
    AncHKU25c_alt: *id001
    AncHKU4a_MAP: *id001
    AncEri_alt: *id001
    AncHKU5b_alt: *id001
    HKU4r_162275: *id001
    HKU4r_JPDB144: *id001
    PDF2180_rpk: *id001
    AncHKU25_MAP: *id001
    AncHKU25b_MAP: *id001
    AncHKU25c_MAP: *id001
    AncHKU25e_MAP: *id001
    AncHKU25f_MAP: *id001
    AncMerbeco3_alt: *id001
    AncHKU25HKU5_alt: *id001
    AncHKU25_alt: *id001
    AncHKU25g_alt: *id001
    AncHKU25h_alt: *id001
    AncHKU4a_alt: *id001
    AncHKU4b_alt: *id001
    AncHKU4d_alt: *id001
    hMERS_AMQ48993: *id001
    cMERS_ANI69835: *id001
    cMERS_AMO03401: *id001
    cMERS_ANI69889: *id001
    HKU4: *id001
    HKU4r_161028: *id001
    HKU5r_24S: *id001
    HKU5r_17S: *id001
    SC2013: *id001
    EriCoV_Ita502651: *id001
    EriCoV_Ita5026512: *id001
    AncMerbeco2_MAP_rpk: *id001
    AncHKU25_MAP_rpk: *id001
    AncHKU25a_MAP_rpk: *id001
    AncHKU25b_MAP_rpk: *id001
    AncMERSHKU4_MAP_rpk: *id001
    AncHKU4_MAP_rpk: *id001
    AncHKU4b_MAP_rpk: *id001
    AncNeoGhEri_alt_rpk: *id001
    AncEri_alt_rpk: *id001
    AncMerbeco2_alt_rpk: *id001
    AncMerbeco3_alt_rpk: *id001
    AncHKU25HKU5_alt_rpk: *id001
    AncHKU25_alt_rpk: *id001
    cMERS_AMO03401_rpk: *id001
    cMERS_ANI69889_rpk: *id001
    HKU4_rpk: *id001
    HKU25_YD131305_rpk: *id001
    NL13892_rpk: *id001
    HKU25_NL140462_rpk: *id001
    NL140455_rpk: *id001
    PaGB01_rpk: *id001
    BatCoV_SWI1_rpk: *id001
    BatCoV_SWI2_rpk: *id001
    EriCoV_Deu174_rpk: *id001
    EriCoV_GB1_rpk: *id001
    EriCoV_Ita5026519_rpk: *id001
    EriCoV_Ita1169881_rpk: *id001
    AncHKU25e_alt: *id001
    EriCoV_Ita5026515_rpk: *id001
    EriCoV_HKU31_rpk: *id001
    Russia_MOW1522_rpk: *id001
    EriCoV_Ita1169881: *id001
    AncHKU25HKU5_MAP: *id001
    AncMerbeco2_alt: *id001
    HKU5r_BY140562: *id001
    HKU5r_20S: *id001
    AncHKU25a_alt: *id001
    HKU4r_panP251T: *id001
    EjCoV3: *id001
    AncNeoEri_MAP_rpk: *id001
    HKU5r_32S: *id001
    AncHKU5a_MAP_rpk: *id001
    AncHKU4c_MAP_rpk: *id001
    HKU5r_BY140568: *id001
    AncHKU25h_alt_rpk: *id001
    HKU4r_panP251T_rpk: *id001
    NeoCoV: *id001
    AncNeoEri_alt_rpk: *id001
    AncNeo_alt_rpk: *id001
    HKU5r_33S_rpk: *id001
    HKU25_NL140462: *id001
    AncHKU4d_MAP_rpk: *id001
    AncHKU25c_alt_rpk: *id001
    AncHKU4c_alt_rpk: *id001
    HKU5r_21S_rpk: *id001
    AncEri_MAP_rpk: *id001
    hMERS_AKM76239_rpk: *id001
    HKU5r_31S_rpk: *id001
    AncMerbeco1_alt_rpk: *id001
    AncHKU25e_alt_rpk: *id001
    HKU5r_30S_rpk: *id001
    HKU5r_BY140535_rpk: *id001
    VsCoVkj15_rpk: *id001
    HKU5r_YD13403_rpk: *id001
    SC2013_rpk: *id001
    HKU5r_20S_rpk: *id001
    HKU4r_GZ160421: *id001
    HsCoV_Ita20664540: *id001
    PkCoV_Ita20664563: *id001
    AncHKU4a_alt_rpk: *id001
    AncMerbeco2_MAP: *id001
    KW2EF93Gha2010_rpk: *id001
    HKU5r_Wenzhou1_rpk: *id001
    cMERS_ALA49341: *id001
    VsCoVa7: *id001
    HKU5_rpk: *id001
    HKU5r_21S: *id001
    HKU5r_30S: *id001
    HKU5r_BY140535: *id001
    VsCoV1_rpk: *id001
    EriCoV_Ita5026517_rpk: *id001
    AncHKU25e_MAP_rpk: *id001
    AncHKU25f_MAP_rpk: *id001
    AncHKU25g_MAP_rpk: *id001
    EjCoV3_rpk: *id001
    EriCoV_Ita502651_rpk: *id001
    AncHKU5a_alt: *id001
    HKU4r_GZ1862: *id001
    EriCoV_Ita5026517: *id001
    EriCoV_Ita5026511: *id001
    EriCoV_Ita5026515: *id001
    EriCoV_HKU31: *id001
    KW2EF93Gha2010: *id001
    Russia_MOW1522: *id001
    AncMerbeco1_MAP_rpk: *id001
    AncHKU25c_MAP_rpk: *id001
    AncHKU25d_MAP_rpk: *id001
    HKU5r_27S: *id001
    HKU25_YD131305: *id001
    PkCoV_Ita20664563_rpk: *id001
    cMERS_ALA49374: *id001
    NL13845: *id001
    HKU5: *id001
    HKU5r_YD13403: *id001
    HKU5r_33S: *id001
    VsCoV1: *id001
    AncHKU4a_MAP_rpk: *id001
    hMERS_EMC2012_rpk: *id001
    HKU4r_BtCoV_1332005_rpk: *id001
    HKU5r_A4342005_rpk: *id001
    AncHKU5b_MAP: *id001
    AncHKU4_MAP: *id001
    AncNeoGhEri_alt: *id001
    MjHKU4r1: *id001
    AncHKU25g_MAP: *id001
    AncHKU25h_MAP: *id001
    AncHKU5_alt: *id001
    EriCoV_Ita5026512_rpk: *id001
    NL13892: *id001
    cMERS_ALA49374_rpk: *id001
    HKU5r_2_rpk: *id001
    VsCoVa7_rpk: *id001
    HKU5r_32S_rpk: *id001
    AncHKU4d_alt_rpk: *id001
    cMERS_ANI69835_rpk: *id001
    HKU4r_162275_rpk: *id001
    HKU5r_17S_rpk: *id001
    EriCoV_GB1: *id001
    AncNeoGhEri_MAP_rpk: *id001
    AncHKU25HKU5_MAP_rpk: *id001
    AncHKU5b_MAP_rpk: *id001
    AncHKU25a_alt_rpk: *id001
    AncHKU25d_alt_rpk: *id001
    AncHKU25g_alt_rpk: *id001
    AncHKU5a_alt_rpk: *id001
    AncMERSHKU4_alt_rpk: *id001
    HKU4r_161028_rpk: *id001
    HKU4r_GZ160421_rpk: *id001
    HKU4r_JPDB144_rpk: *id001
    HKU4r_CZ01_rpk: *id001
    NL13845_rpk: *id001
    NL140422_rpk: *id001
    BtCoV_B12080_rpk: *id001
    HKU5r_BY140568_rpk: *id001
    NeoCoV_rpk: *id001
    AncMerbeco3_MAP: *id001
    NL140422: *id001
    HKU5r_A4342005: *id001
    HKU5r_Wenzhou1: *id001
    BatCoV_SWI1: *id001
    PDF2180: *id001
    EriCoV_Deu174: *id001
    EriCoV_Ita5026519: *id001
    AncNeo_MAP_rpk: *id001
    AncMerbeco3_MAP_rpk: *id001
    AncHKU5_MAP_rpk: *id001
    AncHKU25f_alt_rpk: *id001
    AncHKU5_alt_rpk: *id001
    AncHKU4_alt_rpk: *id001
    cMERS_ALA49341_rpk: *id001
    hMERS_AMQ48993_rpk: *id001
    HKU4r_GZ1862_rpk: *id001
    MjHKU4r1_rpk: *id001
    HKU5r_19S_rpk: *id001
    HKU5r_BY140562_rpk: *id001
    HKU5r_27S_rpk: *id001
    HKU5r_24S_rpk: *id001
    HsCoV_Ita20664540_rpk: *id001
    EriCoV_Ita5026511_rpk: *id001
    BtCoV_NL: *id001
    BtCoV_NL_rpk: *id001
    


## Align CCSs to amplicons
We now align the CCSs to the amplicon and parse features from the resulting alignments using the specs above.

First, we initialize an `alignparse.minimap2.Mapper` to align the reads to SAM files:


```python
mapper = alignparse.minimap2.Mapper(alignparse.minimap2.OPTIONS_CODON_DMS)

print(f"Using `minimap2` {mapper.version} with these options:\n" +
      ' '.join(mapper.options))
```

    Using `minimap2` 2.26-r1175 with these options:
    -A2 -B4 -O12 -E2 --end-bonus=13 --secondary=no --cs


Next, we use `Targets.align_and_parse` to create the alignments and parse them:


```python
readstats, aligned, filtered = targets.align_and_parse(
        df=pacbio_runs,
        mapper=mapper,
        outdir=config['process_ccs_dir'],
        name_col='run',
        group_cols=['name', 'library'],
        queryfile_col='fastq',
        overwrite=True,
        ncpus=config['max_cpus'],
        )
```

First, examine the read stats from the alignment / parsing, both extracting alignment target name and getting stats aggregated by target:


```python
readstats = (
    readstats
    .assign(category_all_targets=lambda x: x['category'].str.split().str[0],
            target=lambda x: x['category'].str.split(None, 1).str[1],
            valid=lambda x: x['category_all_targets'] == 'aligned')
    )
```

Now plot the read stats by run (combining all targets and libraries within a run):


```python
ncol = 2
p = (
    ggplot(readstats
           .groupby(['name', 'category_all_targets', 'valid'])
           .aggregate({'count': 'sum'})
           .reset_index(),
           aes('category_all_targets', 'count', fill='valid')) +
    geom_bar(stat='identity') +
    facet_wrap('~ name', ncol=ncol) +
    theme(axis_text_x=element_text(angle=90),
          figure_size=(1.85 * min(ncol, len(pacbio_runs)),
                       2 * math.ceil(len(pacbio_runs) / ncol)),
          panel_grid_major_x=element_blank(),
          legend_position='none',
          ) +
    scale_fill_manual(values=CBPALETTE)
    )
_ = p.draw()
```


    
![png](process_ccs_panmerbeco_files/process_ccs_panmerbeco_32_0.png)
    


And the read stats by library (combining all targets and runs within a library):


```python
p = (
    ggplot(readstats
           .groupby(['library', 'category_all_targets', 'valid'])
           .aggregate({'count': 'sum'})
           .reset_index(), 
           aes('category_all_targets', 'count', fill='valid')) +
    geom_bar(stat='identity') +
    facet_wrap('~ library', nrow=1) +
    theme(axis_text_x=element_text(angle=90),
          figure_size=(1.5 * pacbio_runs['library'].nunique(), 2),
          panel_grid_major_x=element_blank(),
          legend_position='none',
          ) +
    scale_fill_manual(values=CBPALETTE)
    )
_ = p.draw()
```


    
![png](process_ccs_panmerbeco_files/process_ccs_panmerbeco_34_0.png)
    


And the number of reads by target (combining all libraries and runs for a target):


```python
p = (
    ggplot(readstats
           .groupby(['target'])
           .aggregate({'count': 'sum'})
           .reset_index(), 
           aes('count', 'target')) +
    geom_point(stat='identity', size=3) +
    theme(axis_text_x=element_text(angle=90),
          figure_size=(2, 0.3 * readstats['target'].nunique()),
          panel_grid_major_y=element_blank(),
          ) +
    scale_x_log10(name='number of reads')
    )
_ = p.draw()
```


    
![png](process_ccs_panmerbeco_files/process_ccs_panmerbeco_36_0.png)
    


And read stats by target (combining all libraries and runs for a target):


```python
p = (
    ggplot(readstats
           .groupby(['target', 'valid'])
           .aggregate({'count': 'sum'})
           .reset_index()
           .assign(total=lambda x: x.groupby('target')['count'].transform('sum'),
                   frac=lambda x: x['count'] / x['total'],
                   ), 
           aes('target', 'frac', fill='valid')) +
    geom_bar(stat='identity') +
    theme(axis_text_x=element_text(angle=90),
          figure_size=(0.5 * readstats['target'].nunique(), 2),
          panel_grid_major_x=element_blank(),
          ) +
    scale_fill_manual(values=CBPALETTE)
    )
_ = p.draw()
```


    
![png](process_ccs_panmerbeco_files/process_ccs_panmerbeco_38_0.png)
    


Now let's see **why** we filtered the reads.
First, we do some transformations on the `filtered` dict returned by `Targets.align_and_parse`.
Then we count up the number of CCSs filtered for each reason, and group together "unusual" reasons that represent less than some fraction of all filtering.
For now, we group together all targets to the stats represent all targets combined:


```python
other_cutoff = 0.02  # group as "other" reasons with <= this frac

filtered_df = (
    pd.concat(df.assign(target=target) for target, df in filtered.items())
    .groupby(['library', 'name', 'run', 'filter_reason'])
    .size()
    .rename('count')
    .reset_index()
    .assign(tot_reason_frac=lambda x: (x.groupby('filter_reason')['count']
                                       .transform('sum')) / x['count'].sum(),
            filter_reason=lambda x: numpy.where(x['tot_reason_frac'] > other_cutoff,
                                                x['filter_reason'],
                                                'other')
            )
    )
```

Now plot the filtering reason for all runs:


```python
ncol = 7
nreasons = filtered_df['filter_reason'].nunique()

p = (
    ggplot(filtered_df, aes('filter_reason', 'count')) +
    geom_bar(stat='identity') +
    facet_wrap('~ name', ncol=ncol) +
    theme(axis_text_x=element_text(angle=90),
          figure_size=(0.25 * nreasons * min(ncol, len(pacbio_runs)),
                       2 * math.ceil(len(pacbio_runs) / ncol)),
          panel_grid_major_x=element_blank(),
          )
    )
_ = p.draw()
```


    
![png](process_ccs_panmerbeco_files/process_ccs_panmerbeco_42_0.png)
    


Finally, we take the successfully parsed alignments and read them into a data frame, keeping track of the target that each CCS aligns to.
We also drop the pieces of information we won't use going forward, and rename a few columns:


```python
aligned_df = (
    pd.concat(df.assign(target=target) for target, df in aligned.items())
    .drop(columns=['query_clip5', 'query_clip3', 'run','name'])
    .rename(columns={'barcode_sequence': 'barcode'})
    )

print(f"First few lines of information on the parsed alignments:")
display(HTML(aligned_df.head().to_html(index=False)))
```

    First few lines of information on the parsed alignments:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>library</th>
      <th>query_name</th>
      <th>gene_mutations</th>
      <th>gene_accuracy</th>
      <th>barcode</th>
      <th>barcode_accuracy</th>
      <th>target</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>lib82</td>
      <td>m64296e_240129_011756/262538/ccs</td>
      <td></td>
      <td>1.0</td>
      <td>TTAATAATGTCATATT</td>
      <td>1.0</td>
      <td>AncHKU25d_MAP</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>m64296e_240129_011756/329537/ccs</td>
      <td></td>
      <td>1.0</td>
      <td>TTCCGCCACAGACCTG</td>
      <td>1.0</td>
      <td>AncHKU25d_MAP</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>m64296e_240129_011756/527062/ccs</td>
      <td></td>
      <td>1.0</td>
      <td>ATACTCATTCCCCAAA</td>
      <td>1.0</td>
      <td>AncHKU25d_MAP</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>m64296e_240129_011756/721871/ccs</td>
      <td></td>
      <td>1.0</td>
      <td>CTACGTAGATGGTAGT</td>
      <td>1.0</td>
      <td>AncHKU25d_MAP</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>m64296e_240129_011756/788119/ccs</td>
      <td></td>
      <td>1.0</td>
      <td>AACAATCCAAGAGATG</td>
      <td>1.0</td>
      <td>AncHKU25d_MAP</td>
    </tr>
  </tbody>
</table>


## Write valid CCSs

Write the processed CCSs to a file:


```python
aligned_df.to_csv(config['processed_ccs_file' + '_' + background], index=False)

print("Barcodes and mutations for valid processed CCSs "
      f"have been written to {config['processed_ccs_file' + '_' + background]}.")
```

    Barcodes and mutations for valid processed CCSs have been written to results/process_ccs/processed_ccs_panmerbeco.csv.


Next, we analyze these processed CCSs to build the variants.

# Build barcode variant table
Builds consensus sequences for barcoded variants from the mutations called in the processed PacBio CCSs.
Uses these consensus sequences to build a codon variant table.

Make output directories if needed:


```python
os.makedirs(config['variants_dir'], exist_ok=True)
os.makedirs(config['figs_dir'], exist_ok=True)
```

Read the CSV file with the processed CCSs into a data frame, display first few lines:


```python
processed_ccs = pd.read_csv(config['processed_ccs_file' + '_' + background], na_filter=None)

nlibs = processed_ccs['library'].nunique()  # number of unique libraries

ntargets = processed_ccs['target'].nunique()  # number of unique targets

print(f"Read {len(processed_ccs)} CCSs from {nlibs} libraries and {ntargets} targets.")

display(HTML(processed_ccs.head().to_html(index=False)))
```

    Read 277474 CCSs from 2 libraries and 238 targets.



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>library</th>
      <th>query_name</th>
      <th>gene_mutations</th>
      <th>gene_accuracy</th>
      <th>barcode</th>
      <th>barcode_accuracy</th>
      <th>target</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>lib82</td>
      <td>m64296e_240129_011756/262538/ccs</td>
      <td></td>
      <td>1.0</td>
      <td>TTAATAATGTCATATT</td>
      <td>1.0</td>
      <td>AncHKU25d_MAP</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>m64296e_240129_011756/329537/ccs</td>
      <td></td>
      <td>1.0</td>
      <td>TTCCGCCACAGACCTG</td>
      <td>1.0</td>
      <td>AncHKU25d_MAP</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>m64296e_240129_011756/527062/ccs</td>
      <td></td>
      <td>1.0</td>
      <td>ATACTCATTCCCCAAA</td>
      <td>1.0</td>
      <td>AncHKU25d_MAP</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>m64296e_240129_011756/721871/ccs</td>
      <td></td>
      <td>1.0</td>
      <td>CTACGTAGATGGTAGT</td>
      <td>1.0</td>
      <td>AncHKU25d_MAP</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>m64296e_240129_011756/788119/ccs</td>
      <td></td>
      <td>1.0</td>
      <td>AACAATCCAAGAGATG</td>
      <td>1.0</td>
      <td>AncHKU25d_MAP</td>
    </tr>
  </tbody>
</table>


Optional: discard reads that have a `fwd` or `rev` indicator suggestive of heteroduplexes (or we suspect in these runs, some form of pseudo-duplex).


```python
#start = len(processed_ccs)
#processed_ccs = processed_ccs.loc[~processed_ccs['query_name'].str.contains('fwd')]
#processed_ccs = processed_ccs.loc[~processed_ccs['query_name'].str.contains('rev')]
#end = len(processed_ccs)

#print(f"Went from {start} CCSs to {end} after discarded split/duplexed fwd and rev reads.")

#display(HTML(processed_ccs.head().to_html(index=False)))
```

Overall statistics on number of total CCSs and number of unique barcodes:


```python
display(HTML(
    processed_ccs
    .groupby(['target', 'library'])
    .aggregate(total_CCSs=('barcode', 'size'),
               unique_barcodes=('barcode', 'nunique'))
    .assign(avg_CCSs_per_barcode=lambda x: x['total_CCSs'] / x['unique_barcodes'])
    .round(2)
    .to_html()
    ))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th></th>
      <th>total_CCSs</th>
      <th>unique_barcodes</th>
      <th>avg_CCSs_per_barcode</th>
    </tr>
    <tr>
      <th>target</th>
      <th>library</th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th rowspan="2" valign="top">AncEri_MAP</th>
      <th>lib82</th>
      <td>701</td>
      <td>82</td>
      <td>8.55</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>531</td>
      <td>87</td>
      <td>6.10</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncEri_MAP_rpk</th>
      <th>lib82</th>
      <td>462</td>
      <td>72</td>
      <td>6.42</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>409</td>
      <td>58</td>
      <td>7.05</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncEri_alt</th>
      <th>lib82</th>
      <td>167</td>
      <td>24</td>
      <td>6.96</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>93</td>
      <td>17</td>
      <td>5.47</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncEri_alt_rpk</th>
      <th>lib82</th>
      <td>559</td>
      <td>62</td>
      <td>9.02</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>538</td>
      <td>84</td>
      <td>6.40</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU25HKU5_MAP</th>
      <th>lib82</th>
      <td>944</td>
      <td>129</td>
      <td>7.32</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>655</td>
      <td>93</td>
      <td>7.04</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU25HKU5_MAP_rpk</th>
      <th>lib82</th>
      <td>257</td>
      <td>46</td>
      <td>5.59</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>329</td>
      <td>55</td>
      <td>5.98</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU25HKU5_alt</th>
      <th>lib82</th>
      <td>586</td>
      <td>85</td>
      <td>6.89</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>614</td>
      <td>93</td>
      <td>6.60</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU25HKU5_alt_rpk</th>
      <th>lib82</th>
      <td>459</td>
      <td>73</td>
      <td>6.29</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>664</td>
      <td>87</td>
      <td>7.63</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU25_MAP</th>
      <th>lib82</th>
      <td>448</td>
      <td>70</td>
      <td>6.40</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>293</td>
      <td>46</td>
      <td>6.37</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU25_MAP_rpk</th>
      <th>lib82</th>
      <td>372</td>
      <td>62</td>
      <td>6.00</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>366</td>
      <td>61</td>
      <td>6.00</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU25_alt</th>
      <th>lib82</th>
      <td>594</td>
      <td>107</td>
      <td>5.55</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>583</td>
      <td>103</td>
      <td>5.66</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU25_alt_rpk</th>
      <th>lib82</th>
      <td>1111</td>
      <td>136</td>
      <td>8.17</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>1053</td>
      <td>146</td>
      <td>7.21</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU25a_MAP</th>
      <th>lib82</th>
      <td>570</td>
      <td>82</td>
      <td>6.95</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>582</td>
      <td>103</td>
      <td>5.65</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU25a_MAP_rpk</th>
      <th>lib82</th>
      <td>887</td>
      <td>113</td>
      <td>7.85</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>830</td>
      <td>128</td>
      <td>6.48</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU25a_alt</th>
      <th>lib82</th>
      <td>462</td>
      <td>60</td>
      <td>7.70</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>551</td>
      <td>87</td>
      <td>6.33</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU25a_alt_rpk</th>
      <th>lib82</th>
      <td>603</td>
      <td>100</td>
      <td>6.03</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>650</td>
      <td>122</td>
      <td>5.33</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU25b_MAP</th>
      <th>lib82</th>
      <td>348</td>
      <td>66</td>
      <td>5.27</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>433</td>
      <td>65</td>
      <td>6.66</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU25b_MAP_rpk</th>
      <th>lib82</th>
      <td>899</td>
      <td>124</td>
      <td>7.25</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>790</td>
      <td>134</td>
      <td>5.90</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU25b_alt</th>
      <th>lib82</th>
      <td>754</td>
      <td>117</td>
      <td>6.44</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>932</td>
      <td>136</td>
      <td>6.85</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU25b_alt_rpk</th>
      <th>lib82</th>
      <td>419</td>
      <td>66</td>
      <td>6.35</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>612</td>
      <td>93</td>
      <td>6.58</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU25c_MAP</th>
      <th>lib82</th>
      <td>548</td>
      <td>53</td>
      <td>10.34</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>456</td>
      <td>73</td>
      <td>6.25</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU25c_MAP_rpk</th>
      <th>lib82</th>
      <td>352</td>
      <td>51</td>
      <td>6.90</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>453</td>
      <td>59</td>
      <td>7.68</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU25c_alt</th>
      <th>lib82</th>
      <td>775</td>
      <td>105</td>
      <td>7.38</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>521</td>
      <td>80</td>
      <td>6.51</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU25c_alt_rpk</th>
      <th>lib82</th>
      <td>817</td>
      <td>118</td>
      <td>6.92</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>665</td>
      <td>115</td>
      <td>5.78</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU25d_MAP</th>
      <th>lib82</th>
      <td>627</td>
      <td>85</td>
      <td>7.38</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>564</td>
      <td>88</td>
      <td>6.41</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU25d_MAP_rpk</th>
      <th>lib82</th>
      <td>689</td>
      <td>92</td>
      <td>7.49</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>703</td>
      <td>88</td>
      <td>7.99</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU25d_alt</th>
      <th>lib82</th>
      <td>438</td>
      <td>80</td>
      <td>5.48</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>647</td>
      <td>105</td>
      <td>6.16</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU25d_alt_rpk</th>
      <th>lib82</th>
      <td>766</td>
      <td>109</td>
      <td>7.03</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>592</td>
      <td>84</td>
      <td>7.05</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU25e_MAP</th>
      <th>lib82</th>
      <td>189</td>
      <td>57</td>
      <td>3.32</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>210</td>
      <td>62</td>
      <td>3.39</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU25e_MAP_rpk</th>
      <th>lib82</th>
      <td>445</td>
      <td>59</td>
      <td>7.54</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>427</td>
      <td>69</td>
      <td>6.19</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU25e_alt</th>
      <th>lib82</th>
      <td>937</td>
      <td>118</td>
      <td>7.94</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>752</td>
      <td>114</td>
      <td>6.60</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU25e_alt_rpk</th>
      <th>lib82</th>
      <td>451</td>
      <td>59</td>
      <td>7.64</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>385</td>
      <td>40</td>
      <td>9.62</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU25f_MAP</th>
      <th>lib82</th>
      <td>628</td>
      <td>112</td>
      <td>5.61</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>586</td>
      <td>104</td>
      <td>5.63</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU25f_MAP_rpk</th>
      <th>lib82</th>
      <td>543</td>
      <td>92</td>
      <td>5.90</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>830</td>
      <td>139</td>
      <td>5.97</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU25f_alt</th>
      <th>lib82</th>
      <td>428</td>
      <td>64</td>
      <td>6.69</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>442</td>
      <td>86</td>
      <td>5.14</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU25f_alt_rpk</th>
      <th>lib82</th>
      <td>524</td>
      <td>79</td>
      <td>6.63</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>749</td>
      <td>93</td>
      <td>8.05</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU25g_MAP</th>
      <th>lib82</th>
      <td>686</td>
      <td>109</td>
      <td>6.29</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>603</td>
      <td>82</td>
      <td>7.35</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU25g_MAP_rpk</th>
      <th>lib82</th>
      <td>725</td>
      <td>119</td>
      <td>6.09</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>732</td>
      <td>116</td>
      <td>6.31</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU25g_alt</th>
      <th>lib82</th>
      <td>282</td>
      <td>47</td>
      <td>6.00</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>398</td>
      <td>60</td>
      <td>6.63</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU25g_alt_rpk</th>
      <th>lib82</th>
      <td>496</td>
      <td>70</td>
      <td>7.09</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>597</td>
      <td>80</td>
      <td>7.46</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU25h_MAP</th>
      <th>lib82</th>
      <td>674</td>
      <td>93</td>
      <td>7.25</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>786</td>
      <td>100</td>
      <td>7.86</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU25h_MAP_rpk</th>
      <th>lib82</th>
      <td>356</td>
      <td>63</td>
      <td>5.65</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>434</td>
      <td>68</td>
      <td>6.38</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU25h_alt</th>
      <th>lib82</th>
      <td>736</td>
      <td>124</td>
      <td>5.94</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>623</td>
      <td>97</td>
      <td>6.42</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU25h_alt_rpk</th>
      <th>lib82</th>
      <td>245</td>
      <td>37</td>
      <td>6.62</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>439</td>
      <td>73</td>
      <td>6.01</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU4_MAP</th>
      <th>lib82</th>
      <td>395</td>
      <td>60</td>
      <td>6.58</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>541</td>
      <td>69</td>
      <td>7.84</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU4_MAP_rpk</th>
      <th>lib82</th>
      <td>666</td>
      <td>98</td>
      <td>6.80</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>637</td>
      <td>116</td>
      <td>5.49</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU4_alt</th>
      <th>lib82</th>
      <td>627</td>
      <td>86</td>
      <td>7.29</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>597</td>
      <td>105</td>
      <td>5.69</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU4_alt_rpk</th>
      <th>lib82</th>
      <td>708</td>
      <td>112</td>
      <td>6.32</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>605</td>
      <td>91</td>
      <td>6.65</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU4a_MAP</th>
      <th>lib82</th>
      <td>333</td>
      <td>46</td>
      <td>7.24</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>273</td>
      <td>46</td>
      <td>5.93</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU4a_MAP_rpk</th>
      <th>lib82</th>
      <td>181</td>
      <td>38</td>
      <td>4.76</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>447</td>
      <td>65</td>
      <td>6.88</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU4a_alt</th>
      <th>lib82</th>
      <td>509</td>
      <td>71</td>
      <td>7.17</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>360</td>
      <td>58</td>
      <td>6.21</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU4a_alt_rpk</th>
      <th>lib82</th>
      <td>323</td>
      <td>55</td>
      <td>5.87</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>562</td>
      <td>78</td>
      <td>7.21</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU4b_MAP</th>
      <th>lib82</th>
      <td>511</td>
      <td>70</td>
      <td>7.30</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>785</td>
      <td>104</td>
      <td>7.55</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU4b_MAP_rpk</th>
      <th>lib82</th>
      <td>483</td>
      <td>81</td>
      <td>5.96</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>611</td>
      <td>88</td>
      <td>6.94</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU4b_alt</th>
      <th>lib82</th>
      <td>719</td>
      <td>91</td>
      <td>7.90</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>700</td>
      <td>100</td>
      <td>7.00</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU4b_alt_rpk</th>
      <th>lib82</th>
      <td>825</td>
      <td>109</td>
      <td>7.57</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>688</td>
      <td>91</td>
      <td>7.56</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU4c_MAP</th>
      <th>lib82</th>
      <td>592</td>
      <td>97</td>
      <td>6.10</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>408</td>
      <td>79</td>
      <td>5.16</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU4c_MAP_rpk</th>
      <th>lib82</th>
      <td>644</td>
      <td>108</td>
      <td>5.96</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>716</td>
      <td>101</td>
      <td>7.09</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU4c_alt</th>
      <th>lib82</th>
      <td>804</td>
      <td>120</td>
      <td>6.70</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>763</td>
      <td>112</td>
      <td>6.81</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU4c_alt_rpk</th>
      <th>lib82</th>
      <td>661</td>
      <td>80</td>
      <td>8.26</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>711</td>
      <td>97</td>
      <td>7.33</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU4d_MAP</th>
      <th>lib82</th>
      <td>589</td>
      <td>69</td>
      <td>8.54</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>667</td>
      <td>71</td>
      <td>9.39</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU4d_MAP_rpk</th>
      <th>lib82</th>
      <td>385</td>
      <td>48</td>
      <td>8.02</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>426</td>
      <td>59</td>
      <td>7.22</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU4d_alt</th>
      <th>lib82</th>
      <td>536</td>
      <td>90</td>
      <td>5.96</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>430</td>
      <td>92</td>
      <td>4.67</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU4d_alt_rpk</th>
      <th>lib82</th>
      <td>695</td>
      <td>95</td>
      <td>7.32</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>758</td>
      <td>99</td>
      <td>7.66</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU5_MAP</th>
      <th>lib82</th>
      <td>659</td>
      <td>82</td>
      <td>8.04</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>524</td>
      <td>74</td>
      <td>7.08</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU5_MAP_rpk</th>
      <th>lib82</th>
      <td>607</td>
      <td>84</td>
      <td>7.23</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>701</td>
      <td>93</td>
      <td>7.54</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU5_alt</th>
      <th>lib82</th>
      <td>445</td>
      <td>58</td>
      <td>7.67</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>534</td>
      <td>66</td>
      <td>8.09</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU5_alt_rpk</th>
      <th>lib82</th>
      <td>909</td>
      <td>112</td>
      <td>8.12</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>871</td>
      <td>135</td>
      <td>6.45</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU5a_MAP</th>
      <th>lib82</th>
      <td>734</td>
      <td>119</td>
      <td>6.17</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>749</td>
      <td>150</td>
      <td>4.99</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU5a_MAP_rpk</th>
      <th>lib82</th>
      <td>646</td>
      <td>88</td>
      <td>7.34</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>550</td>
      <td>92</td>
      <td>5.98</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU5a_alt</th>
      <th>lib82</th>
      <td>310</td>
      <td>53</td>
      <td>5.85</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>348</td>
      <td>70</td>
      <td>4.97</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU5a_alt_rpk</th>
      <th>lib82</th>
      <td>533</td>
      <td>74</td>
      <td>7.20</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>612</td>
      <td>78</td>
      <td>7.85</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU5b_MAP</th>
      <th>lib82</th>
      <td>627</td>
      <td>105</td>
      <td>5.97</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>713</td>
      <td>110</td>
      <td>6.48</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU5b_MAP_rpk</th>
      <th>lib82</th>
      <td>446</td>
      <td>74</td>
      <td>6.03</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>357</td>
      <td>57</td>
      <td>6.26</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU5b_alt</th>
      <th>lib82</th>
      <td>818</td>
      <td>98</td>
      <td>8.35</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>604</td>
      <td>105</td>
      <td>5.75</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU5b_alt_rpk</th>
      <th>lib82</th>
      <td>852</td>
      <td>125</td>
      <td>6.82</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>754</td>
      <td>121</td>
      <td>6.23</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncMERSHKU4_MAP</th>
      <th>lib82</th>
      <td>617</td>
      <td>91</td>
      <td>6.78</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>712</td>
      <td>102</td>
      <td>6.98</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncMERSHKU4_MAP_rpk</th>
      <th>lib82</th>
      <td>614</td>
      <td>98</td>
      <td>6.27</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>761</td>
      <td>115</td>
      <td>6.62</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncMERSHKU4_alt</th>
      <th>lib82</th>
      <td>779</td>
      <td>109</td>
      <td>7.15</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>856</td>
      <td>104</td>
      <td>8.23</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncMERSHKU4_alt_rpk</th>
      <th>lib82</th>
      <td>328</td>
      <td>48</td>
      <td>6.83</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>450</td>
      <td>76</td>
      <td>5.92</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncMerbeco1_MAP</th>
      <th>lib82</th>
      <td>480</td>
      <td>85</td>
      <td>5.65</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>461</td>
      <td>64</td>
      <td>7.20</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncMerbeco1_MAP_rpk</th>
      <th>lib82</th>
      <td>620</td>
      <td>86</td>
      <td>7.21</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>600</td>
      <td>87</td>
      <td>6.90</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncMerbeco1_alt</th>
      <th>lib82</th>
      <td>579</td>
      <td>88</td>
      <td>6.58</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>762</td>
      <td>117</td>
      <td>6.51</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncMerbeco1_alt_rpk</th>
      <th>lib82</th>
      <td>482</td>
      <td>60</td>
      <td>8.03</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>525</td>
      <td>74</td>
      <td>7.09</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncMerbeco2_MAP</th>
      <th>lib82</th>
      <td>139</td>
      <td>24</td>
      <td>5.79</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>263</td>
      <td>32</td>
      <td>8.22</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncMerbeco2_MAP_rpk</th>
      <th>lib82</th>
      <td>399</td>
      <td>74</td>
      <td>5.39</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>583</td>
      <td>83</td>
      <td>7.02</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncMerbeco2_alt</th>
      <th>lib82</th>
      <td>374</td>
      <td>67</td>
      <td>5.58</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>679</td>
      <td>87</td>
      <td>7.80</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncMerbeco2_alt_rpk</th>
      <th>lib82</th>
      <td>699</td>
      <td>113</td>
      <td>6.19</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>688</td>
      <td>112</td>
      <td>6.14</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncMerbeco3_MAP</th>
      <th>lib82</th>
      <td>462</td>
      <td>74</td>
      <td>6.24</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>602</td>
      <td>87</td>
      <td>6.92</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncMerbeco3_MAP_rpk</th>
      <th>lib82</th>
      <td>403</td>
      <td>60</td>
      <td>6.72</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>530</td>
      <td>74</td>
      <td>7.16</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncMerbeco3_alt</th>
      <th>lib82</th>
      <td>787</td>
      <td>103</td>
      <td>7.64</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>552</td>
      <td>78</td>
      <td>7.08</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncMerbeco3_alt_rpk</th>
      <th>lib82</th>
      <td>906</td>
      <td>123</td>
      <td>7.37</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>1071</td>
      <td>130</td>
      <td>8.24</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncNeoEri_MAP</th>
      <th>lib82</th>
      <td>416</td>
      <td>46</td>
      <td>9.04</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>309</td>
      <td>44</td>
      <td>7.02</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncNeoEri_MAP_rpk</th>
      <th>lib82</th>
      <td>78</td>
      <td>12</td>
      <td>6.50</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>69</td>
      <td>10</td>
      <td>6.90</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncNeoEri_alt</th>
      <th>lib82</th>
      <td>589</td>
      <td>88</td>
      <td>6.69</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>925</td>
      <td>109</td>
      <td>8.49</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncNeoEri_alt_rpk</th>
      <th>lib82</th>
      <td>706</td>
      <td>107</td>
      <td>6.60</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>647</td>
      <td>104</td>
      <td>6.22</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncNeoGhEri_MAP</th>
      <th>lib82</th>
      <td>868</td>
      <td>112</td>
      <td>7.75</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>647</td>
      <td>96</td>
      <td>6.74</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncNeoGhEri_MAP_rpk</th>
      <th>lib82</th>
      <td>501</td>
      <td>71</td>
      <td>7.06</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>553</td>
      <td>84</td>
      <td>6.58</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncNeoGhEri_alt</th>
      <th>lib82</th>
      <td>383</td>
      <td>51</td>
      <td>7.51</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>492</td>
      <td>62</td>
      <td>7.94</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncNeoGhEri_alt_rpk</th>
      <th>lib82</th>
      <td>337</td>
      <td>57</td>
      <td>5.91</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>314</td>
      <td>43</td>
      <td>7.30</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncNeo_MAP</th>
      <th>lib82</th>
      <td>506</td>
      <td>76</td>
      <td>6.66</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>590</td>
      <td>92</td>
      <td>6.41</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncNeo_MAP_rpk</th>
      <th>lib82</th>
      <td>563</td>
      <td>74</td>
      <td>7.61</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>563</td>
      <td>87</td>
      <td>6.47</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncNeo_alt</th>
      <th>lib82</th>
      <td>887</td>
      <td>110</td>
      <td>8.06</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>682</td>
      <td>92</td>
      <td>7.41</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncNeo_alt_rpk</th>
      <th>lib82</th>
      <td>602</td>
      <td>72</td>
      <td>8.36</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>551</td>
      <td>83</td>
      <td>6.64</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">BatCoV_SWI1</th>
      <th>lib82</th>
      <td>481</td>
      <td>66</td>
      <td>7.29</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>469</td>
      <td>62</td>
      <td>7.56</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">BatCoV_SWI1_rpk</th>
      <th>lib82</th>
      <td>572</td>
      <td>71</td>
      <td>8.06</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>548</td>
      <td>96</td>
      <td>5.71</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">BatCoV_SWI2</th>
      <th>lib82</th>
      <td>765</td>
      <td>104</td>
      <td>7.36</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>530</td>
      <td>91</td>
      <td>5.82</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">BatCoV_SWI2_rpk</th>
      <th>lib82</th>
      <td>634</td>
      <td>91</td>
      <td>6.97</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>741</td>
      <td>106</td>
      <td>6.99</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">BtCoV_B12080</th>
      <th>lib82</th>
      <td>700</td>
      <td>104</td>
      <td>6.73</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>607</td>
      <td>105</td>
      <td>5.78</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">BtCoV_B12080_rpk</th>
      <th>lib82</th>
      <td>466</td>
      <td>59</td>
      <td>7.90</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>661</td>
      <td>93</td>
      <td>7.11</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">EjCoV3</th>
      <th>lib82</th>
      <td>802</td>
      <td>123</td>
      <td>6.52</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>787</td>
      <td>135</td>
      <td>5.83</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">EjCoV3_rpk</th>
      <th>lib82</th>
      <td>573</td>
      <td>82</td>
      <td>6.99</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>625</td>
      <td>90</td>
      <td>6.94</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">EriCoV_Deu174</th>
      <th>lib82</th>
      <td>152</td>
      <td>28</td>
      <td>5.43</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>292</td>
      <td>40</td>
      <td>7.30</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">EriCoV_Deu174_rpk</th>
      <th>lib82</th>
      <td>377</td>
      <td>59</td>
      <td>6.39</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>585</td>
      <td>78</td>
      <td>7.50</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">EriCoV_GB1</th>
      <th>lib82</th>
      <td>615</td>
      <td>113</td>
      <td>5.44</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>603</td>
      <td>102</td>
      <td>5.91</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">EriCoV_GB1_rpk</th>
      <th>lib82</th>
      <td>608</td>
      <td>97</td>
      <td>6.27</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>662</td>
      <td>101</td>
      <td>6.55</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">EriCoV_HKU31</th>
      <th>lib82</th>
      <td>313</td>
      <td>40</td>
      <td>7.82</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>448</td>
      <td>66</td>
      <td>6.79</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">EriCoV_HKU31_rpk</th>
      <th>lib82</th>
      <td>454</td>
      <td>57</td>
      <td>7.96</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>365</td>
      <td>52</td>
      <td>7.02</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">EriCoV_Ita1169881</th>
      <th>lib82</th>
      <td>455</td>
      <td>74</td>
      <td>6.15</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>466</td>
      <td>83</td>
      <td>5.61</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">EriCoV_Ita1169881_rpk</th>
      <th>lib82</th>
      <td>260</td>
      <td>29</td>
      <td>8.97</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>124</td>
      <td>24</td>
      <td>5.17</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">EriCoV_Ita502651</th>
      <th>lib82</th>
      <td>191</td>
      <td>30</td>
      <td>6.37</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>344</td>
      <td>45</td>
      <td>7.64</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">EriCoV_Ita5026511</th>
      <th>lib82</th>
      <td>371</td>
      <td>64</td>
      <td>5.80</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>614</td>
      <td>88</td>
      <td>6.98</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">EriCoV_Ita5026511_rpk</th>
      <th>lib82</th>
      <td>102</td>
      <td>11</td>
      <td>9.27</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>38</td>
      <td>8</td>
      <td>4.75</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">EriCoV_Ita5026512</th>
      <th>lib82</th>
      <td>537</td>
      <td>69</td>
      <td>7.78</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>550</td>
      <td>75</td>
      <td>7.33</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">EriCoV_Ita5026512_rpk</th>
      <th>lib82</th>
      <td>169</td>
      <td>22</td>
      <td>7.68</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>244</td>
      <td>31</td>
      <td>7.87</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">EriCoV_Ita5026515</th>
      <th>lib82</th>
      <td>363</td>
      <td>74</td>
      <td>4.91</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>696</td>
      <td>99</td>
      <td>7.03</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">EriCoV_Ita5026515_rpk</th>
      <th>lib82</th>
      <td>787</td>
      <td>105</td>
      <td>7.50</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>649</td>
      <td>96</td>
      <td>6.76</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">EriCoV_Ita5026517</th>
      <th>lib82</th>
      <td>615</td>
      <td>96</td>
      <td>6.41</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>602</td>
      <td>90</td>
      <td>6.69</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">EriCoV_Ita5026517_rpk</th>
      <th>lib82</th>
      <td>756</td>
      <td>115</td>
      <td>6.57</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>822</td>
      <td>97</td>
      <td>8.47</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">EriCoV_Ita5026519</th>
      <th>lib82</th>
      <td>510</td>
      <td>76</td>
      <td>6.71</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>575</td>
      <td>104</td>
      <td>5.53</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">EriCoV_Ita5026519_rpk</th>
      <th>lib82</th>
      <td>593</td>
      <td>79</td>
      <td>7.51</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>670</td>
      <td>102</td>
      <td>6.57</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">EriCoV_Ita502651_rpk</th>
      <th>lib82</th>
      <td>160</td>
      <td>19</td>
      <td>8.42</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>28</td>
      <td>6</td>
      <td>4.67</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU25_NL140462</th>
      <th>lib82</th>
      <td>689</td>
      <td>110</td>
      <td>6.26</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>891</td>
      <td>153</td>
      <td>5.82</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU25_NL140462_rpk</th>
      <th>lib82</th>
      <td>478</td>
      <td>89</td>
      <td>5.37</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>328</td>
      <td>71</td>
      <td>4.62</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU25_YD131305</th>
      <th>lib82</th>
      <td>647</td>
      <td>103</td>
      <td>6.28</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>481</td>
      <td>86</td>
      <td>5.59</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU25_YD131305_rpk</th>
      <th>lib82</th>
      <td>809</td>
      <td>86</td>
      <td>9.41</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>1119</td>
      <td>107</td>
      <td>10.46</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU4</th>
      <th>lib82</th>
      <td>687</td>
      <td>104</td>
      <td>6.61</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>698</td>
      <td>117</td>
      <td>5.97</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU4_rpk</th>
      <th>lib82</th>
      <td>458</td>
      <td>67</td>
      <td>6.84</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>712</td>
      <td>82</td>
      <td>8.68</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU4r_161028</th>
      <th>lib82</th>
      <td>260</td>
      <td>41</td>
      <td>6.34</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>336</td>
      <td>53</td>
      <td>6.34</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU4r_161028_rpk</th>
      <th>lib82</th>
      <td>746</td>
      <td>107</td>
      <td>6.97</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>757</td>
      <td>108</td>
      <td>7.01</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU4r_162275</th>
      <th>lib82</th>
      <td>461</td>
      <td>60</td>
      <td>7.68</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>560</td>
      <td>86</td>
      <td>6.51</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU4r_162275_rpk</th>
      <th>lib82</th>
      <td>580</td>
      <td>91</td>
      <td>6.37</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>644</td>
      <td>100</td>
      <td>6.44</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU4r_BtCoV_1332005</th>
      <th>lib82</th>
      <td>625</td>
      <td>73</td>
      <td>8.56</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>577</td>
      <td>85</td>
      <td>6.79</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU4r_BtCoV_1332005_rpk</th>
      <th>lib82</th>
      <td>779</td>
      <td>132</td>
      <td>5.90</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>643</td>
      <td>109</td>
      <td>5.90</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU4r_CZ01</th>
      <th>lib82</th>
      <td>263</td>
      <td>50</td>
      <td>5.26</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>500</td>
      <td>84</td>
      <td>5.95</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU4r_CZ01_rpk</th>
      <th>lib82</th>
      <td>562</td>
      <td>91</td>
      <td>6.18</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>447</td>
      <td>70</td>
      <td>6.39</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU4r_GZ160421</th>
      <th>lib82</th>
      <td>646</td>
      <td>97</td>
      <td>6.66</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>476</td>
      <td>75</td>
      <td>6.35</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU4r_GZ160421_rpk</th>
      <th>lib82</th>
      <td>765</td>
      <td>109</td>
      <td>7.02</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>759</td>
      <td>100</td>
      <td>7.59</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU4r_GZ1862</th>
      <th>lib82</th>
      <td>519</td>
      <td>86</td>
      <td>6.03</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>712</td>
      <td>85</td>
      <td>8.38</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU4r_GZ1862_rpk</th>
      <th>lib82</th>
      <td>430</td>
      <td>61</td>
      <td>7.05</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>565</td>
      <td>80</td>
      <td>7.06</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU4r_JPDB144</th>
      <th>lib82</th>
      <td>794</td>
      <td>95</td>
      <td>8.36</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>562</td>
      <td>81</td>
      <td>6.94</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU4r_JPDB144_rpk</th>
      <th>lib82</th>
      <td>753</td>
      <td>97</td>
      <td>7.76</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>952</td>
      <td>135</td>
      <td>7.05</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU4r_panP251T</th>
      <th>lib82</th>
      <td>500</td>
      <td>71</td>
      <td>7.04</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>478</td>
      <td>85</td>
      <td>5.62</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU4r_panP251T_rpk</th>
      <th>lib82</th>
      <td>582</td>
      <td>100</td>
      <td>5.82</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>617</td>
      <td>108</td>
      <td>5.71</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU5</th>
      <th>lib82</th>
      <td>496</td>
      <td>85</td>
      <td>5.84</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>646</td>
      <td>100</td>
      <td>6.46</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU5_rpk</th>
      <th>lib82</th>
      <td>613</td>
      <td>107</td>
      <td>5.73</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>609</td>
      <td>95</td>
      <td>6.41</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU5r_17S</th>
      <th>lib82</th>
      <td>568</td>
      <td>88</td>
      <td>6.45</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>549</td>
      <td>94</td>
      <td>5.84</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU5r_17S_rpk</th>
      <th>lib82</th>
      <td>643</td>
      <td>92</td>
      <td>6.99</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>791</td>
      <td>114</td>
      <td>6.94</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU5r_19S</th>
      <th>lib82</th>
      <td>801</td>
      <td>118</td>
      <td>6.79</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>596</td>
      <td>98</td>
      <td>6.08</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU5r_19S_rpk</th>
      <th>lib82</th>
      <td>606</td>
      <td>99</td>
      <td>6.12</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>584</td>
      <td>104</td>
      <td>5.62</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU5r_2</th>
      <th>lib82</th>
      <td>553</td>
      <td>95</td>
      <td>5.82</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>676</td>
      <td>95</td>
      <td>7.12</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU5r_20S</th>
      <th>lib82</th>
      <td>590</td>
      <td>98</td>
      <td>6.02</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>637</td>
      <td>86</td>
      <td>7.41</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU5r_20S_rpk</th>
      <th>lib82</th>
      <td>848</td>
      <td>85</td>
      <td>9.98</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>919</td>
      <td>116</td>
      <td>7.92</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU5r_21S</th>
      <th>lib82</th>
      <td>364</td>
      <td>55</td>
      <td>6.62</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>424</td>
      <td>65</td>
      <td>6.52</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU5r_21S_rpk</th>
      <th>lib82</th>
      <td>711</td>
      <td>118</td>
      <td>6.03</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>741</td>
      <td>127</td>
      <td>5.83</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU5r_24S</th>
      <th>lib82</th>
      <td>456</td>
      <td>73</td>
      <td>6.25</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>243</td>
      <td>50</td>
      <td>4.86</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU5r_24S_rpk</th>
      <th>lib82</th>
      <td>717</td>
      <td>93</td>
      <td>7.71</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>788</td>
      <td>99</td>
      <td>7.96</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU5r_27S</th>
      <th>lib82</th>
      <td>502</td>
      <td>70</td>
      <td>7.17</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>621</td>
      <td>82</td>
      <td>7.57</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU5r_27S_rpk</th>
      <th>lib82</th>
      <td>474</td>
      <td>81</td>
      <td>5.85</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>505</td>
      <td>67</td>
      <td>7.54</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU5r_2_rpk</th>
      <th>lib82</th>
      <td>982</td>
      <td>120</td>
      <td>8.18</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>825</td>
      <td>129</td>
      <td>6.40</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU5r_30S</th>
      <th>lib82</th>
      <td>670</td>
      <td>94</td>
      <td>7.13</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>617</td>
      <td>106</td>
      <td>5.82</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU5r_30S_rpk</th>
      <th>lib82</th>
      <td>617</td>
      <td>83</td>
      <td>7.43</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>431</td>
      <td>74</td>
      <td>5.82</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU5r_31S</th>
      <th>lib82</th>
      <td>414</td>
      <td>56</td>
      <td>7.39</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>455</td>
      <td>52</td>
      <td>8.75</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU5r_31S_rpk</th>
      <th>lib82</th>
      <td>790</td>
      <td>108</td>
      <td>7.31</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>683</td>
      <td>113</td>
      <td>6.04</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU5r_32S</th>
      <th>lib82</th>
      <td>585</td>
      <td>77</td>
      <td>7.60</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>607</td>
      <td>98</td>
      <td>6.19</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU5r_32S_rpk</th>
      <th>lib82</th>
      <td>529</td>
      <td>72</td>
      <td>7.35</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>634</td>
      <td>83</td>
      <td>7.64</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU5r_33S</th>
      <th>lib82</th>
      <td>812</td>
      <td>116</td>
      <td>7.00</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>787</td>
      <td>98</td>
      <td>8.03</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU5r_33S_rpk</th>
      <th>lib82</th>
      <td>542</td>
      <td>75</td>
      <td>7.23</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>611</td>
      <td>88</td>
      <td>6.94</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU5r_A4342005</th>
      <th>lib82</th>
      <td>623</td>
      <td>109</td>
      <td>5.72</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>554</td>
      <td>96</td>
      <td>5.77</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU5r_A4342005_rpk</th>
      <th>lib82</th>
      <td>659</td>
      <td>108</td>
      <td>6.10</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>806</td>
      <td>130</td>
      <td>6.20</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU5r_BY140535</th>
      <th>lib82</th>
      <td>439</td>
      <td>73</td>
      <td>6.01</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>646</td>
      <td>79</td>
      <td>8.18</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU5r_BY140535_rpk</th>
      <th>lib82</th>
      <td>666</td>
      <td>98</td>
      <td>6.80</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>652</td>
      <td>74</td>
      <td>8.81</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU5r_BY140562</th>
      <th>lib82</th>
      <td>744</td>
      <td>101</td>
      <td>7.37</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>763</td>
      <td>121</td>
      <td>6.31</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU5r_BY140562_rpk</th>
      <th>lib82</th>
      <td>712</td>
      <td>100</td>
      <td>7.12</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>798</td>
      <td>101</td>
      <td>7.90</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU5r_BY140568</th>
      <th>lib82</th>
      <td>475</td>
      <td>79</td>
      <td>6.01</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>853</td>
      <td>107</td>
      <td>7.97</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU5r_BY140568_rpk</th>
      <th>lib82</th>
      <td>525</td>
      <td>75</td>
      <td>7.00</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>446</td>
      <td>81</td>
      <td>5.51</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU5r_Wenzhou1</th>
      <th>lib82</th>
      <td>458</td>
      <td>58</td>
      <td>7.90</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>480</td>
      <td>84</td>
      <td>5.71</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU5r_Wenzhou1_rpk</th>
      <th>lib82</th>
      <td>815</td>
      <td>133</td>
      <td>6.13</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>911</td>
      <td>125</td>
      <td>7.29</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU5r_YD13403</th>
      <th>lib82</th>
      <td>497</td>
      <td>78</td>
      <td>6.37</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>559</td>
      <td>83</td>
      <td>6.73</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU5r_YD13403_rpk</th>
      <th>lib82</th>
      <td>888</td>
      <td>119</td>
      <td>7.46</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>722</td>
      <td>110</td>
      <td>6.56</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HsCoV_Ita20664540</th>
      <th>lib82</th>
      <td>756</td>
      <td>105</td>
      <td>7.20</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>679</td>
      <td>97</td>
      <td>7.00</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HsCoV_Ita20664540_rpk</th>
      <th>lib82</th>
      <td>679</td>
      <td>100</td>
      <td>6.79</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>426</td>
      <td>84</td>
      <td>5.07</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">KW2EF93Gha2010</th>
      <th>lib82</th>
      <td>596</td>
      <td>68</td>
      <td>8.76</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>430</td>
      <td>63</td>
      <td>6.83</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">KW2EF93Gha2010_rpk</th>
      <th>lib82</th>
      <td>457</td>
      <td>64</td>
      <td>7.14</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>378</td>
      <td>65</td>
      <td>5.82</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">MjHKU4r1</th>
      <th>lib82</th>
      <td>724</td>
      <td>93</td>
      <td>7.78</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>644</td>
      <td>99</td>
      <td>6.51</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">MjHKU4r1_rpk</th>
      <th>lib82</th>
      <td>759</td>
      <td>109</td>
      <td>6.96</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>714</td>
      <td>118</td>
      <td>6.05</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">NL13845</th>
      <th>lib82</th>
      <td>582</td>
      <td>95</td>
      <td>6.13</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>637</td>
      <td>89</td>
      <td>7.16</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">NL13845_rpk</th>
      <th>lib82</th>
      <td>503</td>
      <td>64</td>
      <td>7.86</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>295</td>
      <td>60</td>
      <td>4.92</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">NL13892</th>
      <th>lib82</th>
      <td>719</td>
      <td>106</td>
      <td>6.78</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>795</td>
      <td>120</td>
      <td>6.62</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">NL13892_rpk</th>
      <th>lib82</th>
      <td>579</td>
      <td>86</td>
      <td>6.73</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>857</td>
      <td>110</td>
      <td>7.79</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">NL140422</th>
      <th>lib82</th>
      <td>541</td>
      <td>85</td>
      <td>6.36</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>757</td>
      <td>96</td>
      <td>7.89</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">NL140422_rpk</th>
      <th>lib82</th>
      <td>625</td>
      <td>91</td>
      <td>6.87</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>670</td>
      <td>102</td>
      <td>6.57</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">NL140455</th>
      <th>lib82</th>
      <td>1045</td>
      <td>120</td>
      <td>8.71</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>805</td>
      <td>126</td>
      <td>6.39</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">NL140455_rpk</th>
      <th>lib82</th>
      <td>609</td>
      <td>105</td>
      <td>5.80</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>387</td>
      <td>83</td>
      <td>4.66</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">NeoCoV</th>
      <th>lib82</th>
      <td>771</td>
      <td>102</td>
      <td>7.56</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>545</td>
      <td>94</td>
      <td>5.80</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">NeoCoV_rpk</th>
      <th>lib82</th>
      <td>412</td>
      <td>74</td>
      <td>5.57</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>865</td>
      <td>117</td>
      <td>7.39</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">PDF2180</th>
      <th>lib82</th>
      <td>503</td>
      <td>80</td>
      <td>6.29</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>731</td>
      <td>112</td>
      <td>6.53</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">PDF2180_rpk</th>
      <th>lib82</th>
      <td>784</td>
      <td>118</td>
      <td>6.64</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>740</td>
      <td>112</td>
      <td>6.61</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">PaGB01</th>
      <th>lib82</th>
      <td>901</td>
      <td>112</td>
      <td>8.04</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>1075</td>
      <td>146</td>
      <td>7.36</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">PaGB01_rpk</th>
      <th>lib82</th>
      <td>682</td>
      <td>83</td>
      <td>8.22</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>427</td>
      <td>74</td>
      <td>5.77</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">PkCoV_Ita20664563</th>
      <th>lib82</th>
      <td>583</td>
      <td>98</td>
      <td>5.95</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>525</td>
      <td>99</td>
      <td>5.30</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">PkCoV_Ita20664563_rpk</th>
      <th>lib82</th>
      <td>788</td>
      <td>91</td>
      <td>8.66</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>589</td>
      <td>107</td>
      <td>5.50</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">Russia_MOW1522</th>
      <th>lib82</th>
      <td>506</td>
      <td>59</td>
      <td>8.58</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>609</td>
      <td>82</td>
      <td>7.43</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">Russia_MOW1522_rpk</th>
      <th>lib82</th>
      <td>793</td>
      <td>98</td>
      <td>8.09</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>702</td>
      <td>96</td>
      <td>7.31</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">SC2013</th>
      <th>lib82</th>
      <td>562</td>
      <td>82</td>
      <td>6.85</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>420</td>
      <td>71</td>
      <td>5.92</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">SC2013_rpk</th>
      <th>lib82</th>
      <td>374</td>
      <td>49</td>
      <td>7.63</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>371</td>
      <td>52</td>
      <td>7.13</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">VsCoV1</th>
      <th>lib82</th>
      <td>512</td>
      <td>63</td>
      <td>8.13</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>346</td>
      <td>68</td>
      <td>5.09</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">VsCoV1_rpk</th>
      <th>lib82</th>
      <td>246</td>
      <td>47</td>
      <td>5.23</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>379</td>
      <td>58</td>
      <td>6.53</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">VsCoVa7</th>
      <th>lib82</th>
      <td>754</td>
      <td>103</td>
      <td>7.32</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>873</td>
      <td>147</td>
      <td>5.94</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">VsCoVa7_rpk</th>
      <th>lib82</th>
      <td>637</td>
      <td>77</td>
      <td>8.27</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>568</td>
      <td>101</td>
      <td>5.62</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">VsCoVkj15</th>
      <th>lib82</th>
      <td>649</td>
      <td>80</td>
      <td>8.11</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>516</td>
      <td>75</td>
      <td>6.88</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">VsCoVkj15_rpk</th>
      <th>lib82</th>
      <td>537</td>
      <td>80</td>
      <td>6.71</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>473</td>
      <td>81</td>
      <td>5.84</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">cMERS_ALA49341</th>
      <th>lib82</th>
      <td>808</td>
      <td>120</td>
      <td>6.73</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>614</td>
      <td>95</td>
      <td>6.46</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">cMERS_ALA49341_rpk</th>
      <th>lib82</th>
      <td>845</td>
      <td>116</td>
      <td>7.28</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>730</td>
      <td>115</td>
      <td>6.35</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">cMERS_ALA49374</th>
      <th>lib82</th>
      <td>316</td>
      <td>43</td>
      <td>7.35</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>247</td>
      <td>37</td>
      <td>6.68</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">cMERS_ALA49374_rpk</th>
      <th>lib82</th>
      <td>324</td>
      <td>50</td>
      <td>6.48</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>266</td>
      <td>48</td>
      <td>5.54</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">cMERS_AMO03401</th>
      <th>lib82</th>
      <td>458</td>
      <td>69</td>
      <td>6.64</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>588</td>
      <td>78</td>
      <td>7.54</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">cMERS_AMO03401_rpk</th>
      <th>lib82</th>
      <td>541</td>
      <td>110</td>
      <td>4.92</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>826</td>
      <td>112</td>
      <td>7.38</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">cMERS_ANI69835</th>
      <th>lib82</th>
      <td>825</td>
      <td>80</td>
      <td>10.31</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>493</td>
      <td>75</td>
      <td>6.57</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">cMERS_ANI69835_rpk</th>
      <th>lib82</th>
      <td>652</td>
      <td>97</td>
      <td>6.72</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>759</td>
      <td>107</td>
      <td>7.09</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">cMERS_ANI69889</th>
      <th>lib82</th>
      <td>645</td>
      <td>94</td>
      <td>6.86</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>679</td>
      <td>99</td>
      <td>6.86</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">cMERS_ANI69889_rpk</th>
      <th>lib82</th>
      <td>503</td>
      <td>62</td>
      <td>8.11</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>453</td>
      <td>65</td>
      <td>6.97</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">hMERS_AKM76239</th>
      <th>lib82</th>
      <td>703</td>
      <td>120</td>
      <td>5.86</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>999</td>
      <td>131</td>
      <td>7.63</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">hMERS_AKM76239_rpk</th>
      <th>lib82</th>
      <td>387</td>
      <td>62</td>
      <td>6.24</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>416</td>
      <td>62</td>
      <td>6.71</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">hMERS_AMQ48993</th>
      <th>lib82</th>
      <td>371</td>
      <td>57</td>
      <td>6.51</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>405</td>
      <td>51</td>
      <td>7.94</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">hMERS_AMQ48993_rpk</th>
      <th>lib82</th>
      <td>819</td>
      <td>114</td>
      <td>7.18</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>897</td>
      <td>117</td>
      <td>7.67</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">hMERS_EMC2012</th>
      <th>lib82</th>
      <td>405</td>
      <td>57</td>
      <td>7.11</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>306</td>
      <td>51</td>
      <td>6.00</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">hMERS_EMC2012_rpk</th>
      <th>lib82</th>
      <td>506</td>
      <td>73</td>
      <td>6.93</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>790</td>
      <td>86</td>
      <td>9.19</td>
    </tr>
  </tbody>
</table>


### Filter processed CCSs
We have the PacBio `ccs` program's estimated "accuracy" for both the barcode and the gene sequence for each processed CCS.
We will filter the CCSs to only keep ones of sufficiently high accuracy.

First, we want to plot the accuracies.
It is actually visually easier to look at the error rate, which is one minus the accuracy.
Because we want to plot on a log scale (which can't show error rates of zero), we define a *error_rate_floor*, and set all error rates less than this to that value:


```python
error_rate_floor = 1e-9  # error rates < this set to this
if error_rate_floor >= config['max_error_rate']:
    raise ValueError('error_rate_floor must be < max_error_rate')

processed_ccs = (
    processed_ccs
    .assign(barcode_error=lambda x: numpy.clip(1 - x['barcode_accuracy'],
                                               error_rate_floor, None),
            gene_error=lambda x: numpy.clip(1 - x['gene_accuracy'],
                                            error_rate_floor, None)
            )
    )
```

Now plot the error rates, drawing a dashed vertical line at the threshold separating the CCSs we retain for consensus building versus those that we discard:


```python
_ = (
 ggplot(processed_ccs
        .melt(value_vars=['barcode_error', 'gene_error'],
              var_name='feature_type', value_name='error rate'),
        aes('error rate')) +
 geom_histogram(bins=25) +
 geom_vline(xintercept=config['max_error_rate'],
            linetype='dashed',
            color=CBPALETTE[1]) +
 facet_wrap('~ feature_type') +
 theme(figure_size=(4.5, 2)) +
 ylab('number of CCSs') +
 scale_x_log10()
 ).draw()
```


    
![png](process_ccs_panmerbeco_files/process_ccs_panmerbeco_61_0.png)
    


Flag the CCSs to retain, and indicate how many we are retaining and purging due to the accuracy filter:


```python
processed_ccs = (
    processed_ccs
    .assign(retained=lambda x: ((x['gene_error'] < config['max_error_rate']) &
                                (x['barcode_error'] < config['max_error_rate'])))
    )
```

Here are number of retained CCSs:


```python
_ = (
 ggplot(processed_ccs.assign(xlabel=lambda x: x['target'] + ', ' + x['library'])
                     .groupby(['xlabel', 'retained'])
                     .size()
                     .rename('count')
                     .reset_index(),
        aes('count', 'xlabel', color='retained', label='count')) +
 geom_point(size=3) +
 geom_text(va='bottom', size=7, ha='center',format_string='{:.3g}', nudge_y=0.2) +
 theme(figure_size=(3, 0.5 * nlibs * ntargets),
       panel_grid_major_y=element_blank(),
       axis_text_x=element_text(angle=90),
       ) +
 scale_x_log10(name='number of CCSs') +
 ylab('') +
 scale_color_manual(values=CBPALETTE[1:])
 ).draw()
```


    
![png](process_ccs_panmerbeco_files/process_ccs_panmerbeco_65_0.png)
    


### Sequences per barcode
How many times is each barcode sequenced?
This is useful to know for thinking about building the barcode consensus.

Plot the distribution of the number of times each **barcode** is observed among the retained CCSs:


```python
max_count = 8 # in plot, group all barcodes with >= this many counts

p = (
 ggplot(
    processed_ccs
     .query('retained')
     .groupby(['library', 'barcode'])
     .size()
     .rename('nseqs')
     .reset_index()
     .assign(nseqs=lambda x: numpy.clip(x['nseqs'], None, max_count)),
    aes('nseqs')) +
 geom_bar() +
 facet_wrap('~ library', nrow=1) +
 theme(figure_size=(1.75 * nlibs, 2),
       panel_grid_major_x=element_blank(),
       ) +
 ylab('number of barcodes') +
 xlab('CCSs for barcode')
 )

_ = p.draw()
```


    
![png](process_ccs_panmerbeco_files/process_ccs_panmerbeco_67_0.png)
    


### Empirical accuracy of CCSs
We want to directly estimate the accuracy of the gene-barcode link rather than relying on the PacBio `ccs` accuracy, which doesn't include inaccuracies due to things like strand exchange or the same barcode on different sequences.

One way to do this is to examine instances when we have multiple sequences for the same barcode. 
We can calculate the empirical accuracy of the sequences by looking at all instances of multiple sequences of the same barcode and determining how often they are identical.
This calculation is performed by `alignparse.consensus.empirical_accuracy` using the equations described in the docs for that function.

We will do this four for sets of sequences:

 1. All of the CCSs retained above.
 2. CCSs retained by applying a PacBio `ccs` accuracy filter 10-fold more stringent than the one above.
    The rationale is that if this improves the concordance (real accuracy) of the CCSs substantially then maybe we should make the accuracy filter more stringent.
 3. Like (1) but excluding all CCSs with indels.
    the rationale is that we only really care about substitutions, and will exclude sequences with indels anyway.
 4. Like (2) but excluding all CCSs with indels.
 
First, we annotate the sequences with the number of indels and whether they have an indel to enable categorization into the aforementioned sets:


```python
processed_ccs = alignparse.consensus.add_mut_info_cols(processed_ccs,
                                                       mutation_col='gene_mutations',
                                                       n_indel_col='n_indels')

processed_ccs = processed_ccs.assign(has_indel=lambda x: x['n_indels'] > 0)

processed_ccs.head(n=12)
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>library</th>
      <th>query_name</th>
      <th>gene_mutations</th>
      <th>gene_accuracy</th>
      <th>barcode</th>
      <th>barcode_accuracy</th>
      <th>target</th>
      <th>barcode_error</th>
      <th>gene_error</th>
      <th>retained</th>
      <th>n_indels</th>
      <th>has_indel</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>lib82</td>
      <td>m64296e_240129_011756/262538/ccs</td>
      <td></td>
      <td>1.000000</td>
      <td>TTAATAATGTCATATT</td>
      <td>1.000000</td>
      <td>AncHKU25d_MAP</td>
      <td>1.000000e-09</td>
      <td>1.000000e-09</td>
      <td>True</td>
      <td>0</td>
      <td>False</td>
    </tr>
    <tr>
      <th>1</th>
      <td>lib82</td>
      <td>m64296e_240129_011756/329537/ccs</td>
      <td></td>
      <td>1.000000</td>
      <td>TTCCGCCACAGACCTG</td>
      <td>1.000000</td>
      <td>AncHKU25d_MAP</td>
      <td>2.039792e-09</td>
      <td>1.017722e-08</td>
      <td>True</td>
      <td>0</td>
      <td>False</td>
    </tr>
    <tr>
      <th>2</th>
      <td>lib82</td>
      <td>m64296e_240129_011756/527062/ccs</td>
      <td></td>
      <td>1.000000</td>
      <td>ATACTCATTCCCCAAA</td>
      <td>1.000000</td>
      <td>AncHKU25d_MAP</td>
      <td>1.000000e-09</td>
      <td>1.000000e-09</td>
      <td>True</td>
      <td>0</td>
      <td>False</td>
    </tr>
    <tr>
      <th>3</th>
      <td>lib82</td>
      <td>m64296e_240129_011756/721871/ccs</td>
      <td></td>
      <td>1.000000</td>
      <td>CTACGTAGATGGTAGT</td>
      <td>1.000000</td>
      <td>AncHKU25d_MAP</td>
      <td>1.000000e-09</td>
      <td>1.000000e-09</td>
      <td>True</td>
      <td>0</td>
      <td>False</td>
    </tr>
    <tr>
      <th>4</th>
      <td>lib82</td>
      <td>m64296e_240129_011756/788119/ccs</td>
      <td></td>
      <td>1.000000</td>
      <td>AACAATCCAAGAGATG</td>
      <td>1.000000</td>
      <td>AncHKU25d_MAP</td>
      <td>1.000000e-09</td>
      <td>1.000000e-09</td>
      <td>True</td>
      <td>0</td>
      <td>False</td>
    </tr>
    <tr>
      <th>5</th>
      <td>lib82</td>
      <td>m64296e_240129_011756/984844/ccs</td>
      <td></td>
      <td>0.999997</td>
      <td>TTCCACTTTTCCTCAC</td>
      <td>0.999988</td>
      <td>AncHKU25d_MAP</td>
      <td>1.249406e-05</td>
      <td>2.944244e-06</td>
      <td>False</td>
      <td>0</td>
      <td>False</td>
    </tr>
    <tr>
      <th>6</th>
      <td>lib82</td>
      <td>m64296e_240129_011756/1115829/ccs</td>
      <td></td>
      <td>1.000000</td>
      <td>AACAATCCAAGAGATG</td>
      <td>0.999975</td>
      <td>AncHKU25d_MAP</td>
      <td>2.488217e-05</td>
      <td>9.471878e-09</td>
      <td>False</td>
      <td>0</td>
      <td>False</td>
    </tr>
    <tr>
      <th>7</th>
      <td>lib82</td>
      <td>m64296e_240129_011756/1311436/ccs</td>
      <td>del40to40</td>
      <td>0.999188</td>
      <td>CGGTTGCAGTCACGAT</td>
      <td>1.000000</td>
      <td>AncHKU25d_MAP</td>
      <td>1.000000e-09</td>
      <td>8.122975e-04</td>
      <td>False</td>
      <td>1</td>
      <td>True</td>
    </tr>
    <tr>
      <th>8</th>
      <td>lib82</td>
      <td>m64296e_240129_011756/1769708/ccs</td>
      <td></td>
      <td>1.000000</td>
      <td>TGTTCGCGAGAGATCC</td>
      <td>1.000000</td>
      <td>AncHKU25d_MAP</td>
      <td>1.000000e-09</td>
      <td>1.000000e-09</td>
      <td>True</td>
      <td>0</td>
      <td>False</td>
    </tr>
    <tr>
      <th>9</th>
      <td>lib82</td>
      <td>m64296e_240129_011756/2362111/ccs</td>
      <td></td>
      <td>1.000000</td>
      <td>AGTACAGATGTAAATT</td>
      <td>1.000000</td>
      <td>AncHKU25d_MAP</td>
      <td>1.000000e-09</td>
      <td>1.801752e-09</td>
      <td>True</td>
      <td>0</td>
      <td>False</td>
    </tr>
    <tr>
      <th>10</th>
      <td>lib82</td>
      <td>m64296e_240129_011756/2425688/ccs</td>
      <td>T582G G583A T585C A586G A594T A597T C600T A603...</td>
      <td>1.000000</td>
      <td>CATTTCTTATTAAACT</td>
      <td>1.000000</td>
      <td>AncHKU25d_MAP</td>
      <td>1.000000e-09</td>
      <td>1.000000e-09</td>
      <td>True</td>
      <td>0</td>
      <td>False</td>
    </tr>
    <tr>
      <th>11</th>
      <td>lib82</td>
      <td>m64296e_240129_011756/2688099/ccs</td>
      <td></td>
      <td>1.000000</td>
      <td>CCTGTCAGAACAGCAA</td>
      <td>1.000000</td>
      <td>AncHKU25d_MAP</td>
      <td>1.000000e-09</td>
      <td>1.000000e-09</td>
      <td>True</td>
      <td>0</td>
      <td>False</td>
    </tr>
  </tbody>
</table>
</div>



Plot how many sequences have indels:


```python
_ = (
 ggplot(processed_ccs,
        aes('retained', fill='has_indel')) +
 geom_bar(position='dodge') +
 geom_text(aes(label='..count..'), stat='count', va='bottom', size=7,
           position=position_dodge(width=0.9), format_string='{:.2g}') +
 theme(figure_size=(2.5 * nlibs, 3),
       panel_grid_major_x=element_blank(),
       ) +
 ylab('number of CCSs') +
 scale_fill_manual(values=CBPALETTE[1:]) +
 facet_wrap('~ library', nrow=1)
 ).draw()
```


    
![png](process_ccs_panmerbeco_files/process_ccs_panmerbeco_71_0.png)
    


Now get the empirical accuracy for each of the CCS groups mentioned above:


```python
high_acc = config['max_error_rate'] / 10
empirical_acc = []

for desc, query_str in [
        ('retained', 'retained'),
        ('retained, no indel', 'retained and not has_indel'),
        ('10X accuracy',
         f"(gene_error < {high_acc}) and (barcode_error < {high_acc})"),
        ('10X accuracy, no indel',
         f"(gene_error < {high_acc}) and (barcode_error < {high_acc}) and not has_indel")
        ]:
    # get just CCSs in that category
    df = processed_ccs.query(query_str)
    
    # compute empirical accuracy
    empirical_acc.append(
        alignparse.consensus.empirical_accuracy(df,
                                                mutation_col='gene_mutations')
        .assign(description=desc)
        .merge(df
               .groupby('library')
               .size()
               .rename('number_CCSs')
               .reset_index()
               )
        )

# make description categorical to preserve order, and annotate as "actual"
# the category ("retained, no indel") that we will use for building variants.
empirical_acc = (
    pd.concat(empirical_acc, ignore_index=True, sort=False)
    .assign(description=lambda x: pd.Categorical(x['description'],
                                                 x['description'].unique(),
                                                 ordered=True),
            actual=lambda x: numpy.where(x['description'] == 'retained, no indel',
                                         True, False),
            )
    )
```

Display table of the empirical accuracies:


```python
display(HTML(empirical_acc.to_html(index=False)))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>library</th>
      <th>accuracy</th>
      <th>description</th>
      <th>number_CCSs</th>
      <th>actual</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>lib82</td>
      <td>0.993468</td>
      <td>retained</td>
      <td>98864</td>
      <td>False</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>0.992752</td>
      <td>retained</td>
      <td>100851</td>
      <td>False</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>0.998444</td>
      <td>retained, no indel</td>
      <td>95955</td>
      <td>True</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>0.998260</td>
      <td>retained, no indel</td>
      <td>97683</td>
      <td>True</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>0.994350</td>
      <td>10X accuracy</td>
      <td>93286</td>
      <td>False</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>0.993685</td>
      <td>10X accuracy</td>
      <td>95061</td>
      <td>False</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>0.998478</td>
      <td>10X accuracy, no indel</td>
      <td>90622</td>
      <td>False</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>0.998299</td>
      <td>10X accuracy, no indel</td>
      <td>92164</td>
      <td>False</td>
    </tr>
  </tbody>
</table>


Plot the empirical accuracies, using a different color to show the category that we will actually use:


```python
p = (
    ggplot(empirical_acc,
           aes('description', 'accuracy', color='actual', label='accuracy')
           ) +
    geom_point(size=3) +
    geom_text(va='bottom', size=9, format_string='{:.3g}', nudge_y=0.003) +
    facet_wrap('~ library', ncol=8) +
    theme(figure_size=(1.75 * nlibs, 2.25),
          axis_text_x=element_text(angle=90),
          panel_grid_major_x=element_blank(),
          ) +
    xlab('') +
    scale_y_continuous(name='empirical accuracy', limits=(0.9, 1.005)) +
    scale_color_manual(values=CBPALETTE, guide=False)
    )

plotfile = os.path.join(config['figs_dir'], 'empirical_CCS_accuracy.pdf')
print(f"Saving plot to {plotfile}")
_ = p.draw()
```

    Saving plot to results/figures/empirical_CCS_accuracy.pdf



    
![png](process_ccs_panmerbeco_files/process_ccs_panmerbeco_77_1.png)
    


The above analysis shows that if we exclude sequences with indels (which we plan to do among our consensus sequences), then the accuracy of each CCS is around 99%. 
We do **not** get notably higher empirical accuracy by imposing a more stringent filter from the PacBio `ccs` program, indicating that the major sources of error are due to processes that are not modeled in this program's accuracy filter (perhaps strand exchange or barcode sharing).

Note that this empirical accuracy is for a **single** CCS.
When we build the consensus sequences for each barcode below, we will take the consensus of CCSs within a barcode.
So for barcodes with multiple CCSs, the actual accuracy of the consensus sequences will be higher than the empirical accuracy above due to capturing information from multiple CCSs.

### Consensus sequences for barcodes
We call the consensus sequence for each barcode using the simple method implemented in [alignparse.consensus.simple_mutconsensus](https://jbloomlab.github.io/alignparse/alignparse.consensus.html?highlight=simple_mutconsensus#alignparse.consensus.simple_mutconsensus).
The documentation for that function explains the method in detail, but basically it works like this:
 1. When there is just one CCS per barcode, the consensus is just that sequence.
 2. When there are multiple CCSs per barcode, they are used to build a consensus--however, the entire barcode is discarded if there are many differences between CCSs with the barcode, or high-frequency non-consensus mutations. The reason that barcodes are discarded in such cases as many differences between CCSs or high-frequency non-consensus mutations suggest errors such as barcode collisions or strand exchange.
 
First, call the consensus for each barcode including **all** retained sequences, even those with undesirable indels.


```python
consensus, dropped = alignparse.consensus.simple_mutconsensus(
                        processed_ccs.query('retained'),
                        group_cols=('library', 'barcode', 'target'),
                        mutation_col='gene_mutations',
                        )
```

Here are the first few lines of the data frame of consensus sequences for each barcode.
In addition to giving the library, barcode, target, and mutations, it also has a column indicating how many CCSs support the variant call:


```python
display(HTML(consensus.head().to_html(index=False)))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>library</th>
      <th>barcode</th>
      <th>target</th>
      <th>gene_mutations</th>
      <th>variant_call_support</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>lib82</td>
      <td>AAAAAAAAAAATGTTT</td>
      <td>HKU5r_BY140535_rpk</td>
      <td></td>
      <td>1</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AAAAAAAAAAGTACCG</td>
      <td>EriCoV_Ita5026517_rpk</td>
      <td></td>
      <td>6</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AAAAAAAAGTCTTGAG</td>
      <td>PkCoV_Ita20664563</td>
      <td></td>
      <td>2</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AAAAAAACAATGGACC</td>
      <td>EriCoV_Ita5026519</td>
      <td></td>
      <td>4</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AAAAAAACACATGAAA</td>
      <td>cMERS_ANI69835</td>
      <td></td>
      <td>3</td>
    </tr>
  </tbody>
</table>


Since we retain variants with substitutions, add information about substitution mutations, and for troubleshooting, number of indels:


```python
consensus = alignparse.consensus.add_mut_info_cols(
                    consensus,
                    mutation_col='gene_mutations',
                    sub_str_col='substitutions',
                    n_indel_col='number_of_indels',
                    n_sub_col='number_of_substitutions',
                    overwrite_cols=True)

display(HTML(consensus.head().to_html(index=False)))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>library</th>
      <th>barcode</th>
      <th>target</th>
      <th>gene_mutations</th>
      <th>variant_call_support</th>
      <th>substitutions</th>
      <th>number_of_substitutions</th>
      <th>number_of_indels</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>lib82</td>
      <td>AAAAAAAAAAATGTTT</td>
      <td>HKU5r_BY140535_rpk</td>
      <td></td>
      <td>1</td>
      <td></td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AAAAAAAAAAGTACCG</td>
      <td>EriCoV_Ita5026517_rpk</td>
      <td></td>
      <td>6</td>
      <td></td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AAAAAAAAGTCTTGAG</td>
      <td>PkCoV_Ita20664563</td>
      <td></td>
      <td>2</td>
      <td></td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AAAAAAACAATGGACC</td>
      <td>EriCoV_Ita5026519</td>
      <td></td>
      <td>4</td>
      <td></td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AAAAAAACACATGAAA</td>
      <td>cMERS_ANI69835</td>
      <td></td>
      <td>3</td>
      <td></td>
      <td>0</td>
      <td>0</td>
    </tr>
  </tbody>
</table>


Plot distribution of number of CCSs supporting each variant call (consensus), indicating whether or not there is an indel:


```python
max_variant_call_support = 6  # group variants with >= this much support

_ = (
 ggplot(consensus
        .assign(variant_call_support=lambda x: numpy.clip(x['variant_call_support'],
                                                          None,
                                                          max_variant_call_support),
                indel_state=lambda x: numpy.where(x['number_of_indels'] > 0,
                                                  'has indel', 'no indel')
                ),
        aes('variant_call_support')) +
 geom_bar() +
 ylab('number of variants') +
 facet_grid('indel_state ~ library') +
 theme(figure_size=(1.75 * nlibs, 3.5),
       panel_grid_major_x=element_blank(),
       ) 
 ).draw()
```


    
![png](process_ccs_panmerbeco_files/process_ccs_panmerbeco_86_0.png)
    


We see that most variant consensus sequences do **not** have indels, especially if we limit to the more "accurate" ones that have multiple CCSs supporting them.

We will ignore all consensus sequences with indels in the variant-barcode lookup table. 
We do this for two reasons:
 1. When there is just one CCS supporting a consensus, it is less likely to be accurate as indels are the main mode of PacBio error.
 2. For the purposes of our studies, we are interested in point mutations rather than indels anyway.
 
Here are number of valid consensus sequence (no indels) for each library and target:


```python
consensus = consensus.query('number_of_indels < 1')

lib_target_counts = (
    consensus
    .groupby(['library', 'target'])
    .size()
    .rename('consensus sequences')
    .reset_index()
    )

display(HTML(lib_target_counts.to_html(index=False)))

p = (ggplot(lib_target_counts.assign(xlabel=lambda x: x['target'] + ', ' + x['library']),
            aes('xlabel', 'consensus sequences')) +
     geom_point(size=3) +
     theme(figure_size=(0.5 * nlibs * ntargets, 1.75),
           axis_text_x=element_text(angle=90)) +
     xlab('') +
     scale_y_log10()
     )

_ = p.draw()
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>library</th>
      <th>target</th>
      <th>consensus sequences</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>lib82</td>
      <td>AncEri_MAP</td>
      <td>78</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncEri_MAP_rpk</td>
      <td>67</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncEri_alt</td>
      <td>23</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncEri_alt_rpk</td>
      <td>55</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncHKU25HKU5_MAP</td>
      <td>119</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncHKU25HKU5_MAP_rpk</td>
      <td>37</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncHKU25HKU5_alt</td>
      <td>77</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncHKU25HKU5_alt_rpk</td>
      <td>66</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncHKU25_MAP</td>
      <td>64</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncHKU25_MAP_rpk</td>
      <td>57</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncHKU25_alt</td>
      <td>100</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncHKU25_alt_rpk</td>
      <td>126</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncHKU25a_MAP</td>
      <td>73</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncHKU25a_MAP_rpk</td>
      <td>108</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncHKU25a_alt</td>
      <td>54</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncHKU25a_alt_rpk</td>
      <td>87</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncHKU25b_MAP</td>
      <td>60</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncHKU25b_MAP_rpk</td>
      <td>115</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncHKU25b_alt</td>
      <td>100</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncHKU25b_alt_rpk</td>
      <td>58</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncHKU25c_MAP</td>
      <td>50</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncHKU25c_MAP_rpk</td>
      <td>46</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncHKU25c_alt</td>
      <td>93</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncHKU25c_alt_rpk</td>
      <td>103</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncHKU25d_MAP</td>
      <td>77</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncHKU25d_MAP_rpk</td>
      <td>90</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncHKU25d_alt</td>
      <td>71</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncHKU25d_alt_rpk</td>
      <td>106</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncHKU25e_MAP</td>
      <td>49</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncHKU25e_MAP_rpk</td>
      <td>52</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncHKU25e_alt</td>
      <td>97</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncHKU25e_alt_rpk</td>
      <td>56</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncHKU25f_MAP</td>
      <td>100</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncHKU25f_MAP_rpk</td>
      <td>85</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncHKU25f_alt</td>
      <td>61</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncHKU25f_alt_rpk</td>
      <td>72</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncHKU25g_MAP</td>
      <td>101</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncHKU25g_MAP_rpk</td>
      <td>111</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncHKU25g_alt</td>
      <td>35</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncHKU25g_alt_rpk</td>
      <td>63</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncHKU25h_MAP</td>
      <td>90</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncHKU25h_MAP_rpk</td>
      <td>56</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncHKU25h_alt</td>
      <td>111</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncHKU25h_alt_rpk</td>
      <td>34</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncHKU4_MAP</td>
      <td>54</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncHKU4_MAP_rpk</td>
      <td>80</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncHKU4_alt</td>
      <td>77</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncHKU4_alt_rpk</td>
      <td>100</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncHKU4a_MAP</td>
      <td>42</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncHKU4a_MAP_rpk</td>
      <td>34</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncHKU4a_alt</td>
      <td>61</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncHKU4a_alt_rpk</td>
      <td>49</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncHKU4b_MAP</td>
      <td>65</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncHKU4b_MAP_rpk</td>
      <td>67</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncHKU4b_alt</td>
      <td>85</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncHKU4b_alt_rpk</td>
      <td>101</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncHKU4c_MAP</td>
      <td>89</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncHKU4c_MAP_rpk</td>
      <td>96</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncHKU4c_alt</td>
      <td>105</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncHKU4c_alt_rpk</td>
      <td>69</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncHKU4d_MAP</td>
      <td>66</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncHKU4d_MAP_rpk</td>
      <td>47</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncHKU4d_alt</td>
      <td>76</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncHKU4d_alt_rpk</td>
      <td>83</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncHKU5_MAP</td>
      <td>74</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncHKU5_MAP_rpk</td>
      <td>71</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncHKU5_alt</td>
      <td>53</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncHKU5_alt_rpk</td>
      <td>97</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncHKU5a_MAP</td>
      <td>100</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncHKU5a_MAP_rpk</td>
      <td>75</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncHKU5a_alt</td>
      <td>51</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncHKU5a_alt_rpk</td>
      <td>70</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncHKU5b_MAP</td>
      <td>91</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncHKU5b_MAP_rpk</td>
      <td>66</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncHKU5b_alt</td>
      <td>87</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncHKU5b_alt_rpk</td>
      <td>107</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncMERSHKU4_MAP</td>
      <td>82</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncMERSHKU4_MAP_rpk</td>
      <td>87</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncMERSHKU4_alt</td>
      <td>99</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncMERSHKU4_alt_rpk</td>
      <td>43</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncMerbeco1_MAP</td>
      <td>77</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncMerbeco1_MAP_rpk</td>
      <td>73</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncMerbeco1_alt</td>
      <td>76</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncMerbeco1_alt_rpk</td>
      <td>52</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncMerbeco2_MAP</td>
      <td>23</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncMerbeco2_MAP_rpk</td>
      <td>63</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncMerbeco2_alt</td>
      <td>60</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncMerbeco2_alt_rpk</td>
      <td>105</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncMerbeco3_MAP</td>
      <td>62</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncMerbeco3_MAP_rpk</td>
      <td>54</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncMerbeco3_alt</td>
      <td>93</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncMerbeco3_alt_rpk</td>
      <td>106</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncNeoEri_MAP</td>
      <td>38</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncNeoEri_MAP_rpk</td>
      <td>9</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncNeoEri_alt</td>
      <td>80</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncNeoEri_alt_rpk</td>
      <td>94</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncNeoGhEri_MAP</td>
      <td>107</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncNeoGhEri_MAP_rpk</td>
      <td>68</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncNeoGhEri_alt</td>
      <td>43</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncNeoGhEri_alt_rpk</td>
      <td>52</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncNeo_MAP</td>
      <td>71</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncNeo_MAP_rpk</td>
      <td>68</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncNeo_alt</td>
      <td>97</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AncNeo_alt_rpk</td>
      <td>68</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>BatCoV_SWI1</td>
      <td>63</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>BatCoV_SWI1_rpk</td>
      <td>69</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>BatCoV_SWI2</td>
      <td>93</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>BatCoV_SWI2_rpk</td>
      <td>87</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>BtCoV_B12080</td>
      <td>94</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>BtCoV_B12080_rpk</td>
      <td>50</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>EjCoV3</td>
      <td>108</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>EjCoV3_rpk</td>
      <td>76</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>EriCoV_Deu174</td>
      <td>24</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>EriCoV_Deu174_rpk</td>
      <td>55</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>EriCoV_GB1</td>
      <td>99</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>EriCoV_GB1_rpk</td>
      <td>90</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>EriCoV_HKU31</td>
      <td>34</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>EriCoV_HKU31_rpk</td>
      <td>48</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>EriCoV_Ita1169881</td>
      <td>64</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>EriCoV_Ita1169881_rpk</td>
      <td>25</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>EriCoV_Ita502651</td>
      <td>28</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>EriCoV_Ita5026511</td>
      <td>58</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>EriCoV_Ita5026511_rpk</td>
      <td>8</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>EriCoV_Ita5026512</td>
      <td>62</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>EriCoV_Ita5026512_rpk</td>
      <td>22</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>EriCoV_Ita5026515</td>
      <td>66</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>EriCoV_Ita5026515_rpk</td>
      <td>97</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>EriCoV_Ita5026517</td>
      <td>85</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>EriCoV_Ita5026517_rpk</td>
      <td>104</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>EriCoV_Ita5026519</td>
      <td>67</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>EriCoV_Ita5026519_rpk</td>
      <td>76</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>EriCoV_Ita502651_rpk</td>
      <td>19</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>HKU25_NL140462</td>
      <td>104</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>HKU25_NL140462_rpk</td>
      <td>79</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>HKU25_YD131305</td>
      <td>90</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>HKU25_YD131305_rpk</td>
      <td>80</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>HKU4</td>
      <td>94</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>HKU4_rpk</td>
      <td>64</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>HKU4r_161028</td>
      <td>38</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>HKU4r_161028_rpk</td>
      <td>92</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>HKU4r_162275</td>
      <td>57</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>HKU4r_162275_rpk</td>
      <td>79</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>HKU4r_BtCoV_1332005</td>
      <td>65</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>HKU4r_BtCoV_1332005_rpk</td>
      <td>120</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>HKU4r_CZ01</td>
      <td>45</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>HKU4r_CZ01_rpk</td>
      <td>83</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>HKU4r_GZ160421</td>
      <td>88</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>HKU4r_GZ160421_rpk</td>
      <td>97</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>HKU4r_GZ1862</td>
      <td>74</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>HKU4r_GZ1862_rpk</td>
      <td>57</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>HKU4r_JPDB144</td>
      <td>87</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>HKU4r_JPDB144_rpk</td>
      <td>87</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>HKU4r_panP251T</td>
      <td>66</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>HKU4r_panP251T_rpk</td>
      <td>91</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>HKU5</td>
      <td>79</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>HKU5_rpk</td>
      <td>90</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>HKU5r_17S</td>
      <td>79</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>HKU5r_17S_rpk</td>
      <td>81</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>HKU5r_19S</td>
      <td>104</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>HKU5r_19S_rpk</td>
      <td>88</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>HKU5r_2</td>
      <td>88</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>HKU5r_20S</td>
      <td>89</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>HKU5r_20S_rpk</td>
      <td>75</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>HKU5r_21S</td>
      <td>48</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>HKU5r_21S_rpk</td>
      <td>100</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>HKU5r_24S</td>
      <td>66</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>HKU5r_24S_rpk</td>
      <td>83</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>HKU5r_27S</td>
      <td>61</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>HKU5r_27S_rpk</td>
      <td>69</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>HKU5r_2_rpk</td>
      <td>106</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>HKU5r_30S</td>
      <td>83</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>HKU5r_30S_rpk</td>
      <td>79</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>HKU5r_31S</td>
      <td>52</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>HKU5r_31S_rpk</td>
      <td>97</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>HKU5r_32S</td>
      <td>71</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>HKU5r_32S_rpk</td>
      <td>68</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>HKU5r_33S</td>
      <td>100</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>HKU5r_33S_rpk</td>
      <td>66</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>HKU5r_A4342005</td>
      <td>104</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>HKU5r_A4342005_rpk</td>
      <td>101</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>HKU5r_BY140535</td>
      <td>65</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>HKU5r_BY140535_rpk</td>
      <td>88</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>HKU5r_BY140562</td>
      <td>90</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>HKU5r_BY140562_rpk</td>
      <td>90</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>HKU5r_BY140568</td>
      <td>70</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>HKU5r_BY140568_rpk</td>
      <td>69</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>HKU5r_Wenzhou1</td>
      <td>49</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>HKU5r_Wenzhou1_rpk</td>
      <td>116</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>HKU5r_YD13403</td>
      <td>67</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>HKU5r_YD13403_rpk</td>
      <td>114</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>HsCoV_Ita20664540</td>
      <td>98</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>HsCoV_Ita20664540_rpk</td>
      <td>92</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>KW2EF93Gha2010</td>
      <td>58</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>KW2EF93Gha2010_rpk</td>
      <td>58</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>MjHKU4r1</td>
      <td>81</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>MjHKU4r1_rpk</td>
      <td>95</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>NL13845</td>
      <td>87</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>NL13845_rpk</td>
      <td>60</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>NL13892</td>
      <td>91</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>NL13892_rpk</td>
      <td>81</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>NL140422</td>
      <td>79</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>NL140422_rpk</td>
      <td>82</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>NL140455</td>
      <td>109</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>NL140455_rpk</td>
      <td>88</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>NeoCoV</td>
      <td>87</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>NeoCoV_rpk</td>
      <td>65</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>PDF2180</td>
      <td>73</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>PDF2180_rpk</td>
      <td>110</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>PaGB01</td>
      <td>101</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>PaGB01_rpk</td>
      <td>77</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>PkCoV_Ita20664563</td>
      <td>85</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>PkCoV_Ita20664563_rpk</td>
      <td>86</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>Russia_MOW1522</td>
      <td>58</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>Russia_MOW1522_rpk</td>
      <td>92</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>SC2013</td>
      <td>75</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>SC2013_rpk</td>
      <td>43</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>VsCoV1</td>
      <td>60</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>VsCoV1_rpk</td>
      <td>43</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>VsCoVa7</td>
      <td>93</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>VsCoVa7_rpk</td>
      <td>70</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>VsCoVkj15</td>
      <td>76</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>VsCoVkj15_rpk</td>
      <td>74</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>cMERS_ALA49341</td>
      <td>107</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>cMERS_ALA49341_rpk</td>
      <td>101</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>cMERS_ALA49374</td>
      <td>38</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>cMERS_ALA49374_rpk</td>
      <td>47</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>cMERS_AMO03401</td>
      <td>67</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>cMERS_AMO03401_rpk</td>
      <td>101</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>cMERS_ANI69835</td>
      <td>74</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>cMERS_ANI69835_rpk</td>
      <td>90</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>cMERS_ANI69889</td>
      <td>87</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>cMERS_ANI69889_rpk</td>
      <td>53</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>hMERS_AKM76239</td>
      <td>106</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>hMERS_AKM76239_rpk</td>
      <td>53</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>hMERS_AMQ48993</td>
      <td>51</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>hMERS_AMQ48993_rpk</td>
      <td>101</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>hMERS_EMC2012</td>
      <td>48</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>hMERS_EMC2012_rpk</td>
      <td>68</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncEri_MAP</td>
      <td>76</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncEri_MAP_rpk</td>
      <td>52</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncEri_alt</td>
      <td>14</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncEri_alt_rpk</td>
      <td>76</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncHKU25HKU5_MAP</td>
      <td>84</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncHKU25HKU5_MAP_rpk</td>
      <td>48</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncHKU25HKU5_alt</td>
      <td>77</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncHKU25HKU5_alt_rpk</td>
      <td>75</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncHKU25_MAP</td>
      <td>41</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncHKU25_MAP_rpk</td>
      <td>53</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncHKU25_alt</td>
      <td>93</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncHKU25_alt_rpk</td>
      <td>126</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncHKU25a_MAP</td>
      <td>91</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncHKU25a_MAP_rpk</td>
      <td>116</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncHKU25a_alt</td>
      <td>75</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncHKU25a_alt_rpk</td>
      <td>104</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncHKU25b_MAP</td>
      <td>55</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncHKU25b_MAP_rpk</td>
      <td>111</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncHKU25b_alt</td>
      <td>118</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncHKU25b_alt_rpk</td>
      <td>81</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncHKU25c_MAP</td>
      <td>57</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncHKU25c_MAP_rpk</td>
      <td>51</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncHKU25c_alt</td>
      <td>72</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncHKU25c_alt_rpk</td>
      <td>104</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncHKU25d_MAP</td>
      <td>79</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncHKU25d_MAP_rpk</td>
      <td>80</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncHKU25d_alt</td>
      <td>95</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncHKU25d_alt_rpk</td>
      <td>77</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncHKU25e_MAP</td>
      <td>56</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncHKU25e_MAP_rpk</td>
      <td>60</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncHKU25e_alt</td>
      <td>102</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncHKU25e_alt_rpk</td>
      <td>35</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncHKU25f_MAP</td>
      <td>90</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncHKU25f_MAP_rpk</td>
      <td>130</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncHKU25f_alt</td>
      <td>75</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncHKU25f_alt_rpk</td>
      <td>88</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncHKU25g_MAP</td>
      <td>71</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncHKU25g_MAP_rpk</td>
      <td>94</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncHKU25g_alt</td>
      <td>49</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncHKU25g_alt_rpk</td>
      <td>71</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncHKU25h_MAP</td>
      <td>94</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncHKU25h_MAP_rpk</td>
      <td>59</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncHKU25h_alt</td>
      <td>84</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncHKU25h_alt_rpk</td>
      <td>66</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncHKU4_MAP</td>
      <td>59</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncHKU4_MAP_rpk</td>
      <td>102</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncHKU4_alt</td>
      <td>94</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncHKU4_alt_rpk</td>
      <td>86</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncHKU4a_MAP</td>
      <td>37</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncHKU4a_MAP_rpk</td>
      <td>52</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncHKU4a_alt</td>
      <td>52</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncHKU4a_alt_rpk</td>
      <td>70</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncHKU4b_MAP</td>
      <td>90</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncHKU4b_MAP_rpk</td>
      <td>83</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncHKU4b_alt</td>
      <td>83</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncHKU4b_alt_rpk</td>
      <td>81</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncHKU4c_MAP</td>
      <td>64</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncHKU4c_MAP_rpk</td>
      <td>98</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncHKU4c_alt</td>
      <td>99</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncHKU4c_alt_rpk</td>
      <td>88</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncHKU4d_MAP</td>
      <td>64</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncHKU4d_MAP_rpk</td>
      <td>53</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncHKU4d_alt</td>
      <td>79</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncHKU4d_alt_rpk</td>
      <td>85</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncHKU5_MAP</td>
      <td>64</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncHKU5_MAP_rpk</td>
      <td>84</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncHKU5_alt</td>
      <td>60</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncHKU5_alt_rpk</td>
      <td>113</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncHKU5a_MAP</td>
      <td>130</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncHKU5a_MAP_rpk</td>
      <td>81</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncHKU5a_alt</td>
      <td>57</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncHKU5a_alt_rpk</td>
      <td>76</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncHKU5b_MAP</td>
      <td>96</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncHKU5b_MAP_rpk</td>
      <td>45</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncHKU5b_alt</td>
      <td>94</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncHKU5b_alt_rpk</td>
      <td>108</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncMERSHKU4_MAP</td>
      <td>91</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncMERSHKU4_MAP_rpk</td>
      <td>103</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncMERSHKU4_alt</td>
      <td>96</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncMERSHKU4_alt_rpk</td>
      <td>66</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncMerbeco1_MAP</td>
      <td>59</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncMerbeco1_MAP_rpk</td>
      <td>80</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncMerbeco1_alt</td>
      <td>104</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncMerbeco1_alt_rpk</td>
      <td>69</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncMerbeco2_MAP</td>
      <td>27</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncMerbeco2_MAP_rpk</td>
      <td>75</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncMerbeco2_alt</td>
      <td>79</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncMerbeco2_alt_rpk</td>
      <td>95</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncMerbeco3_MAP</td>
      <td>78</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncMerbeco3_MAP_rpk</td>
      <td>66</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncMerbeco3_alt</td>
      <td>70</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncMerbeco3_alt_rpk</td>
      <td>117</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncNeoEri_MAP</td>
      <td>39</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncNeoEri_MAP_rpk</td>
      <td>10</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncNeoEri_alt</td>
      <td>97</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncNeoEri_alt_rpk</td>
      <td>97</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncNeoGhEri_MAP</td>
      <td>85</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncNeoGhEri_MAP_rpk</td>
      <td>76</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncNeoGhEri_alt</td>
      <td>60</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncNeoGhEri_alt_rpk</td>
      <td>37</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncNeo_MAP</td>
      <td>82</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncNeo_MAP_rpk</td>
      <td>77</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncNeo_alt</td>
      <td>86</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>AncNeo_alt_rpk</td>
      <td>77</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>BatCoV_SWI1</td>
      <td>54</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>BatCoV_SWI1_rpk</td>
      <td>84</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>BatCoV_SWI2</td>
      <td>88</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>BatCoV_SWI2_rpk</td>
      <td>95</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>BtCoV_B12080</td>
      <td>84</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>BtCoV_B12080_rpk</td>
      <td>76</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>EjCoV3</td>
      <td>126</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>EjCoV3_rpk</td>
      <td>84</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>EriCoV_Deu174</td>
      <td>36</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>EriCoV_Deu174_rpk</td>
      <td>67</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>EriCoV_GB1</td>
      <td>92</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>EriCoV_GB1_rpk</td>
      <td>88</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>EriCoV_HKU31</td>
      <td>60</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>EriCoV_HKU31_rpk</td>
      <td>48</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>EriCoV_Ita1169881</td>
      <td>76</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>EriCoV_Ita1169881_rpk</td>
      <td>23</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>EriCoV_Ita502651</td>
      <td>41</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>EriCoV_Ita5026511</td>
      <td>80</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>EriCoV_Ita5026511_rpk</td>
      <td>8</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>EriCoV_Ita5026512</td>
      <td>68</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>EriCoV_Ita5026512_rpk</td>
      <td>23</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>EriCoV_Ita5026515</td>
      <td>91</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>EriCoV_Ita5026515_rpk</td>
      <td>80</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>EriCoV_Ita5026517</td>
      <td>84</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>EriCoV_Ita5026517_rpk</td>
      <td>85</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>EriCoV_Ita5026519</td>
      <td>89</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>EriCoV_Ita5026519_rpk</td>
      <td>92</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>EriCoV_Ita502651_rpk</td>
      <td>5</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>HKU25_NL140462</td>
      <td>129</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>HKU25_NL140462_rpk</td>
      <td>62</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>HKU25_YD131305</td>
      <td>82</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>HKU25_YD131305_rpk</td>
      <td>97</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>HKU4</td>
      <td>108</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>HKU4_rpk</td>
      <td>72</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>HKU4r_161028</td>
      <td>46</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>HKU4r_161028_rpk</td>
      <td>101</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>HKU4r_162275</td>
      <td>78</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>HKU4r_162275_rpk</td>
      <td>88</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>HKU4r_BtCoV_1332005</td>
      <td>73</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>HKU4r_BtCoV_1332005_rpk</td>
      <td>96</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>HKU4r_CZ01</td>
      <td>72</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>HKU4r_CZ01_rpk</td>
      <td>63</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>HKU4r_GZ160421</td>
      <td>64</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>HKU4r_GZ160421_rpk</td>
      <td>85</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>HKU4r_GZ1862</td>
      <td>76</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>HKU4r_GZ1862_rpk</td>
      <td>67</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>HKU4r_JPDB144</td>
      <td>78</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>HKU4r_JPDB144_rpk</td>
      <td>117</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>HKU4r_panP251T</td>
      <td>70</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>HKU4r_panP251T_rpk</td>
      <td>93</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>HKU5</td>
      <td>88</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>HKU5_rpk</td>
      <td>83</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>HKU5r_17S</td>
      <td>82</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>HKU5r_17S_rpk</td>
      <td>98</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>HKU5r_19S</td>
      <td>87</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>HKU5r_19S_rpk</td>
      <td>90</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>HKU5r_2</td>
      <td>86</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>HKU5r_20S</td>
      <td>75</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>HKU5r_20S_rpk</td>
      <td>102</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>HKU5r_21S</td>
      <td>54</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>HKU5r_21S_rpk</td>
      <td>115</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>HKU5r_24S</td>
      <td>41</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>HKU5r_24S_rpk</td>
      <td>91</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>HKU5r_27S</td>
      <td>75</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>HKU5r_27S_rpk</td>
      <td>55</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>HKU5r_2_rpk</td>
      <td>114</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>HKU5r_30S</td>
      <td>96</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>HKU5r_30S_rpk</td>
      <td>64</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>HKU5r_31S</td>
      <td>49</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>HKU5r_31S_rpk</td>
      <td>93</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>HKU5r_32S</td>
      <td>90</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>HKU5r_32S_rpk</td>
      <td>71</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>HKU5r_33S</td>
      <td>90</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>HKU5r_33S_rpk</td>
      <td>75</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>HKU5r_A4342005</td>
      <td>87</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>HKU5r_A4342005_rpk</td>
      <td>115</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>HKU5r_BY140535</td>
      <td>70</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>HKU5r_BY140535_rpk</td>
      <td>66</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>HKU5r_BY140562</td>
      <td>109</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>HKU5r_BY140562_rpk</td>
      <td>94</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>HKU5r_BY140568</td>
      <td>94</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>HKU5r_BY140568_rpk</td>
      <td>75</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>HKU5r_Wenzhou1</td>
      <td>68</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>HKU5r_Wenzhou1_rpk</td>
      <td>109</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>HKU5r_YD13403</td>
      <td>75</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>HKU5r_YD13403_rpk</td>
      <td>97</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>HsCoV_Ita20664540</td>
      <td>88</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>HsCoV_Ita20664540_rpk</td>
      <td>71</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>KW2EF93Gha2010</td>
      <td>57</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>KW2EF93Gha2010_rpk</td>
      <td>57</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>MjHKU4r1</td>
      <td>88</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>MjHKU4r1_rpk</td>
      <td>103</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>NL13845</td>
      <td>81</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>NL13845_rpk</td>
      <td>53</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>NL13892</td>
      <td>101</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>NL13892_rpk</td>
      <td>101</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>NL140422</td>
      <td>86</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>NL140422_rpk</td>
      <td>94</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>NL140455</td>
      <td>115</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>NL140455_rpk</td>
      <td>72</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>NeoCoV</td>
      <td>77</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>NeoCoV_rpk</td>
      <td>101</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>PDF2180</td>
      <td>95</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>PDF2180_rpk</td>
      <td>101</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>PaGB01</td>
      <td>135</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>PaGB01_rpk</td>
      <td>69</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>PkCoV_Ita20664563</td>
      <td>88</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>PkCoV_Ita20664563_rpk</td>
      <td>97</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>Russia_MOW1522</td>
      <td>78</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>Russia_MOW1522_rpk</td>
      <td>88</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>SC2013</td>
      <td>64</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>SC2013_rpk</td>
      <td>47</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>VsCoV1</td>
      <td>58</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>VsCoV1_rpk</td>
      <td>49</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>VsCoVa7</td>
      <td>131</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>VsCoVa7_rpk</td>
      <td>90</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>VsCoVkj15</td>
      <td>62</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>VsCoVkj15_rpk</td>
      <td>71</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>cMERS_ALA49341</td>
      <td>85</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>cMERS_ALA49341_rpk</td>
      <td>96</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>cMERS_ALA49374</td>
      <td>34</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>cMERS_ALA49374_rpk</td>
      <td>41</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>cMERS_AMO03401</td>
      <td>69</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>cMERS_AMO03401_rpk</td>
      <td>102</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>cMERS_ANI69835</td>
      <td>62</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>cMERS_ANI69835_rpk</td>
      <td>96</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>cMERS_ANI69889</td>
      <td>90</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>cMERS_ANI69889_rpk</td>
      <td>57</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>hMERS_AKM76239</td>
      <td>114</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>hMERS_AKM76239_rpk</td>
      <td>59</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>hMERS_AMQ48993</td>
      <td>49</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>hMERS_AMQ48993_rpk</td>
      <td>103</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>hMERS_EMC2012</td>
      <td>45</td>
    </tr>
    <tr>
      <td>lib83</td>
      <td>hMERS_EMC2012_rpk</td>
      <td>74</td>
    </tr>
  </tbody>
</table>



    
![png](process_ccs_panmerbeco_files/process_ccs_panmerbeco_88_1.png)
    


We want to drop all barcodes with mutations:


```python
consensus = (
    consensus
    .assign(has_substitutions=lambda x: x['substitutions'].str.len().astype(bool))
    )

has_subs_by_target = (
        consensus
        .groupby(['target', 'library', 'has_substitutions'])
        .aggregate(n_barcodes=pd.NamedAgg('barcode', 'count'))
        .reset_index()
        )

display(HTML(has_subs_by_target
             .pivot_table(index=['target', 'library'],
                          columns='has_substitutions',
                          values='n_barcodes',
                          fill_value=0)
             .to_html()))

p = (ggplot(has_subs_by_target.assign(xlabel=lambda x: x['target'] + ', ' + x['library']),
            aes('xlabel', 'n_barcodes', color='has_substitutions')) +
     geom_point(size=3, alpha=0.7) +
     theme(figure_size=(0.5 * nlibs * ntargets, 1.75),
           axis_text_x=element_text(angle=90)) +
     xlab('') +
     scale_y_log10() +
     scale_color_manual(values=CBPALETTE)
     )

_ = p.draw()

```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>has_substitutions</th>
      <th>False</th>
      <th>True</th>
    </tr>
    <tr>
      <th>target</th>
      <th>library</th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th rowspan="2" valign="top">AncEri_MAP</th>
      <th>lib82</th>
      <td>74</td>
      <td>4</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>74</td>
      <td>2</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncEri_MAP_rpk</th>
      <th>lib82</th>
      <td>65</td>
      <td>2</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>47</td>
      <td>5</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncEri_alt</th>
      <th>lib82</th>
      <td>23</td>
      <td>0</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>11</td>
      <td>3</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncEri_alt_rpk</th>
      <th>lib82</th>
      <td>53</td>
      <td>2</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>72</td>
      <td>4</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU25HKU5_MAP</th>
      <th>lib82</th>
      <td>111</td>
      <td>8</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>79</td>
      <td>5</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU25HKU5_MAP_rpk</th>
      <th>lib82</th>
      <td>37</td>
      <td>0</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>48</td>
      <td>0</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU25HKU5_alt</th>
      <th>lib82</th>
      <td>76</td>
      <td>1</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>71</td>
      <td>6</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU25HKU5_alt_rpk</th>
      <th>lib82</th>
      <td>63</td>
      <td>3</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>70</td>
      <td>5</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU25_MAP</th>
      <th>lib82</th>
      <td>60</td>
      <td>4</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>40</td>
      <td>1</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU25_MAP_rpk</th>
      <th>lib82</th>
      <td>53</td>
      <td>4</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>51</td>
      <td>2</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU25_alt</th>
      <th>lib82</th>
      <td>96</td>
      <td>4</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>88</td>
      <td>5</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU25_alt_rpk</th>
      <th>lib82</th>
      <td>120</td>
      <td>6</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>124</td>
      <td>2</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU25a_MAP</th>
      <th>lib82</th>
      <td>71</td>
      <td>2</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>89</td>
      <td>2</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU25a_MAP_rpk</th>
      <th>lib82</th>
      <td>102</td>
      <td>6</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>109</td>
      <td>7</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU25a_alt</th>
      <th>lib82</th>
      <td>51</td>
      <td>3</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>74</td>
      <td>1</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU25a_alt_rpk</th>
      <th>lib82</th>
      <td>81</td>
      <td>6</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>103</td>
      <td>1</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU25b_MAP</th>
      <th>lib82</th>
      <td>58</td>
      <td>2</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>51</td>
      <td>4</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU25b_MAP_rpk</th>
      <th>lib82</th>
      <td>110</td>
      <td>5</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>103</td>
      <td>8</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU25b_alt</th>
      <th>lib82</th>
      <td>90</td>
      <td>10</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>116</td>
      <td>2</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU25b_alt_rpk</th>
      <th>lib82</th>
      <td>55</td>
      <td>3</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>80</td>
      <td>1</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU25c_MAP</th>
      <th>lib82</th>
      <td>48</td>
      <td>2</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>53</td>
      <td>4</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU25c_MAP_rpk</th>
      <th>lib82</th>
      <td>43</td>
      <td>3</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>47</td>
      <td>4</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU25c_alt</th>
      <th>lib82</th>
      <td>84</td>
      <td>9</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>72</td>
      <td>0</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU25c_alt_rpk</th>
      <th>lib82</th>
      <td>97</td>
      <td>6</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>98</td>
      <td>6</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU25d_MAP</th>
      <th>lib82</th>
      <td>75</td>
      <td>2</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>72</td>
      <td>7</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU25d_MAP_rpk</th>
      <th>lib82</th>
      <td>87</td>
      <td>3</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>79</td>
      <td>1</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU25d_alt</th>
      <th>lib82</th>
      <td>70</td>
      <td>1</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>90</td>
      <td>5</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU25d_alt_rpk</th>
      <th>lib82</th>
      <td>102</td>
      <td>4</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>70</td>
      <td>7</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU25e_MAP</th>
      <th>lib82</th>
      <td>47</td>
      <td>2</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>53</td>
      <td>3</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU25e_MAP_rpk</th>
      <th>lib82</th>
      <td>50</td>
      <td>2</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>59</td>
      <td>1</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU25e_alt</th>
      <th>lib82</th>
      <td>92</td>
      <td>5</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>98</td>
      <td>4</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU25e_alt_rpk</th>
      <th>lib82</th>
      <td>54</td>
      <td>2</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>32</td>
      <td>3</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU25f_MAP</th>
      <th>lib82</th>
      <td>94</td>
      <td>6</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>87</td>
      <td>3</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU25f_MAP_rpk</th>
      <th>lib82</th>
      <td>80</td>
      <td>5</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>121</td>
      <td>9</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU25f_alt</th>
      <th>lib82</th>
      <td>53</td>
      <td>8</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>72</td>
      <td>3</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU25f_alt_rpk</th>
      <th>lib82</th>
      <td>69</td>
      <td>3</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>85</td>
      <td>3</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU25g_MAP</th>
      <th>lib82</th>
      <td>94</td>
      <td>7</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>70</td>
      <td>1</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU25g_MAP_rpk</th>
      <th>lib82</th>
      <td>105</td>
      <td>6</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>93</td>
      <td>1</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU25g_alt</th>
      <th>lib82</th>
      <td>30</td>
      <td>5</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>46</td>
      <td>3</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU25g_alt_rpk</th>
      <th>lib82</th>
      <td>61</td>
      <td>2</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>67</td>
      <td>4</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU25h_MAP</th>
      <th>lib82</th>
      <td>83</td>
      <td>7</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>89</td>
      <td>5</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU25h_MAP_rpk</th>
      <th>lib82</th>
      <td>56</td>
      <td>0</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>57</td>
      <td>2</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU25h_alt</th>
      <th>lib82</th>
      <td>109</td>
      <td>2</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>79</td>
      <td>5</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU25h_alt_rpk</th>
      <th>lib82</th>
      <td>34</td>
      <td>0</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>66</td>
      <td>0</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU4_MAP</th>
      <th>lib82</th>
      <td>51</td>
      <td>3</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>58</td>
      <td>1</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU4_MAP_rpk</th>
      <th>lib82</th>
      <td>72</td>
      <td>8</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>100</td>
      <td>2</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU4_alt</th>
      <th>lib82</th>
      <td>72</td>
      <td>5</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>90</td>
      <td>4</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU4_alt_rpk</th>
      <th>lib82</th>
      <td>96</td>
      <td>4</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>82</td>
      <td>4</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU4a_MAP</th>
      <th>lib82</th>
      <td>42</td>
      <td>0</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>37</td>
      <td>0</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU4a_MAP_rpk</th>
      <th>lib82</th>
      <td>33</td>
      <td>1</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>48</td>
      <td>4</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU4a_alt</th>
      <th>lib82</th>
      <td>55</td>
      <td>6</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>49</td>
      <td>3</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU4a_alt_rpk</th>
      <th>lib82</th>
      <td>48</td>
      <td>1</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>65</td>
      <td>5</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU4b_MAP</th>
      <th>lib82</th>
      <td>63</td>
      <td>2</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>87</td>
      <td>3</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU4b_MAP_rpk</th>
      <th>lib82</th>
      <td>64</td>
      <td>3</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>80</td>
      <td>3</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU4b_alt</th>
      <th>lib82</th>
      <td>84</td>
      <td>1</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>80</td>
      <td>3</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU4b_alt_rpk</th>
      <th>lib82</th>
      <td>96</td>
      <td>5</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>77</td>
      <td>4</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU4c_MAP</th>
      <th>lib82</th>
      <td>85</td>
      <td>4</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>62</td>
      <td>2</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU4c_MAP_rpk</th>
      <th>lib82</th>
      <td>94</td>
      <td>2</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>95</td>
      <td>3</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU4c_alt</th>
      <th>lib82</th>
      <td>102</td>
      <td>3</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>93</td>
      <td>6</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU4c_alt_rpk</th>
      <th>lib82</th>
      <td>65</td>
      <td>4</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>84</td>
      <td>4</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU4d_MAP</th>
      <th>lib82</th>
      <td>63</td>
      <td>3</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>63</td>
      <td>1</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU4d_MAP_rpk</th>
      <th>lib82</th>
      <td>47</td>
      <td>0</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>49</td>
      <td>4</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU4d_alt</th>
      <th>lib82</th>
      <td>72</td>
      <td>4</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>75</td>
      <td>4</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU4d_alt_rpk</th>
      <th>lib82</th>
      <td>80</td>
      <td>3</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>81</td>
      <td>4</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU5_MAP</th>
      <th>lib82</th>
      <td>71</td>
      <td>3</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>60</td>
      <td>4</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU5_MAP_rpk</th>
      <th>lib82</th>
      <td>69</td>
      <td>2</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>84</td>
      <td>0</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU5_alt</th>
      <th>lib82</th>
      <td>50</td>
      <td>3</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>58</td>
      <td>2</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU5_alt_rpk</th>
      <th>lib82</th>
      <td>94</td>
      <td>3</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>109</td>
      <td>4</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU5a_MAP</th>
      <th>lib82</th>
      <td>99</td>
      <td>1</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>125</td>
      <td>5</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU5a_MAP_rpk</th>
      <th>lib82</th>
      <td>74</td>
      <td>1</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>79</td>
      <td>2</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU5a_alt</th>
      <th>lib82</th>
      <td>45</td>
      <td>6</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>55</td>
      <td>2</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU5a_alt_rpk</th>
      <th>lib82</th>
      <td>64</td>
      <td>6</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>73</td>
      <td>3</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU5b_MAP</th>
      <th>lib82</th>
      <td>87</td>
      <td>4</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>92</td>
      <td>4</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU5b_MAP_rpk</th>
      <th>lib82</th>
      <td>62</td>
      <td>4</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>44</td>
      <td>1</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU5b_alt</th>
      <th>lib82</th>
      <td>84</td>
      <td>3</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>91</td>
      <td>3</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHKU5b_alt_rpk</th>
      <th>lib82</th>
      <td>99</td>
      <td>8</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>101</td>
      <td>7</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncMERSHKU4_MAP</th>
      <th>lib82</th>
      <td>81</td>
      <td>1</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>87</td>
      <td>4</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncMERSHKU4_MAP_rpk</th>
      <th>lib82</th>
      <td>84</td>
      <td>3</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>101</td>
      <td>2</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncMERSHKU4_alt</th>
      <th>lib82</th>
      <td>94</td>
      <td>5</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>90</td>
      <td>6</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncMERSHKU4_alt_rpk</th>
      <th>lib82</th>
      <td>43</td>
      <td>0</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>62</td>
      <td>4</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncMerbeco1_MAP</th>
      <th>lib82</th>
      <td>75</td>
      <td>2</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>56</td>
      <td>3</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncMerbeco1_MAP_rpk</th>
      <th>lib82</th>
      <td>68</td>
      <td>5</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>73</td>
      <td>7</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncMerbeco1_alt</th>
      <th>lib82</th>
      <td>70</td>
      <td>6</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>102</td>
      <td>2</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncMerbeco1_alt_rpk</th>
      <th>lib82</th>
      <td>49</td>
      <td>3</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>67</td>
      <td>2</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncMerbeco2_MAP</th>
      <th>lib82</th>
      <td>23</td>
      <td>0</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>26</td>
      <td>1</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncMerbeco2_MAP_rpk</th>
      <th>lib82</th>
      <td>58</td>
      <td>5</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>69</td>
      <td>6</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncMerbeco2_alt</th>
      <th>lib82</th>
      <td>60</td>
      <td>0</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>75</td>
      <td>4</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncMerbeco2_alt_rpk</th>
      <th>lib82</th>
      <td>100</td>
      <td>5</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>90</td>
      <td>5</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncMerbeco3_MAP</th>
      <th>lib82</th>
      <td>60</td>
      <td>2</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>77</td>
      <td>1</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncMerbeco3_MAP_rpk</th>
      <th>lib82</th>
      <td>50</td>
      <td>4</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>65</td>
      <td>1</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncMerbeco3_alt</th>
      <th>lib82</th>
      <td>83</td>
      <td>10</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>65</td>
      <td>5</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncMerbeco3_alt_rpk</th>
      <th>lib82</th>
      <td>101</td>
      <td>5</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>108</td>
      <td>9</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncNeoEri_MAP</th>
      <th>lib82</th>
      <td>36</td>
      <td>2</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>37</td>
      <td>2</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncNeoEri_MAP_rpk</th>
      <th>lib82</th>
      <td>9</td>
      <td>0</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>10</td>
      <td>0</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncNeoEri_alt</th>
      <th>lib82</th>
      <td>74</td>
      <td>6</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>94</td>
      <td>3</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncNeoEri_alt_rpk</th>
      <th>lib82</th>
      <td>89</td>
      <td>5</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>93</td>
      <td>4</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncNeoGhEri_MAP</th>
      <th>lib82</th>
      <td>102</td>
      <td>5</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>81</td>
      <td>4</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncNeoGhEri_MAP_rpk</th>
      <th>lib82</th>
      <td>66</td>
      <td>2</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>72</td>
      <td>4</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncNeoGhEri_alt</th>
      <th>lib82</th>
      <td>41</td>
      <td>2</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>57</td>
      <td>3</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncNeoGhEri_alt_rpk</th>
      <th>lib82</th>
      <td>50</td>
      <td>2</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>36</td>
      <td>1</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncNeo_MAP</th>
      <th>lib82</th>
      <td>68</td>
      <td>3</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>80</td>
      <td>2</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncNeo_MAP_rpk</th>
      <th>lib82</th>
      <td>67</td>
      <td>1</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>76</td>
      <td>1</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncNeo_alt</th>
      <th>lib82</th>
      <td>91</td>
      <td>6</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>83</td>
      <td>3</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncNeo_alt_rpk</th>
      <th>lib82</th>
      <td>66</td>
      <td>2</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>74</td>
      <td>3</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">BatCoV_SWI1</th>
      <th>lib82</th>
      <td>60</td>
      <td>3</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>51</td>
      <td>3</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">BatCoV_SWI1_rpk</th>
      <th>lib82</th>
      <td>69</td>
      <td>0</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>77</td>
      <td>7</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">BatCoV_SWI2</th>
      <th>lib82</th>
      <td>86</td>
      <td>7</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>84</td>
      <td>4</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">BatCoV_SWI2_rpk</th>
      <th>lib82</th>
      <td>81</td>
      <td>6</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>89</td>
      <td>6</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">BtCoV_B12080</th>
      <th>lib82</th>
      <td>91</td>
      <td>3</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>81</td>
      <td>3</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">BtCoV_B12080_rpk</th>
      <th>lib82</th>
      <td>47</td>
      <td>3</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>72</td>
      <td>4</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">EjCoV3</th>
      <th>lib82</th>
      <td>103</td>
      <td>5</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>120</td>
      <td>6</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">EjCoV3_rpk</th>
      <th>lib82</th>
      <td>74</td>
      <td>2</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>78</td>
      <td>6</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">EriCoV_Deu174</th>
      <th>lib82</th>
      <td>22</td>
      <td>2</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>33</td>
      <td>3</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">EriCoV_Deu174_rpk</th>
      <th>lib82</th>
      <td>50</td>
      <td>5</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>61</td>
      <td>6</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">EriCoV_GB1</th>
      <th>lib82</th>
      <td>91</td>
      <td>8</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>87</td>
      <td>5</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">EriCoV_GB1_rpk</th>
      <th>lib82</th>
      <td>86</td>
      <td>4</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>85</td>
      <td>3</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">EriCoV_HKU31</th>
      <th>lib82</th>
      <td>32</td>
      <td>2</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>57</td>
      <td>3</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">EriCoV_HKU31_rpk</th>
      <th>lib82</th>
      <td>47</td>
      <td>1</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>41</td>
      <td>7</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">EriCoV_Ita1169881</th>
      <th>lib82</th>
      <td>63</td>
      <td>1</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>74</td>
      <td>2</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">EriCoV_Ita1169881_rpk</th>
      <th>lib82</th>
      <td>24</td>
      <td>1</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>22</td>
      <td>1</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">EriCoV_Ita502651</th>
      <th>lib82</th>
      <td>27</td>
      <td>1</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>40</td>
      <td>1</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">EriCoV_Ita5026511</th>
      <th>lib82</th>
      <td>55</td>
      <td>3</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>79</td>
      <td>1</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">EriCoV_Ita5026511_rpk</th>
      <th>lib82</th>
      <td>7</td>
      <td>1</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>7</td>
      <td>1</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">EriCoV_Ita5026512</th>
      <th>lib82</th>
      <td>58</td>
      <td>4</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>64</td>
      <td>4</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">EriCoV_Ita5026512_rpk</th>
      <th>lib82</th>
      <td>20</td>
      <td>2</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>21</td>
      <td>2</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">EriCoV_Ita5026515</th>
      <th>lib82</th>
      <td>65</td>
      <td>1</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>91</td>
      <td>0</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">EriCoV_Ita5026515_rpk</th>
      <th>lib82</th>
      <td>93</td>
      <td>4</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>79</td>
      <td>1</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">EriCoV_Ita5026517</th>
      <th>lib82</th>
      <td>82</td>
      <td>3</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>80</td>
      <td>4</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">EriCoV_Ita5026517_rpk</th>
      <th>lib82</th>
      <td>102</td>
      <td>2</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>84</td>
      <td>1</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">EriCoV_Ita5026519</th>
      <th>lib82</th>
      <td>60</td>
      <td>7</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>87</td>
      <td>2</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">EriCoV_Ita5026519_rpk</th>
      <th>lib82</th>
      <td>75</td>
      <td>1</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>91</td>
      <td>1</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">EriCoV_Ita502651_rpk</th>
      <th>lib82</th>
      <td>17</td>
      <td>2</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>4</td>
      <td>1</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU25_NL140462</th>
      <th>lib82</th>
      <td>100</td>
      <td>4</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>122</td>
      <td>7</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU25_NL140462_rpk</th>
      <th>lib82</th>
      <td>76</td>
      <td>3</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>58</td>
      <td>4</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU25_YD131305</th>
      <th>lib82</th>
      <td>86</td>
      <td>4</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>78</td>
      <td>4</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU25_YD131305_rpk</th>
      <th>lib82</th>
      <td>75</td>
      <td>5</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>93</td>
      <td>4</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU4</th>
      <th>lib82</th>
      <td>91</td>
      <td>3</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>104</td>
      <td>4</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU4_rpk</th>
      <th>lib82</th>
      <td>63</td>
      <td>1</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>70</td>
      <td>2</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU4r_161028</th>
      <th>lib82</th>
      <td>34</td>
      <td>4</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>45</td>
      <td>1</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU4r_161028_rpk</th>
      <th>lib82</th>
      <td>90</td>
      <td>2</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>95</td>
      <td>6</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU4r_162275</th>
      <th>lib82</th>
      <td>54</td>
      <td>3</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>76</td>
      <td>2</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU4r_162275_rpk</th>
      <th>lib82</th>
      <td>75</td>
      <td>4</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>84</td>
      <td>4</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU4r_BtCoV_1332005</th>
      <th>lib82</th>
      <td>64</td>
      <td>1</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>70</td>
      <td>3</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU4r_BtCoV_1332005_rpk</th>
      <th>lib82</th>
      <td>116</td>
      <td>4</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>92</td>
      <td>4</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU4r_CZ01</th>
      <th>lib82</th>
      <td>41</td>
      <td>4</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>68</td>
      <td>4</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU4r_CZ01_rpk</th>
      <th>lib82</th>
      <td>79</td>
      <td>4</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>61</td>
      <td>2</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU4r_GZ160421</th>
      <th>lib82</th>
      <td>86</td>
      <td>2</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>60</td>
      <td>4</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU4r_GZ160421_rpk</th>
      <th>lib82</th>
      <td>92</td>
      <td>5</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>80</td>
      <td>5</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU4r_GZ1862</th>
      <th>lib82</th>
      <td>72</td>
      <td>2</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>73</td>
      <td>3</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU4r_GZ1862_rpk</th>
      <th>lib82</th>
      <td>57</td>
      <td>0</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>62</td>
      <td>5</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU4r_JPDB144</th>
      <th>lib82</th>
      <td>84</td>
      <td>3</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>76</td>
      <td>2</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU4r_JPDB144_rpk</th>
      <th>lib82</th>
      <td>85</td>
      <td>2</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>113</td>
      <td>4</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU4r_panP251T</th>
      <th>lib82</th>
      <td>65</td>
      <td>1</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>67</td>
      <td>3</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU4r_panP251T_rpk</th>
      <th>lib82</th>
      <td>85</td>
      <td>6</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>86</td>
      <td>7</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU5</th>
      <th>lib82</th>
      <td>75</td>
      <td>4</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>81</td>
      <td>7</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU5_rpk</th>
      <th>lib82</th>
      <td>88</td>
      <td>2</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>79</td>
      <td>4</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU5r_17S</th>
      <th>lib82</th>
      <td>79</td>
      <td>0</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>77</td>
      <td>5</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU5r_17S_rpk</th>
      <th>lib82</th>
      <td>78</td>
      <td>3</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>93</td>
      <td>5</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU5r_19S</th>
      <th>lib82</th>
      <td>98</td>
      <td>6</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>81</td>
      <td>6</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU5r_19S_rpk</th>
      <th>lib82</th>
      <td>85</td>
      <td>3</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>84</td>
      <td>6</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU5r_2</th>
      <th>lib82</th>
      <td>86</td>
      <td>2</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>84</td>
      <td>2</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU5r_20S</th>
      <th>lib82</th>
      <td>84</td>
      <td>5</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>75</td>
      <td>0</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU5r_20S_rpk</th>
      <th>lib82</th>
      <td>71</td>
      <td>4</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>93</td>
      <td>9</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU5r_21S</th>
      <th>lib82</th>
      <td>45</td>
      <td>3</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>48</td>
      <td>6</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU5r_21S_rpk</th>
      <th>lib82</th>
      <td>96</td>
      <td>4</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>110</td>
      <td>5</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU5r_24S</th>
      <th>lib82</th>
      <td>61</td>
      <td>5</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>37</td>
      <td>4</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU5r_24S_rpk</th>
      <th>lib82</th>
      <td>79</td>
      <td>4</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>90</td>
      <td>1</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU5r_27S</th>
      <th>lib82</th>
      <td>56</td>
      <td>5</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>74</td>
      <td>1</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU5r_27S_rpk</th>
      <th>lib82</th>
      <td>67</td>
      <td>2</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>53</td>
      <td>2</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU5r_2_rpk</th>
      <th>lib82</th>
      <td>101</td>
      <td>5</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>107</td>
      <td>7</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU5r_30S</th>
      <th>lib82</th>
      <td>81</td>
      <td>2</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>94</td>
      <td>2</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU5r_30S_rpk</th>
      <th>lib82</th>
      <td>74</td>
      <td>5</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>62</td>
      <td>2</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU5r_31S</th>
      <th>lib82</th>
      <td>51</td>
      <td>1</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>44</td>
      <td>5</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU5r_31S_rpk</th>
      <th>lib82</th>
      <td>93</td>
      <td>4</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>90</td>
      <td>3</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU5r_32S</th>
      <th>lib82</th>
      <td>68</td>
      <td>3</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>89</td>
      <td>1</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU5r_32S_rpk</th>
      <th>lib82</th>
      <td>62</td>
      <td>6</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>69</td>
      <td>2</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU5r_33S</th>
      <th>lib82</th>
      <td>97</td>
      <td>3</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>82</td>
      <td>8</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU5r_33S_rpk</th>
      <th>lib82</th>
      <td>61</td>
      <td>5</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>71</td>
      <td>4</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU5r_A4342005</th>
      <th>lib82</th>
      <td>96</td>
      <td>8</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>82</td>
      <td>5</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU5r_A4342005_rpk</th>
      <th>lib82</th>
      <td>98</td>
      <td>3</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>107</td>
      <td>8</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU5r_BY140535</th>
      <th>lib82</th>
      <td>62</td>
      <td>3</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>68</td>
      <td>2</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU5r_BY140535_rpk</th>
      <th>lib82</th>
      <td>88</td>
      <td>0</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>64</td>
      <td>2</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU5r_BY140562</th>
      <th>lib82</th>
      <td>84</td>
      <td>6</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>108</td>
      <td>1</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU5r_BY140562_rpk</th>
      <th>lib82</th>
      <td>85</td>
      <td>5</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>91</td>
      <td>3</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU5r_BY140568</th>
      <th>lib82</th>
      <td>66</td>
      <td>4</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>90</td>
      <td>4</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU5r_BY140568_rpk</th>
      <th>lib82</th>
      <td>69</td>
      <td>0</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>72</td>
      <td>3</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU5r_Wenzhou1</th>
      <th>lib82</th>
      <td>49</td>
      <td>0</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>66</td>
      <td>2</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU5r_Wenzhou1_rpk</th>
      <th>lib82</th>
      <td>114</td>
      <td>2</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>106</td>
      <td>3</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU5r_YD13403</th>
      <th>lib82</th>
      <td>64</td>
      <td>3</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>75</td>
      <td>0</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU5r_YD13403_rpk</th>
      <th>lib82</th>
      <td>110</td>
      <td>4</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>91</td>
      <td>6</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HsCoV_Ita20664540</th>
      <th>lib82</th>
      <td>96</td>
      <td>2</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>88</td>
      <td>0</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HsCoV_Ita20664540_rpk</th>
      <th>lib82</th>
      <td>90</td>
      <td>2</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>68</td>
      <td>3</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">KW2EF93Gha2010</th>
      <th>lib82</th>
      <td>58</td>
      <td>0</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>55</td>
      <td>2</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">KW2EF93Gha2010_rpk</th>
      <th>lib82</th>
      <td>58</td>
      <td>0</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>54</td>
      <td>3</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">MjHKU4r1</th>
      <th>lib82</th>
      <td>76</td>
      <td>5</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>84</td>
      <td>4</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">MjHKU4r1_rpk</th>
      <th>lib82</th>
      <td>90</td>
      <td>5</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>98</td>
      <td>5</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">NL13845</th>
      <th>lib82</th>
      <td>83</td>
      <td>4</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>78</td>
      <td>3</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">NL13845_rpk</th>
      <th>lib82</th>
      <td>58</td>
      <td>2</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>52</td>
      <td>1</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">NL13892</th>
      <th>lib82</th>
      <td>85</td>
      <td>6</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>98</td>
      <td>3</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">NL13892_rpk</th>
      <th>lib82</th>
      <td>78</td>
      <td>3</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>95</td>
      <td>6</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">NL140422</th>
      <th>lib82</th>
      <td>79</td>
      <td>0</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>83</td>
      <td>3</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">NL140422_rpk</th>
      <th>lib82</th>
      <td>76</td>
      <td>6</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>92</td>
      <td>2</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">NL140455</th>
      <th>lib82</th>
      <td>105</td>
      <td>4</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>113</td>
      <td>2</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">NL140455_rpk</th>
      <th>lib82</th>
      <td>84</td>
      <td>4</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>71</td>
      <td>1</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">NeoCoV</th>
      <th>lib82</th>
      <td>81</td>
      <td>6</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>74</td>
      <td>3</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">NeoCoV_rpk</th>
      <th>lib82</th>
      <td>63</td>
      <td>2</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>98</td>
      <td>3</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">PDF2180</th>
      <th>lib82</th>
      <td>69</td>
      <td>4</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>86</td>
      <td>9</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">PDF2180_rpk</th>
      <th>lib82</th>
      <td>104</td>
      <td>6</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>101</td>
      <td>0</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">PaGB01</th>
      <th>lib82</th>
      <td>95</td>
      <td>6</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>131</td>
      <td>4</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">PaGB01_rpk</th>
      <th>lib82</th>
      <td>74</td>
      <td>3</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>67</td>
      <td>2</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">PkCoV_Ita20664563</th>
      <th>lib82</th>
      <td>82</td>
      <td>3</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>85</td>
      <td>3</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">PkCoV_Ita20664563_rpk</th>
      <th>lib82</th>
      <td>78</td>
      <td>8</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>89</td>
      <td>8</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">Russia_MOW1522</th>
      <th>lib82</th>
      <td>56</td>
      <td>2</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>75</td>
      <td>3</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">Russia_MOW1522_rpk</th>
      <th>lib82</th>
      <td>90</td>
      <td>2</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>84</td>
      <td>4</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">SC2013</th>
      <th>lib82</th>
      <td>73</td>
      <td>2</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>63</td>
      <td>1</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">SC2013_rpk</th>
      <th>lib82</th>
      <td>42</td>
      <td>1</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>46</td>
      <td>1</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">VsCoV1</th>
      <th>lib82</th>
      <td>58</td>
      <td>2</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>54</td>
      <td>4</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">VsCoV1_rpk</th>
      <th>lib82</th>
      <td>40</td>
      <td>3</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>49</td>
      <td>0</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">VsCoVa7</th>
      <th>lib82</th>
      <td>91</td>
      <td>2</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>125</td>
      <td>6</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">VsCoVa7_rpk</th>
      <th>lib82</th>
      <td>68</td>
      <td>2</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>90</td>
      <td>0</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">VsCoVkj15</th>
      <th>lib82</th>
      <td>70</td>
      <td>6</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>59</td>
      <td>3</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">VsCoVkj15_rpk</th>
      <th>lib82</th>
      <td>71</td>
      <td>3</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>71</td>
      <td>0</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">cMERS_ALA49341</th>
      <th>lib82</th>
      <td>101</td>
      <td>6</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>81</td>
      <td>4</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">cMERS_ALA49341_rpk</th>
      <th>lib82</th>
      <td>95</td>
      <td>6</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>91</td>
      <td>5</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">cMERS_ALA49374</th>
      <th>lib82</th>
      <td>37</td>
      <td>1</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>33</td>
      <td>1</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">cMERS_ALA49374_rpk</th>
      <th>lib82</th>
      <td>44</td>
      <td>3</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>40</td>
      <td>1</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">cMERS_AMO03401</th>
      <th>lib82</th>
      <td>67</td>
      <td>0</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>66</td>
      <td>3</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">cMERS_AMO03401_rpk</th>
      <th>lib82</th>
      <td>97</td>
      <td>4</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>99</td>
      <td>3</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">cMERS_ANI69835</th>
      <th>lib82</th>
      <td>74</td>
      <td>0</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>57</td>
      <td>5</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">cMERS_ANI69835_rpk</th>
      <th>lib82</th>
      <td>87</td>
      <td>3</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>94</td>
      <td>2</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">cMERS_ANI69889</th>
      <th>lib82</th>
      <td>86</td>
      <td>1</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>88</td>
      <td>2</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">cMERS_ANI69889_rpk</th>
      <th>lib82</th>
      <td>50</td>
      <td>3</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>56</td>
      <td>1</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">hMERS_AKM76239</th>
      <th>lib82</th>
      <td>97</td>
      <td>9</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>107</td>
      <td>7</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">hMERS_AKM76239_rpk</th>
      <th>lib82</th>
      <td>50</td>
      <td>3</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>56</td>
      <td>3</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">hMERS_AMQ48993</th>
      <th>lib82</th>
      <td>50</td>
      <td>1</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>48</td>
      <td>1</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">hMERS_AMQ48993_rpk</th>
      <th>lib82</th>
      <td>99</td>
      <td>2</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>99</td>
      <td>4</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">hMERS_EMC2012</th>
      <th>lib82</th>
      <td>48</td>
      <td>0</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>45</td>
      <td>0</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">hMERS_EMC2012_rpk</th>
      <th>lib82</th>
      <td>65</td>
      <td>3</td>
    </tr>
    <tr>
      <th>lib83</th>
      <td>71</td>
      <td>3</td>
    </tr>
  </tbody>
</table>



    
![png](process_ccs_panmerbeco_files/process_ccs_panmerbeco_90_1.png)
    


Are there any barcodes in the same library that are shared across targets? If so, we need to get rid of those as they will be confounded in barcode parsing:


```python
dup_barcodes = (
    consensus
    .groupby(['library', 'barcode'])
    .size()
    .rename('duplicate_count')
    .reset_index()
    .query('duplicate_count > 1')
    )

print('Here are duplicated barcodes:')
display(HTML(dup_barcodes.head().to_html(index=False)))

print(f"\nRemoving the {len(dup_barcodes)} duplicated barcodes."
      f"Started with {len(consensus)} barcodes:")
consensus = (
    consensus
    .merge(dup_barcodes, on=['library', 'barcode'], how='outer')
    .query('duplicate_count.isnull()', engine='python')
    )
print(f"After removing duplicates, there are {len(consensus)} barcodes.")
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
      <td>lib82</td>
      <td>AAAAAAAGGACTTAAA</td>
      <td>2</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AAAAAAATGATGCCGT</td>
      <td>2</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AAAAAACTCAACTATT</td>
      <td>2</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AAAAAGTTCCCCGGCA</td>
      <td>2</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AAAAATCCGACAACAT</td>
      <td>2</td>
    </tr>
  </tbody>
</table>


    
    Removing the 551 duplicated barcodes.Started with 36484 barcodes:
    After removing duplicates, there are 35378 barcodes.


Below we write the retained consensus sequences to a CSV file that links the nucleotide mutations to the barcodes:


```python
print(f"Writing nucleotide variants to {config['nt_variant_table_file' + '_' + background]}")
      
(consensus
 [['target', 'library', 'barcode', 'gene_mutations', 'substitutions','variant_call_support','has_substitutions']]
 .to_csv(config['nt_variant_table_file' + '_' + background], index=False)
 )
      
print('Here are the first few lines of this file:')
display(HTML(
    pd.read_csv(config['nt_variant_table_file' + '_' + background], na_filter=None)
    .head()
    .to_html(index=False)
    ))
```

    Writing nucleotide variants to results/variants/nucleotide_variant_table_panmerbeco.csv
    Here are the first few lines of this file:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>target</th>
      <th>library</th>
      <th>barcode</th>
      <th>gene_mutations</th>
      <th>substitutions</th>
      <th>variant_call_support</th>
      <th>has_substitutions</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>HKU5r_BY140535_rpk</td>
      <td>lib82</td>
      <td>AAAAAAAAAAATGTTT</td>
      <td></td>
      <td></td>
      <td>1</td>
      <td>False</td>
    </tr>
    <tr>
      <td>EriCoV_Ita5026517_rpk</td>
      <td>lib82</td>
      <td>AAAAAAAAAAGTACCG</td>
      <td></td>
      <td></td>
      <td>6</td>
      <td>False</td>
    </tr>
    <tr>
      <td>PkCoV_Ita20664563</td>
      <td>lib82</td>
      <td>AAAAAAAAGTCTTGAG</td>
      <td></td>
      <td></td>
      <td>2</td>
      <td>False</td>
    </tr>
    <tr>
      <td>EriCoV_Ita5026519</td>
      <td>lib82</td>
      <td>AAAAAAACAATGGACC</td>
      <td></td>
      <td></td>
      <td>4</td>
      <td>False</td>
    </tr>
    <tr>
      <td>cMERS_ANI69835</td>
      <td>lib82</td>
      <td>AAAAAAACACATGAAA</td>
      <td></td>
      <td></td>
      <td>3</td>
      <td>False</td>
    </tr>
  </tbody>
</table>


What happened to the barcodes that we "dropped" because we could not construct a reliable consensus?
The `dropped` data frame from [alignparse.consensus.simple_mutconsensus](https://jbloomlab.github.io/alignparse/alignparse.consensus.html?highlight=simple_mutconsensus#alignparse.consensus.simple_mutconsensus) has this information:


```python
display(HTML(dropped.head().to_html(index=False)))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>library</th>
      <th>barcode</th>
      <th>target</th>
      <th>drop_reason</th>
      <th>nseqs</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>lib82</td>
      <td>AAAAGTTCGATAAACC</td>
      <td>AncMerbeco1_alt_rpk</td>
      <td>minor subs too frequent</td>
      <td>7</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AAACAAATATTAGACG</td>
      <td>Russia_MOW1522_rpk</td>
      <td>minor subs too frequent</td>
      <td>7</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AAAGCACCTAACGATG</td>
      <td>AncHKU25_alt_rpk</td>
      <td>subs diff too large</td>
      <td>3</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AACATCACTCTCAGAA</td>
      <td>AncHKU5_alt</td>
      <td>minor subs too frequent</td>
      <td>12</td>
    </tr>
    <tr>
      <td>lib82</td>
      <td>AACGGCGATGTGCAGG</td>
      <td>hMERS_AKM76239_rpk</td>
      <td>subs diff too large</td>
      <td>8</td>
    </tr>
  </tbody>
</table>


Summarize the information in this data frame on dropped barcodes with the plot below.
This plot shows several things.
First, we see that the total number of barcodes dropped is modest (just a few thousand per library) relative to the total number of barcodes per library (seen above to be on the order of hundreds of thousands).
Second, the main reason that barcodes are dropped is that there are CCSs within the same barcode with suspiciously large numbers of mutations relative to the consensus---which we use as a filter to discard the entire barcode as it could indicate strand exchange or some other issue.
In any case, the modest number of dropped barcodes indicates that there probably isn't much of a need to worry: 


```python
max_nseqs = 8  # plot together all barcodes with >= this many sequences

_ = (
 ggplot(
    dropped.assign(nseqs=lambda x: numpy.clip(x['nseqs'], None, max_nseqs)),
    aes('nseqs')) + 
 geom_bar() + 
 scale_x_continuous(limits=(1, None)) +
 xlab('number of sequences for barcode') +
 ylab('number of barcodes') +
 facet_grid('library ~ drop_reason') +
 theme(figure_size=(10, 1.5 * nlibs),
       panel_grid_major_x=element_blank(),
       )
 ).draw()
```


    
![png](process_ccs_panmerbeco_files/process_ccs_panmerbeco_98_0.png)
    

