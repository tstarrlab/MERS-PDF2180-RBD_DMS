"""``snakemake`` file that runs entire analysis."""

# Imports ---------------------------------------------------------------------
import glob
import itertools
import os.path
import os
import textwrap
import urllib.request

import Bio.SeqIO

import dms_variants.codonvarianttable
import dms_variants.illuminabarcodeparser

import pandas as pd

# Configuration  --------------------------------------------------------------
configfile: 'config.yaml'

# run "quick" rules locally:
localrules: make_dag,
            make_summary,
            save_pinned_env,
            get_ccs

# Functions -------------------------------------------------------------------
def nb_markdown(nb):
    """Return path to Markdown results of notebook `nb`."""
    return os.path.join(config['summary_dir'],
                        os.path.basename(os.path.splitext(nb)[0]) + '.md')

# Global variables extracted from config --------------------------------------
pacbio_runs = (pd.read_csv(config['pacbio_runs'], dtype = str)
               .assign(pacbioRun=lambda x: x['library'] + '_' + x['bg'] + '_' + x['run'])
               )
assert len(pacbio_runs['pacbioRun'].unique()) == len(pacbio_runs['pacbioRun'])

# Information on samples and barcode runs -------------------------------------
barcode_runs = pd.read_csv(config['barcode_runs'])

# Rules -----------------------------------------------------------------------

# making this summary is the target rule (in place of `all`) since it
# is first rule listed.
rule make_summary:
    """Create Markdown summary of analysis."""
    input:
        dag=os.path.join(config['summary_dir'], 'dag.svg'),
        env='environment_pinned.yml',
        process_ccs_MERS=nb_markdown('process_ccs_MERS.ipynb'),
        process_ccs_MERS_rpk=nb_markdown('process_ccs_MERS_rpk.ipynb'),
        process_ccs_PDF2180=nb_markdown('process_ccs_PDF2180.ipynb'),
        barcode_variant_table_panmerbeco=config['nt_variant_table_file_panmerbeco'],
       	process_ccs_panmerbeco=nb_markdown('process_ccs_panmerbeco.ipynb'),
        barcode_variant_table_MERS=config['codon_variant_table_file_MERS'],
        barcode_variant_table_MERS_rpk=config['codon_variant_table_file_MERS_rpk'],
        barcode_variant_table_PDF2180=config['codon_variant_table_file_PDF2180'],
        barcode_variant_table_merged=config['codon_variant_table_file_pools'],
        merge_tables='results/summary/merge_pools.md',
        variant_counts_file=config['variant_counts_file'],
        count_variants=nb_markdown('count_variants.ipynb'),
        
        fit_titrations_hDPP4='results/summary/compute_Kd_hDPP4.md',
        hDPP4_Kds_file=config['Titeseq_Kds_file_hDPP4'],
        fit_titrations_ApACE2='results/summary/compute_Kd_ApACE2.md',
        ApACE2_Kds_file=config['Titeseq_Kds_file_ApACE2'],

        fit_AUC_serum='results/summary/compute_serum_AUC.md',
        serum_AUC_file=config['sera_delta_AUC_file'],

        fit_EC50_mAb='results/summary/compute_mAb_EC50.md',
        mAb_EC50_file=config['mAb_EC50_file'],
        
        calculate_expression='results/summary/compute_expression_meanF.md',
        variant_expression_file=config['expression_sortseq_file'],
        
        collapse_scores='results/summary/collapse_scores.md',
        mut_phenos_file_MERS_rpk=config['final_variant_scores_mut_file_MERS_rpk'],
        mut_phenos_file_PDF2180=config['final_variant_scores_mut_file_PDF2180'],
        mut_phenos_file_panmerbeco=config['final_variant_scores_file_panmerbeco'],
        
        heatmap_viz_delta=os.path.join(config['visualization_dir'], "heatmaps_delta.html"),
        heatmap_viz_absolute=os.path.join(config['visualization_dir'], "heatmaps_absolute.html"),
    output:
        summary = os.path.join(config['summary_dir'], 'summary.md')
    run:
        def path(f):
            """Get path relative to `summary_dir`."""
            return os.path.relpath(f, config['summary_dir'])
        with open(output.summary, 'w') as f:
            f.write(textwrap.dedent(f"""
            # Summary

            Analysis run by [Snakefile]({path(workflow.snakefile)})
            using [this config file]({path(workflow.configfiles[0])}).
            See the [README in the top directory]({path('README.md')})
            for details.

            Here is the DAG of the computational workflow:
            ![{path(input.dag)}]({path(input.dag)})

            Here is the Markdown output of each Jupyter notebook in the
            workflow:
            
            1. Process PacBio CCSs for each background: [MERS-CoV]({path(input.process_ccs_MERS)}), [MERS-CoV rpk]({path(input.process_ccs_MERS_rpk)}), [PDF2180]({path(input.process_ccs_PDF2180)}), and [Pan-Merbeco]({path(input.process_ccs_panmerbeco)}). Creates barcode-variant lookup tables for each background: [MERS-CoV]({path(input.barcode_variant_table_MERS)}), [MERS-CoV rpk]({path(input.barcode_variant_table_MERS_rpk)}), [PDF2180]({path(input.barcode_variant_table_PDF2180)}),[Pan-Merbeco]({path(input.barcode_variant_table_panmerbeco)}) .

            2. Merge barcode-variant sublibraries into pooled libraries used for experiments, as done [here]({path(input.merge_tables)}).
            

            3. [Count variants by barcode]({path(input.count_variants)}).
               Creates a [variant counts file]({path(input.variant_counts_file)})
               giving counts of each barcoded variant in each condition.

            4. Fit titration curves for RBD binding to [hDPP4]({path(input.fit_titrations_hDPP4)}) to calculate per-barcode K<sub>D</sub>, recorded in these files for [hDPP4]({path(input.hDPP4_Kds_file)}).
            
            5. Fit titration curves for RBD binding to [ApACE2]({path(input.fit_titrations_ApACE2)}) to calculate per-barcode K<sub>D</sub>, recorded in these files for [ApACE2]({path(input.ApACE2_Kds_file)}).
            
            6. Fit mAb binding curves for RBD binding to  [mAbs]({path(input.fit_EC50_mAb)}) to calculate per-barcode EC<sub>50</sub>, recorded in [this file]({path(input.mAb_EC50_file)}).
            
            7. Fit serum binding curves for RBD binding to  [sera]({path(input.fit_AUC_serum)}) to calculate per-barcode AUC, recorded in [this file]({path(input.serum_AUC_file)}).
            
            8. [Analyze Sort-seq]({path(input.calculate_expression)}) to calculate per-barcode RBD expression, recorded in [this file]({path(input.variant_expression_file)}).
            
            9. [Derive final genotype-level phenotypes from replicate barcoded sequences]({path(input.collapse_scores)}).
               Generates final phenotypes, recorded in these files for  [MERS_rpk]({path(input.mut_phenos_file_MERS_rpk)}), [PDF2180]({path(input.mut_phenos_file_PDF2180)}), and [pan-merbeco]({path(input.mut_phenos_file_panmerbeco)}).
            
            
            10. Make interactive data visualizations, available [here](https://jbloomlab.github.io/MERS-PDF2180-RBD_DMS/)

            """
            ).strip())

rule make_dag:
    # error message, but works: https://github.com/sequana/sequana/issues/115
    input:
        workflow.snakefile
    output:
        os.path.join(config['summary_dir'], 'dag.svg')
    shell:
        "snakemake --forceall --dag | dot -Tsvg > {output}"

rule save_pinned_env:
    input:
        dag=os.path.join(config['summary_dir'], 'dag.svg'),
    log:
    	"environment_pinned.yml"
    shell:
        """
        conda env export > {log}
        """


rule interactive_heatmap_absolute:
    """ Make the interactive heatmaps for absolute expression and binding phenotypes.
    """
    input: 
        scores=config['final_variant_scores_mut_file_MERS_rpk']
    params:
        annotations=config['RBD_sites']
    output:
        html=os.path.join(config['visualization_dir'], "heatmaps_absolute.html")
    notebook: "RBD-Heatmaps-Interactive-Visualization_absolute.ipynb"

rule interactive_heatmap_delta:
    """ Make the interactive heatmaps for delta expression and binding phenotypes.
    """
    input: 
        scores=config['final_variant_scores_mut_file_MERS_rpk']
    params:
        annotations=config['RBD_sites']
    output:
        html=os.path.join(config['visualization_dir'], "heatmaps_delta.html")
    notebook: "RBD-Heatmaps-Interactive-Visualization_delta.ipynb"

rule collapse_scores:
    input:
        config['Titeseq_Kds_file_hDPP4'],
        config['Titeseq_Kds_file_ApACE2'],
        config['expression_sortseq_file'],
        config['sera_delta_AUC_file'],
        config['mAb_EC50_file'],
    output:
        config['final_variant_scores_mut_file_MERS_rpk'],
        config['final_variant_scores_mut_file_PDF2180'],
        config['final_variant_scores_file_panmerbeco'],
        md='results/summary/collapse_scores.md',
        md_files=directory('results/summary/collapse_scores_files')
    envmodules:
        'R/4.1.3'
    params:
        nb='collapse_scores.Rmd',
        md='collapse_scores.md',
        md_files='collapse_scores_files'
    shell:
        """
        R -e \"rmarkdown::render(input=\'{params.nb}\')\";
        mv {params.md} {output.md};
        mv {params.md_files} {output.md_files}
        """

rule fit_EC50_mAbs:
    input:
        config['codon_variant_table_file_pools'],
        config['variant_counts_file']
    output:
        config['mAb_EC50_file'],
        md='results/summary/compute_mAb_EC50.md',
        md_files=directory('results/summary/compute_mAb_EC50_files')
    envmodules:
        'R/4.1.3'
    params:
        nb='compute_mAb_EC50.Rmd',
        md='compute_mAb_EC50.md',
        md_files='compute_mAb_EC50_files'
    shell:
        """
        R -e \"rmarkdown::render(input=\'{params.nb}\')\";
        mv {params.md} {output.md};
        mv {params.md_files} {output.md_files}
        """

rule fit_AUC_sera:
    input:
        config['codon_variant_table_file_pools'],
        config['variant_counts_file']
    output:
        config['sera_delta_AUC_file'],
        md='results/summary/compute_serum_AUC.md',
        md_files=directory('results/summary/compute_serum_AUC_files')
    envmodules:
        'R/4.1.3'
    params:
        nb='compute_serum_AUC.Rmd',
        md='compute_serum_AUC.md',
        md_files='compute_serum_AUC_files'
    shell:
        """
        R -e \"rmarkdown::render(input=\'{params.nb}\')\";
        mv {params.md} {output.md};
        mv {params.md_files} {output.md_files}
        """

rule fit_titrations_ApACE2:
    input:
        config['codon_variant_table_file_pools'],
        config['variant_counts_file']
    output:
        config['Titeseq_Kds_file_ApACE2'],
        md='results/summary/compute_Kd_ApACE2.md',
        md_files=directory('results/summary/compute_Kd_ApACE2_files')
    envmodules:
        'R/4.1.3'
    params:
        nb='compute_Kd_ApACE2.Rmd',
        md='compute_Kd_ApACE2.md',
        md_files='compute_Kd_ApACE2_files'
    shell:
        """
        R -e \"rmarkdown::render(input=\'{params.nb}\')\";
        mv {params.md} {output.md};
        mv {params.md_files} {output.md_files}
        """

rule fit_titrations_hDPP4:
    input:
        config['codon_variant_table_file_pools'],
        config['variant_counts_file']
    output:
        config['Titeseq_Kds_file_hDPP4'],
        md='results/summary/compute_Kd_hDPP4.md',
        md_files=directory('results/summary/compute_Kd_hDPP4_files')
    envmodules:
        'R/4.1.3'
    params:
        nb='compute_Kd_hDPP4.Rmd',
        md='compute_Kd_hDPP4.md',
        md_files='compute_Kd_hDPP4_files'
    shell:
        """
        R -e \"rmarkdown::render(input=\'{params.nb}\')\";
        mv {params.md} {output.md};
        mv {params.md_files} {output.md_files}
        """

rule calculate_expression:
    input:
        config['codon_variant_table_file_pools'],
        config['variant_counts_file']
    output:
        config['expression_sortseq_file'],
        md='results/summary/compute_expression_meanF.md',
        md_files=directory('results/summary/compute_expression_meanF_files')
    envmodules:
        'R/4.1.3'
    params:
        nb='compute_expression_meanF.Rmd',
        md='compute_expression_meanF.md',
        md_files='compute_expression_meanF_files'
    shell:
        """
        R -e \"rmarkdown::render(input=\'{params.nb}\')\";
        mv {params.md} {output.md};
        mv {params.md_files} {output.md_files}
        """




rule count_variants:
    """Count codon variants from Illumina barcode runs."""
    input:
        config['codon_variant_table_file_pools'],
        config['barcode_runs']
    output:
        config['variant_counts_file'],
        nb_markdown=nb_markdown('count_variants.ipynb')
    params:
        nb='count_variants.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"




rule merge_libs_to_pools:
    input:
        config['codon_variant_table_file_MERS'],
        config['codon_variant_table_file_PDF2180'],
        'merge_pools.Rmd'
    output:
        config['codon_variant_table_file_pools'],
        md='results/summary/merge_pools.md'
    envmodules:
        'R/4.1.3'
    params:
        nb='merge_pools.Rmd',
        md='merge_pools.md'
    shell:
        """
        R -e \"rmarkdown::render(input=\'{params.nb}\')\";
        mv {params.md} {output.md};
        """
rule process_ccs_MERS_rpk:
    """Process the PacBio CCSs for MERS_rpk background and build variant table."""
    input:
        expand(os.path.join(config['ccs_dir'], "{pacbioRun}_ccs.fastq.gz"),
               pacbioRun=pacbio_runs['pacbioRun']),
    output:
        config['processed_ccs_file_MERS_rpk'],
    	config['codon_variant_table_file_MERS_rpk'],
        nb_markdown=nb_markdown('process_ccs_MERS_rpk.ipynb')
    params:
        nb='process_ccs_MERS_rpk.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"
        
rule process_ccs_MERS:
    """Process the PacBio CCSs for MERS background and build variant table."""
    input:
        expand(os.path.join(config['ccs_dir'], "{pacbioRun}_ccs.fastq.gz"),
               pacbioRun=pacbio_runs['pacbioRun']),
    output:
        config['processed_ccs_file_MERS'],
    	config['codon_variant_table_file_MERS'],
        nb_markdown=nb_markdown('process_ccs_MERS.ipynb')
    params:
        nb='process_ccs_MERS.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"

rule process_ccs_panmerbeco:
    """Process the PacBio CCSs for panmerbeco background and build variant table."""
    input:
        expand(os.path.join(config['ccs_dir'], "{pacbioRun}_ccs.fastq.gz"),
               pacbioRun=pacbio_runs['pacbioRun']),
    output:
        config['nt_variant_table_file_panmerbeco'],
        nb_markdown=nb_markdown('process_ccs_panmerbeco.ipynb')
    params:
        nb='process_ccs_panmerbeco.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"

rule process_ccs_PDF2180:
    """Process the PacBio CCSs for PDF2180 background and build variant table."""
    input:
        expand(os.path.join(config['ccs_dir'], "{pacbioRun}_ccs.fastq.gz"),
               pacbioRun=pacbio_runs['pacbioRun']),
    output:
        config['processed_ccs_file_PDF2180'],
    	config['codon_variant_table_file_PDF2180'],
        nb_markdown=nb_markdown('process_ccs_PDF2180.ipynb')
    params:
        nb='process_ccs_PDF2180.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"

if config['seqdata_source'] == 'local':

    rule get_ccs:
        """Symbolically link CCS files."""
        input:
            ccs_fastq=lambda wildcards: (pacbio_runs
                                        .set_index('pacbioRun')
                                        .at[wildcards.pacbioRun, 'ccs']
                                        )
        output:
            ccs_fastq=os.path.join(config['ccs_dir'], "{pacbioRun}_ccs.fastq.gz")
        run:
            os.symlink(input.ccs_fastq, output.ccs_fastq)

elif config['seqdata_source'] == 'SRA':
    raise RuntimeError('getting sequence data from SRA not yet implemented')

else:
    raise ValueError(f"invalid `seqdata_source` {config['seqdata_source']}")
