#!/usr/bin/env nextflow


nextflow.enable.dsl=2

project_dir = projectDir



process VIS_MARKERS_DM1 {

    publishDir "${params.output}/DM1", mode: 'copy', overwrite: true

    label 'full_resources'

    input:
    path mae_path
    path r_config

    output:
    path 'visualize_marker_genes_DM1.html'
    path '*.pdf'

    """
    cp -L $project_dir/bin/Analysis/visualize_marker_genes_DM1.Rmd visualize_marker_genes_DM1.Rmd
    cp -L $project_dir/bin/support_functions/Feature_plots.R Feature_plots.R
    Rscript -e  "rmarkdown::render('visualize_marker_genes_DM1.Rmd', output_format = 'html_document', output_file = 'visualize_marker_genes_DM1.html', 
                                    params = list(mae_hdf5_path = '${mae_path}', 
                                                  config_file = '${r_config}'))"
    """
    
}

process DIFF_METH_DM1 {

    publishDir "${params.output}/DM1", mode: 'copy', overwrite: true

    label 'full_resources'
    memory '30 GB'

    input:
    path meth_dir
    path isa_file
    
    output:
    path 'DMRS_DM1_PatientsControls.csv'

    """
    cp -L $project_dir/bin/Analysis/Diff_Meth_DM1_DMPK.Rmd Diff_Meth_DM1_DMPK.Rmd
    Rscript -e  "rmarkdown::render('Diff_Meth_DM1_DMPK.Rmd', output_format = 'html_document', output_file = 'Diff_Meth_DM1_DMPK.html', 
                                    params = list(methrix_qc_dir = '${meth_dir}', isa_path = '${isa_file}' ))"
    """
}

process PLOT_DMPK_REP {

    publishDir "${params.output}/DM1", mode: 'copy', overwrite: true

    label 'full_resources'

    input:
    path meth_dir
    path dmrs_file
    path isa_file
    path mae_object
    
    output:
    path '*.pdf'

    """
    cp -L $project_dir/bin/EMseq_filter_DMPK.Rmd EMseq_filter_DMPK.Rmd
    cp -L $project_dir/bin/support_functions/ReadMethrix2Df.R ReadMethrix2Df.R
    Rscript -e  "rmarkdown::render('EMseq_filter_DMPK.Rmd', output_format = 'html_document', output_file = 'Plots_DM1_DMPK.html', 
                                    params = list(  methrix_dir = '${meth_dir}',
                                                    DMRS_DM1 = '${dmrs_file}',
                                                    isa_file = '${isa_file}',
                                                    mae_hdf5_path = '${mae_object}' ))"
    """
}


process PLOT_DMPK_REP_CLONES {

    publishDir "${params.output}/DM1", mode: 'copy', overwrite: true

    label 'full_resources'

    input:
    path meth_dir
    path isa_file
    path mae_object
    
    output:
    path '*.pdf'

    """
    cp -L $project_dir/bin/EMseq_filter_DMPK_CLONES.Rmd EMseq_filter_DMPK_CLONES.Rmd
    cp -L $project_dir/bin/support_functions/ReadMethrix2Df.R ReadMethrix2Df.R
    Rscript -e  "rmarkdown::render('EMseq_filter_DMPK_CLONES.Rmd', output_format = 'html_document', output_file = 'Diff_Meth_Clones_DM1_DMPK_CLONES.html', 
                                    params = list(  methrix_dir = '${meth_dir}',
                                                    isa_file = '${isa_file}',
                                                    mae_hdf5_path = '${mae_object}' ))"
    """
}
