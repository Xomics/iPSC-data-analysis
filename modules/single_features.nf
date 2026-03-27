#!/usr/bin/env nextflow


nextflow.enable.dsl=2

project_dir = projectDir



process VIS_MARKERS {

    publishDir "${params.output}/MARKERS", mode: 'copy', overwrite: true

    label 'full_resources'

    input:
    path mae_path
    path r_config
    path neuron_markers
    path TMEM_mappings
    
    output:
    path 'visualize_markers.html'

    """
    cp -L $project_dir/bin/Analysis/visualize_marker_genes.Rmd visualize_marker_genes.Rmd
    cp -L $project_dir/bin/support_functions/Feature_plots.R Feature_plots.R
    Rscript -e  "rmarkdown::render('visualize_marker_genes.Rmd', output_format = 'html_document', output_file = 'visualize_markers.html', 
                                    params = list(mae_hdf5_path = '${mae_path}', 
                                                  config_file = '${r_config}',
                                                  neuron_markers = '${neuron_markers}',
                                                  TMEM_mappings = '${TMEM_mappings}' ))"
    """
    
}

process VIS_SYMPATHIC_MARKERS {

    publishDir "${params.output}/MARKERS", mode: 'copy', overwrite: true

    label 'full_resources'

    input:
    path mae_path
    path r_config
    path SYMPATHIC_markers
    
    output:
    path 'visualize_SYMPATHIC_markers.html'

    """
    cp -L $project_dir/bin/Analysis/visualize_SYMPATHIC_marker_genes.Rmd visualize_SYMPATHIC_marker_genes.Rmd
    cp -L $project_dir/bin/support_functions/Feature_plots.R Feature_plots.R
    Rscript -e  "rmarkdown::render('visualize_SYMPATHIC_marker_genes.Rmd', output_format = 'html_document', output_file = 'visualize_SYMPATHIC_markers.html', 
                                    params = list(mae_hdf5_path = '${mae_path}', 
                                                  config_file = '${r_config}',
                                                  SYMPATHIC_markers = '${SYMPATHIC_markers}' ))"
    """
    
}


process DOXYCYCLINE_COR {

    publishDir "${params.output}/MOFA/Factor1", mode: 'copy', overwrite: true

    label 'full_resources'

    input:
    path mae_path
    path r_config
    
    output:
    path 'Doxycycline_cor.html'
    path '*.csv'

    """
    cp -L $project_dir/bin/Analysis/Doxycycline_cor.Rmd Doxycycline_cor.Rmd
    cp -L $project_dir/bin/support_functions/Feature_plots.R Feature_plots.R
    Rscript -e  "rmarkdown::render('Doxycycline_cor.Rmd', output_format = 'html_document', output_file = 'Doxycycline_cor.html', 
                                    params = list(mae_hdf5_path = '${mae_path}', 
                                                  config_file = '${r_config}' ))"
    """
    
}

