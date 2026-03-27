#!/usr/bin/env nextflow


nextflow.enable.dsl=2

project_dir = projectDir

process VIS_MIXED_MODELS {

    publishDir "${params.output}/MixedModels_Plots", mode: 'copy', overwrite: true

    label 'full_resources'

    input:
    path mae_path
    path config_file
    path ModelResults
    
    output:
    path '*.pdf'

    """
    cp -L $project_dir/bin/Analysis/visualize_MixeModels_results.Rmd visualize_MixeModels_results.Rmd
    cp -L $project_dir/bin/support_functions/Feature_plots.R Feature_plots.R
    Rscript -e  "rmarkdown::render('visualize_MixeModels_results.Rmd', output_format = 'html_document', output_file = 'Viz_sign_hits.html', 
                                    params = list(mae_hdf5_path = '${mae_path}', config_file = '${config_file}', MixedModelsResults = '${ModelResults}' ))"
    """
}

process CLUSTER_SIGN_HITS {

    publishDir "${params.output}/RLM", mode: 'copy', overwrite: true

    label 'full_resources'

    input:
    path mae_path
    path config_file
    
    output:
    path '*.csv'

    """
    cp -L $project_dir/bin/Analysis/Cluster_samples_with_sign_hits.Rmd Cluster_samples_with_sign_hits.Rmd
    Rscript -e  "rmarkdown::render('Cluster_samples_with_sign_hits.Rmd', output_format = 'html_document', output_file = 'Cluster_samples_with_sign_hits.html', 
                                    params = list(mae_hdf5_path = '${mae_path}', config_file = '${config_file}' ))"
    """
}


