#!/usr/bin/env nextflow


nextflow.enable.dsl=2

project_dir = projectDir



process RUN_CV_PLOTS {

    publishDir "${params.output}/CV", mode: 'copy', overwrite: true

    label 'full_resources'

    input:
    path mae_path
    path r_config
    
    output:
    path 'run_CV_MAE.Rmd.html'

    """
    cp -L $project_dir/bin/Analysis/run_CV_MAE.Rmd run_CV_MAE.Rmd
    Rscript -e  "rmarkdown::render('run_CV_MAE.Rmd', output_format = 'html_document', output_file = 'run_CV_MAE.Rmd.html', 
                                    params = list(mae_hdf5_path = '${mae_path}', config_file = '${r_config}' ))"
    """
    
}

process RUN_PCA {

    publishDir "${params.output}/PCA", mode: 'copy', overwrite: true

    label 'full_resources'

    input:
    path mae_path
    path r_config
    
    output:
    path 'run_PCA_MAE.html'
    path '*.pdf'

    """
    cp -L $project_dir/bin/Analysis/run_PCA_MAE.Rmd run_PCA_MAE.Rmd
    Rscript -e  "rmarkdown::render('run_PCA_MAE.Rmd', output_format = 'html_document', output_file = 'run_PCA_MAE.html', 
                                    params = list(mae_hdf5_path = '${mae_path}', config_file = '${r_config}' ))"
    """
    
}


process PAIRWISE_CORRELATIONS {

    publishDir "${params.output}/PairCorr", mode: 'copy', overwrite: true

    label 'full_resources'

    input:
    path mae_path
    path r_config
    path isa_file
    path prot_tech
    path mtblmcs_tech
    
    output:
    path 'run_Pairwise_Correlations_MAE.html'
    path '*.pdf'

    """
    cp -L $project_dir/bin/Analysis/run_Pairwise_Correlations_MAE.Rmd run_Pairwise_Correlations_MAE.Rmd
    cp -L $project_dir/bin/support_functions/PairwiseCorrelationsFunctions.R PairwiseCorrelationsFunctions.R
    Rscript -e  "rmarkdown::render('run_Pairwise_Correlations_MAE.Rmd', output_format = 'html_document', output_file = 'run_Pairwise_Correlations_MAE.html', 
                                    params = list(mae_hdf5_path = '${mae_path}', 
                                                  config_file = '${r_config}',
                                                  isa_file  = '${isa_file}',
                                                  prot_tech_reps = '${prot_tech}',
                                                  metab_tech_reps = '${mtblmcs_tech}' ))"
    """
    
}
