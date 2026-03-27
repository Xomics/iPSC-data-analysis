#!/usr/bin/env nextflow


nextflow.enable.dsl=2

project_dir = projectDir


process RUN_MOFA {

    publishDir "${params.output}/MOFA", mode: 'copy', overwrite: true

    label 'full_resources'

    input:
    path support_functions_ch
    path mae_path
    path r_config
    path neuron_markers
    
    output:
    path '*.hdf5'
    path 'SampleScores.csv'

    """
    cp -L $project_dir/bin/Analysis/run_MOFA.Rmd run_MOFA.Rmd
    Rscript -e  "rmarkdown::render('run_MOFA.Rmd', output_format = 'html_document', output_file = 'run_MOFA.html', params = list(mae_hdf5_path = '${mae_path}', config_file = '${r_config}', neuron_markers = '${neuron_markers}' ))"
    """
}


process ANALYZE_MOFA_FACTORS {

    publishDir "${params.output}/MOFA", mode: 'copy', overwrite: true

    label 'full_resources'

    input:
    path mae_path
    path r_config
    path mofa_scores
    
    output:

    """
    cp -L $project_dir/bin/Analysis/AnalyzeMOFAFactors.Rmd AnalyzeMOFAFactors.Rmd
    Rscript -e  "rmarkdown::render('AnalyzeMOFAFactors.Rmd.Rmd', output_format = 'html_document', output_file = 'AnalyzeMOFAFactors.Rmd.html', params = list(mae_hdf5_path = '${mae_path}', config_file = '${r_config}', mofa_scores = '${mofa_scores}' ))"
    """
}


process ANALYZE_MOFA_LIPIDS {

    publishDir "${params.output}/MOFA", mode: 'copy', overwrite: true

    label 'full_resources'

    input:
    path support_functions_ch
    path mae_path
    path r_config
    path mofa_model
    path mofa_scores
    path lipid_class_data
    
    output:
    path "analyze_MOFA_lipids.html"


    """
    cp -L $project_dir/bin/Analysis/run_MOFA_analyzeLipids.Rmd run_MOFA_analyzeLipids.Rmd
    Rscript -e  "rmarkdown::render('run_MOFA_analyzeLipids.Rmd', 
                                    output_format = 'html_document', 
                                    output_file = 'analyze_MOFA_lipids.html', 
                                    params = list(mae_hdf5_path = '${mae_path}', 
                                                  config_file = '${r_config}', 
                                                  mofa_model = '${mofa_model}', 
                                                  mofa_scores = '${mofa_scores}',  
                                                  lipid_class_data = '${lipid_class_data}'))"
    """
}

process RUN_MOFA_NANS {

    publishDir "${params.output}/MOFA_NANS", mode: 'copy', overwrite: true

    label 'full_resources'

    input:
    path mae_path
    path r_config
    
    output:
    path '*.pdf'
    path 'run_MOFA_NANS.html'


    """
    cp -L $project_dir/bin/Analysis/run_MOFA_NANS.Rmd run_MOFA_NANS.Rmd
    Rscript -e  "rmarkdown::render('run_MOFA_NANS.Rmd', output_format = 'html_document', output_file = 'run_MOFA_NANS.html', params = list(mae_hdf5_path = '${mae_path}', config_file = '${r_config}'))"
    """
}

process RUN_MOFA_DM1 {

    publishDir "${params.output}/MOFA_DM1", mode: 'copy', overwrite: true

    label 'full_resources'

    input:
    path mae_path
    path r_config
    
    output:
    path '*.pdf'
    path 'run_MOFA_DM1.html'

    """
    cp -L $project_dir/bin/Analysis/run_MOFA_DM1.Rmd run_MOFA_DM1.Rmd
    Rscript -e  "rmarkdown::render('run_MOFA_DM1.Rmd', output_format = 'html_document', output_file = 'run_MOFA_DM1.html', params = list(mae_hdf5_path = '${mae_path}', config_file = '${r_config}'))"
    """
}

process RUN_MOFA_CHD2 {

    publishDir "${params.output}/MOFA_CHD2", mode: 'copy', overwrite: true

    label 'full_resources'

    input:
    path mae_path
    path r_config
    
    output:
    path '*.pdf'
    path 'run_MOFA_CHD2.html'

    """
    cp -L $project_dir/bin/Analysis/run_MOFA_CHD2.Rmd run_MOFA_CHD2.Rmd
    Rscript -e  "rmarkdown::render('run_MOFA_CHD2.Rmd', output_format = 'html_document', output_file = 'run_MOFA_CHD2.html', params = list(mae_hdf5_path = '${mae_path}', config_file = '${r_config}'))"
    """
}
