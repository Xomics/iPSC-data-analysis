#!/usr/bin/env nextflow


nextflow.enable.dsl=2

project_dir = projectDir



process RUN_RLM {

    publishDir "${params.output}/RLM", mode: 'copy', overwrite: true

    label 'full_resources'

    input:
    path mae_path
    path config_file
    
    output:
    path '*.csv'

    """
    cp -L $project_dir/bin/Analysis/run_RLM_MAE.Rmd run_RLM_MAE.Rmd
    Rscript -e  "rmarkdown::render('run_RLM_MAE.Rmd', output_format = 'html_document', output_file = 'run_RLM_MAE.html', 
                                    params = list(mae_hdf5_path = '${mae_path}', config_file = '${config_file}' ))"
    """
}


process RUN_MIXED_MODELS {

    publishDir "${params.output}/Mixed_LM", mode: 'copy', overwrite: true

    label 'full_resources'

    input:
    path mae_path
    path config_file
    
    output:
    path 'MixedModelsResults.csv'
    path '*.xlsx'

    """
    cp -L $project_dir/bin/Analysis/run_linear_mixed_model_MAE.Rmd run_linear_mixed_model_MAE.Rmd
    cp -L $project_dir/bin/Analysis/Mixed_Model_functions.R Mixed_Model_functions.R
    Rscript -e  "rmarkdown::render('run_linear_mixed_model_MAE.Rmd', output_format = 'html_document', output_file = 'run_linear_mixed_model_MAE.html', 
                                    params = list(mae_hdf5_path = '${mae_path}', config_file = '${config_file}' ))"
    """
}

process RUN_MIXED_MODELS_CORRECTED {

    publishDir "${params.output}/Mixed_LM_CORRECTED", mode: 'copy', overwrite: true

    label 'full_resources'

    input:
    path mae_path
    path config_file
    path mofa_scores
    
    output:
    path 'MixedModelsResults_Corrected.csv'
    path '*.xlsx'

    """
    cp -L $project_dir/bin/Analysis/run_linear_mixed_model_MAE_CORRECTED.Rmd run_linear_mixed_model_MAE_CORRECTED.Rmd
    cp -L $project_dir/bin/Analysis/Mixed_Model_functions.R Mixed_Model_functions.R
    Rscript -e  "rmarkdown::render('run_linear_mixed_model_MAE_CORRECTED.Rmd', output_format = 'html_document', output_file = 'run_linear_mixed_model_MAE_CORRECTED.html', 
                                    params = list(mae_hdf5_path = '${mae_path}', config_file = '${config_file}', mofa_scores = '${mofa_scores}' ))"
    """
}


process PREP_CHUNKS_MIXED_MODELS {

    //publishDir "${params.output}/Mixed_LM", mode: 'copy', overwrite: true

    label 'full_resources'

    input:
    path mae_path
    path config_file
    
    output:
    path '*.csv'

    """
    cp -L $project_dir/bin/Analysis/prepare_chunks_linear_models.Rmd prepare_chunks_linear_models.Rmd
    cp -L $project_dir/bin/Analysis/Mixed_Model_functions.R Mixed_Model_functions.R
    Rscript -e  "source('Mixed_Model_functions.R'); \
                 rmarkdown::render('prepare_chunks_linear_models.Rmd', output_format = 'html_document', output_file = 'prepare_chunks_linear_models.html', 
                                    params = list(mae_hdf5_path = '${mae_path}', config_file = '${config_file}' ))"
    """
}

process PREP_CORRECTED_CHUNKS_MIXED_MODELS {

    //publishDir "${params.output}/Mixed_LM", mode: 'copy', overwrite: true

    label 'full_resources'

    input:
    path mae_path
    path config_file
    path mofa_scores
    
    output:
    path '*.csv'

    """
    cp -L $project_dir/bin/Analysis/prepare_corrected_chunks_linear_models.Rmd prepare_corrected_chunks_linear_models.Rmd
    cp -L $project_dir/bin/Analysis/Mixed_Model_functions.R Mixed_Model_functions.R
    Rscript -e  "source('Mixed_Model_functions.R'); \
                 rmarkdown::render('prepare_corrected_chunks_linear_models.Rmd', output_format = 'html_document', output_file = 'prepare_corrected_chunks_linear_models.html', 
                                    params = list(mae_hdf5_path = '${mae_path}', config_file = '${config_file}', mofa_scores = '${mofa_scores}' ))"
    """
}


process MIXED_LINEAR_MODELS {

    //publishDir "${params.output}/Mixed_LM", mode: 'copy', overwrite: true

    label 'full_resources'

    input:
    each path(mae_path)
    path chunk_path
    
    output:
    path 'Model_*.csv'

    """
    cp -L $project_dir/bin/Analysis/run_linear_mixed_model_CHUNK.Rmd run_linear_mixed_model_CHUNK.Rmd
    cp -L $project_dir/bin/Analysis/Mixed_Model_functions.R Mixed_Model_functions.R
    Rscript -e  "source('Mixed_Model_functions.R'); \
                 rmarkdown::render('run_linear_mixed_model_CHUNK.Rmd', output_format = 'html_document', output_file = 'run_linear_mixed_model_CHUNK.html', 
                                    params = list(chunk = '${chunk_path}', mae_object = '${mae_path}' ))"
    """
    
}


process MERGE_CSV_FILES {

    publishDir "${params.output}/Mixed_LM_CHUNKS", mode: 'copy', overwrite: true

    input:
    path csv_files

    output:
    path 'Model_Results_merged.csv'

    script:
    """
    first_file=\$(ls ${csv_files} | head -n 1)
    echo "Mergin files..."
    head -n 1 \$first_file > Model_Results_merged.csv
    for file in ${csv_files}; do
	tail -n +2 \$file >> Model_Results_merged.csv
    done
    """

}

process MERGE_EMSEQ_MAE_RESULTS {

    publishDir "${params.output}/Mixed_LM", mode: 'copy', overwrite: true

    input:
    path chunks_results
    path mae_results

    output:
    path 'Model_Results_All_merged.csv'

    script:
    """
    cp -L $project_dir/bin/Analysis/Merge_linear_models_chunks_mae.Rmd Merge_linear_models_chunks_mae.Rmd
    Rscript -e  "rmarkdown::render('Merge_linear_models_chunks_mae.Rmd.Rmd', output_format = 'html_document',  
                                    params = list(merged_chunks = '${chunks_results}', merged_mae_results = '${mae_results}' ))"
    """

}
