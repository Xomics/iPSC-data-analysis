#!/usr/bin/env nextflow


nextflow.enable.dsl=2

project_dir = projectDir



process DIFF_METH_CLONES {

    publishDir "${params.output}/Meth", mode: 'copy', overwrite: true

    label 'full_resources'

    input:
    path meth_dir
    path isa_file
    
    output:
    path '*.csv'

    """
    cp -L $project_dir/bin/Analysis/Diff_Meth_Clones.Rmd Diff_Meth_Clones.Rmd
    Rscript -e  "rmarkdown::render('Diff_Meth_Clones.Rmd', output_format = 'html_document', output_file = 'Diff_Meth_Clones.html', 
                                    params = list(methrix_qc_dir = '${meth_dir}', isa_path = '${isa_file}' ))"
    """
}

process METH_COR_IPSC_NEURONS {

    publishDir "${params.output}/Meth", mode: 'copy', overwrite: true

    label 'full_resources'

    input:
    path meth_ipsc
    path meth_neuron
    
    output:
    path '*.pdf'

    """
    cp -L $project_dir/bin/Analysis/Meth_correlations_iPSC_Neurons.Rmd Meth_correlations_iPSC_Neurons.Rmd
    Rscript -e  "rmarkdown::render('Meth_correlations_iPSC_Neurons.Rmd', output_format = 'html_document', output_file = 'Meth_correlations_iPSC_Neurons.html', 
                                    params = list(meth_iPSC = '${meth_ipsc}', meth_iNeuron = '${meth_neuron}' ))"
    """
}
