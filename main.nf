#!/usr/bin/env nextflow

nextflow.enable.dsl=2

project_dir = projectDir


////////////////////////////////////////////////////
/*    --               Functions               -- */
////////////////////////////////////////////////////+


def helpMessage() {
  log.info """
        Usage:
        The typical command for running the pipeline is as follows:
        nextflow run main.nf 
            --output dir/of/choice
            --mae_object path/to/mae_object
	    	--mudata path/to/mudata.hd5	
	    	--r_config path/to/r_config.yml
			--container_dir	Dir where Singularity (.sif files) images are stored
			--markers_dir path/to/markers Mapping files for molecular markers
		Optional:
			--use_corrected Perform correction for confounding MOFA factors before running linear models
			--use_chunks Run linear mixed models for Emseq data in chunks to process faster
			--mofa_scores path/to/mofa_scores If provided
			--alt_omics_dir Directory with additional omics data files not in the MAE object
			--meth_dir Directory with additonal DNA methylation data files not in the MAE object
	"""
}



////////////////////////////////////////////////////
/* --          Main input data file            -- */
////////////////////////////////////////////////////+


mae_object_ch = Channel.fromPath("${params.mae_object}", type: 'dir')
mudata = Channel.fromPath("${params.mudata}")  
r_config_ch = Channel.fromPath("${params.r_config}")
support_functions_ch = Channel.fromPath("$project_dir/bin/support_functions", type: 'dir')
isa_file_ch = Channel.fromPath("${params.isa_file}")


////////////////////////////////////////////////////
/* --      Input files from directories         -- */
////////////////////////////////////////////////////+

neuron_markers_ch = Channel.fromPath("${params.markers_dir}/Neuron_markers.csv")
tmem_markers_ch = Channel.fromPath("${params.markers_dir}/TMEM_mappings.csv")
sympathic_markers_ch = Channel.fromPath("${params.markers_dir}/Sympathic_marker_genes.csv")

lipid_class_data_ch = Channel.fromPath("${params.alt_omics_dir}/lipidomics_T_Class_data.csv")
mtblmcs_tech_ch = Channel.fromPath("${params.alt_omics_dir}/metabolomics_tml_data_tech_rep_complete.csv")
prot_tech_ch = Channel.fromPath("${params.alt_omics_dir}/proteomics_data_tech_reps_complete.csv")

meth_ipsc_raw   	= file("${params.meth_dir}/iPSC_h5_raw")   
meth_ipsc_qc  		= file("${params.meth_dir}/iPSC_h5_QC_SNP_filtered") 
meth_neuron_raw   	= file("${params.meth_dir}/h5_raw")   
meth_neuron_qc    	= file("${params.meth_dir}/h5_QC_SNP_filtered") 
meth_neuron_sub    	= file("${params.meth_dir}/h5_QC_SNP_filtered_sub_SD")   



// Temp files for debugging
mofa_scores = Channel.fromPath("${params.mofa_scores}")
temp_files_ch = Channel.fromPath("${params.temp_file}/*.csv")
temp_file_ch = Channel.fromPath("${params.temp_file}")



////////////////////////////////////////////////////
/* --                  Modules                 -- */
////////////////////////////////////////////////////+


include { PREP_CHUNKS_MIXED_MODELS as PREP_CHUNKS} from './modules/linear_models.nf'
include { PREP_CORRECTED_CHUNKS_MIXED_MODELS as PREP_CORRECTED_CHUNKS } from './modules/linear_models.nf'
include { MIXED_LINEAR_MODELS; MERGE_CSV_FILES; RUN_RLM; RUN_MIXED_MODELS; RUN_MIXED_MODELS_CORRECTED } from './modules/linear_models.nf'
include { MERGE_EMSEQ_MAE_RESULTS as MERGE_EMSEQ_MAE_RESULTS_CORRECTED; MERGE_EMSEQ_MAE_RESULTS as MERGE_EMSEQ_MAE_RESULTS_NORMAL } from './modules/linear_models.nf'
include { VIS_MIXED_MODELS as VIS_MIXED_MODELS_CORRECTED; VIS_MIXED_MODELS as VIS_MIXED_MODELS_NORMAL } from './modules/viz_linear_models.nf'
include { RUN_CV_PLOTS; RUN_PCA; PAIRWISE_CORRELATIONS } from './modules/single_omics.nf'
include { VIS_MARKERS; VIS_SYMPATHIC_MARKERS; DOXYCYCLINE_COR } from './modules/single_features.nf'
include { RUN_MOFA; ANALYZE_MOFA_FACTORS; ANALYZE_MOFA_LIPIDS; RUN_MOFA_NANS; RUN_MOFA_CHD2; RUN_MOFA_DM1 } from './modules/multi_omics.nf'
include { METH_COR_IPSC_NEURONS; DIFF_METH_CLONES} from  './modules/methylation.nf'
include { VIS_MARKERS_DM1; DIFF_METH_DM1; PLOT_DMPK_REP; PLOT_DMPK_REP_CLONES } from  './modules/DM1.nf'


////////////////////////////////////////////////////
/* --                 Workflow                 -- */
////////////////////////////////////////////////////+


workflow {

	// Show help message
	if (params.help) {
	helpMessage()
	exit 0
	}


	//
	// ────────────────────────────────
	// Single omics analyses
	// ────────────────────────────────
	//
	//RUN_CV_PLOTS(mae_object_ch, r_config_ch)
	//RUN_PCA(mae_object_ch, r_config_ch)
	//PAIRWISE_CORRELATIONS(mae_object_ch, r_config_ch, isa_file_ch, prot_tech_ch, mtblmcs_tech_ch)


	//
	// ────────────────────────────────
	// Single features analyses
	// ────────────────────────────────
	//
	//VIS_MARKERS(mae_object_ch, r_config_ch, neuron_markers_ch, tmem_markers_ch)
	//VIS_SYMPATHIC_MARKERS(mae_object_ch, r_config_ch, sympathic_markers_ch)
	//DOXYCYCLINE_COR(mae_object_ch, r_config_ch)


	//
	// ────────────────────────────────
	// Linear models
	// ────────────────────────────────
	//
	//RUN_RLM(mae_object_ch, r_config_ch)
	//VIS_RLM()

	// Choose prep process based on a parameter
	if (params.use_chunks) {
		// Prep chunks (corrected or standard)
		if (params.use_corrected) {
			PREP_CORRECTED_CHUNKS(params.mae_object, r_config_ch, mofa_scores)
			prep_out = PREP_CORRECTED_CHUNKS.out.flatten()
		} else {
			PREP_CHUNKS(params.mae_object, r_config_ch)
			prep_out = PREP_CHUNKS.out.flatten()
		}

		// Run MLM on chunks
		MIXED_LINEAR_MODELS(mae_object_ch, prep_out)
		MERGE_CSV_FILES(MIXED_LINEAR_MODELS.out.toList())

		// Merge chunk results with mixed models, then visualize
		if (params.use_corrected) {
			RUN_MIXED_MODELS_CORRECTED(mae_object_ch, r_config_ch, mofa_scores)
			MERGE_EMSEQ_MAE_RESULTS_CORRECTED(MERGE_CSV_FILES.out, RUN_MIXED_MODELS_CORRECTED.out[0])
			VIS_MIXED_MODELS_CORRECTED(params.mae_object, r_config_ch, MERGE_EMSEQ_MAE_RESULTS_CORRECTED.out)
		} else {
			RUN_MIXED_MODELS(mae_object_ch, r_config_ch)
			MERGE_EMSEQ_MAE_RESULTS_NORMAL(MERGE_CSV_FILES.out, RUN_MIXED_MODELS.out[0])
			VIS_MIXED_MODELS_NORMAL(params.mae_object, r_config_ch, MERGE_EMSEQ_MAE_RESULTS_NORMAL.out)
		}

	} else {
	// No chunks — run directly and visualize
	if (params.use_corrected) {
		RUN_MIXED_MODELS_CORRECTED(mae_object_ch, r_config_ch, mofa_scores)
		VIS_MIXED_MODELS_CORRECTED(params.mae_object, r_config_ch, RUN_MIXED_MODELS_CORRECTED.out[0])
	} else {
		//RUN_MIXED_MODELS(mae_object_ch, r_config_ch)
		//VIS_MIXED_MODELS_NORMAL(params.mae_object, r_config_ch, RUN_MIXED_MODELS.out[0])
	}
	}
	//CLUSTER_SIGN_HITS(params.mae_object, r_config_ch)


	//
	// ────────────────────────────────
	// Multi-omics analyses
	// ────────────────────────────────
	//
	//RUN_MOFA(support_functions_ch, params.mae_object, r_config_ch, neuron_markers_ch)
	//ANALYZE_MOFA_FACTORS(params.mae_object, r_config_ch, RUN_MOFA.out[1])
	//ANALYZE_MOFA_FACTORS(params.mae_object, r_config_ch, mofa_scores)
	//ANALYZE_MOFA_LIPIDS(support_functions_ch, params.mae_object, r_config_ch, temp_file_ch, mofa_scores, lipid_class_data_ch)
	//ANALYZE_MOFA_LIPIDS(support_functions_ch, params.mae_object, r_config_ch, RUN_MOFA.out[1], mofa_scores, lipid_class_data_ch)
	//RUN_MOFA_NANS(params.mae_object, r_config_ch)
	//RUN_MOFA_CHD2(params.mae_object, r_config_ch)
	//RUN_MOFA_DM1(params.mae_object, r_config_ch)


	//
	// ────────────────────────────────
	// DNA methylation analyses
	// ────────────────────────────────
	//
	//METH_COR_IPSC_NEURONS(meth_ipsc_qc, meth_neuron_sub)
	//DIFF_METH_CLONES(meth_neuron_qc, isa_file_ch)


	//
	// ────────────────────────────────
	// 			DM1 analyses
	// ────────────────────────────────
	//
	//VIS_MARKERS_DM1(mae_object_ch, r_config_ch)
	DIFF_METH_DM1(meth_neuron_qc , isa_file_ch)
	PLOT_DMPK_REP(meth_neuron_raw, DIFF_METH_DM1.out, isa_file_ch, mae_object_ch)
	//PLOT_DMPK_REP_CLONES(meth_ipsc_raw, isa_file_ch, mae_object_ch)
}
