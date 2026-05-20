#!/usr/bin/env nextflow

nextflow.enable.dsl=2

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
		System.exit(0)
	}

	// Initialise required channels
	if (params.mae_object) {
		mae_object_ch = Channel.fromPath(params.mae_object, type: 'dir', checkIfExists: true)
	} else {
		log.info "mae_object path needs to exist!"
		System.exit(0)
	}

	if (params.support_functions) {
		support_functions_ch = Channel.fromPath(params.support_functions, type: 'dir', checkIfExists: true)
	} else {
		log.info "support_functions path needs to exist!"
		System.exit(0)
	}

	// if (params.temp_file) {
	// 	temp_files_ch = Channel.fromPath("${params.temp_file}/*.csv", checkIfExists: true)
	// 	temp_file_ch = Channel.fromPath("${params.temp_file}", type: 'dir', checkIfExists: true)
	// } else {
	// 	log.info "temp_file path needs to exist!"
	// 	System.exit(0)
	// }

	
	//
	// ────────────────────────────────
	// Single omics analyses
	// ────────────────────────────────
	//

	// RUN_CV_PLOTS(mae_object_ch, params.r_config)
	RUN_PCA(mae_object_ch, params.r_config)
	//PAIRWISE_CORRELATIONS(mae_object_ch, params.r_config, params.isa_file, params.prot_tech, params.mtblmcs_tech)


	//
	// ────────────────────────────────
	// Single features analyses
	// ────────────────────────────────
	//
	//VIS_MARKERS(mae_object_ch, params.r_config, params.neuron_markers, params.tmem_markers)
	//VIS_SYMPATHIC_MARKERS(mae_object_ch, params.r_config, params.sympathic_markers)
	//DOXYCYCLINE_COR(mae_object_ch, params.r_config)


	//
	// ────────────────────────────────
	// Linear models
	// ────────────────────────────────
	//
	// RUN_RLM(mae_object_ch, params.r_config)

	// Choose prep process based on a parameter
	if (params.use_chunks) {
		// Prep chunks (corrected or standard)
		if (params.use_corrected) {
			PREP_CORRECTED_CHUNKS(mae_object_ch, params.r_config, params.mofa_scores)
			prep_out = PREP_CORRECTED_CHUNKS.out.flatten()
		} else {
			PREP_CHUNKS(mae_object_ch, params.r_config)
			prep_out = PREP_CHUNKS.out.flatten()
		}

		// Run MLM on chunks
		MIXED_LINEAR_MODELS(mae_object_ch, prep_out)
		MERGE_CSV_FILES(MIXED_LINEAR_MODELS.out.toList())

		// Merge chunk results with mixed models, then visualize
		if (params.use_corrected) {
			RUN_MIXED_MODELS_CORRECTED(mae_object_ch, params.r_config, params.mofa_scores)
			MERGE_EMSEQ_MAE_RESULTS_CORRECTED(MERGE_CSV_FILES.out, RUN_MIXED_MODELS_CORRECTED.out[0])
			VIS_MIXED_MODELS_CORRECTED(mae_object_ch, params.r_config, MERGE_EMSEQ_MAE_RESULTS_CORRECTED.out)
		} else {
			// RUN_MIXED_MODELS(mae_object_ch, params.r_config)
			MERGE_EMSEQ_MAE_RESULTS_NORMAL(MERGE_CSV_FILES.out, RUN_MIXED_MODELS.out[0])
			VIS_MIXED_MODELS_NORMAL(params.mae_object, params.r_config, MERGE_EMSEQ_MAE_RESULTS_NORMAL.out)
		}

	} else {
	// No chunks — run directly and visualize
	if (params.use_corrected) {
		RUN_MIXED_MODELS_CORRECTED(mae_object_ch, params.r_config, params.mofa_scores)
		VIS_MIXED_MODELS_CORRECTED(mae_object_ch, params.r_config, RUN_MIXED_MODELS_CORRECTED.out[0])
	} else {
		// RUN_MIXED_MODELS(mae_object_ch, params.r_config)
		//VIS_MIXED_MODELS_NORMAL(mae_object_ch, params.r_config, RUN_MIXED_MODELS.out[0])
	}
	}
	//CLUSTER_SIGN_HITS(mae_object_ch, params.r_config)


	//
	// ────────────────────────────────
	// Multi-omics analyses
	// ────────────────────────────────
	//
	RUN_MOFA(support_functions_ch, mae_object_ch, params.r_config, params.neuron_markers)
	//ANALYZE_MOFA_FACTORS(mae_object_ch, params.r_config, RUN_MOFA.out[1])
	//ANALYZE_MOFA_FACTORS(mae_object_ch, params.r_config, params.mofa_scores)
	//ANALYZE_MOFA_LIPIDS(support_functions_ch, mae_object_ch, params.r_config, temp_file_ch, params.mofa_scores, params.lipid_class_data)
	//ANALYZE_MOFA_LIPIDS(support_functions_ch, mae_object_ch, params.r_config, RUN_MOFA.out[1], params.mofa_scores, params.lipid_class_data)
	//RUN_MOFA_NANS(mae_object_ch, params.r_config)
	//RUN_MOFA_CHD2(mae_object_ch, params.r_config)
	//RUN_MOFA_DM1(mae_object_ch, params.r_config)


	//
	// ────────────────────────────────
	// DNA methylation analyses
	// ────────────────────────────────
	//
	//METH_COR_IPSC_NEURONS(meth_ipsc_qc, meth_neuron_sub)
	//DIFF_METH_CLONES(meth_neuron_qc, params.isa_file)


	//
	// ────────────────────────────────
	// 			DM1 analyses
	// ────────────────────────────────
	//
	//VIS_MARKERS_DM1(mae_object_ch, params.r_config)
	//DIFF_METH_DM1(meth_neuron_qc , params.isa_file)
	//PLOT_DMPK_REP(meth_neuron_raw, DIFF_METH_DM1.out, params.isa_file, mae_object_ch)
	//PLOT_DMPK_REP_CLONES(meth_ipsc_raw, params.isa_file, mae_object_ch)
}
