// This file is used to generate config/configfile.yaml
// To regenerate YAML, run
//     make config/configfile.yaml

conda_environment: "workflow/envs/nextstrain.yaml"
genesforglycosylation: ["G", "F"]
builds_to_run: ["genome", "G", "F"]
resolutions_to_run: ["all-time", "6y", "3y"]
description: "config/description.md"
strain_id_field: "accession"
display_strain_field: "strain"
subtypes: ["a", "b"]

subsample: {
	for build in builds_to_run
	for resolution in resolutions_to_run {
		"\(build)/\(resolution)": {
			samples: {
				for sampleName, sampleConfig in _resolutionSamples[resolution] {
					"\(sampleName)": _buildConfigs[build] & sampleConfig
				}
			}
		}
	}
}

_buildConfigs: {
	genome: {
		min_length: 10000
		query: "genome_coverage>0.3 & missing_data<1000"
	}
	G: {
		min_length: 600
		query: "G_coverage>0.3 & missing_data<1000"
	}
	F: {
		min_length: 1200
		query: "F_coverage>0.3 & missing_data<1000"
	}
}

_resolutionSamples: {
	"all-time": {
		"all-time": {
			min_date: "1975-01-01"
			max_sequences: 3000
			exclude_where: [
				"qc.overallStatus=bad",
			]
		} & _defaultSubsampleConfig
	}

	for resolutionName, cutoff in _cutoffs {
		"\(resolutionName)": {
			recent: {
				min_date: cutoff
				max_sequences: 3000
				exclude_where: [
					"qc.overallStatus=bad",
				]
			} & _defaultSubsampleConfig
			background: {
				max_date: cutoff
				max_sequences: 300
				exclude_where: [
					"qc.overallStatus=bad",
					"qc.overallStatus=mediocre",
				]
			} & _defaultSubsampleConfig
		}
	}
}

_cutoffs: {
	"6y": "6Y"
	"3y": "3Y"
}

_defaultSubsampleConfig: {
	exclude: "config/outliers_ppx.txt"
	group_by: ["year", "country"]
}

files: {
	auspice_config: "config/auspice_config.json"
	auspice_config_additional_colorings: "config/auspice_config_additional_colorings.json"
}

refine: {
	coalescent: "opt"
	date_inference: "marginal"
	clock_filter_iqd: 4
}

ancestral: {
	inference: "joint"
}

cds: {
	F: "F"
	G: "G"
	genome: "F"
}

traits: {
	columns: "country"
}

frequencies: {
	resolutions: {
		"all-time": {
			min_date: "1975-01-01"
		}
		"6y": {
			min_date: "6Y"
		}
		"3y": {
			min_date: "3Y"
		}
	}
}

nextclade_attributes: {
	a: {
		name: "RSV-A NextClade using real-time tree"
		reference_name: "hRSV/A/England/397/2017"
		accession: "EPI_ISL_412866"
	}
	b: {
		name: "RSV-B NextClade using real-time tree"
		reference_name: "hRSV/B/Australia/VIC-RCH056/2019"
		accession: "EPI_ISL_1653999"
	}
}
