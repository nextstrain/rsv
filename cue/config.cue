package cue

import (
	"strings"

	augurSubsample "github.com/nextstrain/augur/schemas/subsample"
)


// Helpers (not exported)

let BuildInfo = {
	coverage_col:    string
	min_coverage:    number & >=0 & <=1
	min_length:      int & >0
	recent_max_seqs: int & >0
}

let Builds = {
	[string]: BuildInfo

	genome: {
		coverage_col:    "genome_coverage"
		min_coverage:    0.3
		min_length:      10000
		recent_max_seqs: 3000
	}
	G: {
		coverage_col:    "G_coverage"
		min_coverage:    0.3
		min_length:      600
		recent_max_seqs: 3000
	}
	F: {
		coverage_col:    "F_coverage"
		min_coverage:    0.3
		min_length:      1200
		recent_max_seqs: 3000
	}
	"F-antibody-escape": {
		coverage_col:    "F_coverage"
		min_coverage:    0.75
		min_length:      1200
		recent_max_seqs: 2000
	}
}

let Subtypes = ["a", "b"]
let GroupBy = ["year", "country"]
let ExcludeFile = "config/outliers_ppx.txt"
let RecentExclude = ["qc.overallStatus=bad"]
let BackgroundExclude = ["qc.overallStatus=bad", "qc.overallStatus=mediocre"]


// Reusable references (#defs)

#InputItem: {
	name: string
	metadata?: string
	sequences?: string
} & ({metadata: string} | {sequences: string})


// Schema (for config.schema.yaml)

#Schema: {
	conda_environment?:     string
	genesforglycosylation?: [...string]
	builds_to_run?: [...string]
	resolutions_to_run?: [...string]
	subtypes?: [...string]
	inputs?: [...#InputItem]
	additional_inputs?: [...#InputItem]
	exclude?:              string
	description?:          string
	strain_id_field?:      string
	display_strain_field?: string

	// Subsampling configuration. When using --configfile, it is recommended to use 'custom_subsample' instead to ignore default subsampling configuration.
	subsample?: {
		// <build_name>
		[string]: augurSubsample.#SubsampleConfigUnaligned
	}

	// Custom subsampling configuration. When using --configfile, this is recommended over 'subsample' to ignore default subsampling configuration.
	custom_subsample?: {
		// <build_name>
		[string]: augurSubsample.#SubsampleConfigUnaligned
	}

	files?: {
		auspice_config?:                      string
		auspice_config_additional_colorings?: string
		auspice_config_f_antibody_escape?:    string
		"auspice_config_non-genome_builds"?:  string
	}

	refine?: {
		coalescent?:       string
		date_inference?:   string
		clock_filter_iqd?: number
		divergence_units?: string
	}

	ancestral?: {
		inference?: string
	}

	cds?: {
		// <build_name>
		[string]: string
	}

	traits?: {
		columns?: string | [...string]
	}

	frequencies?: {
		resolutions?: {
			// <resolution>
			[string]: {
				min_date?: string
			}
		}
	}

	nextclade_attributes?: {
		// <subtype>
		[string]: {
			name?: string
			reference_name?: string
			accession?: string
		}
	}

	filter_for_f_antibody_escape?: {
		group_by?: string
		min_length?: {
			// <build_name>
			[string]: int & >0
		}
		min_coverage?: {
			// <build_name>
			[string]: number & >=0 & <=1
		}

		resolutions?: {
			// <resolution>
			[string]: {
				min_date?:            string
				background_min_date?: string
			}
		}
	}

	f_dms_data?: string

	f_dms_antibodies?: [...string]

	dms_only_positive_escape?: bool

	enrich_antibody_escape?: {
		// <build_name>
		[string]: {
			nseqs_per_antibody_scoretype?: int & >0
			group_by?: [...string]
			max_identical_f_prot_muts?: int & >=0
			max_identical_max_escape_mut?: int & >=0
		}
	}

	// Custom Snakemake rule files to include. If used, this will disable config schema validation.
	custom_rules?: [...string]
}

// Default values (for config/configfile.yaml)

// Set the "type" of the default config
#Schema

conda_environment:     "workflow/envs/nextstrain.yaml"
genesforglycosylation: ["G", "F"]
builds_to_run:         ["genome", "G", "F", "F-antibody-escape"]
resolutions_to_run:    ["all-time", "6y", "3y"]
subtypes:              Subtypes

inputs: [
	{
		name:      "ppx_open"
		metadata:  "https://data.nextstrain.org/files/workflows/rsv/{a_or_b}/metadata.tsv.gz"
		sequences: "https://data.nextstrain.org/files/workflows/rsv/{a_or_b}/sequences.fasta.xz"
	},
	{
		name:      "ppx_restricted"
		metadata:  "https://data.nextstrain.org/files/workflows/rsv/{a_or_b}/metadata_restricted.tsv.gz"
		sequences: "https://data.nextstrain.org/files/workflows/rsv/{a_or_b}/sequences_restricted.fasta.xz"
	},
]

exclude:              ExcludeFile
description:          "config/description.md"
strain_id_field:      "accession"
display_strain_field: "strain"

subsample: {
	for a_or_b in Subtypes {
		for build_name, build_info in Builds {
			let recentQuery = "\(build_info.coverage_col)>\(build_info.min_coverage) & missing_data<1000"
			let backgroundQuery = "\(build_info.coverage_col)>\(build_info.min_coverage) & missing_data<1000 & clade.str.startswith(\"\(strings.ToUpper(a_or_b)).D\", na=False)"

			"\(a_or_b)/\(build_name)/all-time": samples: sample: {
				exclude:       ExcludeFile
				exclude_where: RecentExclude
				group_by:      GroupBy
				max_sequences: build_info.recent_max_seqs
				min_date:      "1975-01-01"
				min_length:    build_info.min_length
				query:         recentQuery
			}

			"\(a_or_b)/\(build_name)/6y": samples: {
				recent: {
					exclude:       ExcludeFile
					exclude_where: RecentExclude
					group_by:      GroupBy
					max_sequences: build_info.recent_max_seqs
					min_date:      "6Y"
					min_length:    build_info.min_length
					query:         recentQuery
				}
				background: {
					include:       "config/include_\(a_or_b).txt"
					exclude:       ExcludeFile
					exclude_where: BackgroundExclude
					group_by:      GroupBy
					max_sequences: div(build_info.recent_max_seqs, 10)
					min_date:      "12Y"
					max_date:      "6Y"
					min_length:    build_info.min_length
					query:         backgroundQuery
				}
			}

			"\(a_or_b)/\(build_name)/3y": samples: {
				recent: {
					exclude:       ExcludeFile
					exclude_where: RecentExclude
					group_by:      GroupBy
					max_sequences: build_info.recent_max_seqs
					min_date:      "3Y"
					min_length:    build_info.min_length
					query:         recentQuery
				}
				background: {
					include:       "config/include_\(a_or_b).txt"
					exclude:       ExcludeFile
					exclude_where: BackgroundExclude
					group_by:      GroupBy
					max_sequences: div(build_info.recent_max_seqs, 10)
					min_date:      "12Y"
					max_date:      "3Y"
					min_length:    build_info.min_length
					query:         backgroundQuery
				}
			}
		}
	}
}

files: {
	auspice_config:                      "config/auspice_config.json"
	auspice_config_additional_colorings: "config/auspice_config_additional_colorings.json"
	auspice_config_f_antibody_escape:    "config/auspice_config_dms-defaults.json"
	"auspice_config_non-genome_builds":  "config/auspice_config_non-genome.json"
}

refine: {
	coalescent:       "opt"
	date_inference:   "marginal"
	clock_filter_iqd: 4
	divergence_units: "mutations-per-site"
}

ancestral: inference: "joint"

cds: {
	F:                   "F"
	G:                   "G"
	genome:              "F"
	"F-antibody-escape": "F"
}

traits: columns: "country region"

frequencies: resolutions: {
	"all-time": min_date: "1975-01-01"
	"6y": min_date:       "6Y"
	"3y": min_date:       "3Y"
}

nextclade_attributes: {
	a: {
		name:           "RSV-A NextClade using real-time tree"
		reference_name: "hRSV/A/England/397/2017"
		accession:      "EPI_ISL_412866"
	}
	b: {
		name:           "RSV-B NextClade using real-time tree"
		reference_name: "hRSV/B/Australia/VIC-RCH056/2019"
		accession:      "EPI_ISL_1653999"
	}
}

filter_for_f_antibody_escape: {
	min_length: {
		genome:              10000
		G:                   600
		F:                   1200
		"F-antibody-escape": 1200
	}
	min_coverage: {
		genome:              0.3
		G:                   0.3
		F:                   0.3
		"F-antibody-escape": 0.75
	}
	resolutions: {
		"all-time": min_date: "1975-01-01"
		"6y": min_date:       "6Y"
		"3y": min_date:       "3Y"
	}
}

f_dms_data: "dms-data/all_antibodies.csv"
f_dms_antibodies: [
	"Clesrovimab-Fab",
	"Clesrovimab-IgG",
	"Nirsevimab-Fab",
	"Nirsevimab-IgG",
]
dms_only_positive_escape: true
enrich_antibody_escape: "F-antibody-escape": {
	nseqs_per_antibody_scoretype: 500
	group_by: ["country", "year"]
	max_identical_f_prot_muts:    2
	max_identical_max_escape_mut: 6
}
