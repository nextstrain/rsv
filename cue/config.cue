package cue

import "strings"

// ---------------------------------------------------------------------------
// Schemas & Definitions
// ---------------------------------------------------------------------------

// Input dataset specification for metadata and sequences
#InputItem: {
	// Unique identifier name for this input dataset
	name: string

	// Path or URL to metadata TSV file (supports {a_or_b} template)
	metadata?: string

	// Path or URL to sequence FASTA file (supports {a_or_b} template)
	sequences?: string
}

// RSV Phylogenetic Workflow Configuration Schema
#Config: {
	// Path to conda environment definition file used by Snakemake
	conda_environment?: string

	// Target builds to run
	builds_to_run?: [...string]

	// Temporal resolutions to build trees for
	resolutions_to_run?: [...string]

	// RSV subtypes to analyze (e.g. "a", "b")
	subtypes?: [...string]

	// Genes to analyze for N-linked glycosylation site changes
	genesforglycosylation?: [...string]

	// Mapping of build names to their corresponding CDS gene name
	cds?: {
		// <build_name>
		[string]: string
	}

	// Default input datasets to fetch and merge
	inputs?: [...#InputItem]

	// Additional user-provided input datasets to merge with defaults
	additional_inputs?: [...#InputItem]

	// Path to file with outlier sequence accessions to exclude
	exclude?: string

	// Path to markdown description file displayed in Auspice
	description?: string

	// Metadata column identifying the unique sequence/strain ID
	strain_id_field?: string

	// Metadata column used for display names on tree tips
	display_strain_field?: string

	// Subsampling configuration. When using --configfile, it is recommended to use 'custom_subsample' instead to ignore default subsampling configuration.
	subsample?: {
		// <build_name>
		[string]: #SubsampleConfigUnaligned
	}

	// Custom subsampling configuration. When using --configfile, this is recommended over 'subsample' to ignore default subsampling configuration.
	custom_subsample?: {
		// <build_name>
		[string]: #SubsampleConfigUnaligned
	}

	// Paths to Auspice configuration JSON files
	files?: {
		// Main Auspice dataset configuration JSON
		auspice_config?: string

		// Additional colorings configuration JSON
		auspice_config_additional_colorings?: string

		// Auspice configuration for F antibody escape DMS defaults
		auspice_config_f_antibody_escape?: string

		// Auspice configuration for non-genome (gene-specific) builds
		"auspice_config_non-genome_builds"?: string
	}

	// TreeTime refine parameters for molecular clock and date inference
	refine?: {
		// Coalescent tree prior ('opt', 'skyline', or float)
		coalescent?: string

		// Date inference mode ('marginal' or 'joint')
		date_inference?: string

		// Number of interquartile distances for clock outlier filtering
		clock_filter_iqd?: number

		// Divergence units on phylogenetic tree ('mutations-per-site' or 'mutations')
		divergence_units?: string
	}

	// Ancestral sequence reconstruction settings
	ancestral?: {
		// Inference mode for ancestral sequences ('joint' or 'marginal')
		inference?: string
	}

	// Metadata traits for ancestral state reconstruction
	traits?: {
		// Column name(s) in metadata to reconstruct on the tree
		columns?: string | [...string]
	}

	// Clade/tip frequency estimation parameters per resolution
	frequencies?: {
		resolutions?: {
			// <resolution>
			[string]: {
				// Minimum/start date for frequency calculation window
				min_date?: string
			}
		}
	}

	// Nextclade dataset reference attributes per subtype
	nextclade_attributes?: {
		// <subtype>
		[string]: {
			// Human-readable Nextclade dataset display name
			name?: string

			// Reference strain identifier
			reference_name?: string

			// Reference sequence accession number
			accession?: string
		}
	}

	// Sequence length and coverage filtering thresholds for antibody escape builds
	filter_for_f_antibody_escape?: {
		// Metadata column name for grouping during escape sequence filtering
		group_by?: string

		// Minimum sequence length per build
		min_length?: {
			// <build_name>
			[string]: int & >0
		}

		// Minimum fraction coverage required per build
		min_coverage?: {
			// <build_name>
			[string]: number & >=0 & <=1
		}

		// Date bounds per temporal resolution
		resolutions?: {
			// <resolution>
			[string]: {
				min_date?:             string
				background_min_date?: string
			}
		}
	}

	// Path to deep mutational scanning (DMS) per-mutation escape CSV
	f_dms_data?: string

	// Monoclonal antibodies to evaluate from DMS dataset
	f_dms_antibodies?: [...string]

	// Whether to clamp negative DMS escape scores to zero
	dms_only_positive_escape?: bool

	// Criteria for enriching phylogenetic trees with high escape mutants
	enrich_antibody_escape?: {
		// <build_name>
		[string]: {
			// Number of high-escape sequences to select per antibody score
			nseqs_per_antibody_scoretype?: int & >0

			// Metadata variables to group by when selecting escape sequences
			group_by?: [...string]

			// Maximum sequences per group with identical F protein mutations
			max_identical_f_prot_muts?: int & >=0

			// Maximum sequences per group with the same top escape mutation
			max_identical_max_escape_mut?: int & >=0
		}
	}

	// Custom Snakemake rule files to include (bypasses schema validation)
	custom_rules?: [...string]
}

// ---------------------------------------------------------------------------
// Build Specification Helper Schema (used internally for subsampling matrix)
// ---------------------------------------------------------------------------

#BuildSpec: {
	coverage_col:    string
	min_coverage:    number & >=0 & <=1
	min_length:      int & >0
	recent_max_seqs: int & >0
}

// ---------------------------------------------------------------------------
// RSV Concrete Configuration & Matrix Generation
// ---------------------------------------------------------------------------

#Subtypes: ["a", "b"]

#Builds: [string]: #BuildSpec
#Builds: {
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

let GroupBy = ["year", "country"]
let ExcludeFile = "config/outliers_ppx.txt"
let RecentExclude = ["qc.overallStatus=bad"]
let BackgroundExclude = ["qc.overallStatus=bad", "qc.overallStatus=mediocre"]

// Enforce that the concrete output conforms to #Config
#Config

conda_environment:     "workflow/envs/nextstrain.yaml"
genesforglycosylation: ["G", "F"]
builds_to_run:         ["genome", "G", "F", "F-antibody-escape"]
resolutions_to_run:    ["all-time", "6y", "3y"]
subtypes:              #Subtypes

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
	for a_or_b in #Subtypes {
		for build_name, build_info in #Builds {
			let recentQuery = "\(build_info.coverage_col)>\(build_info.min_coverage) & missing_data<1000"
			let backgroundQuery = "\(build_info.coverage_col)>\(build_info.min_coverage) & missing_data<1000 & clade.str.startswith(\"\(strings.ToUpper(a_or_b)).D\", na=False)"

			"\(a_or_b)/\(build_name)/all-time": samples: sample: {
				exclude:        ExcludeFile
				exclude_where:  RecentExclude
				group_by:       GroupBy
				max_sequences:  build_info.recent_max_seqs
				min_date:       "1975-01-01"
				min_length:     build_info.min_length
				query:          recentQuery
			}

			"\(a_or_b)/\(build_name)/6y": samples: {
				recent: {
					exclude:        ExcludeFile
					exclude_where:  RecentExclude
					group_by:       GroupBy
					max_sequences:  build_info.recent_max_seqs
					min_date:       "6Y"
					min_length:     build_info.min_length
					query:          recentQuery
				}
				background: {
					include:        "config/include_\(a_or_b).txt"
					exclude:        ExcludeFile
					exclude_where:  BackgroundExclude
					group_by:       GroupBy
					max_sequences:  div(build_info.recent_max_seqs, 10)
					min_date:       "12Y"
					max_date:       "6Y"
					min_length:     build_info.min_length
					query:          backgroundQuery
				}
			}

			"\(a_or_b)/\(build_name)/3y": samples: {
				recent: {
					exclude:        ExcludeFile
					exclude_where:  RecentExclude
					group_by:       GroupBy
					max_sequences:  build_info.recent_max_seqs
					min_date:       "3Y"
					min_length:     build_info.min_length
					query:          recentQuery
				}
				background: {
					include:        "config/include_\(a_or_b).txt"
					exclude:        ExcludeFile
					exclude_where:  BackgroundExclude
					group_by:       GroupBy
					max_sequences:  div(build_info.recent_max_seqs, 10)
					min_date:       "12Y"
					max_date:       "3Y"
					min_length:     build_info.min_length
					query:          backgroundQuery
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
