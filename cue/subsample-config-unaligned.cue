// Configuration file to be supplied to `augur subsample --config` (intended for
// inputs where sequences are unaligned)
package cue

import "struct"

@jsonschema(schema="http://json-schema.org/draft-07/schema#")
@jsonschema(id="https://nextstrain.org/schemas/augur/subsample-config-unaligned/v1")

#SubsampleConfigUnaligned: close({
	// Default properties applied across all sample partitions
	defaults?: #defaultProperties

	// Map of sample partition names (e.g. 'sample', 'recent', 'background') to their subsampling parameters
	samples!: struct.MinFields(1) & {
		// <sample_name>
		[string]: #filterSampleProperties
		...
	}
})

#defaultProperties: close({
	// File(s) with list of strains to exclude.
	exclude?: string | [...string]

	// Exclude all strains by default. Use this with the include arguments to
	// select a specific subset of strains.
	exclude_all?: bool

	// Exclude ambiguous dates by day (e.g., 2020-09-XX), month (e.g.,
	// 2020-XX-XX), year (e.g., 200X-10-01), or any date fields.
	exclude_ambiguous_dates_by?: "any" | "day" | "month" | "year"

	// Exclude sequences matching these conditions. Ex: "host=rat" or
	// "host!=rat". Multiple values are processed as OR (matching any of those
	// specified will be excluded), not AND.
	exclude_where?: string | [...string]

	// File(s) with list of strains to include regardless of priorities,
	// subsampling, or absence of an entry in sequences.
	include?: string | [...string]

	// Include sequences with these values. ex: host=rat. Multiple values are
	// processed as OR (having any of those specified will be included), not
	// AND.
	include_where?: string | [...string]

	// Minimal cutoff for date (inclusive).
	min_date?: int | string

	// Maximal cutoff for date (inclusive).
	max_date?: int | string

	// Minimal length of the sequences, only counting valid characters.
	min_length?: int

	// Maximum length of the sequences, only counting valid characters.
	max_length?: int

	// Exclude sequences that contain invalid characters.
	exclude_invalid?: bool

	// Deprecated, please use 'exclude_invalid' instead.
	non_nucleotide?: bool

	// Filter sequences by attribute. Uses Pandas DataFrame query syntax.
	query?: string

	// Use alongside query to specify columns and data types in the format 'column:type'.
	query_columns?: string | [...string]
})

#filterSampleProperties: close({
	// File(s) with list of strains to exclude.
	exclude?: string | [...string]

	// Exclude all strains by default. Use this with the include arguments to
	// select a specific subset of strains.
	exclude_all?: bool

	// Exclude ambiguous dates by day (e.g., 2020-09-XX), month (e.g.,
	// 2020-XX-XX), year (e.g., 200X-10-01), or any date fields.
	exclude_ambiguous_dates_by?: "any" | "day" | "month" | "year"

	// Exclude sequences matching these conditions. Ex: "host=rat" or
	// "host!=rat". Multiple values are processed as OR (matching any of those
	// specified will be excluded), not AND.
	exclude_where?: string | [...string]

	// File(s) with list of strains to include regardless of priorities,
	// subsampling, or absence of an entry in sequences.
	include?: string | [...string]

	// Include sequences with these values. ex: host=rat. Multiple values are
	// processed as OR (having any of those specified will be included), not
	// AND.
	include_where?: string | [...string]

	// Minimal cutoff for date (inclusive).
	min_date?: int | string

	// Maximal cutoff for date (inclusive).
	max_date?: int | string

	// Minimal length of the sequences, only counting valid characters.
	min_length?: int

	// Maximum length of the sequences, only counting valid characters.
	max_length?: int

	// Exclude sequences that contain invalid characters.
	exclude_invalid?: bool

	// Deprecated, please use 'exclude_invalid' instead.
	non_nucleotide?: bool

	// Filter sequences by attribute. Uses Pandas DataFrame query syntax.
	query?: string

	// Use alongside query to specify columns and data types in the format 'column:type'.
	query_columns?: string | [...string]

	// Use the outputs from another sample as the inputs for this sample. Value must be a sample name.
	context_sample?: string

	// Drop this sample from the final output
	drop_sample?: bool

	// Grouping columns for subsampling.
	group_by?: string | [...string]

	// TSV file defining weights for grouping.
	group_by_weights?: string

	// Allow probabilistic sampling during subsampling.
	probabilistic_sampling?: bool

	// Select no more than this number of sequences per category.
	sequences_per_group?: int

	// Select no more than this number of sequences (i.e. total sample size).
	max_sequences?: int
})
