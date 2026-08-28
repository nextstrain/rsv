<!-- [DO NOT EDIT] This file was generated automatically from cue/ definitions. -->
# Configuration Reference

This reference is automatically generated from [`cue/`](../cue/).

## Top-Level Settings

RSV Phylogenetic Workflow Configuration Schema

| Parameter Path | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `builds_to_run` | `list[string]` | `["genome", "G", "F", "F-antibody-escape"]` | Target builds to run |
| `cds` | `object` | `{"F": "F", "G": "G", "genome": "F", "F-antibody-escape": "F"}` | Mapping of build names to their corresponding CDS gene name |
| `conda_environment` | `string` | `workflow/envs/nextstrain.yaml` | Path to conda environment definition file used by Snakemake |
| `custom_rules` | `list[string]` | - | Custom Snakemake rule files to include (bypasses schema validation) |
| `description` | `string` | `config/description.md` | Path to markdown description file displayed in Auspice |
| `display_strain_field` | `string` | `strain` | Metadata column used for display names on tree tips |
| `dms_only_positive_escape` | `boolean` | `true` | Whether to clamp negative DMS escape scores to zero |
| `exclude` | `string` | `config/outliers_ppx.txt` | Path to file with outlier sequence accessions to exclude |
| `f_dms_antibodies` | `list[string]` | `["Clesrovimab-Fab", "Clesrovimab-IgG", "Nirsevimab-Fab", "Nirsevimab-IgG"]` | Monoclonal antibodies to evaluate from DMS dataset |
| `f_dms_data` | `string` | `dms-data/all_antibodies.csv` | Path to deep mutational scanning (DMS) per-mutation escape CSV |
| `genesforglycosylation` | `list[string]` | `["G", "F"]` | Genes to analyze for N-linked glycosylation site changes |
| `resolutions_to_run` | `list[string]` | `["all-time", "6y", "3y"]` | Temporal resolutions to build trees for |
| `strain_id_field` | `string` | `accession` | Metadata column identifying the unique sequence/strain ID |
| `subtypes` | `list[string]` | `["a", "b"]` | RSV subtypes to analyze (e.g. "a", "b") |

### `additional_inputs[]`

Additional user-provided input datasets to merge with defaults

| Parameter Path | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `additional_inputs[].metadata` | `string` | - | Path or URL to metadata TSV file (supports {a_or_b} template) |
| `additional_inputs[].name` | `string` | - | Unique identifier name for this input dataset |
| `additional_inputs[].sequences` | `string` | - | Path or URL to sequence FASTA file (supports {a_or_b} template) |

### `ancestral`

Ancestral sequence reconstruction settings

| Parameter Path | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `ancestral.inference` | `string` | `joint` | Inference mode for ancestral sequences ('joint' or 'marginal') |

#### `custom_subsample.<build_name>.defaults`

Default properties applied across all sample partitions

| Parameter Path | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `custom_subsample.<build_name>.defaults.exclude` | `string \| list[string]` | - | File(s) with list of strains to exclude. |
| `custom_subsample.<build_name>.defaults.exclude_all` | `boolean` | - | Exclude all strains by default. Use this with the include arguments to select a specific subset of strains. |
| `custom_subsample.<build_name>.defaults.exclude_ambiguous_dates_by` | `"any" \| "day" \| "month" \| "year"` | - | Exclude ambiguous dates by day (e.g., 2020-09-XX), month (e.g., 2020-XX-XX), year (e.g., 200X-10-01), or any date fields. |
| `custom_subsample.<build_name>.defaults.exclude_invalid` | `boolean` | - | Exclude sequences that contain invalid characters. |
| `custom_subsample.<build_name>.defaults.exclude_where` | `string \| list[string]` | - | Exclude sequences matching these conditions. Ex: "host=rat" or "host!=rat". Multiple values are processed as OR (matching any of those specified will be excluded), not AND. |
| `custom_subsample.<build_name>.defaults.include` | `string \| list[string]` | - | File(s) with list of strains to include regardless of priorities, subsampling, or absence of an entry in sequences. |
| `custom_subsample.<build_name>.defaults.include_where` | `string \| list[string]` | - | Include sequences with these values. ex: host=rat. Multiple values are processed as OR (having any of those specified will be included), not AND. |
| `custom_subsample.<build_name>.defaults.max_date` | `integer \| string` | - | Maximal cutoff for date (inclusive). |
| `custom_subsample.<build_name>.defaults.max_length` | `integer` | - | Maximum length of the sequences, only counting valid characters. |
| `custom_subsample.<build_name>.defaults.min_date` | `integer \| string` | - | Minimal cutoff for date (inclusive). |
| `custom_subsample.<build_name>.defaults.min_length` | `integer` | - | Minimal length of the sequences, only counting valid characters. |
| `custom_subsample.<build_name>.defaults.non_nucleotide` | `boolean` | - | Deprecated, please use 'exclude_invalid' instead. |
| `custom_subsample.<build_name>.defaults.query` | `string` | - | Filter sequences by attribute. Uses Pandas DataFrame query syntax. |
| `custom_subsample.<build_name>.defaults.query_columns` | `string \| list[string]` | - | Use alongside query to specify columns and data types in the format 'column:type'. |

#### `custom_subsample.<build_name>.samples.<sample_name>`

Map of sample partition names (e.g. 'sample', 'recent', 'background') to their subsampling parameters

| Parameter Path | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `custom_subsample.<build_name>.samples.<sample_name>.context_sample` | `string` | - | Use the outputs from another sample as the inputs for this sample. Value must be a sample name. |
| `custom_subsample.<build_name>.samples.<sample_name>.drop_sample` | `boolean` | - | Drop this sample from the final output |
| `custom_subsample.<build_name>.samples.<sample_name>.exclude` | `string \| list[string]` | - | File(s) with list of strains to exclude. |
| `custom_subsample.<build_name>.samples.<sample_name>.exclude_all` | `boolean` | - | Exclude all strains by default. Use this with the include arguments to select a specific subset of strains. |
| `custom_subsample.<build_name>.samples.<sample_name>.exclude_ambiguous_dates_by` | `"any" \| "day" \| "month" \| "year"` | - | Exclude ambiguous dates by day (e.g., 2020-09-XX), month (e.g., 2020-XX-XX), year (e.g., 200X-10-01), or any date fields. |
| `custom_subsample.<build_name>.samples.<sample_name>.exclude_invalid` | `boolean` | - | Exclude sequences that contain invalid characters. |
| `custom_subsample.<build_name>.samples.<sample_name>.exclude_where` | `string \| list[string]` | - | Exclude sequences matching these conditions. Ex: "host=rat" or "host!=rat". Multiple values are processed as OR (matching any of those specified will be excluded), not AND. |
| `custom_subsample.<build_name>.samples.<sample_name>.group_by` | `string \| list[string]` | - | Grouping columns for subsampling. |
| `custom_subsample.<build_name>.samples.<sample_name>.group_by_weights` | `string` | - | TSV file defining weights for grouping. |
| `custom_subsample.<build_name>.samples.<sample_name>.include` | `string \| list[string]` | - | File(s) with list of strains to include regardless of priorities, subsampling, or absence of an entry in sequences. |
| `custom_subsample.<build_name>.samples.<sample_name>.include_where` | `string \| list[string]` | - | Include sequences with these values. ex: host=rat. Multiple values are processed as OR (having any of those specified will be included), not AND. |
| `custom_subsample.<build_name>.samples.<sample_name>.max_date` | `integer \| string` | - | Maximal cutoff for date (inclusive). |
| `custom_subsample.<build_name>.samples.<sample_name>.max_length` | `integer` | - | Maximum length of the sequences, only counting valid characters. |
| `custom_subsample.<build_name>.samples.<sample_name>.max_sequences` | `integer` | - | Select no more than this number of sequences (i.e. total sample size). |
| `custom_subsample.<build_name>.samples.<sample_name>.min_date` | `integer \| string` | - | Minimal cutoff for date (inclusive). |
| `custom_subsample.<build_name>.samples.<sample_name>.min_length` | `integer` | - | Minimal length of the sequences, only counting valid characters. |
| `custom_subsample.<build_name>.samples.<sample_name>.non_nucleotide` | `boolean` | - | Deprecated, please use 'exclude_invalid' instead. |
| `custom_subsample.<build_name>.samples.<sample_name>.probabilistic_sampling` | `boolean` | - | Allow probabilistic sampling during subsampling. |
| `custom_subsample.<build_name>.samples.<sample_name>.query` | `string` | - | Filter sequences by attribute. Uses Pandas DataFrame query syntax. |
| `custom_subsample.<build_name>.samples.<sample_name>.query_columns` | `string \| list[string]` | - | Use alongside query to specify columns and data types in the format 'column:type'. |
| `custom_subsample.<build_name>.samples.<sample_name>.sequences_per_group` | `integer` | - | Select no more than this number of sequences per category. |

### `enrich_antibody_escape.<build_name>`

Criteria for enriching phylogenetic trees with high escape mutants

| Parameter Path | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `enrich_antibody_escape.<build_name>.group_by` | `list[string]` | - | Metadata variables to group by when selecting escape sequences |
| `enrich_antibody_escape.<build_name>.max_identical_f_prot_muts` | `object` | - | Maximum sequences per group with identical F protein mutations |
| `enrich_antibody_escape.<build_name>.max_identical_max_escape_mut` | `object` | - | Maximum sequences per group with the same top escape mutation |
| `enrich_antibody_escape.<build_name>.nseqs_per_antibody_scoretype` | `object` | - | Number of high-escape sequences to select per antibody score |

### `files`

Paths to Auspice configuration JSON files

| Parameter Path | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `files.auspice_config` | `string` | `config/auspice_config.json` | Main Auspice dataset configuration JSON |
| `files.auspice_config_additional_colorings` | `string` | `config/auspice_config_additional_colorings.json` | Additional colorings configuration JSON |
| `files.auspice_config_f_antibody_escape` | `string` | `config/auspice_config_dms-defaults.json` | Auspice configuration for F antibody escape DMS defaults |
| `files.auspice_config_non-genome_builds` | `string` | `config/auspice_config_non-genome.json` | Auspice configuration for non-genome (gene-specific) builds |

### `filter_for_f_antibody_escape`

Sequence length and coverage filtering thresholds for antibody escape builds

| Parameter Path | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `filter_for_f_antibody_escape.group_by` | `string` | - | Metadata column name for grouping during escape sequence filtering |
| `filter_for_f_antibody_escape.min_coverage` | `object` | `{"genome": 0.3, "G": 0.3, "F": 0.3, "F-antibody-escape": 0.75}` | Minimum fraction coverage required per build |
| `filter_for_f_antibody_escape.min_length` | `object` | `{"genome": 10000, "G": 600, "F": 1200, "F-antibody-escape": 1200}` | Minimum sequence length per build |

#### `filter_for_f_antibody_escape.resolutions.<resolution>`

Date bounds per temporal resolution

| Parameter Path | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `filter_for_f_antibody_escape.resolutions.<resolution>.background_min_date` | `string` | - | - |
| `filter_for_f_antibody_escape.resolutions.<resolution>.min_date` | `string` | - | - |

#### `frequencies.resolutions.<resolution>`

<resolution>

| Parameter Path | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `frequencies.resolutions.<resolution>.min_date` | `string` | - | Minimum/start date for frequency calculation window |

### `inputs[]`

Default input datasets to fetch and merge

| Parameter Path | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `inputs[].metadata` | `string` | - | Path or URL to metadata TSV file (supports {a_or_b} template) |
| `inputs[].name` | `string` | - | Unique identifier name for this input dataset |
| `inputs[].sequences` | `string` | - | Path or URL to sequence FASTA file (supports {a_or_b} template) |

### `nextclade_attributes.<subtype>`

Nextclade dataset reference attributes per subtype

| Parameter Path | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `nextclade_attributes.<subtype>.accession` | `string` | - | Reference sequence accession number |
| `nextclade_attributes.<subtype>.name` | `string` | - | Human-readable Nextclade dataset display name |
| `nextclade_attributes.<subtype>.reference_name` | `string` | - | Reference strain identifier |

### `refine`

TreeTime refine parameters for molecular clock and date inference

| Parameter Path | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `refine.clock_filter_iqd` | `number` | `4` | Number of interquartile distances for clock outlier filtering |
| `refine.coalescent` | `string` | `opt` | Coalescent tree prior ('opt', 'skyline', or float) |
| `refine.date_inference` | `string` | `marginal` | Date inference mode ('marginal' or 'joint') |
| `refine.divergence_units` | `string` | `mutations-per-site` | Divergence units on phylogenetic tree ('mutations-per-site' or 'mutations') |

#### `subsample.<build_name>.defaults`

Default properties applied across all sample partitions

| Parameter Path | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `subsample.<build_name>.defaults.exclude` | `string \| list[string]` | - | File(s) with list of strains to exclude. |
| `subsample.<build_name>.defaults.exclude_all` | `boolean` | - | Exclude all strains by default. Use this with the include arguments to select a specific subset of strains. |
| `subsample.<build_name>.defaults.exclude_ambiguous_dates_by` | `"any" \| "day" \| "month" \| "year"` | - | Exclude ambiguous dates by day (e.g., 2020-09-XX), month (e.g., 2020-XX-XX), year (e.g., 200X-10-01), or any date fields. |
| `subsample.<build_name>.defaults.exclude_invalid` | `boolean` | - | Exclude sequences that contain invalid characters. |
| `subsample.<build_name>.defaults.exclude_where` | `string \| list[string]` | - | Exclude sequences matching these conditions. Ex: "host=rat" or "host!=rat". Multiple values are processed as OR (matching any of those specified will be excluded), not AND. |
| `subsample.<build_name>.defaults.include` | `string \| list[string]` | - | File(s) with list of strains to include regardless of priorities, subsampling, or absence of an entry in sequences. |
| `subsample.<build_name>.defaults.include_where` | `string \| list[string]` | - | Include sequences with these values. ex: host=rat. Multiple values are processed as OR (having any of those specified will be included), not AND. |
| `subsample.<build_name>.defaults.max_date` | `integer \| string` | - | Maximal cutoff for date (inclusive). |
| `subsample.<build_name>.defaults.max_length` | `integer` | - | Maximum length of the sequences, only counting valid characters. |
| `subsample.<build_name>.defaults.min_date` | `integer \| string` | - | Minimal cutoff for date (inclusive). |
| `subsample.<build_name>.defaults.min_length` | `integer` | - | Minimal length of the sequences, only counting valid characters. |
| `subsample.<build_name>.defaults.non_nucleotide` | `boolean` | - | Deprecated, please use 'exclude_invalid' instead. |
| `subsample.<build_name>.defaults.query` | `string` | - | Filter sequences by attribute. Uses Pandas DataFrame query syntax. |
| `subsample.<build_name>.defaults.query_columns` | `string \| list[string]` | - | Use alongside query to specify columns and data types in the format 'column:type'. |

#### `subsample.<build_name>.samples.<sample_name>`

Map of sample partition names (e.g. 'sample', 'recent', 'background') to their subsampling parameters

| Parameter Path | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `subsample.<build_name>.samples.<sample_name>.context_sample` | `string` | - | Use the outputs from another sample as the inputs for this sample. Value must be a sample name. |
| `subsample.<build_name>.samples.<sample_name>.drop_sample` | `boolean` | - | Drop this sample from the final output |
| `subsample.<build_name>.samples.<sample_name>.exclude` | `string \| list[string]` | - | File(s) with list of strains to exclude. |
| `subsample.<build_name>.samples.<sample_name>.exclude_all` | `boolean` | - | Exclude all strains by default. Use this with the include arguments to select a specific subset of strains. |
| `subsample.<build_name>.samples.<sample_name>.exclude_ambiguous_dates_by` | `"any" \| "day" \| "month" \| "year"` | - | Exclude ambiguous dates by day (e.g., 2020-09-XX), month (e.g., 2020-XX-XX), year (e.g., 200X-10-01), or any date fields. |
| `subsample.<build_name>.samples.<sample_name>.exclude_invalid` | `boolean` | - | Exclude sequences that contain invalid characters. |
| `subsample.<build_name>.samples.<sample_name>.exclude_where` | `string \| list[string]` | - | Exclude sequences matching these conditions. Ex: "host=rat" or "host!=rat". Multiple values are processed as OR (matching any of those specified will be excluded), not AND. |
| `subsample.<build_name>.samples.<sample_name>.group_by` | `string \| list[string]` | - | Grouping columns for subsampling. |
| `subsample.<build_name>.samples.<sample_name>.group_by_weights` | `string` | - | TSV file defining weights for grouping. |
| `subsample.<build_name>.samples.<sample_name>.include` | `string \| list[string]` | - | File(s) with list of strains to include regardless of priorities, subsampling, or absence of an entry in sequences. |
| `subsample.<build_name>.samples.<sample_name>.include_where` | `string \| list[string]` | - | Include sequences with these values. ex: host=rat. Multiple values are processed as OR (having any of those specified will be included), not AND. |
| `subsample.<build_name>.samples.<sample_name>.max_date` | `integer \| string` | - | Maximal cutoff for date (inclusive). |
| `subsample.<build_name>.samples.<sample_name>.max_length` | `integer` | - | Maximum length of the sequences, only counting valid characters. |
| `subsample.<build_name>.samples.<sample_name>.max_sequences` | `integer` | - | Select no more than this number of sequences (i.e. total sample size). |
| `subsample.<build_name>.samples.<sample_name>.min_date` | `integer \| string` | - | Minimal cutoff for date (inclusive). |
| `subsample.<build_name>.samples.<sample_name>.min_length` | `integer` | - | Minimal length of the sequences, only counting valid characters. |
| `subsample.<build_name>.samples.<sample_name>.non_nucleotide` | `boolean` | - | Deprecated, please use 'exclude_invalid' instead. |
| `subsample.<build_name>.samples.<sample_name>.probabilistic_sampling` | `boolean` | - | Allow probabilistic sampling during subsampling. |
| `subsample.<build_name>.samples.<sample_name>.query` | `string` | - | Filter sequences by attribute. Uses Pandas DataFrame query syntax. |
| `subsample.<build_name>.samples.<sample_name>.query_columns` | `string \| list[string]` | - | Use alongside query to specify columns and data types in the format 'column:type'. |
| `subsample.<build_name>.samples.<sample_name>.sequences_per_group` | `integer` | - | Select no more than this number of sequences per category. |

### `traits`

Metadata traits for ancestral state reconstruction

| Parameter Path | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `traits.columns` | `string \| list[string]` | `country region` | Column name(s) in metadata to reconstruct on the tree |

