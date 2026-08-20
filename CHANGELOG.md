# CHANGELOG

We use this CHANGELOG to document breaking changes, new features, bug fixes, and config value changes that may affect both the usage of the workflows and the outputs of the workflows.

## 2026

* TBD: The `filter` section in phylogenetic workflow configuration has been replaced by `subsample`/`custom_subsample` for subsampling, and `filter_for_f_antibody_escape` for initial quality filtering. **This is a breaking change**.
    * NOTE: The workflow does not yet support proximal samples.
* 11 August 2026: Phylogenetic workflow configuration is now validated against a strict schema. The workflow will error if your configuration has extraneous entries that were previously ignored.
