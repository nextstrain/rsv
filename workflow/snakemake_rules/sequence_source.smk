"""
This part of the workflow adds a `sequence_source` column to the metadata to
distinguish sequences from different input sources.

If `additional_inputs` is defined in the config, sequences from those inputs
are labeled with their input name, and all other sequences are labeled as
"Pathoplexus". If no additional inputs are defined, the metadata is passed
through unchanged.

Additional inputs with `keep_all: True` bypass subsampling (force-included
via augur filter --include). All other additional inputs are subsampled
normally alongside the default inputs.

INPUTS:
    metadata_merged = results/{a_or_b}/metadata_merged.tsv

OUTPUTS:
    metadata           = results/{a_or_b}/metadata.tsv
    additional_include = results/{a_or_b}/additional_include.txt
"""

_additional_inputs = config.get("additional_inputs", [])

# Validate that no additional input uses the reserved name "Pathoplexus"
for _ai in _additional_inputs:
    if _ai["name"].lower() == "pathoplexus":
        raise ValueError(
            f"Additional input name '{_ai['name']}' conflicts with the reserved "
            "name used for background sequences. Please choose a different name."
        )

# Validate keep_all is boolean if present
for _ai in _additional_inputs:
    if "keep_all" in _ai and not isinstance(_ai["keep_all"], bool):
        raise ValueError(
            f"Additional input '{_ai['name']}' has keep_all={_ai['keep_all']!r} "
            "but it must be True or False (a boolean, not a string)."
        )


if _additional_inputs:

    def _get_additional_sequence_files(wildcards):
        """Get all sequence files from additional_inputs for this subtype."""
        files = []
        for ai in _additional_inputs:
            if "sequences" in ai:
                seq_path = ai["sequences"].replace("{a_or_b}", wildcards.a_or_b)
                files.append(seq_path)
        return files

    def _get_keep_all_sequence_files(wildcards):
        """Get sequence files from additional_inputs with keep_all: True."""
        files = []
        for ai in _additional_inputs:
            if ai.get("keep_all", False) and "sequences" in ai:
                seq_path = ai["sequences"].replace("{a_or_b}", wildcards.a_or_b)
                files.append(seq_path)
        return files

    rule add_sequence_source:
        """
        Add sequence_source column to metadata based on additional_inputs.
        Sequences found in additional_inputs FASTAs are labeled with their
        input name; all others are labeled 'Pathoplexus'.
        """
        input:
            metadata="results/{a_or_b}/metadata_merged.tsv",
            additional_sequences=_get_additional_sequence_files,
        output:
            metadata="results/{a_or_b}/metadata.tsv",
        log:
            "logs/add_sequence_source_{a_or_b}.txt",
        benchmark:
            "benchmarks/add_sequence_source_{a_or_b}.txt"
        run:
            import csv
            from Bio import SeqIO

            with open(log[0], "w") as log_file:
                # Build mapping: accession -> source name from additional FASTAs
                accession_to_source = {}
                for ai in _additional_inputs:
                    if "sequences" not in ai:
                        continue
                    seq_path = ai["sequences"].replace("{a_or_b}", wildcards.a_or_b)
                    count = 0
                    for record in SeqIO.parse(seq_path, "fasta"):
                        accession_to_source[record.id] = ai["name"]
                        count += 1
                    print(
                        f"Found {count} sequences from additional input '{ai['name']}' "
                        f"in {seq_path}",
                        file=log_file,
                    )

                print(
                    f"Total additional accessions: {len(accession_to_source)}",
                    file=log_file,
                )

                # Read merged metadata, add sequence_source column, write out
                with open(input.metadata) as f_in, \
                     open(output.metadata, "w", newline="") as f_out:
                    reader = csv.DictReader(f_in, delimiter="\t")
                    fieldnames = list(reader.fieldnames)
                    if "sequence_source" not in fieldnames:
                        fieldnames.append("sequence_source")
                    writer = csv.DictWriter(
                        f_out,
                        fieldnames=fieldnames,
                        delimiter="\t",
                        lineterminator="\n",
                    )
                    writer.writeheader()
                    n_additional = 0
                    n_pathoplexus = 0
                    for row in reader:
                        accession = row.get("accession", "")
                        source = accession_to_source.get(accession, "Pathoplexus")
                        row["sequence_source"] = source
                        if source != "Pathoplexus":
                            n_additional += 1
                        else:
                            n_pathoplexus += 1
                        writer.writerow(row)

                print(
                    f"Labeled {n_additional} additional and {n_pathoplexus} Pathoplexus sequences",
                    file=log_file,
                )


    rule generate_additional_include_list:
        """
        Generate a list of accessions from additional_inputs with keep_all: True
        to force-include in subsampling filters.
        """
        input:
            keep_all_sequences=_get_keep_all_sequence_files,
        output:
            include_list="results/{a_or_b}/additional_include.txt",
        log:
            "logs/generate_additional_include_list_{a_or_b}.txt",
        benchmark:
            "benchmarks/generate_additional_include_list_{a_or_b}.txt"
        run:
            from Bio import SeqIO

            accessions = []
            for seq_file in input.keep_all_sequences:
                for record in SeqIO.parse(seq_file, "fasta"):
                    accessions.append(record.id)

            with open(output.include_list, "w") as f:
                for acc in accessions:
                    f.write(acc + "\n")

            with open(log[0], "w") as log_file:
                print(
                    f"Wrote {len(accessions)} keep_all accessions to {output.include_list}",
                    file=log_file,
                )


else:

    rule passthrough_metadata:
        """
        No additional inputs defined; pass metadata through unchanged.
        """
        input:
            metadata="results/{a_or_b}/metadata_merged.tsv",
        output:
            metadata="results/{a_or_b}/metadata.tsv",
        log:
            "logs/passthrough_metadata_{a_or_b}.txt",
        benchmark:
            "benchmarks/passthrough_metadata_{a_or_b}.txt"
        shell:
            r"""
            exec &> >(tee {log:q})

            cp {input.metadata} {output.metadata}
            """


    rule generate_empty_additional_include_list:
        """
        No additional inputs; create empty include list (no-op for augur filter).
        """
        output:
            include_list="results/{a_or_b}/additional_include.txt",
        log:
            "logs/generate_empty_additional_include_list_{a_or_b}.txt",
        benchmark:
            "benchmarks/generate_empty_additional_include_list_{a_or_b}.txt"
        shell:
            r"""
            exec &> >(tee {log:q})

            touch {output.include_list}
            """
