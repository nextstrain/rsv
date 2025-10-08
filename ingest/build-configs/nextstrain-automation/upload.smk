"""
This part of the workflow handles uploading files to a specified destination.

Uses predefined wildcard `file_to_upload` determine input and predefined
wildcard `remote_file_name` as the remote file name in the specified destination.

Produces output files as `data/upload/{upload_target_name}/{file_to_upload}-to-{remote_file_name}.done`.

Currently only supports uploads to AWS S3, but additional upload rules can
be easily added as long as they follow the output pattern described above.
"""
import os

rule upload_to_s3:
    input:
        file_to_upload = "data/{file_to_upload}"
    output:
        touch("data/upload/s3/{file_to_upload}-to-{remote_file_name}.done")
    params:
        s3_dst = config["s3_dst"],
        cloudfront_domain = config["upload"].get("s3", {}).get("cloudfront_domain", ""),
        current_basedir = str(workflow.current_basedir),
    shell:
        """
        {params.current_basedir}/../../../shared/vendored/scripts/upload-to-s3 \
            --quiet \
            {input:q} \
            {params.s3_dst:q}/{wildcards.remote_file_name:q} \
            {params.cloudfront_domain}
        """



def _get_all_targets(wildcards):
    # Default targets are the metadata TSV and sequences FASTA files
    all_targets = []

    # Add additional targets based on upload config
    upload_config = config.get("upload", {})

    for target, params in upload_config.items():
        files_to_upload = params.get("files_to_upload", [])
        remote_file_names = params.get("remote_file_names", [])

        if len(files_to_upload) != len(remote_file_names):
            print(
                f"Skipping file upload for {target!r} because the number of",
                "files to upload does not match the number of remote file names."
            )
        elif len(remote_file_names) != len(set(remote_file_names)):
            print(f"Skipping file upload for {target!r} because there are duplicate remote file names.")
        elif "s3_dst" not in config:
            print(f"Skipping file upload for {target!r} because the destintion was not defined.")
        else:
            all_targets.extend(
                expand(
                    [f"data/upload/{target}/{{file_to_upload}}-to-{{remote_file_name}}.done"],
                    zip,
                    file_to_upload=files_to_upload,
                    remote_file_name=remote_file_names
                )
            )
    return all_targets


rule upload_all:
    input:
        _get_all_targets
