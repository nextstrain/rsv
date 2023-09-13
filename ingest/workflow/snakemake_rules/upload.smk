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
        cloudfront_domain = config["upload"].get("s3", {}).get("cloudfront_domain", "")
    shell:
        """
        ./vendored/upload-to-s3 --quiet {input:q} {params.s3_dst:q}/{wildcards.remote_file_name:q} {params.cloudfront_domain}
        """
