"""
This part of the workflow handles the bells and whistles for automated
Nextstrain builds. Designed to be used internally by the Nextstrain team.

This part of the workflow continues after the main workflow, so the first
rule `deploy` expects input files from `rules.all.input`.

Requires `build_dir` to be defined upstream.
"""
DEPLOY_URL = config["deploy_url"]


rule deploy:
    input:
        *rules.all.input,
    output:
        touch(build_dir + f"/deploy.done"),
    params:
        deploy_url=DEPLOY_URL,
    shell:
        """
        nextstrain remote upload {params.deploy_url} {input}
        """
