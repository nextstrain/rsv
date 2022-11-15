"""
This part of the workflow handles the bells and whistles for automated
Nextstrain builds. Designed to be used internally by the Nextstrain team.

This part of the workflow continues after the main workflow, so the first
rule `deploy` expects input files from `rules.all.input`.

Requires `build_dir` to be defined upstream.
"""
DEPLOY_URL = config["deploy_url"]
SLACK_TS_FILE = build_dir + f"/slack_thread_ts.txt"


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


rule notify_on_deploy:
    input:
        deploy_flag=build_dir + f"/deploy.done",
    output:
        touch(build_dir + f"/notify_on_deploy.done"),
    params:
        deploy_url=DEPLOY_URL,
        slack_ts=SLACK_TS_FILE,
    shell:
        """
        ./bin/notify-on-deploy {params.deploy_url} {params.slack_ts}
        """


onstart:
    # Saves onstart Slack message thread timestamp to file SLACK_TS_FILE
    shell(
        f"./bin/notify-on-start {config.get('build_name', 'unknown')} {SLACK_TS_FILE}"
    )


onsuccess:
    shell(f"./bin/notify-on-success {SLACK_TS_FILE}")


onerror:
    shell(f"./bin/notify-on-error {SLACK_TS_FILE}")
