"""
This part of the workflow deals with configuration.

OUTPUTS:

    results/run_config.yaml
"""
def main():
    write_config("results/run_config.yaml")
    write_subsample_config()


def write_subsample_config():
    # TODO: Support custom build names in the workflow and infer from
    # config["builds"].
    for a_or_b in ["a", "b"]:
        for build_name in ["genome", "G", "F", "F-antibody-escape"]:
            for resolution in ["all-time", "6y", "3y"]:
                build = f"{a_or_b}/{build_name}/{resolution}"
                if "custom_subsample" in config:
                    section = ["custom_subsample", build]
                else:
                    section = ["subsample", build]
                write_config(f"results/{build}/subsample_config.yaml", section=section)


main()
