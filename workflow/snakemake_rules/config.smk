"""
This part of the workflow deals with configuration.

OUTPUTS:

    results/run_config.yaml
"""
def get_subtype(build: str) -> str:
    """Extract subtype from build name (e.g. 'a/genome/all-time' -> 'a')."""
    return build.split("/")[0]


def get_gene_build(build: str) -> str:
    """Extract gene build from build name (e.g. 'a/F-antibody-escape/all-time' -> 'F-antibody-escape')."""
    return build.split("/")[1]


def get_gene(build: str) -> str:
    """Extract actual CDS gene from build name (e.g. 'a/F-antibody-escape/all-time' -> 'F')."""
    gene_build = get_gene_build(build)
    if gene_build == "genome":
        return "F"
    else:
        return gene_build.split("-")[0]


def get_resolution(build: str) -> str:
    """Extract resolution from build name (e.g. 'a/genome/all-time' -> 'all-time')."""
    return build.split("/")[2]


write_config("results/run_config.yaml")
