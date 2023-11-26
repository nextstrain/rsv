import json, argparse

def get_clade_configs(name):
    return {
    "G_clade": {
        "name": "G_clade",
        "displayName": "G clades",
        "description": "Legacy clade system by Goya et al, superseeded by the new Consortium Nomenclature."
    },
    }.get(name, {'name':name, "displayName":name, "description":""})


if __name__=="__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--lineage", required=True, help="attribute info")
    parser.add_argument("--reference", required=True, help="attribute info")
    parser.add_argument("--reference-name", required=True, help="attribute info")
    parser.add_argument("--auspice-config", required=True, help="Auspice config JSON with coloring entry to have scale added to")
    parser.add_argument("--pathogen-jsons", nargs='+', required=True, help="name of the coloring field in the Auspice config JSON")
    parser.add_argument("--clades", nargs="+", required=True, help="list of values to assign colors to")
    parser.add_argument("--output-auspice", required=True, help="Auspice config JSON with scale added to the requested coloring")
    parser.add_argument("--output-pathogen", required=True, help="Auspice config JSON with scale added to the requested coloring")
    args = parser.parse_args()

    pathogen_json = {}
    for p in args.pathogen_jsons:
        with open(p) as fh:
            pathogen_json.update(json.load(fh))

    with open(args.auspice_config) as fh:
        auspice_json = json.load(fh)

    RSV_type = {'a':'RSV-A', 'b':'RSV-B'}[args.lineage]

    pathogen_json['attributes'] = {"name": f"{RSV_type}",
                                   "reference accession": args.reference,
                                   "reference name": args.reference_name}

    pathogen_json['geneOrderPreference'] = ["F", "G", "L"]

    if len(args.clades):
        auspice_json['extensions']['nextclade']["clade_node_attrs"] =  [
            get_clade_configs(c) for c in args.clades if c!='default'
        ]

    with open(args.output_pathogen, 'w') as fh:
        json.dump(pathogen_json, fh, indent=2)

    with open(args.output_auspice, 'w') as fh:
        json.dump(auspice_json, fh, indent=2)

