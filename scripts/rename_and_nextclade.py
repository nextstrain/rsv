import json, argparse

def replace_clade_recursive(node):
    if "genome_clade_annotation" in node["node_attrs"]:
        if "labels" not in node["branch_attrs"]:
            node["branch_attrs"]["labels"] = {}
        node["branch_attrs"]["labels"]["genome_clade"] = node["node_attrs"]["genome_clade_annotation"]["value"]
        node["node_attrs"].pop("genome_clade_annotation")
    if "children" in node:
        for child in node["children"]:
            replace_clade_recursive(child)

if __name__=="__main__":
    parser = argparse.ArgumentParser(
        description="fix genome clade info",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--input-auspice-json', type=str, required=True, help="input auspice_json")
    parser.add_argument('--pathogen-json', type=str, required=True, help="pathogen json")
    parser.add_argument('--reference', type=str, required=True, help="reference")
    parser.add_argument('--build-name', type=str, required=True, help="nextclade build name")
    parser.add_argument('--reference-accession', type=str, required=True, help="reference accession")
    parser.add_argument('--output', type=str, metavar="JSON", required=True, help="output Auspice JSON")
    args = parser.parse_args()

    # read pathogen json
    with open(args.pathogen_json, 'r') as fh:
        pathogen_data = json.load(fh)

    with open(args.input_auspice_json, 'r') as fh:
        data = json.load(fh)

    data["meta"]["colorings"] = [x for x in data["meta"]["colorings"]
                                 if x["key"] != "genome_clade_annotation"]
    replace_clade_recursive(data['tree'])

    # remove unneeded files structure
    pathogen_data.pop("files")

    pathogen_data["attributes"] = {"reference accession": args.reference_accession, "reference name": args.reference, "name": args.build_name}
    pathogen_data["experimental"] = True
    data["meta"]["extensions"] = {'nextclade': {'pathogen': pathogen_data}}

    with open(args.output, 'w') as fh:
        json.dump(data, fh, indent=0)
