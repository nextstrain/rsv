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
    parser.add_argument('--output', type=str, metavar="JSON", required=True, help="output Auspice JSON")
    args = parser.parse_args()

    with open(args.input_auspice_json, 'r') as fh:
        data = json.load(fh)

    data["meta"]["colorings"] = [x for x in data["meta"]["colorings"]
                                 if x["key"] != "genome_clade_annotation"]
    replace_clade_recursive(data['tree'])

    with open(args.output, 'w') as fh:
        json.dump(data, fh, indent=0)
