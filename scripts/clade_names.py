import json, argparse

def replace_clade_recursive(node):
    if 'genome_clade' in node:
        node["branch_attrs"]["labels"]= node["node_attrs"]["genome_clade_annotation"]
        node['node_attrs'].pop('genome_clade')
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

    replace_clade_recursive(data['tree'])

    with open(args.output, 'w') as fh:
        json.dump(data, fh, indent=2)
