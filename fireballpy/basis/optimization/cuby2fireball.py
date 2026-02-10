import argparse
import os
import json
from io import StringIO
import yaml

def interaction_energies(definition_file: str, destination_file: str) -> None:
    with open(definition_file, 'r') as fp:
        raw = fp.read().splitlines()

    # Item section variables
    itemnum = 0
    items = []

    # Alternative reference variables
    alternative_reference = []

    # Iterate to reformat
    section = ''
    sections = ['description', 'items', 'plots', 'alternative_reference']
    for line in raw:
        for s in sections:
            if line.startswith(f'{s}:'):
                section = s
                break

        if section == 'items':
            if line.startswith('- !ruby'):
                line = f'- item: {itemnum}'
                itemnum += 1
            items.append(line)
        elif section == 'alternative_reference':
            alternative_reference.append(line)
    aryaml = yaml.safe_load(StringIO('\n'.join(alternative_reference)))

    # Read dictionaries from strings
    aryaml = yaml.safe_load(StringIO('\n'.join(alternative_reference)))
    itemyaml = yaml.safe_load(StringIO('\n'.join(items)))

    # Build final dictionary
    idat = 0
    units = 'eV'
    curvedict = {}
    mappings = {}
    for item in itemyaml['items']:
        key = item['curve_name']
        folder, geom = item['geometry'].split(':')
        if folder not in mappings:
            with open(folder + os.sep + folder + '.yaml', 'r') as fp:
                mappings[folder] = yaml.safe_load(fp)
        filename = os.path.join(folder, mappings[folder][geom])

        if key not in curvedict:
            curvedict[key] = {'x': [], 'files': [], 'blocks': {}}
            curvedict[key]['references'] = {dataset: [] for dataset in aryaml['alternative_reference']}
            curvedict[key]['references']['reference_value'] = []
            with open(filename, 'r') as fp:
                info = map(lambda x: x.replace(' ', ''), fp.read().splitlines()[1].split())
                info = {x.split('=')[0]: x.split('=')[1] for x in info}
                for k in info:
                    if k.startswith('selection_'):
                        _, b = k.split('_')
                        curvedict[key]['blocks'][b] = {'selection': info[k], 'charge': int(info[f'charge_{b}'])}
                    elif k == 'reference_unit':
                        if info[k] == 'kcalmol':
                            units = 'kcal/mol'
                    elif k == 'charge':
                        curvedict[key]['charge'] = int(info['charge'])

        curvedict[key]['x'].append(item['curve_x'])
        curvedict[key]['files'].append(filename)
        curvedict[key]['references']['reference_value'].append(item['reference_value'])
        for dataset in aryaml['alternative_reference']:
            curvedict[key]['references'][dataset].append(aryaml['alternative_reference'][dataset][idat])
        idat += 1

    savedict = {'name': os.path.splitext(definition_file)[0],
                'method': 'interaction_energies',
                'units': units,
                'curves': curvedict}
    with open(destination_file, 'w', encoding='utf-8') as fp:
        json.dump(savedict, fp, ensure_ascii=False, indent=4)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='cuby2fireball',
                                     description='Converts cuby4 yaml description files into valid Fireball json optimization files')
    parser.add_argument('input', type=str, help='cuby4 yaml dataset description file')
    parser.add_argument('method', type=str, help='type of dataset')
    parser.add_argument('-o', '--output', type=str, help='custom save path for the new Fireball json optimization file')
    args = parser.parse_args()
    output = args.output if args.output else os.path.splitext(args.input)[0] + '.json'

    if args.method == 'interaction_energies':
        interaction_energies(args.input, output)
    else:
        raise ValueError("Invalid method")

