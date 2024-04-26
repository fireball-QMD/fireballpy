import re
import os

retype = re.compile(
    r'integer|real\*8|complex\*16|character\*?[0-9]*\(len=[0-9]+\)|logical')
refixed = re.compile(r'dimension\((?:\d(?::\d+)?[,|\)])+')
realloc = re.compile(r'dimension\((?::[,|\)])+')


def get_lines(fpath):
    with open(fpath, 'r') as fp:
        for line in fp:
            yield line


def get_vars_from_module(mod):
    variables = {}
    mod_iter = get_lines(mod)
    for i, line in enumerate(mod_iter):
        try:
            decls, vars = line.split('!')[0] \
                .strip() \
                .lower() \
                .split('::')
        except ValueError:
            continue
        decls = decls.replace(' ', '')
        while vars[-1] == '&':
            nline = next(mod_iter).split('!')[0] \
                .strip() \
                .lower()
            vars = vars + os.linesep + '  ' + nline
        if '=' in vars:
            try:
                var, val = vars.split('=')
            except ValueError:
                raise ValueError("Two values assigned in one line"
                                 f"in file {mod}, line {i+1}: '{line}'")
            vars = [var.strip() + ' = ' + val.strip()]
        else:
            vars = vars.replace(' ', '').split(',')
        for var in vars:
            type_match = retype.findall(decls)
            if len(type_match) != 1:
                raise ValueError("Invalid variable type in file"
                                 f"{mod}, line {i+1}: '{line}'")
            fixed_match = refixed.findall(decls)
            alloc_match = realloc.findall(decls)
            dims_match = fixed_match + alloc_match

            variables[var] = {
                'TYPE': type_match[0],
                'DIMS': dims_match[0] if dims_match else '',
                'ALLOC': 'allocatable' if alloc_match else '',
                'PARAM': 'parameter' if 'parameter' in decls else ''
            }
    return variables


def get_calls(fpath):
    pass


MODULES = ["src/M_fdata.f90"]
vars = get_vars_from_module(MODULES[0])

for var in vars:
    pvar = vars[var]['TYPE']
    if vars[var]['DIMS']:
        pvar += ', ' + vars[var]['DIMS']
    if vars[var]['ALLOC']:
        pvar += ', ' + vars[var]['ALLOC']
    if vars[var]['PARAM']:
        pvar += ', ' + vars[var]['PARAM']
    pvar += ' :: ' + var
    print(pvar)
