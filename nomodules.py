import copy
import re
import os

retype = re.compile(
    r'(integer|real\*8|complex\*16|character\*?[\d]*\(len=[\d]+\)|logical)')
refixeddim = re.compile(r'dimension(\([\w\d\s,-:+\*]+\))')
recall = re.compile(r'call\s+([\w\d]+)')
realloc = re.compile(
    r'allocate\s*\(\s*([\w\d]+)\s*(\([\w\d\s,-:+\*]+\))\s*\)')
resubr = re.compile(r'subroutine\s+([\w\d]+)\s*\(([\w\d,\s]*)\s*\)')
reuse = re.compile(r'use\s+(m_\w+)\s*,\s*only\s*:([\w\d,\s]*)')
reintent = re.compile(r'(implicit)\s+(none)')
resize = re.compile(r'([\w\d]+)')


def get_lines(fpath):
    with open(fpath, 'r') as fp:
        for line in fp:
            yield line


def get_vars_from_module(mod):
    variables = {}
    inblock = False
    mod_iter = enumerate(get_lines(mod))
    for i, line in mod_iter:
        iold = i
        comment = '' if '!' not in line else line.split('!')[1].strip()
        inblock = True if comment == 'f2py start inblock' else inblock
        inblock = False if comment == 'f2py end inblock' else inblock
        uline = line.split('!')[0].strip().lower()
        if not uline or 'module' in uline:
            continue
        try:
            decls, vars = line.split('!')[0] \
                .strip() \
                .lower() \
                .split('::')
        except ValueError:
            raise ValueError(f"Invalid declaration in file {mod} "
                             f", line {iold+1}: '{uline}'. "
                             "Remember '::' is mandatory")
        decls = decls.replace(' ', '')
        while vars[-1] == '&':
            nline = next(mod_iter)[1].split('!')[0] \
                .strip() \
                .lower()
            vars = vars + os.linesep + '  ' + nline
        if '=' in vars:
            try:
                var, val = vars.split('=')
            except ValueError:
                raise ValueError("Two values assigned in one line "
                                 f"in file {mod}, line {iold+1}: '{uline}'")
            vars = [var.strip() + ' = ' + val.strip()]
        else:
            vars = vars.replace(' ', '').split(',')
        for var in vars:
            try:
                type_match = retype.search(decls).group(1)
            except IndexError:
                raise ValueError("Invalid variable type in file "
                                 f"{mod}, line {iold+1}: '{uline}'")
            if 'dimension' in decls:
                dims = True
                if 'allocatable' in decls:
                    alloc = True
                    dims_match = '(' + (decls.count(':')*':,')[:-1] + ')'
                else:
                    alloc = False
                    dims_match = refixeddim.search(decls).group(1)
            else:
                dims = False
                alloc = False

            varname = var.split('=')[0].strip()
            variables[varname] = {
                'VALUE': var if not inblock else varname,
                'TYPE': type_match,
                'DIMS': dims_match if dims else dims,
                'ALLOC': alloc,
                'PARAM': 'parameter' in decls,
                'INVAR': inblock
            }
    return variables


def _prep_line(line, iter):
    oline = line.split('!')[0].rstrip()
    uuline = oline.lower()
    if not uuline:
        return uuline, oline, line.split('!')[0]
    uline = uuline.lstrip()
    oline = oline.lstrip()
    indent = (len(uuline) - len(uline))*' '
    while uline[-1] == '&':
        nline = next(iter)[1].split('!')[0].strip()
        uline = uline[:-1] + nline[1:].lower()
        oline = oline[:-1] + nline[1:]
    return uline, oline, indent


def _check_subname(uline, subname, fpath):
    subr_match = resubr.search(uline)
    if subr_match is not None:
        if subname != subr_match.group(1):
            raise ValueError(
                f"Could not identify subroutine name in {fpath}")


def _get_self_depends(uline, depends):
    subr_match = resubr.search(uline)
    use_match = reuse.search(uline)
    if subr_match is not None:
        vars = subr_match.group(2).replace(' ', '').split(',')
        varss = list(filter(lambda x: x, vars))
        depends['SELF'] = varss
        return True
    if use_match is not None:
        mod = use_match.group(1)
        vars = use_match.group(2).replace(' ', '').split(',')
        varss = list(filter(lambda x: x, vars))
        if not varss:
            return True
        depends['MODS'].append(mod)
        depends[mod] = varss
        depends['TOTAL'][mod] = varss
        return True
    return False


def _assert_call_file(uline, callname, fpath, iold):
    callpath = os.path.join(os.path.dirname(fpath), f"{callname}.f90")
    if not os.path.isfile(callpath):
        raise ValueError(f"File with subroutine {callname} not found "
                         f"in {callpath}. Please take in account that "
                         "subroutine files should have full lowercase "
                         "names and should be named as the subroutine "
                         "they contain. Error in "
                         f"file {fpath}, line {iold+1}: '{uline}'")
    return callpath


def _update_depends(depends, subname, callname):
    for mod in depends[callname]['MODS']:
        if mod not in depends[subname]['MODS']:
            depends[subname]['MODS'].append(mod)
            depends[subname]['TOTAL'][mod] = []
        depends[subname]['TOTAL'][mod] += depends[callname]['TOTAL'][mod]


def _order_depends(depends, mod_vars):
    for mod in depends['MODS']:
        params, arrays, scalar = [], [], []
        for var in set(depends[mod]):
            if mod_vars[mod][var]['PARAM']:
                params.append(var)
            elif mod_vars[mod][var]['DIMS']:
                arrays.append(var)
            else:
                scalar.append(var)
        params.sort()
        arrays.sort()
        scalar.sort()
        depends[mod] = params + scalar + arrays
    for mod in depends['MODS']:
        params, arrays, scalar = [], [], []
        for var in set(depends['TOTAL'][mod]):
            if mod_vars[mod][var]['PARAM']:
                params.append(var)
            elif mod_vars[mod][var]['DIMS']:
                arrays.append(var)
            else:
                scalar.append(var)
            params.sort()
            arrays.sort()
            scalar.sort()
            depends['TOTAL'][mod] = params + scalar + arrays


def _long_line_split(s, max_len, indent=None):
    indent = indent if indent is not None else ''
    line_len = len(indent)
    pieces = s.split(',')
    splitted = indent
    for i, piece in enumerate(pieces):
        line_len += len(piece)
        if line_len > max_len - 2:
            if i == 0:
                raise ValueError(f"Cannot split line '{s}'. Achived "
                                 f"max number of characters ({max_len}) "
                                 "before the first variable")
            splitted += ' &' + os.linesep + indent + '  &'
            line_len = len(indent) + 3
        splitted += piece + ','
    if splitted[-1] == ',':
        splitted = splitted[:-1]
    return splitted


def _alloc_var(uline, mod_vars, depends):
    alloc_match = realloc.search(uline)
    if (alloc_match is None) and ('allocate' in uline):
        if 'deallocate' in uline:
            return
        raise ValueError(f"More than one allocate in line '{uline}'")
    if alloc_match is not None:
        varname = alloc_match.group(1)
        vardims = alloc_match.group(2).replace(' ', '')
        for mod in mod_vars:
            if varname in mod_vars[mod]:
                mod_vars[mod][varname]['DIMS'] = vardims
                mod_vars[mod][varname]['ALLOC'] = False


def _print_var(var, intent=None):
    pvar = var['TYPE']
    if var['DIMS']:
        pvar += f", dimension{var['DIMS']}"
    if var['ALLOC']:
        pvar += ', allocatable'
    if var['PARAM']:
        pvar += ', parameter'
    else:
        pvar += f", intent({intent})"
    return f"{pvar} :: {var['VALUE']}" + os.linesep


def _print_vars(vars, parms, mod_vars, intent):
    cline = ''
    for mod in mod_vars:
        for parm in parms:
            if parm not in mod_vars[mod]:
                continue
            cline += _print_var(mod_vars[mod][parm])
        for var in vars:
            if var not in mod_vars[mod]:
                continue
            cline += _print_var(
                mod_vars[mod][var],
                'in' if mod_vars[mod][var]['INVAR'] else intent)
    return cline


def _mod_code_mods(subname, depends, mod_vars,
                   parms, scalars, arrays, own):
    for mod in depends['MODS']:
        for var in depends['TOTAL'][mod]:
            if mod_vars[mod][var]['PARAM']:
                parms.append(var)
    for mod in depends['MODS']:
        for var in depends['TOTAL'][mod]:
            if mod_vars[mod][var]['DIMS']:
                arrays.append(var)
                vardims = mod_vars[mod][var]['DIMS'][1:-1]
                dims_match = resize.findall(vardims)
                for dims in dims_match:
                    if dims.isnumeric():
                        continue
                    if dims in mod_vars[mod]:
                        if mod_vars[mod][dims]['PARAM']:
                            parms.append(dims)
                        else:
                            scalars.append(dims)
                        continue
                    scalars.append(dims)
                    mod_vars[mod][dims] = {
                        'VALUE': dims,
                        'TYPE': 'integer',
                        'DIMS': False,
                        'ALLOC': False,
                        'PARAM': False,
                        'INVAR': False
                    }
            elif not mod_vars[mod][var]['PARAM']:
                scalars.append(var)


def _mod_code_self(subname, depends, mod_vars,
                   parms, scalars, arrays, own):
    for var in depends['SELF']:
        for mod in mod_vars:
            if (var not in mod_vars[mod]) or (not mod_vars[mod][var]['DIMS']):
                continue
            vardims = mod_vars[mod][var]['DIMS'][1:-1]
            dims_match = resize.findall(vardims)
            for dims in dims_match:
                if dims.isnumeric():
                    continue
                if dims in mod_vars[mod]:
                    if mod_vars[mod][dims]['PARAM']:
                        parms.append(dims)
                    else:
                        scalars.append(dims)
                    continue
                mod_vars[mod][dims] = {
                    'VALUE': dims,
                    'TYPE': 'integer',
                    'DIMS': False,
                    'ALLOC': False,
                    'PARAM': False
                }
    for var in depends['SELF']:
        if (var not in scalars) and (var not in arrays):
            own.append(var)


def _mod_code(subname, depends, mod_vars):
    parms, scalars, arrays, own = [], [], [], []
    _mod_code_mods(subname, depends, mod_vars,
                   parms, scalars, arrays, own)
    _mod_code_self(subname, depends, mod_vars,
                   parms, scalars, arrays, own)
    parms = list(set(parms))
    scalars = list(set(scalars))
    arrays = list(set(arrays))
    scalars.sort()
    arrays.sort()
    parms.sort()
    varss = list(filter(lambda x: x, scalars + arrays))
    owns = list(filter(lambda x: x, own))
    vars = ', '.join(varss + owns)
    cline = f"subroutine {subname}({vars})"
    return _long_line_split(cline, 128) + os.linesep, varss, parms, owns


def _join_code(code, subname, depends, mod_vars, intent):
    code_lines = code.splitlines()
    for i, line in enumerate(code_lines):
        intent_match = reintent.search(line)
        if (intent_match is not None) and (len(intent_match.groups()) == 2):
            intent_end = i + 1
            break
    try:
        intent_end = intent_end
    except NameError:
        raise ValueError(f"'implicit none' not found in subroutine {subname}")

    for i, line in enumerate(code_lines[::-1]):
        if 'end' in line:
            sub_end = len(code_lines) - i - 1
            break

    try:
        sub_end = sub_end
    except NameError:
        raise ValueError(f"Subroutine {subname} was not closed")

    mod_code, vars, parms, owns = _mod_code(
        subname, depends[subname], mod_vars)
    code = mod_code + \
        os.linesep.join(code_lines[:intent_end]) + \
        os.linesep + '  ' + \
        _print_vars(vars, parms, mod_vars, intent) \
        .replace(os.linesep, os.linesep + '  ')[:-2] + \
        os.linesep.join(code_lines[intent_end:sub_end]) + \
        os.linesep + \
        f"end subroutine {subname}"
    return code, vars


def edit_file(fpath, mod_vars, depends=None):
    subname = os.path.splitext(os.path.basename(fpath))[0]
    intent = 'out' if depends is None else 'inout'
    depends = {} if depends is None else depends

    depends[subname] = {'SELF': [], 'MODS': [], 'TOTAL': {}}
    cmod_vars = copy.deepcopy(mod_vars)
    code = ''
    file_iter = enumerate(get_lines(fpath))
    for i, line in file_iter:
        iold = i
        uline, oline, indent = _prep_line(line, file_iter)
        if not uline:
            code += line
            continue

        _check_subname(uline, subname, fpath)
        if _get_self_depends(uline, depends[subname]):
            continue
        _alloc_var(uline, mod_vars, depends)
        if 'call' in uline:
            callname = recall.search(uline).group(1)
            if callname in depends:
                continue
            callpath = _assert_call_file(uline, callname, fpath, iold)
            vars = edit_file(callpath, mod_vars, depends)
            _update_depends(depends, subname, callname)
            cline = f"call {callname}({', '.join(vars)})"
            code += _long_line_split(cline, 128, indent) + os.linesep
            continue
        code += _long_line_split(oline, 128, indent) + os.linesep

    _order_depends(depends[subname], mod_vars)
    code, vars = _join_code(code, subname, depends, cmod_vars, intent)
    wrtpath = os.path.join(os.path.dirname(fpath), f"f2py_{subname}.f90")
    with open(wrtpath, 'w') as fp:
        fp.write(code)

    return vars


def main(modules, entries):
    mod_vars = {}
    for mod in modules:
        modname = os.path.splitext(os.path.basename(mod))[0].lower()
        mod_vars[modname] = get_vars_from_module(mod)

    for entry in entries:
        edit_file(entry, mod_vars)


if __name__ == '__main__':
    MODULES = ["parallel/src/m_fdata.f90"]
    ENTRIES = ["parallel/src/LOADFDATA/load_fdata.f90",
               "parallel/src/LOADFDATA/calc_me_max.f90"]

    main(MODULES, ENTRIES)
