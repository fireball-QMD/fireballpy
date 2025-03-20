import sys

print(sys.argv[1])
with open(sys.argv[1], 'r') as fp:
    lines = fp.read().splitlines()

inone = -1
precisionin = False
newlines = []
for i, l in enumerate(lines):
    ls = l.lstrip()
    indent = ' '*(len(l) - len(ls))
    #if 'write' in l.lower() and '(*' in l:
    #    newlines.append('!' + l)
    #    continue
    if len(ls) > 1 and ls[0] == '&':
        for j in range(len(newlines[-1])):
            if newlines[-1][len(newlines[-1])-j-1] == '&':
                break
            if newlines[-1][len(newlines[-1])-j-1] == ' ':
                continue
            if newlines[-1][len(newlines[-1])-j-1] != '&':
                newlines[-1] += ' &'
                break
    if len(ls) > 1 and ls[0].isnumeric() and 'format' not in l.lower() and 'continue' not in l.lower() and 'write' not in l.lower():
        newlines[-1] += ' &'
        l = indent + '&' + ls[1:]
    if ls.lower().startswith('real'):
        td = ls.split(' ')[0]
        rep = 'real(kind=long)'
        if td[-1] == ',':
            rep += ','
        newl = l.replace(td, rep)
        newlines.append(newl)
        continue
    if ls.lower().startswith('complex'):
        td = ls.split(' ')[0]
        rep = 'complex(kind=long)'
        if td[-1] == ',':
            rep += ','
        newl = l.replace(td, rep)
        newlines.append(newl)
        continue
    if ls.lower().startswith('character'):
        td = ls.split(' ')[0]
        if 'len' not in ls.lower() and '*' in td:
            tdl = td.split('*')[1]
            rep = f'character(len={tdl})'
            if td[-1] == ',':
                rep += ','
            newl = l.replace(td, rep)
            newlines.append(newl)
            continue
    if 'use' in l.lower() and 'precision' in l.lower():
        precisionin = True
    if not precisionin and 'implicit none' in l.lower():
        newlines.append(indent + 'use precision')
        inone = i - 1
    for err in ['. or.', '. or .', '.or .']:
        if err in l.lower():
            l = l.replace(err, ' .or. ')
    for err in ['. and.', '. and .', '.and .']:
        if err in l.lower():
            l = l.replace(err, ' .and. ')

    for err, funname in [('dble', 'real'), ('dfloat', 'real'), ('dcmplx', 'cmplx')]:
        if err in l.lower():
            for cc in range(l.lower().count(err)):
                for j in range(len(err), len(l)):
                    if l.lower()[(j-len(err)):j] == err:
                        count = 0
                        for k in range(j, len(l)):
                            if l[k] == '(':
                                count += 1
                            if l[k] == ')':
                                count -= 1
                            if count == 0:
                                l = l[:(j-len(err))] + funname + l[j:k] + ', kind=long' + l[k:]
                                break
                        break

    for err in ['dsin', 'dcos', 'dsqrt']:
        l.replace(err, err[1:])
    newlines.append(l)

if inone == -1 and not precisionin:
    print('NO IMPLICIT NONE. CHANGE AND RUN AGAIN')
    exit()

if sys.argv[1].endswith('.f90'):
    with open(sys.argv[1], 'w') as fp:
        fp.write('\n'.join(newlines))
else:
    with open(sys.argv[1] + '90', 'w') as fp:
        fp.write('\n'.join(newlines))
