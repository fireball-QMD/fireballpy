import sys

with open(sys.argv[1], 'r') as fp:
    lines = fp.read().splitlines()

inone = -1
precisionin = False
newlines = []
for i, l in enumerate(lines):
    ls = l.lstrip()
    indent = ' '*(len(l) - len(ls))
    if 'write' in l.lower() and '(*' in l:
        newlines.append('!' + l)
        continue
    if len(ls) > 1 and ls[0].isnumeric() and 'format' not in l.lower():
        newlines[-1] += ' &'
        newl = indent + '&' + ls[1:]
        newlines.append(newl)
        continue
    if ls.lower().startswith('real'):
        td = ls.split(' ')[0]
        newl = l.replace(td, 'real(kind=long)')
        newlines.append(newl)
        continue
    if ls.lower().startswith('character'):
        td = ls.split(' ')[0]
        tdl = '*' if '*' not in td else td.split('*')[1]
        newl = l.replace(td, f'character(len={tdl})')
        newlines.append(newl)
        continue
    if 'use' in l.lower() and 'precision' in l.lower():
        precisionin = True
    if not precisionin and 'implicit none' in l.lower():
        newlines.append(indent + 'use precision')
        inone = i - 1
    newlines.append(l)

if inone == -1 and not precisionin:
    print('NO IMPLICIT NONE. CHANGE AND RUN AGAIN')
    exit()

with open(sys.argv[1], 'w') as fp:
    fp.write('\n'.join(newlines))
