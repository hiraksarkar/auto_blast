import sys

cont_dict = {}
with open(sys.argv[1]) as fp:
    for line in fp:
        if line[0] == '#':
            continue
        else:
            if line.strip():
                elems = list(filter(lambda x: len(x) != 0, line.strip().split('\t')))
                sp, seq = elems[0], elems[1]
                cont_dict[sp] = seq

with open(sys.argv[2], 'w') as fp:
    for sp in cont_dict:
        fp.write('>{}\n'.format(sp))
        fp.write('{}\n'.format(cont_dict[sp]))
