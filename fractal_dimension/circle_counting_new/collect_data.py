import glob
import re
import sys

precs = glob.glob('output*')
types = list(map(lambda x: x.replace('data/', ''), glob.glob('data/**/*.txt', recursive=True)))

dims = {}
for prec in precs:
    dims[prec] = {}
    for type in types:
        try:
            with open('{}/{}'.format(prec, type)) as f:
                lines = f.readlines()
                dimension = [line[19:-1] for line in lines if line[0:19] == "Fractal dimension: "]
                dims[prec][type] = dimension[0]
        except:
            pass

with open('dimension.csv', 'w') as f:
    for type in types:
        f.write(',{}'.format(re.sub(r'^.*/', '', type).replace('.txt', '').replace('_', ' ')))
    f.write('\n')
    
    for prec in precs:
        f.write('{}'.format(prec.replace('output', '')))
        for type in types:
            try:
                f.write(',{}'.format(dims[prec][type]))
            except:
                f.write(',MISSING')
        f.write('\n')
