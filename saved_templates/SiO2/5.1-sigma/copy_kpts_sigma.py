
import pyparsing as pp
import re
import numpy as np

# Create the pattern. 
pattern = ... + pp.Literal('begin kpoints') + ... + pp.Literal('end') + pp.Regex(r'.*', re.DOTALL)

# Parse epsilon.inp.
parsed_text = pattern.parse_file('sigma.inp').as_list()

# Create text to substitute. 
subst_text = ''
data = np.loadtxt('../2.1-wfn/kgrid.out', skiprows=2)
shape = data.shape
if len(shape) == 1:
    data = data.reshape(1, shape[0])
np.savetxt('sigma_kpoints.csv', data, delimiter=' ', fmt='%10.6f')
with open('sigma_kpoints.csv', 'r') as f: subst_text = f.read()

# Create the new text. 
new_text = ''.join([
    parsed_text[0],
    parsed_text[1],
    '\n',
    subst_text,
    parsed_text[3],
    parsed_text[4]
])

with open('sigma.inp', 'w') as f: f.write(new_text)
