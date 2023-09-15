
import pyparsing as pp
import re
import numpy as np

# Create the pattern. 
pattern = ... + pp.Literal('begin qpoints\n') + ... + pp.Literal('end') + pp.Regex(r'.*', re.DOTALL)

# Parse epsilon.inp.
parsed_text = pattern.parse_file('epsilon.inp').as_list()

# Create text to substitute. 
subst_text = ''
data = np.loadtxt('../2.1-wfn/kgrid.out', skiprows=2)
shape = data.shape
if len(shape) == 1:
    data = data.reshape(1, shape[0])
max_rows = data.shape[0]
new_data = np.zeros((max_rows, 5), dtype='f8')
new_data[:, 0:4] = data
new_data[0, 0] = 0.0  # For the qshifted point.
new_data[0, 1] = 0.0   # For the qshifted point.
new_data[0, 2] = 0.001   # For the qshifted point.
new_data[0, 3] = 1.0  # For the qshifted point.
new_data[0, 4] = 1    # For the qshifted point.

with open('epsilon_qpoints.csv', 'w') as f:
    for row in new_data:
        f.write(f'{row[0]:10.6f} {row[1]:10.6f} {row[2]:10.6f} {row[3]:10.6f} {int(row[4])}\n')

with open('epsilon_qpoints.csv', 'r') as f: subst_text = f.read()

# Create the new text. 
new_text = ''.join([
    parsed_text[0],
    parsed_text[1],
    subst_text,
    parsed_text[3],
    parsed_text[4]
])

with open('epsilon.inp', 'w') as f: f.write(new_text)
