
import pyparsing as pp
import re

pattern = ... + pp.Literal('K_POINTS') + ... + pp.Char(pp.srange('[A-Z]')) + pp.Regex(r'.*', re.DOTALL)

parsed_text = pattern.parse_file('./pwbands.in').as_list()

subst_text = ''
with open('kgrid.out', 'r') as f: 
    subst_text = f.read()

new_text = ''.join([
    parsed_text[0],
    subst_text,
    parsed_text[3],
    parsed_text[4]
])

with open('pwbands.in', 'w') as f: 
    f.write(new_text)
