
import pyparsing as pp
import re


# Copy to pwbands.in
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


# Copy to ph.in and matdyn.in. 
pattern = ... + pp.Literal('/\n')  + pp.Regex(r'.*', re.DOTALL)

ph_parsed_text = pattern.parse_file('ph.in').as_list()
matdyn_parsed_text = pattern.parse_file('matdyn.in').as_list()

subst_text = ''
with open('kgrid.out', 'r') as f: subst_text = ''.join(f.readlines()[1:])

ph_new_text = ''.join([
    ph_parsed_text[0],
    ph_parsed_text[1],
    subst_text
])

matdyn_new_text = ''.join([
    matdyn_parsed_text[0],
    matdyn_parsed_text[1],
    subst_text
])

# with open('ph.in', 'w') as f: f.write(ph_new_text)        # Commented out. Not writing to ph.in for now. 
with open('matdyn.in', 'w') as f: f.write(matdyn_new_text)

# Copy number of irreducible k-points to create_epw_save.py. 
pattern = ...+ pp.Literal('nqpt') + '=' + pp.Word(pp.nums) + pp.Regex(r'.*', re.DOTALL)
parsed_text = pattern.parse_file('./create_epw_save.py')

nirk = len(open('kgrid.out', 'r').readlines()) - 2
new_text = ''.join([
    parsed_text[0],
    parsed_text[1],
    parsed_text[2],
    str(nirk),
    parsed_text[4]
])
open('create_epw_save.py', 'w').write(new_text)
