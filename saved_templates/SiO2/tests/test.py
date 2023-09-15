import pyparsing as pp
import re 

epw_band_pattern = \
... \
+ pp.Literal('nbndsub') \
+ pp.Literal('=') \
+ ... \
+ pp.Literal('!') \
+ pp.Regex('.*', re.DOTALL)

num_epw_bands = int(epw_band_pattern.parse_file('../1.3-epw/epw.in')[3])

print(num_epw_bands)