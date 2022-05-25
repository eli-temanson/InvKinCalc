# Python script to create a json file of nuclear masses in MeV/c^2
#
# Author: Eli Temanson
# Date: May 25, 2022

import urllib.request
import pandas as pd

# the service URL
livechart = "https://nds.iaea.org/relnsd/v0/data?"

def lc_read_csv(url):
    req = urllib.request.Request(url)
    req.add_header('User-Agent', 'Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:77.0) Gecko/20100101 Firefox/77.0')
    return pd.read_csv(urllib.request.urlopen(req))

# Get the ground state data for all nuclei
df=lc_read_csv(livechart+"fields=ground_states&nuclides=all")

# Create another column with the number of nucleons
df['A'] = df['z'] + df['n']

# Make a unique ID, i.e. combine the number of nucleons and symbol ex. "25Al"
df['isotope'] = df['A'].astype(str) + df['symbol']

# Convert the atomic mass from a string to a double, empty strings get value "null"
## AND convert from micro-AMU to AMU
df['atomic_mass'] = pd.to_numeric(df['atomic_mass'], errors='coerce')*1.0e-6

amu_to_MeV = 931.4940954
electron_mass = 0.000548579909
df['isotope_mass'] = (df['atomic_mass'] - df['z']*electron_mass)*amu_to_MeV

mass_table = df[['isotope','isotope_mass']]

# Set the ID (key) as the isotope symbol name
mass_table.set_index("isotope", drop=True, inplace=True)

# Convert the dataframe to a json file, 
js = mass_table.to_json("nuclear_masses.json", orient="index")
# put the final json file into a https://jsonformatter.org/
