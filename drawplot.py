import os, sys

# read all *.slha files in the current directory
slha_files = [f for f in os.listdir('.') if f.endswith('.slha')]

for slha_file in slha_files:
    os.system(f'slhaplot --maxmass=1800 --format=pdf --br=0% ' + slha_file)