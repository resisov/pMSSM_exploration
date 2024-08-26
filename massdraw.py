# -*- coding: utf-8 -*-
import os, sys
import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np
import argparse
from adjustText import adjust_text # type: ignore

# Define the parser
parser = argparse.ArgumentParser(description='Draw the mass spectrum from SLHA file')
parser.add_argument('--maxmass', type=float, default=10000, help='Maximum mass to draw')
args = parser.parse_args()

def readslha(file):
    global run_number, iter_number
    run_number = file.split('_')[2].split('.')[0]
    iter_number = file.split('_')[3].split('.')[0]
    print('Run number:', run_number, 'Iteration number:', iter_number)

    with open(file, 'r') as f:
        lines = f.readlines()
    masses = {}
    susy = []
    decays = []
    for i, line in enumerate(lines):
        ### Block MASS
        #print(line)
        if 'Block MASS' in line:
            for j in range(i+2, len(lines)):
                if 'Block' in lines[j]:
                    break
                else:
                    pdgid = lines[j].split()[0]
                    mass = lines[j].split()[1]
                    name = lines[j].split()[-1]
                    if float(mass) > args.maxmass or float(pdgid) < 25:
                        continue
                    elif pdgid == '25':
                        masses[pdgid] = {
                            'column': 1,
                            'mass': abs(float(mass)),
                            'color': 'blue',
                            'name': r'$h^0$'
                        }
                        susy.append(pdgid)
                    elif pdgid == '35':
                        masses[pdgid] = {
                            'column': 1,
                            'mass': abs(float(mass)),
                            'color': 'blue',
                            'name': r'$H^0$'
                        }
                        susy.append(pdgid)
                    elif pdgid == '36':
                        masses[pdgid] = {
                            'column': 1,
                            'mass': abs(float(mass)),
                            'color': 'blue',
                            'name': r'$A^0$'
                        }
                        susy.append(pdgid)
                    elif pdgid == '37':
                        masses[pdgid] = {
                            'column': 1,
                            'mass': abs(float(mass)),
                            'color': 'red',
                            'name': r'$H^{\pm}$'
                        }
                        susy.append(pdgid)
                    elif pdgid == '1000001':
                        masses[pdgid] = {
                            'column': 2,
                            'mass': abs(float(mass)),
                            'color': 'blue',
                            'name': r'$\~{q}_L$'
                        }
                        susy.append(pdgid)
                    elif pdgid == '2000001':
                        masses[pdgid] = {
                            'column': 2,
                            'mass': abs(float(mass)),
                            'color': 'blue',
                            'name': r'$\~{q}_R$'
                        }
                        susy.append(pdgid)
                    elif pdgid == '1000005':
                        masses[pdgid] = {
                            'column': 2,
                            'mass': abs(float(mass)),
                            'color': 'green',
                            'name': r'$\~{b}_1$'
                        }
                        susy.append(pdgid)
                    elif pdgid == '2000005':
                        masses[pdgid] = {
                            'column': 2,
                            'mass': abs(float(mass)),
                            'color': 'green',
                            'name': r'$\~{b}_2$'
                        }
                        susy.append(pdgid)
                    elif pdgid == '1000006':
                        masses[pdgid] = {
                            'column': 2,
                            'mass': abs(float(mass)),
                            'color': 'red',
                            'name': r'$\~{t}_1$'
                        }
                        susy.append(pdgid)
                    elif pdgid == '2000006':
                        masses[pdgid] = {
                            'column': 2,
                            'mass': abs(float(mass)),
                            'color': 'red',
                            'name': r'$\~{t}_2$'
                        }
                        susy.append(pdgid)
                    elif pdgid == '1000011':
                        masses[pdgid] = {
                            'column': 4,
                            'mass': abs(float(mass)),
                            'color': 'blue',
                            'name': r'$\~{\ell}_L$'
                        }
                        susy.append(pdgid)
                    elif pdgid == '2000011':
                        masses[pdgid] = {
                            'column': 4,
                            'mass': abs(float(mass)),
                            'color': 'blue',
                            'name': r'$\~{\ell}_R$'
                        }
                        susy.append(pdgid)
                    elif pdgid == '1000015':
                        masses[pdgid] = {
                            'column': 4,
                            'mass': abs(float(mass)),
                            'color': 'blue',
                            'name': r'$\~{\tau}_1$'
                        }
                        susy.append(pdgid)
                    elif pdgid == '2000015':
                        masses[pdgid] = {
                            'column': 4,
                            'mass': abs(float(mass)),
                            'color': 'blue',
                            'name': r'$\~{\tau}_2$'
                        }
                        susy.append(pdgid)
                    elif pdgid == '1000016':
                        masses[pdgid] = {
                            'column': 4,
                            'mass': abs(float(mass)),
                            'color': 'red',
                            'name': r'$\~{\nu}_{\tau}$'
                        }
                        susy.append(pdgid)
                    elif pdgid == '1000021':
                        masses[pdgid] = {
                            'column': 2,
                            'mass': abs(float(mass)),
                            'color': 'green',
                            'name': r'$\~{g}$'
                        }
                        susy.append(pdgid)
                    elif pdgid == '1000022':
                        masses[pdgid] = {
                            'column': 3,
                            'mass': abs(float(mass)),
                            'color': 'blue',
                            'name': r'$\~{\chi}_1^0$'
                        }
                        susy.append(pdgid)
                    elif pdgid == '1000023':
                        masses[pdgid] = {
                            'column': 3,
                            'mass': abs(float(mass)),
                            'color': 'blue',
                            'name': r'$\~{\chi}_2^0$'
                        }
                        susy.append(pdgid)
                    elif pdgid == '1000025':
                        masses[pdgid] = {
                            'column': 3,
                            'mass': abs(float(mass)),
                            'color': 'blue',
                            'name': r'$\~{\chi}_3^0$'
                        }
                        susy.append(pdgid)
                    elif pdgid == '1000035':
                        masses[pdgid] = {
                            'column': 3,
                            'mass': abs(float(mass)),
                            'color': 'blue',
                            'name': r'$\~{\chi}_4^0$'
                        }
                        susy.append(pdgid)
                    elif pdgid == '1000024':
                        masses[pdgid] = {
                            'column': 3,
                            'mass': abs(float(mass)),
                            'color': 'red',
                            'name': r'$\~{\chi}_1^{\pm}$'
                        }
                        susy.append(pdgid)
                    elif pdgid == '1000037':
                        masses[pdgid] = {
                            'column': 3,
                            'mass': abs(float(mass)),
                            'color': 'red',
                            'name': r'$\~{\chi}_2^{\pm}$'
                        }
                        susy.append(pdgid)

        ### Decays
        if 'DECAY' in line:
            pdgid = line.split()[1]
            if pdgid in susy:
                total_width = float(line.split()[2])
                print('Total width of', pdgid, 'is', total_width, 'GeV')
                if total_width == 0.0:
                    print('This particle is LSP', pdgid)
            
                for j in range(i+2, len(lines)):
                    br = lines[j].split()[0]
                    prong = lines[j].split()[1]
                    if prong == '2':
                        d1 = lines[j].split()[2]
                        d2 = lines[j].split()[3]
                        if d1 in susy or d2 in susy:
                            if pdgid != d1 and pdgid != d2:
                                if d1 in masses.keys():
                                    #print('Decay of', pdgid, '->', d1, d2, 'with BR:', float(br))
                                    decay_set = [pdgid, d1, float(br)]
                                    decays.append(decay_set)
                                if d2 in masses.keys():
                                    #print('Decay of', pdgid, '->', d1, d2, 'with BR:', float(br))
                                    decay_set = [pdgid, d2, float(br)]
                                    decays.append(decay_set)

                    elif prong == '3':
                        d1 = lines[j].split()[2]
                        d2 = lines[j].split()[3]
                        d3 = lines[j].split()[4]
                        if d1 in susy or d2 in susy or d3 in susy:
                            if pdgid != d1 and pdgid != d2 and pdgid != d3:
                                if d1 in masses.keys():
                                    #print('Decay of', pdgid, '->', d1, d2, 'with BR:', float(br))
                                    decay_set = [pdgid, d1, float(br)]
                                    decays.append(decay_set)
                                if d2 in masses.keys():
                                    #print('Decay of', pdgid, '->', d1, d2, 'with BR:', float(br))
                                    decay_set = [pdgid, d2, float(br)]
                                    decays.append(decay_set)
                                if d3 in masses.keys():
                                    #print('Decay of', pdgid, '->', d1, d2, 'with BR:', float(br))
                                    decay_set = [pdgid, d3, float(br)]
                                    decays.append(decay_set)                                    
                    if 'DECAY' in lines[j]:
                        break

    return masses, decays

def deleteDuplication(decays):
    import pandas as pd
    df = pd.DataFrame(decays, columns=['parent', 'daughter', 'width'])
    # if parent and daughter are the same, sum up the width and delete the duplication
    df = df.groupby(['parent', 'daughter']).sum().reset_index()
    decays = df.values.tolist()
    return decays

def addArrow(parent,daughter,dwidth,ax):
    if parent['color'] == 'blue':
        start = (parent['column'], parent['mass'])
    else:
        start = (parent['column'] + 0.1, parent['mass'])

    if daughter['color'] == 'blue':
        end = (daughter['column'], daughter['mass'])
    else:
        end = (daughter['column'] + 0.1, daughter['mass'])
    ax.annotate("", xy=end, xytext=start, arrowprops=dict(arrowstyle="->", color='black', linestyle='dashed', alpha=dwidth))

def drawplot(masses, decays):
    hep.style.use(hep.style.CMS)
    fig, ax = plt.subplots(figsize=(9, 8))
    texts = []
    for pdgid, properties in masses.items():
        if properties['color'] == 'blue':
            ax.errorbar(properties['column'], properties['mass'], xerr=0.3, fmt='', color=properties['color'], label=properties['name'], linewidth=2)
            texts.append(ax.text(properties['column'] - 0.55, properties['mass']-(args.maxmass*0.02), properties['name'], fontsize=20, fontdict={'family': 'Baskerville'}))
        else:
            ax.errorbar(properties['column'] + 0.1, properties['mass'], xerr=0.3, fmt='', color=properties['color'], label=properties['name'], linewidth=2)
            texts.append(ax.text(properties['column'] + 0.45, properties['mass']-(args.maxmass*0.02), properties['name'], fontsize=20, fontdict={'family': 'Baskerville'}))
    adjust_text(texts)
    deleteDuplication(decays)
    for l in decays:
        parent = masses[l[0]]
        daughter = masses[l[1]]
        dwidth = l[2]
        addArrow(parent, daughter, dwidth, ax)
    ax.set_ylabel('mass (GeV)')
    ax.set_ylim(0, args.maxmass)
    ax.set_xlim(-0.1,5.1)
    # delete the x-axis label and ticks
    ax.set_xticks([])
    ax.set_xlabel('')
    fig.tight_layout()
    fig.savefig('./xsecplots/mass_spectrum_'+run_number+'_'+iter_number+'.png')
    plt.close()

if __name__ == '__main__':
    # Read the SLHA file
    filelist = os.listdir('./slhafiles/')
    for file in filelist:
        if file.endswith('.slha'):
            masses, decays = readslha(file)
            # Draw the mass spectrum
            drawplot(masses, decays)