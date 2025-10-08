#! /usr/bin/python3
'''
When given a trial directory in the outputFiles subdirectory (ie, with defaults something like ubq_1),
will generate a stacked plot displaying various changes in droplet composition as the simulation progresses.
This includes number of each molecule in the droplet, protein charge, and ratio of droplet charge to Rayleigh
limit. Requires protein mass be inputted to the --mass arg in kDa (ie., 8.56 for ubiquitin).

Example call: python composition.py -d ubq_1 -mass 8.56
'''
import matplotlib.pyplot as plt
from matplotlib import gridspec
import os
import math
import argparse
import numpy as np

''' ----- Constants ----- '''
Na = 6.0221408*10**23   # Avagrados Number
e  = 1.60217663*10**-19 # e- charge (C)
eo = 8.854187*10**-12   # Vacuum permitivity (F/m)
st = 0.0727             # Water surface tension (N m^-1)
water_d   = 1000        # Water density (kg/m^3)
prot_d    = 1220        # Protein density (kg/m^3)
water_vol = 0.01801 / (Na * water_d) # Single water molecule volume (m^3/molec)

''' ----- Functions ----- '''
def rayleigh_limit(num_water):
    rl_const = (8*math.pi / e) * ((3 * eo * st) / (4 * math.pi))**(1/2)
    Vprot = (args.m / Na) / prot_d
    return rl_const * (Vprot + water_vol * num_water)**(1/2)
#
def openfile(fname, dpath):
    data = {}
    with open(os.path.join(dpath, fname), 'r') as f:
        header = f.readline().strip()
        header = header.replace('#', '').split()
        values = np.loadtxt(f, skiprows=1).T
        for k, v in zip(header, values):
            data[k] = v
    return data

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', type=str,   required=True, help="Name of trial dir in outputFiles")
    parser.add_argument('-m', type=float, required=True, help="Mass of protein in kDa") 
    parser.add_argument('-s', action="store_true",       help="Auto-save plot") 
    args = parser.parse_args()

    # Get droplet composition, system/protein charge, and simulation time
    dpath = os.path.join('../outputFiles', args.d, 'data')
    composition = openfile("composition.txt", dpath)
    charge = openfile("charge.txt", dpath)
    time = np.array([i*0.004 for i in range(len(composition['SOL']))])

    # Get ratio of system charge to Rayleigh limit
    charge['RL'] = rayleigh_limit(composition['SOL'])
    charge['QZR'] = charge['System'] / charge['RL']

    print(f'Final protein charge is {charge["Protein"][-1]}.')
    print(f'Final system charge (including adducts) is {charge["System"][-1]}.')

    ''' ----- Plotting ----- '''
    fig = plt.figure()
    fig.set_size_inches(8, 7)
    plt.rcParams["font.size"] = 15
    plt.rcParams['lines.linewidth'] = 0.75
    gs = gridspec.GridSpec(5, 1)

    # Waters
    ax0 = plt.subplot(gs[0])
    ax0.set_xlim(0, time[-1])
    ax0.set_ylim(0, np.max(composition['SOL']))
    ax0.set_ylabel('H$_{2}$O')
    ax0.get_xaxis().set_ticks([])
    ax0.plot(time, composition['SOL'], color='black')

    # Acetate/acetic acid
    ax1 = plt.subplot(gs[1])
    ax1.set_xlim(0, time[-1])
    buff_max = max(np.max(composition['ATX']), np.max(composition['NXH']))
    ax1.set_ylim(0, buff_max)
    ax1.set_ylabel('Ac')
    ax1.get_xaxis().set_ticks([])
    ax1.plot(time, composition['AHX'], label='CH$_{3}$CO$_{2}$H', color='black')
    ax1.plot(time, composition['ATX'], label='CH$_{3}$CO$_{2}^{-}$', color='red')
    ax1.legend(loc='upper right', frameon=False)

    # Ammonium/ammonia
    ax2 = plt.subplot(gs[2])
    ax2.set_xlim(0, time[-1])
    ax2.set_ylim(0, buff_max)
    ax2.set_ylabel('NH$_{x}$')
    ax2.get_xaxis().set_ticks([])
    ax2.plot(time, composition['NXH'],  label='NH$_{4}^{+}$', color='black')
    ax2.plot(time, composition['NXX'],  label='NH$_{3}$', color='red')
    ax2.legend(loc='upper right', frameon=False)

    # Protein charge
    ax3 = plt.subplot(gs[3])
    ax3.set_xlim(0, time[-1])
    ax3.set_ylabel('Protein $z$')
    ax3.get_xaxis().set_ticks([])
    ax3.plot(time, charge['Protein'], color='black')

    # Ratio of droplet charge to Rayleigh ratio
    ax4 = plt.subplot(gs[4])
    ax4.set_xlim(0, time[-1])
    ax4.set_ylim(0, max(1.2, np.max(charge['QZR'])))
    ax4.set_xlabel('Time (ns)')
    ax4.set_ylabel('$z/z_{r}$')
    ax4.hlines(y=1.0, xmin=0, xmax=time[-1], linestyle='dashed', color='black', label='Rayleigh Limit')
    ax4.plot(time, charge['QZR'], color='black')
    ax4.legend(loc='lower right', frameon=False)

    fig.subplots_adjust(hspace=0.15)
    fig.align_labels()

    if args.s:
        if not os.path.exists('pics'):
            os.mkdir('pics')
        plt.draw()
        plt.savefig(f'pics/composition_{args.d}.png', dpi=300, bbox_inches='tight')
    else:
        plt.show()