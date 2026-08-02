"""Functions for scanning Guth VCP files"""
import sys, os
import datetime
import re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tck


def ask_YorN(s_desc):
# Part of APUF (Andy's Python Useful Functions).
# Ask yes or no question and return a logical.
# Argument is a string containing the question, the function automatically adds syntax guidance onto the end of this.
    while True:
        s_YorN=input(s_desc+' Enter "y" or "n", or "e" to exit: ')
        if s_YorN=='y' or s_YorN=='Y':
            return True
        elif s_YorN=='n' or s_YorN=='N':
            return False
        elif s_YorN=='e' or s_YorN=='E':
            sys.exit('Exiting.')
        elif s_YorN=='42':
            print ('Nice reference you hoopy frood, but try again.')
        else:
            print ('Unrecognised input, please try again.')

def _read_file(filepath):
    """
    Reads in tab separated Guth VCP data
    All comments (#) are stripped and each line is an element in the returned
    list. Each line element is stripped of whitespace at either end, and is
    partitioned at the first whitespace character.
    Further splitting of elements will be required based on what file type is
    being read.
    """
    file = []
    with open(filepath, "r") as fp:
        for line in fp:
            # .partition returns a tuple: everything before the hash, the hash
            # and everything after the hash.
            # By indexing with [0] takes just the part before the hash.
            line = line.partition("#")[0]
            line = re.split('[ \t]', line)
            # print(line)

            # if line not empty append to file list
            if line[0]:
                file.append(line)
    return file


def _get_var(ifile, find_str):
    """
    Given (all) lines of a file and a specific string to search for
    return the item which follows (as a string).
    """
    y = [x[1] for x in ifile if x[0] == find_str]
    if y:
        var = y[0]
    else:
        var = None

    return var

# Got this from www.geeksforgeeks.org/adding-value-labels-on-a-matplotlib-bar-chart/
def addlabels(x,y):
    for i in range(len(x)):
        if larger_fonts:
            plt.text(i, y[i], y[i], ha = 'center', fontsize=11)
        else:
            plt.text(i, y[i], y[i], ha = 'center')

# Parse command line.
# s_vcpFile='glr.vcp'
s_vcpFile=sys.argv[1]
if not os.path.isfile(s_vcpFile): sys.exit('Error: input file "'+s_vcpFile+'" does not exist. Exiting.')

# Main program.
# look for line beginning  -60.000000
guth_data=False
found_guth=False
dft= []

introduction_a = 'This script scans a Radiance Guth VCP report and creates a polar plot.'
introduction_b = 'Give it the name of the vcp file and it should create a *.png file.'
print (introduction_a)
print (introduction_b)
message = 'Proceed scanning ' + s_vcpFile + '(Y/N)?'
doit = ask_YorN(message)
if not doit:
    sys.exit('Exiting...')

# Scan the vcp report into an overall entity holding all the lines to work with.
vcp_rep = _read_file(s_vcpFile)

larger_fonts = ask_YorN('Use larger fonts in graphs?')

# Get user opinion on verbosity.
verbose = ask_YorN('Display data as it is scanned?')

# Loop through each line, find the number of tokens in the current line.
for s_line in vcp_rep:
    n_tokens = len(s_line)
    if verbose:
        print (s_line)

# Parse skipping lines until -60.00000 then append to dft.
    if s_line[0] == 'glarendx':
        found_guth=True
        continue
    if s_line[0] == 'findglare':
        continue
    if s_line[0] == 'VIEW=':
        continue
    if s_line[0] == 'FORMAT=ascii':
        continue
    if s_line[0] == '-60.000000':
        dft.append(float(s_line[1]))
        continue
    if s_line[0] == '-50.000000':
        dft.append(float(s_line[1]))
        continue
    if s_line[0] == '-40.000000':
        dft.append(float(s_line[1]))
        continue
    if s_line[0] == '-30.000000':
        dft.append(float(s_line[1]))
        continue
    if s_line[0] == '-20.000000':
        dft.append(float(s_line[1]))
        continue
    if s_line[0] == '-10.000000':
        dft.append(float(s_line[1]))
        continue
    if s_line[0] == '0.000000':
        dft.append(float(s_line[1]))
        continue
    if s_line[0] == '10.000000':
        dft.append(float(s_line[1]))
        continue
    if s_line[0] == '20.000000':
        dft.append(float(s_line[1]))
        continue
    if s_line[0] == '30.000000':
        dft.append(float(s_line[1]))
        continue
    if s_line[0] == '40.000000':
        dft.append(float(s_line[1]))
        continue
    if s_line[0] == '50.000000':
        dft.append(float(s_line[1]))
        continue
    if s_line[0] == '60.000000':
        dft.append(float(s_line[1]))
        guth_data=True
        found_guth=True

# A guth plot uses a polar plot constrained between -60 and 60 degrees. Convert the
# angles into radians into the array theta. Note: adjust the '80' below to better
# match the scale found from running findglare and glarendx on the Radiance model.
    if found_guth and guth_data:
        found_guth=False
        guth_data = False
        fig, axs = plt.subplots(subplot_kw={'projection': 'polar'})    # instanciate a new polar plot.
        axs.set_thetamin(-60)
        axs.set_thetamax(60)
        gx = np.arange(-60.0,70.0,10)
        theta = (np.pi/180.0 )*gx    # in radians
        l1 = np.array(dft)    # convert list to np array.
        axs.plot(theta, l1, c='blue', lw=2.5)
        axs.set_ylim(0,80)          # ADJUST as required
        axs.set_yticks(np.arange(0,80,10))
        if larger_fonts:
            plt.xlabel('Guth visual comfort probability', fontsize=12)
            axs.tick_params(axis='x', labelsize=11)
            axs.tick_params(axis='y', labelsize=11)
        else:
            plt.xlabel('Guth visual comfort probability')
        figname = 'guth_vcp' + '.png'
        plt.savefig(figname, format='png', bbox_inches='tight')
        plt.close()
        dft = []
        continue
        
