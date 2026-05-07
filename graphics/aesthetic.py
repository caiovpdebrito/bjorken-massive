#===================================================#
#======== GLOBAL AESTHETIC OPTIONS FOR PLOTS =======#
#===================================================#

import matplotlib.pyplot as plt

main_font = {'family': 'STIXGeneral', 'color': 'black', 'size': 22} # font parameters to be used in labels and text
tick_font = {'family': 'STIXGeneral', 'color': 'black', 'size': 16} # font parameters to be used in ticks

# General parameters for plots
plt.rcParams["figure.figsize"] = [10, 14]  # set the default figure size to the tuple 'figsize'
plt.rcParams['axes.labelsize'] = 30 # set the default font size of the axis
plt.rcParams['legend.fontsize'] = 35
plt.rcParams['xtick.labelsize'] = 30
plt.rcParams['ytick.labelsize'] = 30
plt.rcParams['lines.linewidth'] = 3
plt.rcParams.update({'font.size': 30})

plt.rcParams["legend.framealpha"] = 0 # perhaps exlude this option
plt.rcParams["legend.handlelength"] = 1.0

# LaTeX style font
plt.rcParams['mathtext.fontset'] = 'stix' 
plt.rcParams['font.family'] = 'STIXGeneral'

plt.rcParams['figure.titlesize'] = 27


FONTE = 20  # fontsize of the text inside the plots