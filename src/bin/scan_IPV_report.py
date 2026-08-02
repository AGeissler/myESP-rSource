"""Functions for importing and reading IPV report files"""
import sys, os
import datetime
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
    Reads in comma separated IPV report files.
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

            # Split line at commas and remove all whitespace
            line = [x.strip() for x in line.strip().split(",")]

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
s_ipvFile=sys.argv[1]
# s_ipvFile='ipv_3tweeks_4tsph.rep'
# s_ipvFile='Mid_floor_1965_tpl-5ftnipv.rep'
if not os.path.isfile(s_ipvFile): sys.exit('Error: input file "'+s_ipvFile+'" does not exist. Exiting.')

# Main program.

# looking for '*report  1 distribution' specifically '*report'  and  'distribution'
# AC logic is (I think) scan lines looking for specific key words, when found set a logic
# variable so that as subsequent lines are read the info needed can be decoded.
datacount=0     # which data line active
figurecount=0   # for numbering of figures generated
histo_count=0
demand_count=0
capacity_count=0
end_count=0
graph_count=0
comfort_count=0
dbt_count=0
solar_count=0
number_of_simulations=0
simulation_count=0
assessment_count=0
infiltration_count=0
casual_count=0
rh_count=0
aggregate_count=0
stats_row =0
seasonal_row=0
summary_energy_count=0
diversified_count=0
multidemand_count=0
emission_count=0
read_periods=False        # toggle which topic is active
focus_days=False
found_distribution=False
found_summary_distribution=False
found_data=False
found_diversified=False
found_demand=False
found_dispersed=False
found_daylight=False
found_guth=False
daylight_data=False
guth_data=False
found_graph=False
graph_data=False
found_demand_per_unit_time=False
found_delivered_per_unit_time=False
found_emissions=False
found_thermal_comfort=False
found_dbT=False
found_infiltration=False
found_solar=False
found_casual=False
found_rh=False
found_title=False
found_stats=False
found_seasonal=False
found_summary=False
found_power=False
found_energy=False
start_date  = []
finish_date = []
period_name = []
focus_date  = []
hits = []           # clear arrays holding data
hits_s1 = []        # clear arrays holding 1st season data
hits_s2 = []        # clear arrays holding 2nd season data
hits_s3 = []        # clear arrays holding 3rd season data
hits_s4 = []        # clear arrays holding 4th season data
hits_s5 = []        # clear arrays holding 5th season data
percents = []
emission_ht =[]
emission_cl =[]
emission_lt =[]
emission_fn =[]
emission_sp =[]
emission_hw =[]
power_cap   =[]
linelbl = []        # labels array 
line_1 = []         # supports up to 12 lines in a graph
line_2 = []
line_3 = []
line_4 = []
line_5 = []
line_6 = []
line_7 = []
line_8 = []
line_9 = []
line_10= []
line_11= []
line_12= []
line_13= []
line_14= []
line_15= []
t1 = []              # array for time in graphs
dft= []

# Define the cells of capacity or demand tables. 
cell_text = [[' ',' '],[' ',' '],[' ',' '],[' ',' '],[' ',' '],[' ',' '],[' ',' '],[' ',' ']]

# Allow up to 15 rows for the summary table.
overview_text = [[' ',' '],[' ',' '],[' ',' '],[' ',' '],
                 [' ',' '],[' ',' '],[' ',' '],[' ',' '],
                 [' ',' '],[' ',' '],[' ',' '],[' ',' ']]

# Define the cells of multiple season capacity tables. 
capacity_text = [['Season','within','Heating','Cooling','Lighting\nuncnctld','Lighting\ncontrolled','Fans','Small\nPower','DHW'],
                 [' ',' ',' ',' ',' ',' ',' ',' ',' '],
                 [' ',' ',' ',' ',' ',' ',' ',' ',' '],
                 [' ',' ',' ',' ',' ',' ',' ',' ',' '],
                 [' ',' ',' ',' ',' ',' ',' ',' ',' '],
                 [' ',' ',' ',' ',' ',' ',' ',' ',' '],
                 [' ',' ',' ',' ',' ',' ',' ',' ',' '],
                 [' ',' ',' ',' ',' ',' ',' ',' ',' '],
                 [' ',' ',' ',' ',' ',' ',' ',' ',' '],
                 [' ',' ',' ',' ',' ',' ',' ',' ',' '],
                 [' ',' ',' ',' ',' ',' ',' ',' ',' ']]

# Define the cells of multiple season demand tables. 
demand_text = [['Season','within','Heating','Cooling','Lighting\nuncnctld','Lighting\ncontrolled','Fans','Small\nPower','DHW'],
                 [' ',' ',' ',' ',' ',' ',' ',' ',' '],
                 [' ',' ',' ',' ',' ',' ',' ',' ',' '],
                 [' ',' ',' ',' ',' ',' ',' ',' ',' '],
                 [' ',' ',' ',' ',' ',' ',' ',' ',' '],
                 [' ',' ',' ',' ',' ',' ',' ',' ',' '],
                 [' ',' ',' ',' ',' ',' ',' ',' ',' '],
                 [' ',' ',' ',' ',' ',' ',' ',' ',' '],
                 [' ',' ',' ',' ',' ',' ',' ',' ',' '],
                 [' ',' ',' ',' ',' ',' ',' ',' ',' '],
                 [' ',' ',' ',' ',' ',' ',' ',' ',' '],
                 [' ',' ',' ',' ',' ',' ',' ',' ',' ']]

# Allow up to 32 rows for the stats table.
stats_text = [['Season','within','Topic','Maximum','Minum','Average'],
              [' ',' ',' ',' ',' ',' '],[' ',' ',' ',' ',' ',' '],
              [' ',' ',' ',' ',' ',' '],[' ',' ',' ',' ',' ',' '],
              [' ',' ',' ',' ',' ',' '],[' ',' ',' ',' ',' ',' '],
              [' ',' ',' ',' ',' ',' '],[' ',' ',' ',' ',' ',' '],
              [' ',' ',' ',' ',' ',' '],[' ',' ',' ',' ',' ',' '],
              [' ',' ',' ',' ',' ',' '],[' ',' ',' ',' ',' ',' '],
              [' ',' ',' ',' ',' ',' '],[' ',' ',' ',' ',' ',' '],
              [' ',' ',' ',' ',' ',' '],[' ',' ',' ',' ',' ',' '],
              [' ',' ',' ',' ',' ',' '],[' ',' ',' ',' ',' ',' '],
              [' ',' ',' ',' ',' ',' '],[' ',' ',' ',' ',' ',' '],
              [' ',' ',' ',' ',' ',' '],[' ',' ',' ',' ',' ',' '],
              [' ',' ',' ',' ',' ',' '],[' ',' ',' ',' ',' ',' '],
              [' ',' ',' ',' ',' ',' '],[' ',' ',' ',' ',' ',' '],
              [' ',' ',' ',' ',' ',' '],[' ',' ',' ',' ',' ',' '],
              [' ',' ',' ',' ',' ',' '],[' ',' ',' ',' ',' ',' '],
              [' ',' ',' ',' ',' ',' '],[' ',' ',' ',' ',' ',' ']]

# Allow up to 18 rows for the seasonal summary table.
seasonal_text = [['Season','Integrated\nDemand','Heating','Cooling','Lighting','Fans\netc','Small\nPower','DHW'],
              [' ',' ',' ',' ',' ',' ',' ',' '],
              [' ',' ',' ',' ',' ',' ',' ',' '],
              [' ',' ',' ',' ',' ',' ',' ',' '],
              [' ',' ',' ',' ',' ',' ',' ',' '],
              [' ',' ',' ',' ',' ',' ',' ',' '],
              [' ',' ',' ',' ',' ',' ',' ',' '],
              [' ',' ',' ',' ',' ',' ',' ',' '],
              [' ',' ',' ',' ',' ',' ',' ',' '],
              [' ',' ',' ',' ',' ',' ',' ',' '],
              [' ',' ',' ',' ',' ',' ',' ',' '],
              [' ',' ',' ',' ',' ',' ',' ',' '],
              [' ',' ',' ',' ',' ',' ',' ',' '],
              [' ',' ',' ',' ',' ',' ',' ',' ']]

# Sized for a 3 season assessment.
seasonal_text3 = [['Season','Integrated\nDemand','Heating','Cooling','Lighting','Fans\netc','Small\nPower','DHW'],
              [' ',' ',' ',' ',' ',' ',' ',' '],
              [' ',' ',' ',' ',' ',' ',' ',' '],
              [' ',' ',' ',' ',' ',' ',' ',' '],
              [' ',' ',' ',' ',' ',' ',' ',' '],
              [' ',' ',' ',' ',' ',' ',' ',' '],
              [' ',' ',' ',' ',' ',' ',' ',' '],
              [' ',' ',' ',' ',' ',' ',' ',' '],
              [' ',' ',' ',' ',' ',' ',' ',' '],
              [' ',' ',' ',' ',' ',' ',' ',' ']]


# Define the cells of annual demand tables. 
annual_demand_text = [['Topic','kWh/m^2.a','kWh.a'],
            ['Heating',' ',' '],
	    ['Cooling',' ',' '],
	    ['Lighting',' ',' '],
            ['Fans',' ',' '],
	    ['Small power',' ',' '],
	    ['DHW',' ',' ']]

introduction_a = 'This script scans an ESP-r IPV report file (generated from\ndirectives in a IPV definition via one or more assessments.'
introduction_b = 'It creates a series of time-value graphs and histograms and\ntables which can be manually merged into user reports.'
print (introduction_a)
print (introduction_b)
message = 'Proceed scanning ' + s_ipvFile + '(Y/N)?'
doit = ask_YorN(message)
if not doit:
    sys.exit('Exiting...')

# Scan the IPV report into an overall entity holding all the lines to work with.
IPV_rep = _read_file(s_ipvFile)

# Get user opinion as to frequency bin vertical axis.
as_hits = ask_YorN('Display frequency graphs as hits (Y) or percent (N) ?')

larger_fonts = ask_YorN('Use larger fonts in graphs?')

# Get user opinion on verbosity.
verbose = ask_YorN('Display data as it is scanned?')

# Loop through each line, find the number of tokens in the current line.
for s_line in IPV_rep:
    n_tokens = len(s_line)
    if verbose:
        print (s_line)

# Parse the heading of report e.g. the date weather site lat & long and synopsis phrase.
    if s_line[0] == '*title' and not found_title:
        project_title = s_line[1]
        found_title = True
        continue

    if s_line[0] == '*date':
        project_date = s_line[1]
        continue

    if s_line[0] == '*climate':
        project_weather = s_line[1]
        continue

    if s_line[0] == '*latitude':
        project_latitude = s_line[1]
        continue

    if s_line[0] == '*longitude':
        project_longitude = s_line[1]
        continue

    if s_line[0] == '*synopsis':
        project_synopsis = s_line[1]
        continue

# Remember the number of assessments set focus on periods for subsequent line reads.
    if s_line[0] == '*simulations':
        number_of_simulations= int(s_line[1])
        read_periods = True
        continue

# Set to scan for each assessment focus day & unfocus read of assessment periods.
    if s_line[0] == '*days':
        focus_days = True
        read_periods = False
        simulation_count=0
        continue

# If current focus is simulation periods.
    if read_periods:
        if simulation_count < number_of_simulations:
            start_date.append(s_line[0])
            finish_date.append(s_line[1])
            period_name.append(s_line[2])
            simulation_count=simulation_count +1
            continue
        else:
            read_periods = False
            simulation_count=0
            continue

# If current focus is focus days.
    if focus_days:
        if simulation_count < number_of_simulations:
            focus_date.append(s_line[0])
            simulation_count=simulation_count +1
            continue
        else:

# Sufficient information to generate an overview table.
            focus_days = False
            simulation_count=0
            overview_text[0][0] = 'Project '
            overview_text[0][1] =  project_title
            overview_text[1][0] = 'Date'
            overview_text[1][1] =  project_date
            overview_text[2][0] = 'Weather'
            overview_text[2][1] =  project_weather
            overview_text[3][0] = 'Location'
            overview_text[3][1] = 'Latitude ' + project_latitude + ' Longitude ' + project_longitude
            if number_of_simulations == 1:
                overview_text[4][0] =start_date[0] 
                overview_text[4][1] ='to ' + finish_date[0] + ' focus ' + focus_date[0]
            elif number_of_simulations == 2:
                overview_text[4][0] =start_date[0] 
                overview_text[4][1] ='to ' + finish_date[0] + ' focus ' + focus_date[0]
                overview_text[5][0] =start_date[1] 
                overview_text[5][1] ='to ' + finish_date[1] + ' focus ' + focus_date[1]
            elif number_of_simulations == 3:
                overview_text[4][0] =start_date[0] 
                overview_text[4][1] ='to ' + finish_date[0] + ' focus ' + focus_date[0]
                overview_text[5][0] =start_date[1] 
                overview_text[5][1] ='to ' + finish_date[1] + ' focus ' + focus_date[1]
                overview_text[6][0] =start_date[2] 
                overview_text[6][1] ='to ' + finish_date[2] + ' focus ' + focus_date[2]
            elif number_of_simulations == 4:
                overview_text[4][0] =start_date[0] 
                overview_text[4][1] ='to ' + finish_date[0] + ' focus ' + focus_date[0]
                overview_text[5][0] =start_date[1] 
                overview_text[5][1] ='to ' + finish_date[1] + ' focus ' + focus_date[1] 
                overview_text[6][0] =start_date[2] 
                overview_text[6][1] ='to ' + finish_date[2] + ' focus ' + focus_date[2]
                overview_text[7][0] =start_date[3] 
                overview_text[7][1] ='to ' + finish_date[3] + ' focus ' + focus_date[3]
            elif number_of_simulations == 5:
                overview_text[4][0] =start_date[0] 
                overview_text[4][1] ='to ' + finish_date[0] + ' focus ' + focus_date[0]
                overview_text[5][0] =start_date[1] 
                overview_text[5][1] ='to ' + finish_date[1] + ' focus ' + focus_date[1]
                overview_text[6][0] =start_date[2] 
                overview_text[6][1] ='to ' + finish_date[2] + ' focus ' + focus_date[2]
                overview_text[7][0] =start_date[3] 
                overview_text[7][1] ='to ' + finish_date[3] + ' focus ' + focus_date[3]
                overview_text[8][0] =start_date[4] 
                overview_text[8][1] ='to ' + finish_date[4] + ' focus ' + focus_date[4]
            if verbose:
                print (overview_text)
            if number_of_simulations == 1:
                fig, axs = plt.subplots(figsize=(5.0,2.5), layout='constrained')
            elif number_of_simulations == 3:
                fig, axs = plt.subplots(figsize=(5.0,3.0), layout='constrained')
            elif number_of_simulations == 5:
                fig, axs = plt.subplots(figsize=(5.2,3.5), layout='constrained')

            the_table = plt.table(cellText=overview_text, loc='center', edges='open')
            the_table.auto_set_column_width(0)
            the_table.auto_set_column_width(1)
            the_table.scale(1.3,1.4)
            if number_of_simulations == 1:
                plt.text(0.5, 0.6, project_synopsis, ha='center', va='top', wrap=True)
            elif number_of_simulations == 2:
                plt.text(0.5, 0.55, project_synopsis, ha='center', va='top', wrap=True)
            elif number_of_simulations == 3:
                plt.text(0.5, 0.35, project_synopsis, ha='center', va='top', wrap=True)
            elif number_of_simulations == 4:
                plt.text(0.5, 0.30, project_synopsis, ha='center', va='top', wrap=True)
            elif number_of_simulations == 5:
                plt.text(0.5, 0.20, project_synopsis, ha='center', va='top', wrap=True)
            plt.axis('off')
            figname = 'overview' + '.png'
#             plt.savefig(figname, format='png', bbox_inches='tight')
            plt.savefig(figname, format='png')
            plt.close()            
            print ('created ' + figname)
            continue


# Graph of daylight factors.
    if found_daylight:
        if s_line[0] == '*title':
            graph_title = s_line[1]     # remember title for graph
            units = s_line[2]
            continue
        elif s_line[0] == '*format':
            steps = int(s_line[2])
            lines = int(s_line[3])
            continue
        elif s_line[0] == '*fields':
            if lines == 2:
              linelbl.append(s_line[2])
            elif lines == 3:
              linelbl.append(s_line[2])
              linelbl.append(s_line[3])
        elif s_line[0] == '*data':
            daylight_data=True
            continue
        elif s_line[0] == '*end_report': # have all information needed
            found_daylight=False
            daylight_data = False
            fig, axs = plt.subplots()    # instanciate a new plot.
            if lines == 2:
                l1 = np.array(line_1)    # convert list to np array.
                plt.plot(dft, l1, c='red', label=linelbl[0])
            elif lines == 3:
                l1 = np.array(line_1)
                l2 = np.array(line_2)
                plt.plot(dft, l1, c='red', label=linelbl[0])
                plt.plot(dft, l2, c='blue', label=linelbl[1])
            if larger_fonts:
                plt.ylabel(units, fontsize=12)
                axs.tick_params(axis='x', labelsize=11)
                axs.tick_params(axis='y', labelsize=11)
                plt.xlabel('Daylight factors - distance from facade (m)', fontsize=12)
            else:
                plt.ylabel(units)
                plt.xlabel('Daylight factors - distance from facade (m)')
            figname = 'daylight_factors' + str(demand_count) + '.png'
            if larger_fonts:
                plt.legend(loc='best', fontsize=11)
            else:
                plt.legend(loc='best')
            plt.grid(axis='x', color='0.6')
            plt.grid(axis='y', color='0.6')
            plt.savefig(figname, format='png', bbox_inches='tight')
            plt.close()
            linelbl = []  # Clear arrays in case of subsequent graph.
            line_1 = []
            line_2 = []
            dft = []
            continue


# Extract daylight foctor data for graphing.
    if found_daylight and daylight_data:
        if lines == 2:
            dft.append(float(s_line[0]))
            line_1.append(float(s_line[1]))     # one data to plot
            continue
        elif lines == 3:
            dft.append(float(s_line[0]))
            line_1.append(float(s_line[1]))     # two data to plot
            line_2.append(float(s_line[2]))
            continue

# Graph of Guth glare indicators.
    if found_guth:
        if s_line[0] == '*title':
            graph_title = s_line[1]     # remember title for graph
            units = s_line[2]
            continue
        elif s_line[0] == '*format':
            steps = int(s_line[2])
            lines = int(s_line[3])
            continue
        elif s_line[0] == '*fields':
            linelbl.append(s_line[2])
        elif s_line[0] == '*data':
            guth_data=True
            continue
        elif s_line[0] == '*end_report': # have all information needed

# A guth plot uses a polar plot constrained between -60 and 60 degrees. Convert the
# angles into radians into the array theta. Note: adjust the '80' below to better
# match the scale found from running findglare and glarendx on the Radiance model.
            found_guth=False
            guth_data = False
            fig, axs = plt.subplots(subplot_kw={'projection': 'polar'})    # instanciate a new polar plot.
            axs.set_thetamin(-60)
            axs.set_thetamax(60)
            gx = np.arange(-60.0,70.0,10)
            theta = (np.pi/180.0 )*gx    # in radians
            l1 = np.array(line_1)    # convert list to np array.
            axs.plot(theta, l1, c='blue', lw=2.5)
            axs.set_ylim(0,80)          # ADJUST as required
            axs.set_yticks(np.arange(0,80,10))
            if larger_fonts:
                plt.xlabel('Guth visual comfort probability', fontsize=12)
                axs.tick_params(axis='x', labelsize=11)
                axs.tick_params(axis='y', labelsize=11)
            else:
                plt.xlabel('Guth visual comfort probability')
            figname = 'guth_factors' + str(demand_count) + '.png'
#            plt.legend(loc='best')
#            plt.grid(axis='x', color='0.6')
#            plt.grid(axis='y', color='0.6')
            plt.savefig(figname, format='png', bbox_inches='tight')
            plt.close()
            linelbl = []  # Clear arrays in case of subsequent graph.
            line_1 = []
            dft = []
            continue


# Extract daylight foctor data for graphing.
    if found_guth and guth_data:
        if lines == 2:
            dft.append(float(s_line[0]))
            line_1.append(float(s_line[1]))     # one data to plot
            continue


# Graph of energy demands, thermal comfort, solar or zone dbT.
    if found_graph:
        if s_line[0] == '*title':
            graph_title = s_line[1]     # remember title for graph
            units = s_line[2]
            continue
        elif s_line[0] == '*format':
            steps = int(s_line[2])
            lines = int(s_line[3])
            continue
        elif s_line[0] == '*fields':
            if lines == 2:
              linelbl.append(s_line[2])
            elif lines == 3:
              linelbl.append(s_line[2])
              linelbl.append(s_line[3])
            elif lines == 4:
              linelbl.append(s_line[2])
              linelbl.append(s_line[3])
              linelbl.append(s_line[4])
            elif lines == 5:
              linelbl.append(s_line[2])
              linelbl.append(s_line[3])
              linelbl.append(s_line[4])
              linelbl.append(s_line[5])
            elif lines == 6:
              linelbl.append(s_line[2])
              linelbl.append(s_line[3])
              linelbl.append(s_line[4])
              linelbl.append(s_line[5])
              linelbl.append(s_line[6])
            elif lines == 7:
              linelbl.append(s_line[2])
              linelbl.append(s_line[3])
              linelbl.append(s_line[4])
              linelbl.append(s_line[5])
              linelbl.append(s_line[6])
              linelbl.append(s_line[7])
            elif lines == 8:
              linelbl.append(s_line[2])
              linelbl.append(s_line[3])
              linelbl.append(s_line[4])
              linelbl.append(s_line[5])
              linelbl.append(s_line[6])
              linelbl.append(s_line[7])
              linelbl.append(s_line[8])
            elif lines == 9:
              linelbl.append(s_line[2])
              linelbl.append(s_line[3])
              linelbl.append(s_line[4])
              linelbl.append(s_line[5])
              linelbl.append(s_line[6])
              linelbl.append(s_line[7])
              linelbl.append(s_line[8])
              linelbl.append(s_line[9])
            elif lines == 10:
              linelbl.append(s_line[2])
              linelbl.append(s_line[3])
              linelbl.append(s_line[4])
              linelbl.append(s_line[5])
              linelbl.append(s_line[6])
              linelbl.append(s_line[7])
              linelbl.append(s_line[8])
              linelbl.append(s_line[9])
              linelbl.append(s_line[10])
            elif lines == 11:
              linelbl.append(s_line[2])
              linelbl.append(s_line[3])
              linelbl.append(s_line[4])
              linelbl.append(s_line[5])
              linelbl.append(s_line[6])
              linelbl.append(s_line[7])
              linelbl.append(s_line[8])
              linelbl.append(s_line[9])
              linelbl.append(s_line[10])
              linelbl.append(s_line[11])
            elif lines == 12:
              linelbl.append(s_line[2])
              linelbl.append(s_line[3])
              linelbl.append(s_line[4])
              linelbl.append(s_line[5])
              linelbl.append(s_line[6])
              linelbl.append(s_line[7])
              linelbl.append(s_line[8])
              linelbl.append(s_line[9])
              linelbl.append(s_line[10])
              linelbl.append(s_line[11])
              linelbl.append(s_line[12])
            elif lines == 13:
              linelbl.append(s_line[2])
              linelbl.append(s_line[3])
              linelbl.append(s_line[4])
              linelbl.append(s_line[5])
              linelbl.append(s_line[6])
              linelbl.append(s_line[7])
              linelbl.append(s_line[8])
              linelbl.append(s_line[9])
              linelbl.append(s_line[10])
              linelbl.append(s_line[11])
              linelbl.append(s_line[12])
              linelbl.append(s_line[13])
            elif lines == 14:
              linelbl.append(s_line[2])
              linelbl.append(s_line[3])
              linelbl.append(s_line[4])
              linelbl.append(s_line[5])
              linelbl.append(s_line[6])
              linelbl.append(s_line[7])
              linelbl.append(s_line[8])
              linelbl.append(s_line[9])
              linelbl.append(s_line[10])
              linelbl.append(s_line[11])
              linelbl.append(s_line[12])
              linelbl.append(s_line[13])
              linelbl.append(s_line[14])
            if verbose:
                print (linelbl)
            continue
        elif s_line[0] == '*data':
            graph_data=True
            continue
        elif s_line[0] == '*end_report': # have all information needed
            if verbose:
                print (linelbl)
            found_graph=False
            graph_data = False
            start = t1[0]                # get time of X start and finish
            finish = t1[steps-1]
            fig, axs = plt.subplots()    # instanciate a new plot.
            if lines == 2:
                l1 = np.array(line_1)    # convert list to np array.
                plt.plot(t1, l1, c='red', label=linelbl[0])
            elif lines == 3:
                l1 = np.array(line_1)
                l2 = np.array(line_2)
                plt.plot(t1, l1, c='red', label=linelbl[0])
                plt.plot(t1, l2, c='blue', label=linelbl[1])
            elif lines == 4:
                l1 = np.array(line_1)
                l2 = np.array(line_2)
                l3 = np.array(line_3)
                plt.plot(t1, l1, c='red', label=linelbl[0])
                plt.plot(t1, l2, c='blue', label=linelbl[1])
                plt.plot(t1, l3, c='orange', label=linelbl[2])
            elif lines == 5:
                l1 = np.array(line_1)
                l2 = np.array(line_2)
                l3 = np.array(line_3)
                l4 = np.array(line_4)
                plt.plot(t1, l1, c='red', label=linelbl[0])
                plt.plot(t1, l2, c='blue', label=linelbl[1])
                plt.plot(t1, l3, c='orange', label=linelbl[2])
                plt.plot(t1, l4, c='olivedrab', label=linelbl[3])
            elif lines == 6:
                l1 = np.array(line_1)
                l2 = np.array(line_2)
                l3 = np.array(line_3)
                l4 = np.array(line_4)
                l5 = np.array(line_5)
                plt.plot(t1, l1, c='red', label=linelbl[0])
                plt.plot(t1, l2, c='blue', label=linelbl[1])
                plt.plot(t1, l3, c='orange', label=linelbl[2])
                plt.plot(t1, l4, c='olivedrab', label=linelbl[3])
                plt.plot(t1, l5, c='aquamarine', label=linelbl[4])
            elif lines == 7:
                if found_demand_per_unit_time:  # Use a stack plot for this topic
                    l1 = np.array(line_1)
                    l2 = np.array(line_2)
                    l3 = np.array(line_3)
                    l4 = np.array(line_4)
                    l5 = np.array(line_5)
                    l6 = np.array(line_6)
                    plt.plot([], [], c='red', label=linelbl[0])
                    plt.plot([], [], c='blue', label=linelbl[1])
                    plt.plot([], [], c='orange', label=linelbl[2])
                    plt.plot([], [], c='olivedrab', label=linelbl[3])
                    plt.plot([], [], c='aquamarine', label=linelbl[4])
                    plt.plot([], [], c='silver', label=linelbl[5])
                    plt.stackplot(t1, l1, l2, l3, l4, l5, l6, linewidth=0.7, colors=['red', 'blue', 'orange', 'olivedrab', 'aquamarine', 'silver'])

# There may be a sequence of focus days in different seasons. If you want to set a
# consistent Y axis uncomment and edit the following line:
#                    axs.set(ylim=(0,1800))
                elif found_delivered_per_unit_time:
                    l1 = np.array(line_1)
                    l2 = np.array(line_2)
                    l3 = np.array(line_3)
                    l4 = np.array(line_4)
                    l5 = np.array(line_5)
                    l6 = np.array(line_6)
                    plt.plot([], [], c='red', label=linelbl[0])
                    plt.plot([], [], c='blue', label=linelbl[1])
                    plt.plot([], [], c='orange', label=linelbl[2])
                    plt.plot([], [], c='olivedrab', label=linelbl[3])
                    plt.plot([], [], c='aquamarine', label=linelbl[4])
                    plt.plot([], [], c='silver', label=linelbl[5])
                    plt.stackplot(t1, l1, l2, l3, l4, l5, l6, linewidth=0.7, colors=['red', 'blue', 'orange', 'olivedrab', 'aquamarine', 'silver'])

# There may be a sequence of focus days in different seasons. If you want to set a
# consistent Y axis uncomment and edit the following line:
#                    axs.set(ylim=(0,90))
                else:
                    l1 = np.array(line_1)
                    l2 = np.array(line_2)
                    l3 = np.array(line_3)
                    l4 = np.array(line_4)
                    l5 = np.array(line_5)
                    l6 = np.array(line_6)
                    plt.plot(t1, l1, c='red', label=linelbl[0])
                    plt.plot(t1, l2, c='blue', label=linelbl[1])
                    plt.plot(t1, l3, c='orange', label=linelbl[2])
                    plt.plot(t1, l4, c='olivedrab', label=linelbl[3])
                    plt.plot(t1, l5, c='aquamarine', label=linelbl[4])
                    plt.plot(t1, l6, c='silver', label=linelbl[5])
            elif lines == 8:
                l1 = np.array(line_1)
                l2 = np.array(line_2)
                l3 = np.array(line_3)
                l4 = np.array(line_4)
                l5 = np.array(line_5)
                l6 = np.array(line_6)
                l7 = np.array(line_7)
                plt.plot(t1, l1, c='red', label=linelbl[0])
                plt.plot(t1, l2, c='blue', label=linelbl[1])
                plt.plot(t1, l3, c='orange', label=linelbl[2])
                plt.plot(t1, l4, c='olivedrab', label=linelbl[3])
                plt.plot(t1, l5, c='aquamarine', label=linelbl[4])
                plt.plot(t1, l6, c='silver', label=linelbl[5])
                plt.plot(t1, l7, c='chocolate', label=linelbl[6])
            elif lines == 9:
                l1 = np.array(line_1)
                l2 = np.array(line_2)
                l3 = np.array(line_3)
                l4 = np.array(line_4)
                l5 = np.array(line_5)
                l6 = np.array(line_6)
                l7 = np.array(line_7)
                l8 = np.array(line_8)
                plt.plot(t1, l1, c='red', label=linelbl[0])
                plt.plot(t1, l2, c='blue', label=linelbl[1])
                plt.plot(t1, l3, c='orange', label=linelbl[2])
                plt.plot(t1, l4, c='olivedrab', label=linelbl[3])
                plt.plot(t1, l5, c='aquamarine', label=linelbl[4])
                plt.plot(t1, l6, c='silver', label=linelbl[5])
                plt.plot(t1, l7, c='chocolate', label=linelbl[6])
                plt.plot(t1, l8, c='goldenrod', label=linelbl[7])
            elif lines == 10:
                l1 = np.array(line_1)
                l2 = np.array(line_2)
                l3 = np.array(line_3)
                l4 = np.array(line_4)
                l5 = np.array(line_5)
                l6 = np.array(line_6)
                l7 = np.array(line_7)
                l8 = np.array(line_8)
                l9 = np.array(line_9)
                plt.plot(t1, l1, c='red', label=linelbl[0])
                plt.plot(t1, l2, c='blue', label=linelbl[1])
                plt.plot(t1, l3, c='orange', label=linelbl[2])
                plt.plot(t1, l4, c='olivedrab', label=linelbl[3])
                plt.plot(t1, l5, c='aquamarine', label=linelbl[4])
                plt.plot(t1, l6, c='silver', label=linelbl[5])
                plt.plot(t1, l7, c='chocolate', label=linelbl[6])
                plt.plot(t1, l8, c='goldenrod', label=linelbl[7])
                plt.plot(t1, l9, c='salmon', label=linelbl[8])
            elif lines == 11:
                l1 = np.array(line_1)
                l2 = np.array(line_2)
                l3 = np.array(line_3)
                l4 = np.array(line_4)
                l5 = np.array(line_5)
                l6 = np.array(line_6)
                l7 = np.array(line_7)
                l8 = np.array(line_8)
                l9 = np.array(line_9)
                l10= np.array(line_10)
                plt.plot(t1, l1, c='red', label=linelbl[0])
                plt.plot(t1, l2, c='blue', label=linelbl[1])
                plt.plot(t1, l3, c='orange', label=linelbl[2])
                plt.plot(t1, l4, c='olivedrab', label=linelbl[3])
                plt.plot(t1, l5, c='aquamarine', label=linelbl[4])
                plt.plot(t1, l6, c='silver', label=linelbl[5])
                plt.plot(t1, l7, c='chocolate', label=linelbl[6])
                plt.plot(t1, l8, c='goldenrod', label=linelbl[7])
                plt.plot(t1, l9, c='salmon', label=linelbl[8])
                plt.plot(t1, l10, c='limegreen', label=linelbl[9])
            elif lines == 12:
                l1 = np.array(line_1)
                l2 = np.array(line_2)
                l3 = np.array(line_3)
                l4 = np.array(line_4)
                l5 = np.array(line_5)
                l6 = np.array(line_6)
                l7 = np.array(line_7)
                l8 = np.array(line_8)
                l9 = np.array(line_9)
                l10= np.array(line_10)
                l11= np.array(line_11)
                plt.plot(t1, l1, c='red', label=linelbl[0])
                plt.plot(t1, l2, c='blue', label=linelbl[1])
                plt.plot(t1, l3, c='orange', label=linelbl[2])
                plt.plot(t1, l4, c='olivedrab', label=linelbl[3])
                plt.plot(t1, l5, c='aquamarine', label=linelbl[4])
                plt.plot(t1, l6, c='silver', label=linelbl[5])
                plt.plot(t1, l7, c='chocolate', label=linelbl[6])
                plt.plot(t1, l8, c='goldenrod', label=linelbl[7])
                plt.plot(t1, l9, c='salmon', label=linelbl[8])
                plt.plot(t1, l10, c='limegreen', label=linelbl[9])
                plt.plot(t1, l11, c='orchid', label=linelbl[10])
            elif lines == 13:
                l1 = np.array(line_1)
                l2 = np.array(line_2)
                l3 = np.array(line_3)
                l4 = np.array(line_4)
                l5 = np.array(line_5)
                l6 = np.array(line_6)
                l7 = np.array(line_7)
                l8 = np.array(line_8)
                l9 = np.array(line_9)
                l10= np.array(line_10)
                l11= np.array(line_11)
                l12= np.array(line_12)
                plt.plot(t1, l1, c='red', label=linelbl[0])
                plt.plot(t1, l2, c='blue', label=linelbl[1])
                plt.plot(t1, l3, c='orange', label=linelbl[2])
                plt.plot(t1, l4, c='olivedrab', label=linelbl[3])
                plt.plot(t1, l5, c='aquamarine', label=linelbl[4])
                plt.plot(t1, l6, c='silver', label=linelbl[5])
                plt.plot(t1, l7, c='chocolate', label=linelbl[6])
                plt.plot(t1, l8, c='goldenrod', label=linelbl[7])
                plt.plot(t1, l9, c='salmon', label=linelbl[8])
                plt.plot(t1, l10, c='limegreen', label=linelbl[9])
                plt.plot(t1, l11, c='orchid', label=linelbl[10])
                plt.plot(t1, l12, c='chocolate', label=linelbl[11])
            elif lines == 14:
                l1 = np.array(line_1)
                l2 = np.array(line_2)
                l3 = np.array(line_3)
                l4 = np.array(line_4)
                l5 = np.array(line_5)
                l6 = np.array(line_6)
                l7 = np.array(line_7)
                l8 = np.array(line_8)
                l9 = np.array(line_9)
                l10= np.array(line_10)
                l11= np.array(line_11)
                l12= np.array(line_12)
                l13= np.array(line_13)
                plt.plot(t1, l1, c='red', label=linelbl[0])
                plt.plot(t1, l2, c='blue', label=linelbl[1])
                plt.plot(t1, l3, c='orange', label=linelbl[2])
                plt.plot(t1, l4, c='olivedrab', label=linelbl[3])
                plt.plot(t1, l5, c='aquamarine', label=linelbl[4])
                plt.plot(t1, l6, c='silver', label=linelbl[5])
                plt.plot(t1, l7, c='chocolate', label=linelbl[6])
                plt.plot(t1, l8, c='goldenrod', label=linelbl[7])
                plt.plot(t1, l9, c='salmon', label=linelbl[8])
                plt.plot(t1, l10, c='limegreen', label=linelbl[9])
                plt.plot(t1, l11, c='orchid', label=linelbl[10])
                plt.plot(t1, l12, c='chocolate', label=linelbl[11])
                plt.plot(t1, l13, c='orchid', label=linelbl[12])
            else:      # catch-all in case of more than 14 lines.
                l1 = np.array(line_1)
                l2 = np.array(line_2)
                l3 = np.array(line_3)
                l4 = np.array(line_4)
                l5 = np.array(line_5)
                l6 = np.array(line_6)
                l7 = np.array(line_7)
                l8 = np.array(line_8)
                l9 = np.array(line_9)
                l10= np.array(line_10)
                l11= np.array(line_11)
                l12= np.array(line_12)
                l13= np.array(line_13)
                plt.plot(t1, l1, c='red', label=linelbl[0])
                plt.plot(t1, l2, c='blue', label=linelbl[1])
                plt.plot(t1, l3, c='orange', label=linelbl[2])
                plt.plot(t1, l4, c='olivedrab', label=linelbl[3])
                plt.plot(t1, l5, c='aquamarine', label=linelbl[4])
                plt.plot(t1, l6, c='silver', label=linelbl[5])
                plt.plot(t1, l7, c='chocolate', label=linelbl[6])
                plt.plot(t1, l8, c='cyan', label=linelbl[7])
                plt.plot(t1, l9, c='teal', label=linelbl[8])
                plt.plot(t1, l10, c='limegreen', label=linelbl[9])
                plt.plot(t1, l11, c='orchid', label=linelbl[10])
                plt.plot(t1, l12, c='darkgrey', label=linelbl[11])
                plt.plot(t1, l13, c='orchid', label=linelbl[12])
 
# further colours: dimgrey firebrick salmon sandybrown bisque khaki
# palegreen limegreen lightblue cornflowerblue darkviolet
# orchid chocolate
# Based on current focus set the X axis label.
            if int(start) == int(finish):
                axs.set_xlim(int(start),int(start+1))
                if larger_fonts:
                    plt.xlabel( context + ' ' + tbl_title + ' focus day', fontsize=12)
                    axs.tick_params(axis='x', labelsize=11)
                    axs.tick_params(axis='y', labelsize=11)
                else:
                    plt.xlabel( context + ' ' + tbl_title + ' focus day')
            else:
                axs.set_xlim(int(start),finish)
                if larger_fonts:
                    plt.xlabel( context + ' ' + tbl_title + ' day(s)', fontsize=12)
                    axs.tick_params(axis='x', labelsize=11)
                    axs.tick_params(axis='y', labelsize=11)
                else:
                    plt.xlabel( context + ' ' + tbl_title + ' day(s)')

            if found_demand_per_unit_time:
                if larger_fonts:
                    plt.ylabel(graph_title + ' ' + units, fontsize=12)
                    axs.tick_params(axis='x', labelsize=11)
                    axs.tick_params(axis='y', labelsize=11)
                else:
                    plt.ylabel(graph_title + ' ' + units)
                figname = 'demand_gr' + str(demand_count) + '.png'

            if found_delivered_per_unit_time:
                if larger_fonts:
                    plt.xlabel( tbl_title + ' day(s)', fontsize=12)
                    plt.ylabel(graph_title + ' ' + units, fontsize=12)
                    axs.tick_params(axis='x', labelsize=11)
                    axs.tick_params(axis='y', labelsize=11)
                else:
                    plt.xlabel( tbl_title + ' day(s)')
                    plt.ylabel(graph_title + ' ' + units)
                figname = 'delivered_time_gr' + str(aggregate_count) + '.png'

            if found_thermal_comfort:
                if larger_fonts:
                    plt.ylabel(graph_title + ' ' + units, fontsize=12)
                    axs.tick_params(axis='x', labelsize=11)
                    axs.tick_params(axis='y', labelsize=11)
                else:
                    plt.ylabel(graph_title + ' ' + units)
                figname = 'comfort_gr' + str(comfort_count) + '.png'

            if found_dbT:
                if larger_fonts:
                    plt.ylabel(graph_title + ' ' + units, fontsize=12)
                    axs.tick_params(axis='x', labelsize=11)
                    axs.tick_params(axis='y', labelsize=11)
                else:
                    plt.ylabel(graph_title + ' ' + units)
                figname = 'dbt_gr' + str(dbt_count) + '.png'

            if found_infiltration:
                if larger_fonts:
                    plt.ylabel(graph_title + ' ' + units, fontsize=12)
                    axs.tick_params(axis='x', labelsize=11)
                    axs.tick_params(axis='y', labelsize=11)
                else:
                    plt.ylabel(graph_title + ' ' + units)
                figname = 'infiltration_gr' + str(infiltration_count) + '.png'

            if found_casual:
                if larger_fonts:
                    plt.ylabel(graph_title + ' ' + units, fontsize=12)
                    axs.tick_params(axis='x', labelsize=11)
                    axs.tick_params(axis='y', labelsize=11)
                else:
                    plt.ylabel(graph_title + ' ' + units)
                figname = 'casual_gr' + str(casual_count) + '.png'

            if found_rh:
                if larger_fonts:
                    plt.ylabel(graph_title + ' ' + units, fontsize=12)
                    axs.tick_params(axis='x', labelsize=11)
                    axs.tick_params(axis='y', labelsize=11)
                else:
                    plt.ylabel(graph_title + ' ' + units)
                figname = 'rh_gr' + str(rh_count) + '.png'

            if found_solar:
                if larger_fonts:
                    plt.ylabel(graph_title + ' ' + units, fontsize=12)
                    axs.tick_params(axis='x', labelsize=11)
                    axs.tick_params(axis='y', labelsize=11)
                else:
                    plt.ylabel(graph_title + ' ' + units)
                figname = 'solar_gr' + str(solar_count) + '.png'


# Auto placement of the legend.
            if larger_fonts:
                plt.legend(loc='best', fontsize=11)
            else:
                plt.legend(loc='best')
            plt.grid(axis='x', color='0.6')
            plt.grid(axis='y', color='0.6')

            if found_demand_per_unit_time:
                graph_count = graph_count + 1
                found_demand_per_unit_time=False

# Increment counts for later use in deciding which tables to generate at
# the end of the scanning.
            if found_thermal_comfort:
                comfort_count = comfort_count + 1
                found_thermal_comfort=False

            if found_dbT:
                dbt_count = dbt_count + 1
                found_dbT=False

            if found_infiltration:
                infiltration_count = infiltration_count + 1
                found_infiltration=False

            if found_delivered_per_unit_time:
                aggregate_count = aggregate_count + 1
                found_delivered_per_unit_time=False

            if found_casual:
                casual_count = casual_count + 1
                found_casual=False

            if found_rh:
                rh_count = rh_count + 1
                found_rh=False

            if found_solar:
                solar_count = solar_count + 1
                found_solar=False
                
            plt.savefig(figname, format='png', bbox_inches='tight')
            plt.close()
#             plt.show
            print ('created ' + figname)
            graph_demand=False
            linelbl = []  # Clear arrays in case of subsequent graph.
            line_1 = []
            line_2 = []
            line_3 = []
            line_4 = []
            line_5 = []
            line_6 = []
            line_7 = []
            line_8 = []
            line_9 = []
            line_10= []
            line_11= []
            line_12= []
            line_13= []
            line_14= []
            line_15= []
            t1 = []
            continue

# Extract timestep data for energy demands.
    if found_graph and graph_data:
        if lines == 2:
            t1.append(float(s_line[0]))
            line_1.append(float(s_line[1]))     # one data to plot
            continue
        elif lines == 3:
            t1.append(float(s_line[0]))
            line_1.append(float(s_line[1]))     # two data to plot
            line_2.append(float(s_line[2]))
            continue
        elif lines == 4:
            t1.append(float(s_line[0]))
            line_1.append(float(s_line[1]))     # 3 data to plot
            line_2.append(float(s_line[2]))
            line_3.append(float(s_line[3]))
            continue
        elif lines == 5:
            t1.append(float(s_line[0]))
            line_1.append(float(s_line[1]))     # 4 data to plot
            line_2.append(float(s_line[2]))
            line_3.append(float(s_line[3]))
            line_4.append(float(s_line[4]))
            continue
        elif lines == 6:
            t1.append(float(s_line[0]))
            line_1.append(float(s_line[1]))     # 5 data to plot
            line_2.append(float(s_line[2]))
            line_3.append(float(s_line[3]))
            line_4.append(float(s_line[4]))
            line_5.append(float(s_line[5]))
            continue
        elif lines == 7:
            t1.append(float(s_line[0]))
            line_1.append(float(s_line[1]))     # 6 data to plot
            line_2.append(float(s_line[2]))
            line_3.append(float(s_line[3]))
            line_4.append(float(s_line[4]))
            line_5.append(float(s_line[5]))
            line_6.append(float(s_line[6]))
            continue
        elif lines == 8:
            t1.append(float(s_line[0]))
            line_1.append(float(s_line[1]))     # 7 data to plot
            line_2.append(float(s_line[2]))
            line_3.append(float(s_line[3]))
            line_4.append(float(s_line[4]))
            line_5.append(float(s_line[5]))
            line_6.append(float(s_line[6]))
            line_7.append(float(s_line[7]))
            continue
        elif lines == 9:
            t1.append(float(s_line[0]))
            line_1.append(float(s_line[1]))     # 8 data to plot
            line_2.append(float(s_line[2]))
            line_3.append(float(s_line[3]))
            line_4.append(float(s_line[4]))
            line_5.append(float(s_line[5]))
            line_6.append(float(s_line[6]))
            line_7.append(float(s_line[7]))
            line_8.append(float(s_line[8]))
            continue
        elif lines == 10:
            t1.append(float(s_line[0]))
            line_1.append(float(s_line[1]))     # 9 data to plot
            line_2.append(float(s_line[2]))
            line_3.append(float(s_line[3]))
            line_4.append(float(s_line[4]))
            line_5.append(float(s_line[5]))
            line_6.append(float(s_line[6]))
            line_7.append(float(s_line[7]))
            line_8.append(float(s_line[8]))
            line_9.append(float(s_line[9]))
            continue
        elif lines == 11:
            t1.append(float(s_line[0]))
            line_1.append(float(s_line[1]))     # 10 data to plot
            line_2.append(float(s_line[2]))
            line_3.append(float(s_line[3]))
            line_4.append(float(s_line[4]))
            line_5.append(float(s_line[5]))
            line_6.append(float(s_line[6]))
            line_7.append(float(s_line[7]))
            line_8.append(float(s_line[8]))
            line_9.append(float(s_line[9]))
            line_10.append(float(s_line[10]))
            continue
        elif lines == 12:
            t1.append(float(s_line[0]))
            line_1.append(float(s_line[1]))     # 11 data to plot
            line_2.append(float(s_line[2]))
            line_3.append(float(s_line[3]))
            line_4.append(float(s_line[4]))
            line_5.append(float(s_line[5]))
            line_6.append(float(s_line[6]))
            line_7.append(float(s_line[7]))
            line_8.append(float(s_line[8]))
            line_9.append(float(s_line[9]))
            line_10.append(float(s_line[10]))
            line_11.append(float(s_line[11]))
            continue
        elif lines == 13:
            t1.append(float(s_line[0]))
            line_1.append(float(s_line[1]))     # 12 data to plot
            line_2.append(float(s_line[2]))
            line_3.append(float(s_line[3]))
            line_4.append(float(s_line[4]))
            line_5.append(float(s_line[5]))
            line_6.append(float(s_line[6]))
            line_7.append(float(s_line[7]))
            line_8.append(float(s_line[8]))
            line_9.append(float(s_line[9]))
            line_10.append(float(s_line[10]))
            line_11.append(float(s_line[11]))
            line_12.append(float(s_line[12]))
            continue
        else:
            t1.append(float(s_line[0]))
            line_1.append(float(s_line[1]))     # catch-out for more lines.
            line_2.append(float(s_line[2]))
            line_3.append(float(s_line[3]))
            line_4.append(float(s_line[4]))
            line_5.append(float(s_line[5]))
            line_6.append(float(s_line[6]))
            line_7.append(float(s_line[7]))
            line_8.append(float(s_line[8]))
            line_9.append(float(s_line[9]))
            line_10.append(float(s_line[10]))
            line_11.append(float(s_line[11]))
            line_12.append(float(s_line[12]))
            continue

# Table of diversified capacity. Based on 1st token on the line grab
# relevant information.             
    if found_diversified:
        if s_line[0] == '*title':
            capacity_title = s_line[1]     # remember title for graph
            units = s_line[2]
        elif s_line[0] == '*format':       # ignore
            continue
        elif s_line[0] == '*fields':
            cell_text[0][0] = 'Within ' + context + '\n' + tbl_title
            cell_text[1][0] = s_line[1]
            cell_text[2][0] = s_line[2]
            cell_text[3][0] = s_line[3]
            cell_text[4][0] = s_line[4]
            cell_text[5][0] = s_line[5]
            cell_text[6][0] = s_line[6]
            cell_text[7][0] = s_line[7]
            continue
        elif s_line[0] == '*data':

# Put values in cell_text column and also fill capacity_text array for use later.
            cell_text[0][1] = 'Capacity W'
            cell_text[1][1] = s_line[1]
            cell_text[2][1] = s_line[2]
            cell_text[3][1] = s_line[3]
            cell_text[4][1] = s_line[4]
            cell_text[5][1] = s_line[5]
            cell_text[6][1] = s_line[6]
            cell_text[7][1] = s_line[7]
            capacity_text[diversified_count][0] = season_title
            capacity_text[diversified_count][1] = context
            capacity_text[diversified_count][2] = s_line[1]
            capacity_text[diversified_count][3] = s_line[2]
            capacity_text[diversified_count][4] = s_line[3]
            capacity_text[diversified_count][5] = s_line[4]
            capacity_text[diversified_count][6] = s_line[5]
            capacity_text[diversified_count][7] = s_line[6]
            capacity_text[diversified_count][8] = s_line[7]
            continue
        elif s_line[0] == '*end_report': # have all information needed
            if verbose:
                print (cell_text)
            capacity_count = capacity_count + 1
            fig, axs = plt.subplots()
            the_table = plt.table(cellText=cell_text, loc='center', edges='open')
            the_table.auto_set_column_width(0)
            the_table.auto_set_column_width(1)
            the_table.scale(1.5,1.4)
            plt.axis('off')
            figname = fig_topic + str(capacity_count) + '.png'
# Uncomment next line to enable separate assessment & zone grouping tables
#             plt.savefig(figname, format='png', bbox_inches='tight')
            plt.close()
            found_diversified=False
            continue

# Table of energy energy demands. Based on 1st token on the line grab
# relevant information.          
    if found_demand:
        if s_line[0] == '*title':
            demand_title = s_line[1]     # remember title for graph
        elif s_line[0] == '*format':
            continue
        elif s_line[0] == '*fields':
            cell_text[0][0] = 'Within ' + context + '\n' + tbl_title
            cell_text[1][0] = s_line[1]
            cell_text[2][0] = s_line[2]
            cell_text[3][0] = s_line[3]
            cell_text[4][0] = s_line[4]
            cell_text[5][0] = s_line[5]
            cell_text[6][0] = s_line[6]
            cell_text[7][0] = s_line[7]
            continue
        elif s_line[0] == '*data':
            cell_text[0][1] = 'Demand\nkWhrs'
            cell_text[1][1] = s_line[1]
            cell_text[2][1] = s_line[2]
            cell_text[3][1] = s_line[3]
            cell_text[4][1] = s_line[4]
            cell_text[5][1] = s_line[5]
            cell_text[6][1] = s_line[6]
            cell_text[7][1] = s_line[7]
            demand_text[multidemand_count][0] = season_title
            demand_text[multidemand_count][1] = context
            demand_text[multidemand_count][2] = s_line[1]
            demand_text[multidemand_count][3] = s_line[2]
            demand_text[multidemand_count][4] = s_line[3]
            demand_text[multidemand_count][5] = s_line[4]
            demand_text[multidemand_count][6] = s_line[5]
            demand_text[multidemand_count][7] = s_line[6]
            demand_text[multidemand_count][8] = s_line[7]
            continue
        elif s_line[0] == '*end_report': # have all information needed
            if verbose:
                print (cell_text)
            demand_count = demand_count + 1
            fig, axs = plt.subplots()
            the_table = plt.table(cellText=cell_text, loc='center', edges='open')
            the_table.auto_set_column_width(0)
            the_table.auto_set_column_width(1)
            the_table.scale(1.5,1.4)
            plt.axis('off')
            figname = 'demand_' + str(demand_count) + '.png'
# Uncomment next line to enable separate assessment & zone grouping tables
#             plt.savefig(figname, format='png', bbox_inches='tight')
            plt.close()
            found_demand=False
            continue

# Table of energy dispersed demands. Based on 1st token on the line grab
# relevant information. Currently no figure generated.          
    if found_dispersed:
        if s_line[0] == '*title':
            demand_title = s_line[1]     # remember title for graph
        elif s_line[0] == '*format':
            continue
        elif s_line[0] == '*fields':
            cell_text[0][0] = 'Within ' + context + '\n' + tbl_title
            cell_text[1][0] = s_line[1]
            cell_text[2][0] = s_line[2]
            cell_text[3][0] = s_line[3]
            cell_text[4][0] = s_line[4]
            cell_text[5][0] = s_line[5]
            cell_text[6][0] = s_line[6]
            continue
        elif s_line[0] == '*data':
            cell_text[0][1] = 'Demand\nkWhrs'
            cell_text[1][1] = s_line[1]
            cell_text[2][1] = s_line[2]
            cell_text[3][1] = s_line[3]
            cell_text[4][1] = s_line[4]
            cell_text[5][1] = s_line[5]
            cell_text[6][1] = s_line[6]
        elif s_line[0] == '*end_report': # have all information needed
            if verbose:
                print (cell_text)
            found_dispersed=False
            continue

# Histogram for a specific topic. Based on 1st token on the line grab
# relevant information.          
    if found_distribution:
        if s_line[0] == '*title':
            hist_title = s_line[1]     # remember title for graph
            units = s_line[2]
            continue
        elif s_line[0] == '*format' and s_line[1] == 'frequency':
            nb_bins = int(s_line[2])   # how many bins
            bstart  = float(s_line[4]) # from
            bincr   = float(s_line[5]) # increment
            bend    = float(s_line[6]) # to
            bstartlbl = bstart - bincr # adjust labels to deal with 1st less than data
            bendlbl = bend + bincr     # adjust labels to deal with greater than data
            continue
        elif s_line[0] == '*fields':
            continue
        elif s_line[0] == '*data':
            found_data=True
            continue
        elif s_line[0] == '*end_report': # have all information needed
            if verbose:
                print (hits)
                print (percents)
            frq_hits = np.array(hits)
            frq_pc = np.array(percents)
            found_distribution=False
            found_data=False
            datacount=0

# Adapt whether vertical axis is hits or percent. The critical array size
# is the number of frq_hits or prq_pc. Adjust the creation of labels to
# match otherwise the graph will fail.
            if as_hits:
                histo_count = histo_count + 1
                fig, axs = plt.subplots()
                x = np.arange(len(frq_hits))
                if verbose:
                    print (x)
                if float(bincr) < 1.0:
                    if len(x) == 12:
                      labels = ['0.2','0.4','0.6','0.8','1.0','1.2','1.4','1.6','1.8','2.0','2.2','2.4']
                    if len(x) == 13:
                      labels = ['0.2','0.4','0.6','0.8','1.0','1.2','1.4','1.6','1.8','2.0','2.2','2.4','2.6']
                if float(bincr) >= 1.0:
                    labels = np.arange(int(bstartlbl),int(bendlbl),int(bincr))
                if len(labels) < len(x):
                    delta = len(x) - len(labels)
                    if float(bincr) < 1.0:
                        if len(x) == 12:
                            labels = ['0.2','0.4','0.6','0.8','1.0','1.2','1.4','1.6','1.8','2.0','2.2','2.4']
                        if len(x) == 13:
                            labels = ['0.2','0.4','0.6','0.8','1.0','1.2','1.4','1.6','1.8','2.0','2.2','2.4','2.6']
                    if float(bincr) >= 1.0:
                        labels = np.arange(int(bstartlbl),int(bendlbl)+delta,int(bincr))
                if len(labels) > len(x):
                    delta = len(labels) - len(x)
                    if float(bincr) < 1.0:
                        if len(x) == 12:
                            labels = ['0.2','0.4','0.6','0.8','1.0','1.2','1.4','1.6','1.8','2.0','2.2','2.4']
                        if len(x) == 13:
                            labels = ['0.2','0.4','0.6','0.8','1.0','1.2','1.4','1.6','1.8','2.0','2.2','2.4','2.6']
                    if float(bincr) >= 1.0:
                        labels = np.arange(int(bstartlbl),int(bendlbl)-delta,int(bincr))
                axs.bar(x, frq_hits, width=1, color='0.7', edgecolor="white", linewidth=0.7)
                axs.set_xticks(x)
                if larger_fonts:
                    axs.set_xticklabels(labels, fontsize=11)
                    axs.tick_params(axis='y', labelsize=11)
                    plt.xlabel(hist_title + '  ' + units + ' : ' + context + ' : ' + tbl_title, fontsize=12)
                    plt.ylabel('Hits', fontsize=12)
                else:
                    axs.set_xticklabels(labels)
                    plt.xlabel(hist_title + '  ' + units +' : ' + context + ' : ' + tbl_title)
                    plt.ylabel('Hits')
                plt.grid(axis='y', color='0.6')
                figname = fig_topic + 'frq_' + str(histo_count) + '.png'
                plt.savefig(figname, format='png', bbox_inches='tight')
                plt.close()
                print ('created ' + figname)
                hits = []      # clear the arrays
                percents = []
            else:
                histo_count = histo_count + 1
                fig, axs = plt.subplots()
                x = np.arange(len(frq_pc))
                labels = np.arange(int(bstartlbl),int(bendlbl),int(bincr))
                if len(labels) < len(x):
                    delta = len(x) - len(labels)
                    labels = np.arange(int(bstartlbl),int(bendlbl)+delta,int(bincr))
                if len(labels) > len(x):
                    delta = len(labels) - len(x)
                    labels = np.arange(int(bstartlbl),int(bendlbl)-delta,int(bincr))
                axs.bar(x, frq_pc, width=1, color='0.7', edgecolor="white", linewidth=0.7)
                axs.set_xticks(x)
                if larger_fonts:
                    axs.set_xticklabels(labels, fontsize=11)
                    axs.tick_params(axis='y', labelsize=11)
                    plt.xlabel(hist_title + ' : ' + context + ' : ' + tbl_title, fontsize=12)
                    plt.ylabel('%', fontsize=12)
                else:
                    axs.set_xticklabels(labels)
                    plt.xlabel(hist_title + ' : ' + context + ' : ' + tbl_title)
                    plt.ylabel('%')
                plt.grid(axis='y', color='0.6')
                figname = fig_topic + 'frq_' + str(histo_count) + '.png'
                plt.savefig(figname, format='png', bbox_inches='tight')
                plt.close()
                print ('created ' + figname)
                hits = []      # clear the arrays
                percents = []
            
            continue

# Histogram for a specific topic. Based on 1st token on the line grab
# relevant information.          
    if found_summary_distribution:
        if s_line[0] == '*title':
            hist_title = s_line[1]     # remember title for graph
            units = s_line[2]
            continue
        elif s_line[0] == '*format' and s_line[1] == 'frequency':
            nb_bins = int(s_line[2])   # how many bins
            bstart  = float(s_line[4]) # from
            bincr   = float(s_line[5]) # increment
            bend    = float(s_line[6]) # to
            bstartlbl = bstart - bincr # adjust labels to deal with 1st less than data
            bendlbl = bend + bincr     # adjust labels to deal with greater than data
            continue
        elif s_line[0] == '*fields':
            continue
        elif s_line[0] == '*data':
            found_data=True
            continue
        elif s_line[0] == '*end_report': # have all information needed
            if verbose:
                print (hits_s1)
                print (hits_s2)
                print (hits_s3)
                print (hits_s4)
                print (hits_s5)
            frq_hits_s1 = np.array(hits_s1)
            frq_hits_s2 = np.array(hits_s2)
            frq_hits_s3 = np.array(hits_s3)
            frq_hits_s4 = np.array(hits_s4)
            frq_hits_s5 = np.array(hits_s5)
            found_summary_distribution=False
            found_data=False
            datacount=0

            histo_count = histo_count + 1
            fig, axs = plt.subplots()
            x = np.arange(len(frq_hits_s1))
            labels = np.arange(int(bstartlbl),int(bendlbl),int(bincr))
            if verbose:
                 print (len(frq_hits))
                 print (len(x))
                 print (len(labels))
                 print (x)
                 print (labels)
            if len(labels) < len(x):
                delta = len(x) - len(labels)
                labels = np.arange(int(bstartlbl),int(bendlbl)+delta,int(bincr))
            if len(labels) > len(x):
                delta = len(labels) - len(x)
                labels = np.arange(int(bstartlbl),int(bendlbl)-delta,int(bincr))
            axs.bar(x - 0.3, frq_hits_s1, width=0.2, label='1st winter', color='lightblue', edgecolor="white", linewidth=0.7)
            axs.bar(x - 0.1, frq_hits_s2, width=0.2, label='spring', color='gold', edgecolor="white", linewidth=0.7)
            axs.bar(x + 0.1, frq_hits_s3, width=0.2, label='summer', color='salmon', edgecolor="white", linewidth=0.7)
            if assessment_count > 3:
                axs.bar(x + 0.3, frq_hits_s4, width=0.2, label='autumn', color='orange', edgecolor="white", linewidth=0.7)
                axs.bar(x + 0.5, frq_hits_s5, width=0.2, label='2nd winter', color='silver', edgecolor="white", linewidth=0.7)
            axs.set_xticks(x)
            if larger_fonts:
                axs.set_xticklabels(labels, fontsize=11)
                axs.tick_params(axis='y', labelsize=11)
                plt.xlabel(hist_title + ' : seasonal frequency', fontsize=12)
                plt.ylabel('Hits', fontsize=12)
                plt.legend(loc='best', fontsize=11)
            else:
                axs.set_xticklabels(labels)
                plt.xlabel(hist_title + ' : seasonal frequency')
                plt.ylabel('Hits')
                plt.legend(loc='best')
            plt.grid(axis='y', color='0.6')
            figname = 'summary_comfort_distribution' + '.png'
            plt.savefig(figname, format='png', bbox_inches='tight')
            plt.close()
            print ('created ' + figname)
            hits_s1 = []      # clear the arrays
            hits_s2 = []
            hits_s3 = []
            hits_s4 = []
            hits_s5 = [] 
            continue

# If focus is on stats then slot data based on the current assessment.
    if found_stats:
        if s_line[0] == '*title':
            stats_text[stats_row][0] = season_title    # season
            stats_text[stats_row][1] = context         # which zones
            stats_text[stats_row][2] = s_line[1]       # topic
            continue
        elif s_line[0] == '*format':
            continue
        elif s_line[0] == '*fields':
            continue
        elif s_line[0] == '*data':
            stats_text[stats_row][3] = float(s_line[1])     # max
            stats_text[stats_row][4] = float(s_line[2])     # min
            stats_text[stats_row][5] = float(s_line[3])     # avg
            continue
        elif s_line[0] == '*end_report':
            found_stats=False
            continue
        continue

# If focus is on emissions build up multi bar chart.
    if found_emissions:
        if s_line[0] == '*title':
            hist_title = s_line[1]   # title
            units = s_line[2]        # units
            continue
        elif s_line[0] == '*format':
            continue
        elif s_line[0] == '*fields':
            continue
        elif s_line[0] == '*data':
            if emission_count == 0:
                emission_ht.append(float(s_line[1]))
                emission_cl.append(float(s_line[2]))
                emission_lt.append(float(s_line[3]))
                emission_fn.append(float(s_line[4]))
                emission_sp.append(float(s_line[5]))
                emission_hw.append(float(s_line[6]))
            elif emission_count ==1:
                emission_ht.append(float(s_line[1]))
                emission_cl.append(float(s_line[2]))
                emission_lt.append(float(s_line[3]))
                emission_fn.append(float(s_line[4]))
                emission_sp.append(float(s_line[5]))
                emission_hw.append(float(s_line[6]))
            elif emission_count == 2:
                emission_ht.append(float(s_line[1]))
                emission_cl.append(float(s_line[2]))
                emission_lt.append(float(s_line[3]))
                emission_fn.append(float(s_line[4]))
                emission_sp.append(float(s_line[5]))
                emission_hw.append(float(s_line[6]))
            emission_count = emission_count + 1
            continue
        elif s_line[0] == '*end_report':
            frq_ht = np.array(emission_ht)
            frq_cl = np.array(emission_cl)
            frq_lt = np.array(emission_lt)
            frq_fn = np.array(emission_fn)
            frq_sp = np.array(emission_sp)
            frq_hw = np.array(emission_hw)
            datacount=0

            histo_count = histo_count + 1
            fig, axs = plt.subplots()
            x = np.arange(3)
            labels = ('CO2','NOx','SOx')
            axs.bar(x - 0.3, frq_ht, width=0.15, label='heating', color='red', edgecolor="white", linewidth=0.7)
            axs.bar(x - 0.15, frq_cl, width=0.15, label='cooling', color='blue', edgecolor="white", linewidth=0.7)
            axs.bar(x, frq_lt, width=0.15, label='lighting', color='orange', edgecolor="white", linewidth=0.7)
            axs.bar(x + 0.15, frq_fn, width=0.15, label='fans', color='olivedrab', edgecolor="white", linewidth=0.7)
            axs.bar(x + 0.3, frq_sp, width=0.15, label='small power', color='silver', edgecolor="white", linewidth=0.7)
            axs.bar(x + 0.45, frq_hw, width=0.15, label='DHW', color='aquamarine', edgecolor="white", linewidth=0.7)
            axs.set_xticks(x)
            axs.set_yscale('log')
            if larger_fonts:
                axs.set_xticklabels(labels, fontsize=11)
                axs.tick_params(axis='y', labelsize=11)
                plt.xlabel(hist_title + ' : annual', fontsize=12)
                plt.ylabel(units, fontsize=12)
                plt.legend(loc='best', fontsize=11)
            else:
                axs.set_xticklabels(labels)
                plt.xlabel(hist_title + ' : annual')
                plt.ylabel(units)
                plt.legend(loc='best')
            plt.grid(axis='y', color='0.6')
            figname = 'emission_freq' + str(histo_count) + '.png'
            plt.savefig(figname, format='png', bbox_inches='tight')
            plt.close()
            print ('created ' + figname)
            emission_ht =[]
            emission_cl =[]
            emission_lt =[]
            emission_fn =[]
            emission_sp =[]
            emission_hw =[]
            emission_count = 0
            found_emissions=False
            continue

# If focus is on summary capacity build bar chart.
    if found_power:
        if s_line[0] == '*title':
            hist_title = s_line[1]   # title
            units = s_line[2]        # units
            continue
        elif s_line[0] == '*format':
            continue
        elif s_line[0] == '*fields':
            continue
        elif s_line[0] == '*data':
            power_cap.append(float(s_line[1]))
            power_cap.append(float(s_line[2]))
            power_cap.append(float(s_line[3]))
            power_cap.append(float(s_line[4]))
            power_cap.append(float(s_line[5]))
            power_cap.append(float(s_line[6]))
            continue
        elif s_line[0] == '*end_report':
            power_bar = np.array(power_cap)
            topics = ['Heating', 'Cooling', 'Lighting', 'Fans', 'SPL', 'DHW']
            topic_colours = ['red', 'blue', 'orange', 'olivedrab', 'silver', 'aquamarine']
            fig, axs = plt.subplots()
            axs.bar(topics, power_bar, width=0.5, color=topic_colours, edgecolor="white", linewidth=0.7)
            if larger_fonts:
                plt.xlabel( 'Maximum capacity ' + ' annual', fontsize=12)
                plt.ylabel(units, fontsize=12)
                axs.tick_params(axis='y', labelsize=11)
            else:
                plt.xlabel( 'Maximum capacity ' + ' annual')
                plt.ylabel(units, fontsize=12)
            plt.grid(axis='y', color='0.6')
            addlabels(topics,power_bar)  # Call local function.
            figname = 'max_capacity' + str(histo_count) + '.png'
            plt.savefig(figname, format='png', bbox_inches='tight')
            plt.close()
            print ('created ' + figname)
            power_cap =[]
            found_power=False
            continue


# Process summary energy add to the seasonal table.
    if found_energy:
        if s_line[0] == '*title':
            seasonal_row = seasonal_row +1
            seasonal_text[seasonal_row][0] = 'annual'    # season
            seasonal_text[seasonal_row][1] = s_line[2]   # unit
            if number_of_simulations == 3:
                seasonal_text3[seasonal_row][0] = 'annual'    # season
                seasonal_text3[seasonal_row][1] = s_line[2]   # unit
            continue
        elif s_line[0] == '*format':
            continue
        elif s_line[0] == '*fields':
            continue
        elif s_line[0] == '*data':
            seasonal_text[seasonal_row][2] = float(s_line[1])     # heating
            seasonal_text[seasonal_row][3] = float(s_line[2])     # cooling
            seasonal_text[seasonal_row][4] = float(s_line[3])     # lighting
            seasonal_text[seasonal_row][5] = float(s_line[4])     # fans
            seasonal_text[seasonal_row][6] = float(s_line[5])     # small power
            seasonal_text[seasonal_row][7] = float(s_line[6])     # dhw
            if number_of_simulations == 3:
                seasonal_text3[seasonal_row][2] = float(s_line[1])     # heating
                seasonal_text3[seasonal_row][3] = float(s_line[2])     # cooling
                seasonal_text3[seasonal_row][4] = float(s_line[3])     # lighting
                seasonal_text3[seasonal_row][5] = float(s_line[4])     # fans
                seasonal_text3[seasonal_row][6] = float(s_line[5])     # small power
                seasonal_text3[seasonal_row][7] = float(s_line[6])     # dhw

# Also fill in the vertical summary of energy demand cells. 1st will be kWh/m^2.a.
            summary_energy_count = summary_energy_count +1
            if summary_energy_count == 1:
                annual_demand_text[1][1] = float(s_line[1])     # heating
                annual_demand_text[2][1] = float(s_line[2])     # cooling
                annual_demand_text[3][1] = float(s_line[3])     # lighting
                annual_demand_text[4][1] = float(s_line[4])     # fans
                annual_demand_text[5][1] = float(s_line[5])     # small power
                annual_demand_text[6][1] = float(s_line[6])     # dhw
            elif summary_energy_count == 2:
                annual_demand_text[1][2] = float(s_line[1])     # heating
                annual_demand_text[2][2] = float(s_line[2])     # cooling
                annual_demand_text[3][2] = float(s_line[3])     # lighting
                annual_demand_text[4][2] = float(s_line[4])     # fans
                annual_demand_text[5][2] = float(s_line[5])     # small power
                annual_demand_text[6][2] = float(s_line[6])     # dhw
            continue
        elif s_line[0] == '*end_report':
            continue

        
# Gather hits and percentages for the histogram.
    if found_distribution and found_data:
        hits.append(float(s_line[2]))     # the number of hits
        percents.append(float(s_line[3])) # the percents for this bin
        datacount = datacount +1
        continue
        
# Gather hits and percentages for the histogram.
    if found_summary_distribution and found_data:
        hits_s1.append(float(s_line[1]))
        hits_s2.append(float(s_line[2]))
        hits_s3.append(float(s_line[3]))
        if assessment_count > 3:
            hits_s4.append(float(s_line[4]))
            hits_s5.append(float(s_line[5]))
        datacount = datacount +1
        continue

# If there are multiple assessments reflect this in the title.
# If status table partly filled generate figure.
    if s_line[0] == '*assessment':
        assessment_count = int(s_line[1])
        run_title = 'Period: ' + s_line[2]
        tbl_title =  s_line[2]
        season_title =  s_line[2]            
        continue

# If end of initial report reached and there are stats table entries, generate.
# the first *end should be followed by seasonal reports.
    if s_line[0] == '*end':
        if stats_row > 1:
            fig, axs = plt.subplots()
            the_table = plt.table(cellText=stats_text, loc='center', edges='open')
            the_table.auto_set_column_width(0)  # adjust width of each coloum
            the_table.auto_set_column_width(1)
            the_table.auto_set_column_width(2)
            the_table.auto_set_column_width(3)
            the_table.auto_set_column_width(4)
            the_table.auto_set_column_width(5)
            the_table.scale(1.5,1.4)
            plt.axis('off')
            figname = 'summary_stats' + '.png'
            plt.savefig(figname, format='png', bbox_inches='tight')
            plt.close()
            print ('created ' + figname)
            stats_row = 0
            found_stats = False

# Option to also include in this the seasonal summaries as well as the annual.
# Depending on the number of seasons fine tune the size of the figure to limit white space.
        if seasonal_row > 1 and found_summary:
            if number_of_simulations == 1:
                fig, axs = plt.subplots(figsize=(8.0,2.5), layout='constrained')
                the_table = plt.table(cellText=seasonal_text, loc='center', edges='open')
            elif number_of_simulations == 3:
                fig, axs = plt.subplots(figsize=(6.5,2.6), layout='constrained')
                the_table = plt.table(cellText=seasonal_text3, loc='center', edges='open')
            elif number_of_simulations == 5:
                fig, axs = plt.subplots(figsize=(8.0,4.0), layout='constrained')
                the_table = plt.table(cellText=seasonal_text, loc='center', edges='open')
            the_table.auto_set_column_width(0)  # adjust width of each coloum
            the_table.auto_set_column_width(1)
            the_table.auto_set_column_width(2)
            the_table.auto_set_column_width(3)
            the_table.auto_set_column_width(4)
            the_table.auto_set_column_width(5)
            the_table.auto_set_column_width(6)
            the_table.auto_set_column_width(7)
            the_table.scale(1.5,1.4)
            plt.axis('off')
            figname = 'summary_seasonal' + '.png'
#            plt.savefig(figname, format='png', bbox_inches='tight')
            plt.savefig(figname, format='png', bbox_inches='tight')
            plt.close()
            print ('created ' + figname)
            seasonal_row = 0

# Option to also create table of annual building demands.
        if summary_energy_count > 1 and found_summary:
            fig, axs = plt.subplots(figsize=(3.2,2.6), layout='constrained')
            the_table = plt.table(cellText=annual_demand_text, loc='center', edges='open')
            the_table.auto_set_column_width(0)  # adjust width of each coloum
            the_table.auto_set_column_width(1)
            the_table.auto_set_column_width(2)
            the_table.scale(1.5,1.4)
            plt.axis('off')
            plt.text(0.5, 0.05, 'Annual energy demands', ha='center', va='top', wrap=True)
            figname = 'annual_demands_v' + '.png'
            plt.savefig(figname, format='png', bbox_inches='tight')
            plt.close()
            print ('created ' + figname)
            seasonal_row = 0

        if diversified_count > 1:
            fig, axs = plt.subplots()
            the_table = plt.table(cellText=capacity_text, loc='center', edges='open')
            the_table.auto_set_column_width(0)  # adjust width of each coloum
            the_table.auto_set_column_width(1)
            the_table.auto_set_column_width(2)
            the_table.auto_set_column_width(3)
            the_table.auto_set_column_width(4)
            the_table.auto_set_column_width(5)
            the_table.auto_set_column_width(6)
            the_table.auto_set_column_width(7)
            the_table.auto_set_column_width(8)
            the_table.auto_set_column_width(9)
            the_table.scale(1.5,1.4)

# Adjust position of where the table lable is placed at the bottom of the entries.
            if diversified_count <= 3:
                plt.text(0.5, 0.55, capacity_title + ' W', ha='center', va='top', wrap=True)
            elif diversified_count > 3 and diversified_count <= 7:
                plt.text(0.5, 0.45, capacity_title + ' W', ha='center', va='top', wrap=True)
            elif diversified_count > 7 and diversified_count <= 9:
                plt.text(0.5, 0.4, capacity_title + ' W', ha='center', va='top', wrap=True)
            elif diversified_count > 9 and diversified_count <= 12:
                plt.text(0.5, 0.35, capacity_title + ' W', ha='center', va='top', wrap=True)
            elif diversified_count > 12 and diversified_count <= 16:
                plt.text(0.5, 0.25, capacity_title + ' W', ha='center', va='top', wrap=True)
            elif diversified_count > 16 and diversified_count <= 20:
                plt.text(0.5, -0.1, capacity_title + ' W', ha='center', va='top', wrap=True)
            else:
                plt.text(0.5, -0.2, capacity_title + ' W', ha='center', va='top', wrap=True)
	    
#             plt.title(capacity_title + ' ' + units )
            plt.axis('off')
            figname = 'summary_diversified_capacity' + '.png'
            plt.savefig(figname, format='png', bbox_inches='tight')
            plt.close()
            print ('created ' + figname)
            diversified_count = 0

# Generate a summary table of the various demands collected during the scan.
        if multidemand_count > 1:
            fig, axs = plt.subplots()
            the_table = plt.table(cellText=demand_text, loc='center', edges='open')
            the_table.auto_set_column_width(0)
            the_table.auto_set_column_width(1)
            the_table.auto_set_column_width(2)
            the_table.auto_set_column_width(3)
            the_table.auto_set_column_width(4)
            the_table.auto_set_column_width(5)
            the_table.auto_set_column_width(6)
            the_table.auto_set_column_width(7)
            the_table.auto_set_column_width(8)
            the_table.auto_set_column_width(9)
            the_table.scale(1.5,1.4)

# Adjust position of where the table lable is placed at the bottom of the entries.
            if multidemand_count <= 3:
                plt.text(0.5, 0.55, demand_title + ' kWh', ha='center', va='top', wrap=True)
            elif multidemand_count > 3 and multidemand_count <= 7:
                plt.text(0.5, 0.45, demand_title + ' kWh', ha='center', va='top', wrap=True)
            elif multidemand_count > 7 and multidemand_count <= 9:
                plt.text(0.5, 0.4, demand_title + ' kWh', ha='center', va='top', wrap=True)
            elif multidemand_count > 9 and multidemand_count <= 12:
                plt.text(0.5, 0.35, demand_title + ' kWh', ha='center', va='top', wrap=True)
            elif multidemand_count > 12 and multidemand_count <= 16:
                plt.text(0.5, 0.25, demand_title + ' kWh', ha='center', va='top', wrap=True)
            elif multidemand_count > 16 and multidemand_count <= 20:
                plt.text(0.5, -0.1, demand_title + ' kWh', ha='center', va='top', wrap=True)
            else:
                plt.text(0.5, -0.2, demand_title + ' kWh', ha='center', va='top', wrap=True)
            plt.axis('off')
            figname = 'summary_demands' + '.png'
            plt.savefig(figname, format='png', bbox_inches='tight')
            plt.close()
            print ('created ' + figname)
            multidemand_count = 0
            continue
        continue

# Process seasonal assessments. In this case ignore *end_report lines so as
# keep reading the various annual reports.
    if found_seasonal:
        if s_line[0] == '*title':
            seasonal_row = seasonal_row +1
            seasonal_text[seasonal_row][0] = current_season    # season
            seasonal_text[seasonal_row][1] = s_line[2]         # unit
            if number_of_simulations == 3:
                seasonal_text3[seasonal_row][0] = current_season    # season
                seasonal_text3[seasonal_row][1] = s_line[2]         # unit
            continue
        elif s_line[0] == '*format':
            continue
        elif s_line[0] == '*fields':
            continue
        elif s_line[0] == '*data':
            seasonal_text[seasonal_row][2] = float(s_line[1])     # heating
            seasonal_text[seasonal_row][3] = float(s_line[2])     # cooling
            seasonal_text[seasonal_row][4] = float(s_line[3])     # lighting
            seasonal_text[seasonal_row][5] = float(s_line[4])     # fans
            seasonal_text[seasonal_row][6] = float(s_line[5])     # small power
            seasonal_text[seasonal_row][7] = float(s_line[6])     # dhw
            if number_of_simulations == 3:
                seasonal_text3[seasonal_row][2] = float(s_line[1])     # heating
                seasonal_text3[seasonal_row][3] = float(s_line[2])     # cooling
                seasonal_text3[seasonal_row][4] = float(s_line[3])     # lighting
                seasonal_text3[seasonal_row][5] = float(s_line[4])     # fans
                seasonal_text3[seasonal_row][6] = float(s_line[5])     # small power
                seasonal_text3[seasonal_row][7] = float(s_line[6])     # dhw
            continue
        elif s_line[0] == '*end_report':
            continue
        elif s_line[0] == '*Seasonal_summary':
            current_season = s_line[1]
            continue
        elif s_line[0] == '*Summary':
            current_season = 'Annual'
            found_seasonal=False
            found_summary=True
            continue
        elif s_line[0] == '*End_seasonal_summary':
            continue
        elif s_line[0] == '*report' and s_line[2] == 'emissions':
            found_emissions=True     # getting towards the end
            context = s_line[3]      # remember context of report
            continue
        elif s_line[0] == '*report':
            continue


# Process summary section. ?? In this case ignore *end_report lines so as
# keep reading the various annual reports.
    if found_summary:
        if s_line[0] == '*report' and s_line[2] == 'emissions':
            found_emissions=True     # getting towards the end
            context = s_line[3]      # remember context of report
            continue
        elif s_line[0] == '*report' and s_line[2] == 'power':
            found_power=True     # getting towards the end
            context = s_line[3]      # remember context of report
            continue
        elif s_line[0] == '*report' and s_line[2] == 'energy':
            found_energy=True        # getting towards the end 
            context = s_line[3]      # remember context of report
            continue
        elif s_line[0] == '*report' and s_line[2] == 'distribution':
            found_summary_distribution=True
            fig_topic = s_line[3]
            context = s_line[3]     # the overall seasonal
            continue
        elif s_line[0] == '*title':
            seasonal_row = seasonal_row +1
            seasonal_text[seasonal_row][0] = 'annual'    # season
            seasonal_text[seasonal_row][1] = s_line[1] + ' ' + s_line[2]    # topic
            seasonal_text3[seasonal_row][0] = 'annual'    # season
            seasonal_text3[seasonal_row][1] = s_line[1] + ' ' + s_line[2]    # topic
            continue
        elif s_line[0] == '*format':
            continue
        elif s_line[0] == '*fields':
            continue
        elif s_line[0] == '*data':  # hopefully not used
            seasonal_text[seasonal_row][2] = float(s_line[1])     # heating
            seasonal_text[seasonal_row][3] = float(s_line[2])     # cooling
            seasonal_text[seasonal_row][4] = float(s_line[3])     # lighting
            seasonal_text[seasonal_row][5] = float(s_line[4])     # fans
            seasonal_text[seasonal_row][6] = float(s_line[5])     # small power
            seasonal_text[seasonal_row][7] = float(s_line[6])     # dhw
            seasonal_text3[seasonal_row][2] = float(s_line[1])     # heating
            seasonal_text3[seasonal_row][3] = float(s_line[2])     # cooling
            seasonal_text3[seasonal_row][4] = float(s_line[3])     # lighting
            seasonal_text3[seasonal_row][5] = float(s_line[4])     # fans
            seasonal_text3[seasonal_row][6] = float(s_line[5])     # small power
            seasonal_text3[seasonal_row][7] = float(s_line[6])     # dhw
            continue
        elif s_line[0] == '*end_report':
            continue
        elif s_line[0] == '*Seasonal_summary':
            current_season = s_line[1]
            continue
        elif s_line[0] == '*Summary':
            current_season = 'Annual'
            found_seasonal=False
            found_summary=True
            continue
        elif s_line[0] == '*End_seasonal_summary':
            continue
        elif s_line[0] == '*report' and s_line[2] == 'emissions':
            found_emissions=True     # getting towards the end
            context = s_line[3]      # remember context of report
            continue
        elif s_line[0] == '*report':
            continue

        
# If '*report' & 'distribution' set found_distribution.
    if s_line[0] == '*report' and s_line[2] == 'distribution':
        if s_line[3] == 'thermal_comfort':
            found_summary_distribution=True
            fig_topic = s_line[3]
            context = s_line[4]     # the overall seasonal
        else:
            fig_topic = s_line[3]
            found_distribution=True
            context = s_line[4]     # remember context of report
        continue

# If '*report' & 'diversified' set capacity table. Increment diversified_count
# each time it is encountered.
    if s_line[0] == '*report' and s_line[2] == 'diversified':
        found_diversified=True
        diversified_count=diversified_count+1  
        fig_topic = s_line[3]
        context = s_line[4]     # remember context of report
        continue

# If '*report' & 'demand' & 'integrated' set demand table.
    if s_line[0] == '*report' and s_line[2] == 'demand' and s_line[3] == 'integrated':
        if s_line[4] == 'Dispersed':
            fig_topic = s_line[3]
            found_dispersed=True
            context = s_line[4]     # remember context of report
        else:
            fig_topic = s_line[3]
            found_demand=True
            multidemand_count = multidemand_count +1
            context = s_line[4]     # remember context of report
        continue

# If '*report' & 'demand' & 'per_unit_time' set for energy graph.
    if s_line[0] == '*report' and s_line[2] == 'demand' and s_line[3] == 'per_unit_time':
        if s_line[4] == 'Dispersed':
            context = s_line[4]     # remember context of report
        else:
            found_graph=True
            found_demand_per_unit_time=True
            context = s_line[4]     # remember context of report
        continue

# If '*report' & '*delivered_per_unit_time' set for energy graph.
    if s_line[0] == '*report' and s_line[2] == '*delivered_per_unit_time':
        found_graph=True
        found_delivered_per_unit_time=True
        context = s_line[3]     # remember context of report
        continue

# If '*report' & 'thermal_comfort' set for temperature graph.
    if s_line[0] == '*report' and s_line[2] == 'thermal_comfort':
        found_graph=True
        found_thermal_comfort=True
        context = s_line[4]     # remember context of report
        continue

# If '*report' & 'zone_dbt' set for temperature graph.
    if s_line[0] == '*report' and s_line[2] == 'zone_dbt':
        found_graph=True
        found_dbT=True
        context = s_line[4]     # remember context of report
        continue

# If '*report' & 'relative humidity' set for resultant temperature graph.
    if s_line[0] == '*report' and s_line[2] == 'relative humidity':
        found_graph=True
        found_rh=True
        context = s_line[4]     # remember context of report
        continue

# If '*report' & 'total casual gains' set for casual gains graph.
    if s_line[0] == '*report' and s_line[2] == 'total casual gains':
        found_graph=True
        found_casual=True
        context = s_line[4]     # remember context of report
        continue

# If '*report' & 'infiltration load' set for zone flux graph.
    if s_line[0] == '*report' and s_line[2] == 'infiltration load':
        found_graph=True
        found_infiltration=True
        context = s_line[4]     # remember context of report
        continue

# If '*report' & 'solar entering from outside' set for solar graph.
    if s_line[0] == '*report' and s_line[2] == 'solar entering from outside':
        found_graph=True
        found_solar=True
        context = s_line[4]     # remember context of report
        continue

# If '*report' & 'emissions' set for emissions graph.
    if s_line[0] == '*report' and s_line[2] == 'emissions':
        found_emissions=True
        context = s_line[3]     # remember context of report
        continue

# If '*report' & 'stats' set for stats table.
    if s_line[0] == '*report' and s_line[2] == 'stats':
        found_stats=True
        stats_row = stats_row + 1  # increment row in multstats report.
        context = s_line[4]        # remember context of report
        continue

# If '*report' & 'Daylight' set for DF Graph.
    if s_line[0] == '*report' and s_line[2] == 'Daylight':
        found_daylight=True
        context = s_line[4]        # remember context of report
        continue

# If '*report' & 'Guth' set for Guth Graph.
    if s_line[0] == '*report' and s_line[2] == 'Guth':
        found_guth=True
        context = s_line[4]        # remember context of report
        continue

# If '*Seasonal_summary set for summary table.
    if s_line[0] == '*Seasonal_summary':
        found_seasonal=True
        current_season = s_line[1]
        continue

# If other reports types need to be noted this is where they would be added.
        
