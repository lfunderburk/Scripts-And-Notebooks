#Authors: Laura Gutierrez Funderburk, Bryce Hayley
# Date: February 8 - 12 2018
# COMAP 2018

# Import libraries
import csv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from numpy.random import rand
from numpy import arange

"""Part 0: Constants"""

## Profile States: Selected Categories

energy = ['Coal','Jet Fuel','Natural Gas','Nuclear power',\
                    'Biomass','Electricity','Geothermal','Hydro',\
                    'Motor Gas','Photovalic and Solar Thermal',\
                    'Wind']

msn_array=['CLTCB','JFTCB','NNTCB','NUETB','BMTCB','ESTCB',\
           'GETCB','HYTCB','MMTCB','SOTCB','WYTCB']


## Profiles: Consumption and Production

energy_con = ['Coal','Natural Gas w/supp. fuels', 'Motor Gasoline exld. Ethenol',\
         'Distillate fuel', 'Jet Fuel', 'LGP', 'Residual Fuels', 'Other Petroleum',\
         'Biomass','Total Electricity']

msn_array_con=['CLTCB','NGTCB','MMTCB', 'DFTCB','JFTCB','LGTCB','RFTCB','POTCB', \
        'BMTCB', 'ESTCB']
        
energy_pro = ['Coal','Natural Gas - Marketed','Crude Oil','Electricity from Nuclear', 'Wind Electricity',\
              'Electricity from Solar Energy', 'Hydroelectricity', 'Geothermal Electricity']

msn_array_pro=['CLPRB','NGMPB','PAPRB','NUETB','WYTCB','SOEGB','HYTCB','GEEGB']

## Scoring Scheme: Selected Categories

badCodes=['CLTCB','NGTCB','MMTCB', 'DFTCB','JFTCB','LGTCB','RFTCB']
badNames=['Coal','Natural Gas', 'Motor Gasoline','Distillate','Jet Fuel','Propane/LGP','Residual']
badNumbers = np.array([95.35, 53.07,71.30, 73.16,70.9, 63.07,68.79])*1000

"""Part I: File reading"""

# COMPLETESEDS DATA
# Read complete csv file
complete_seeds_file = "./completeseds19602009.csv"
# Read file using pandas
comp_seed_df = pd.read_csv(complete_seeds_file)
# Read file using pandas, with a bit more specific parameters
data = pd.read_csv(complete_seeds_file,sep=',',na_values='Nothing')
# Dataframes!
df = pd.DataFrame(data)


# PROBLEM C DATA
# Read problemC excel dataset
problemCData = "./ProblemCData.xlsx"
prob_c_data = pd.ExcelFile(problemCData)

# Extract labels from Excel file

labels = prob_c_data.parse('msncodes',parse_cols=[0,1,2])
label_dictionary = {labels['MSN'][i]:(labels['Description'][i],labels['Unit'][i]) for i in range(len(labels))}

"""Part II: Identifying states"""

# Identify states of interest
california = df.loc[df['StateCode'] == "CA"]
arizona = df.loc[df['StateCode']=="AZ"]
newmexico = df.loc[df['StateCode']=="NM"]
texas = df.loc[df['StateCode']=='TX']

""" Part III: Function definition"""
########
# Simple histogram for a given state and two given years.
########
def plot_state_label(state,start_p,end_p):
    st = state[start_p:end_p]    
    
    code_p = st['MSN'].iloc[0]
    st_code = st['StateCode'].iloc[0]
    desc = label_dictionary[code_p][0]
    unit = label_dictionary[code_p][1]

    
    y = st['Data']
    x = st['Year']
    plt.plot(x,y,label=code_p)
    plt.xlabel('Year')
    plt.ylabel(unit)
    plt.title(desc + st_code + " 1960 - 2009")
    plt.show()
    
########
# Simple histogram for a given state and MSN.
########
def plot_state_msn(state, msn):
    dfstate = pd.DataFrame(state)
    access_msn = state.loc[state['MSN']==msn]
    
    x = access_msn['Year']
    y = access_msn['Data']
    
    desc = label_dictionary[msn][0]
    unit = label_dictionary[msn][1]
    
    plt.plot(x,y,'bs',label=msn)
    plt.xlabel('Year')
    plt.ylabel(unit)
    plt.title(desc  + str(x.iloc[0]) + "-" + str(x.iloc[-1]))
    plt.show()
    
########
# This function will take as input a state, and an array with MSN values and output an array with 'Data' values from 2009
########    
def get_data_array(state,msn_array,year):
    msn = [state.loc[state['MSN']==item] for item in msn_array]
    twentynine_data = [item.iloc[year]['Data'] for item in msn]
    
    return twentynine_data

########
# This function will take as input an array with data from 2009 for a given state (see function above), and a state name. This function plots a histogram with what will later become the state profile
######## 

def plot_profile(state_name_con,state_con_array,msn_con,\
                state_name_pro,state_pro_array,msn_pro):
    pos_con = np.arange(len(msn_con))
    #print(pos)
    plt.rcdefaults()
    fix,(ax1,ax2) = plt.subplots(2,1,sharex=True)
    ax1.grid(True)
    ax2.grid(True)

    ax1.barh(pos_con,np.array(state_con_array), align='center',height=0.5)
    ax1.set_yticks(pos_con)
    for tick in ax1.get_xticklabels():
        tick.set_rotation(45)
    #ax1.xticks(rotation=90)
    ax1.set_yticklabels(msn_con)
    ax1.invert_yaxis()
    ax1.set_xlabel("Billion Btu")
    ax1.set_title(str(state_name_con))
    
    #print("\n")
    pos_pro = np.arange(len(msn_pro))
    ax2.barh(pos_pro,np.array(state_pro_array), align='center',height=0.5)
    ax2.set_yticks(pos_pro)
    for tick in ax2.get_xticklabels():
        tick.set_rotation(45)
    ax2.set_yticklabels(msn_pro)
    ax2.invert_yaxis()
    ax2.set_xlabel("Billion Btu")
    ax2.set_title(str(state_name_pro))
    plt.show()
    
####
# This takes data from a given state and MSC code and stores the data in an array
####
def dataframe_to_array(state,msn):
    dfstate = pd.DataFrame(state)
    access_msn = state.loc[state['MSN']==msn]
    y = access_msn['Data']
    data_values = [y.iloc[i] for i in range(len(y))]
    return data_values



#### 
# This function calculates states scores
####
def state_grade(state,year):
    bad_state = np.array(get_data_array(state,badCodes,year))
    state_score = bad_state*badNumbers
    state_grade = np.sum(state_score)
    return state_grade


####
# This function will plot energy consumption for a given state over time
####

def plot_consumption_per_state(state,state_name):
    dfstate = pd.DataFrame(state)
    access_CLTCB = state.loc[state['MSN']=='CLTCB']
    access_NGTCB = state.loc[state['MSN']=='NGTCB']
    access_JFTCB = state.loc[state['MSN']=='JFTCB']
    access_ESTCB = state.loc[state['MSN']=='ESTCB']
    access_LGTCB = state.loc[state['MSN']=='LGTCB']

    x = access_CLTCB['Year']
    y = access_CLTCB['Data']

    a = access_NGTCB['Year']
    b = access_NGTCB['Data']

    c = access_JFTCB['Year']
    d = access_JFTCB['Data']

    e = access_ESTCB["Year"]
    f = access_ESTCB['Data']
    
    g = access_LGTCB["Year"]
    h = access_LGTCB['Data']
    
    desc = label_dictionary['CLTCB'][0]
    unit = label_dictionary['CLTCB'][1]
    
    plt.plot(x,y,'b--',label='Coal')
    plt.plot(a,b,'g--',label='Natural Gas')
    plt.plot(c,d,'r--',label='Jet Fuel')
    plt.plot(e,f,'y--',label='Total electricity')
    plt.plot(g,h,'m--',label='LPG')
    plt.xlabel('Year')
    plt.ylabel(unit)
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.title(state_name+ " Energy Consumption"  + str(x.iloc[0]) + "-" + str(x.iloc[-1]))
    plt.show()

    
####
# This function will plot energy production for a given state over time
####

def plot_electricity_production_per_state(state,state_name):
    dfstate = pd.DataFrame(state)
    access_CLTCB = state.loc[state['MSN']=='CLPRB']
    access_NGTCB = state.loc[state['MSN']=='NGMPB']
    access_JFTCB = state.loc[state['MSN']=='NUETB']
    access_ESTCB = state.loc[state['MSN']=='WYTCB']
    acces_SOEGB = state.loc[state['MSN']=='SOEGB']
    access_HYTCB = state.loc[state['MSN']=='HYTCB']
    access_GEEGB = state.loc[state['MSN']=='GEEGB']

    x = access_CLTCB['Year']
    y = access_CLTCB['Data']

    a = access_NGTCB['Year']
    b = access_NGTCB['Data']

    c = access_JFTCB['Year']
    d = access_JFTCB['Data']

    e = access_ESTCB["Year"]
    f = access_ESTCB['Data']
    
    g = acces_SOEGB['Year']
    h = acces_SOEGB['Data']
    
    i = access_HYTCB['Year']
    j = access_HYTCB['Data']
    
    l = access_GEEGB['Year']
    m = access_GEEGB['Data']
    
    desc = label_dictionary['CLTCB'][0]
    unit = label_dictionary['CLTCB'][1]
    
    plt.plot(x,y,'b--',label='Coal')
    plt.plot(a,b,'g--',label='Natural Gas')
    plt.plot(c,d,'r--',label='Nuclear Power')
    plt.plot(e,f,'y--',label='Eolic')
    plt.plot(g,h,'c--',label='Solar')
    plt.plot(i,j,'m--',label='Hydroelectricity')
    plt.plot(l,m,'k--',label='Geothermal')
    plt.xlabel('Year')
    plt.ylabel(unit)
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.title(state_name+ " Energy Production"  + str(x.iloc[0]) + "-" + str(x.iloc[-1]))
    plt.show()
    
""" Part IV: Score Scheme"""

def LeastS(x,y):
    A = np.vstack([x, np.ones(len(x))]).T
    m, c = np.linalg.lstsq(A, y)[0]
    return m, c

## Calculate data for all states over 50 years
california_grades_50_years = [state_grade(california,i) for i in range(50)]
arizona_grades_50_years = [state_grade(arizona,i) for i in range(50)]
newmexico_grades_50_years =  [state_grade(newmexico,i) for i in range(50)]
texas_grades_50_years = [state_grade(texas,i) for i in range(50)]
years = np.array(range(1960,2010))

### Population arrays per state 1960-2009
## Population in California
cal_pop = np.array([ 15870, 16497, 17072, 17668, 18151, 18585, 18858, 19176, 19394, 19711,
    20007, 20346, 20585, 20869, 21174, 21538, 21936, 22352, 22836, 23257,
    23801, 24286, 24820, 25360, 25844, 26441, 27102, 27777, 28464, 29218,
     29960, 30471, 30975, 31275, 31484, 31697, 32019, 32486, 32988, 33499,
        33988, 34479, 34872, 35253, 35575, 35828, 36021, 36250, 36604, 36961])*1000
## Population in Arizona 1960 - 2009
ari_pop = np.array([1321, 1407, 1471, 1521, 1556, 1584, 1614, 1646, 1682,1737, 
 1792, 1896, 2008, 2124, 2223, 2285, 2346, 2425, 2515, 2636,
     2738, 2810, 2890, 2969, 3067, 3184, 3308, 3437, 3535, 3622,
         3684, 3789, 3916, 4065, 4245, 4432, 4587, 4737, 4883, 5024,
             5161, 5273, 5396, 5510, 5652, 5839, 6029, 6168, 6280, 6343])*1000
## Population in New Mexico 1960 - 2009
newm_pop = np.array([954, 965, 979, 989, 1006, 1012, 1007, 1000, 994, 1011,
1023, 1054, 1079, 1106, 1131, 1160, 1189, 1216, 1238, 1285,
1309, 1333, 1364, 1394, 1417, 1438, 1463, 1479, 1490, 1504,
1522, 1555, 1595, 1636, 1682, 1720, 1752, 1775, 1793, 1808,
1821, 1832, 1855, 1878, 1904, 1932, 1962, 1990, 2011, 2037])*1000
## Population in Texas 1960-2009
tex_pop = np.array([9624, 9820, 10053, 10159, 10270, 10378, 10492, 10599, 10819, 11045,
11236, 11510, 11759, 12020, 12269, 12569, 12904, 13193, 13500, 13888,
14338, 14746, 15331, 15752, 16007, 16273, 16561, 16622, 16667, 16807,
17057, 17398, 17760, 18162, 18564, 18959, 19340, 19740, 20158, 20558,
 20944, 21320, 21690, 22031, 22394, 22778, 23360, 23832, 24309, 24802])*1000

# Divide grades by population
end_grades_cal = california_grades_50_years/cal_pop
end_grades_ari = arizona_grades_50_years/ari_pop
end_grades_nm = newmexico_grades_50_years/newm_pop
end_grades_tx = texas_grades_50_years/tex_pop

def plot_scores():
    plt.plot(years,end_grades_cal,'bs',label='California Score')
    plt.plot(years,end_grades_ari,'gs',label='Arizona Score')
    plt.plot(years,end_grades_nm,'rs',label='New Mexico Score')
    plt.plot(years,end_grades_tx,'ys',label = 'Texas Score')
    plt.xlabel("Year")
    plt.ylabel("Scores")
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.title("State Scores between 1960 and 2009")
    
    plt.show()