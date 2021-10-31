#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 27 15:47:18 2021

@author: anon
"""

#!/usr/bin/env python
# coding: utf-8

# In[202]:


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
#import numpy as np
import statsmodels.stats.api as sms
#from scipy.optimize import curve_fit
#from scipy import stats
sns.set(color_codes=True)
sns.set(font_scale=3) 
plt.rcParams["figure.figsize"] = (30,10)

# In[203]:


df = pd.read_csv (r'~/Documents/MATLAB/Model_identifiability/HaptoSAADataTable.csv')

dfSAT1SAA = df.loc[df['SAT'] == 'SAT1SAA']
dfSAT2SAA = df.loc[df['SAT'] == 'SAT2SAA']
dfSAT3SAA = df.loc[df['SAT'] == 'SAT3SAA']
dfSAT1Hapto = df.loc[df['SAT'] == 'SAT1Hapto']
dfSAT2Hapto = df.loc[df['SAT'] == 'SAT2Hapto']
dfSAT3Hapto = df.loc[df['SAT'] == 'SAT3Hapto']
dfHapto = df.loc[df['InnateMeasure'] == 'Hapto']
dfSAA = df.loc[df['InnateMeasure' ] == 'SAA']
# In[ ]:





# In[207]:



h = sns.boxplot(x="ID", y=df['r'],
            hue="SAT",
            data=df, palette="Paired")
h.set(
    xlabel='Host', 
    ylabel='Range'
)
xticks = [.2, 2.2, 4.2, 6.2, 8.6,10.6,12.6,14.8,16.8,18.8]
h.set_xticks(xticks)
h.set_xticklabels(["ID 2", "ID 4", "ID 19", "ID 33", "ID 5", 'ID 9', 'ID 22','ID 12','ID 16','ID 17'])
sns.despine(offset=10, trim=True)
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.title('Viral Growth Rate')
plt.show()

h = sns.boxplot(x="SAT", y=df['r'],
            hue="SAT",
            data=df, palette="Paired")
h.set(
    xlabel='Host', 
    ylabel='Range'
)
xticks = [.2, 2.5, 4.8]
h.set_xticks(xticks)
h.set_xticklabels(["SAT 1", "SAT 2", "SAT 3"])
sns.despine(offset=10, trim=True)
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.title('Viral Growth Rate')
#plt.rcParams["figure.figsize"] = (30,10)
#plt.rcParams["figure.figsize"] = (30,10)
plt.show()


# In[208]:


S1, S2, S3 = dfSAT1SAA['r'], dfSAT2SAA['r'],dfSAT3SAA['r']
H1, H2, H3 = dfSAT1Hapto['r'], dfSAT2Hapto['r'],dfSAT3Hapto['r']

print('########################################################################')
print('\nviral growth rate\n')
print('difference across SATs, same measure (SAA)')
cm = sms.CompareMeans(sms.DescrStatsW(S1), sms.DescrStatsW(S2))
print(cm.tconfint_diff(usevar='unequal'))
cm = sms.CompareMeans(sms.DescrStatsW(S1), sms.DescrStatsW(S3))
print(cm.tconfint_diff(usevar='unequal'))
cm = sms.CompareMeans(sms.DescrStatsW(S2), sms.DescrStatsW(S3))
print(cm.tconfint_diff(usevar='unequal'))
print('\ndifference across SATs, same measure (Hapto)')
cm = sms.CompareMeans(sms.DescrStatsW(H1), sms.DescrStatsW(H2))
print(cm.tconfint_diff(usevar='unequal'))
cm = sms.CompareMeans(sms.DescrStatsW(H1), sms.DescrStatsW(H3))
print(cm.tconfint_diff(usevar='unequal'))
cm = sms.CompareMeans(sms.DescrStatsW(H2), sms.DescrStatsW(H3))
print(cm.tconfint_diff(usevar='unequal'))

print('\n\ndifference across measure, same SATs')
cm = sms.CompareMeans(sms.DescrStatsW(S1), sms.DescrStatsW(H1))
print(cm.tconfint_diff(usevar='unequal'))
cm = sms.CompareMeans(sms.DescrStatsW(S2), sms.DescrStatsW(H2))
print(cm.tconfint_diff(usevar='unequal'))
cm = sms.CompareMeans(sms.DescrStatsW(S3), sms.DescrStatsW(H3))
print(cm.tconfint_diff(usevar='unequal'))


# In[ ]:





# In[209]:


h = sns.boxplot(x="ID", y=df['P0'],
            hue="SAT",
            data=df, palette="Paired")
sns.despine(offset=10, trim=True)
h.set(
    xlabel='Host', 
    ylabel='Range'
)
xticks = [.2, 2.2, 4.2, 6.2, 8.6,10.6,12.6,14.8,16.8,18.8]
h.set_xticks(xticks)
h.set_xticklabels(["ID 2", "ID 4", "ID 19", "ID 33", "ID 5", 'ID 9', 'ID 22','ID 12','ID 16','ID 17'])
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.title('Initial Viral Load')
#plt.rcParams["figure.figsize"] = (30,10)
plt.show()

h = sns.boxplot(x="SAT", y=df['P0'],
            hue="SAT",
            data=df, palette="Paired")
h.set(
    xlabel='Host', 
    ylabel='Range'
)
xticks = [.2, 2.5, 4.8]
h.set_xticks(xticks)
h.set_xticklabels(["SAT 1", "SAT 2", "SAT 3"])
sns.despine(offset=10, trim=True)
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.title('Initial Viral Load')
#plt.rcParams["figure.figsize"] = (30,10)
#plt.rcParams["figure.figsize"] = (30,10)
plt.show()


# In[ ]:


S1, S2, S3 = dfSAT1SAA['P0'], dfSAT2SAA['P0'],dfSAT3SAA['P0']
H1, H2, H3 = dfSAT1Hapto['P0'], dfSAT2Hapto['P0'],dfSAT3Hapto['P0']

print('########################################################################')
print('\nInitial Viral Load\n')
print('difference across SATs, same measure (SAA)')
cm = sms.CompareMeans(sms.DescrStatsW(S1), sms.DescrStatsW(S2))
print(cm.tconfint_diff(usevar='unequal'))
cm = sms.CompareMeans(sms.DescrStatsW(S1), sms.DescrStatsW(S3))
print(cm.tconfint_diff(usevar='unequal'))
cm = sms.CompareMeans(sms.DescrStatsW(S2), sms.DescrStatsW(S3))
print(cm.tconfint_diff(usevar='unequal'))
print('\ndifference across SATs, same measure (Hapto)')
cm = sms.CompareMeans(sms.DescrStatsW(H1), sms.DescrStatsW(H2))
print(cm.tconfint_diff(usevar='unequal'))
cm = sms.CompareMeans(sms.DescrStatsW(H1), sms.DescrStatsW(H3))
print(cm.tconfint_diff(usevar='unequal'))
cm = sms.CompareMeans(sms.DescrStatsW(H2), sms.DescrStatsW(H3))
print(cm.tconfint_diff(usevar='unequal'))

print('\n\ndifference across measure, same SATs')
cm = sms.CompareMeans(sms.DescrStatsW(S1), sms.DescrStatsW(H1))
print(cm.tconfint_diff(usevar='unequal'))
cm = sms.CompareMeans(sms.DescrStatsW(S2), sms.DescrStatsW(H2))
print(cm.tconfint_diff(usevar='unequal'))
cm = sms.CompareMeans(sms.DescrStatsW(S3), sms.DescrStatsW(H3))
print(cm.tconfint_diff(usevar='unequal'))


# In[ ]:


h = sns.boxplot(x="ID", y=df['k'],
            hue="SAT",
            data=df, palette="Paired")
h.set(
    xlabel='Host', 
    ylabel='Range'
)
xticks = [.2, 2.2, 4.2, 6.2, 8.6,10.6,12.6,14.8,16.8,18.8]
h.set_xticks(xticks)
h.set_xticklabels(["ID 2", "ID 4", "ID 19", "ID 33", "ID 5", 'ID 9', 'ID 22','ID 12','ID 16','ID 17'])
sns.despine(offset=10, trim=True)
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.title('Activation Rate of Innate Response in Presence of Pathogen')
#plt.rcParams["figure.figsize"] = (30,10)
plt.show()

h = sns.boxplot(x="SAT", y=df['k'],
            hue="SAT",
            data=df, palette="Paired")
h.set(
    xlabel='Host', 
    ylabel='Range'
)
xticks = [.2, 2.5, 4.8]
h.set_xticks(xticks)
h.set_xticklabels(["SAT 1", "SAT 2", "SAT 3"])
sns.despine(offset=10, trim=True)
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.title('Activation Rate of Innate Response in Presence of Pathogen')
#plt.rcParams["figure.figsize"] = (30,10)
#plt.rcParams["figure.figsize"] = (30,10)
plt.show()


# In[ ]:


S1, S2, S3 = dfSAT1SAA['k'], dfSAT2SAA['k'],dfSAT3SAA['k']
H1, H2, H3 = dfSAT1Hapto['k'], dfSAT2Hapto['k'],dfSAT3Hapto['k']

print('########################################################################')
print('\nActivation rate of Innate Response in Presence of Pathogen\n')
print('difference across SATs, same measure (SAA)')
cm = sms.CompareMeans(sms.DescrStatsW(S1), sms.DescrStatsW(S2))
print(cm.tconfint_diff(usevar='unequal'))
cm = sms.CompareMeans(sms.DescrStatsW(S1), sms.DescrStatsW(S3))
print(cm.tconfint_diff(usevar='unequal'))
cm = sms.CompareMeans(sms.DescrStatsW(S2), sms.DescrStatsW(S3))
print(cm.tconfint_diff(usevar='unequal'))
print('\ndifference across SATs, same measure (Hapto)')
cm = sms.CompareMeans(sms.DescrStatsW(H1), sms.DescrStatsW(H2))
print(cm.tconfint_diff(usevar='unequal'))
cm = sms.CompareMeans(sms.DescrStatsW(H1), sms.DescrStatsW(H3))
print(cm.tconfint_diff(usevar='unequal'))
cm = sms.CompareMeans(sms.DescrStatsW(H2), sms.DescrStatsW(H3))
print(cm.tconfint_diff(usevar='unequal'))

print('\n\ndifference across measure, same SATs')
cm = sms.CompareMeans(sms.DescrStatsW(S1), sms.DescrStatsW(H1))
print(cm.tconfint_diff(usevar='unequal'))
cm = sms.CompareMeans(sms.DescrStatsW(S2), sms.DescrStatsW(H2))
print(cm.tconfint_diff(usevar='unequal'))
cm = sms.CompareMeans(sms.DescrStatsW(S3), sms.DescrStatsW(H3))
print(cm.tconfint_diff(usevar='unequal'))


# 

# In[ ]:


h = sns.boxplot(x="ID", y=df['K'],
            hue="SAT",
            data=df, palette="Paired")
h.set(
    xlabel='Host', 
    ylabel='Range'
)
xticks = [.2, 2.2, 4.2, 6.2, 8.6,10.6,12.6,14.8,16.8,18.8]
h.set_xticks(xticks)
h.set_xticklabels(["ID 2", "ID 4", "ID 19", "ID 33", "ID 5", 'ID 9', 'ID 22','ID 12','ID 16','ID 17'])
sns.despine(offset=10, trim=True)
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.title('Viral Carrying Capacity')
#plt.rcParams["figure.figsize"] = (30,10)
plt.show()

h = sns.boxplot(x="SAT", y=df['K'],
            hue="SAT",
            data=df, palette="Paired")
h.set(
    xlabel='Host', 
    ylabel='Range'
)
xticks = [.2, 2.5, 4.8]
h.set_xticks(xticks)
h.set_xticklabels(["SAT 1", "SAT 2", "SAT 3"])
sns.despine(offset=10, trim=True)
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.title('Viral Carrying Capacity')
#plt.rcParams["figure.figsize"] = (30,10)
#plt.rcParams["figure.figsize"] = (30,10)
plt.show()


# In[ ]:


S1, S2, S3 = dfSAT1SAA['K'], dfSAT2SAA['K'],dfSAT3SAA['K']
H1, H2, H3 = dfSAT1Hapto['K'], dfSAT2Hapto['K'],dfSAT3Hapto['K']

print('########################################################################')
print('\nViral Carrying Capacity\n')

print('difference across SATs, same measure (SAA)')
cm = sms.CompareMeans(sms.DescrStatsW(S1), sms.DescrStatsW(S2))
print(cm.tconfint_diff(usevar='unequal'))
cm = sms.CompareMeans(sms.DescrStatsW(S1), sms.DescrStatsW(S3))
print(cm.tconfint_diff(usevar='unequal'))
cm = sms.CompareMeans(sms.DescrStatsW(S2), sms.DescrStatsW(S3))
print(cm.tconfint_diff(usevar='unequal'))
print('\ndifference across SATs, same measure (Hapto)')
cm = sms.CompareMeans(sms.DescrStatsW(H1), sms.DescrStatsW(H2))
print(cm.tconfint_diff(usevar='unequal'))
cm = sms.CompareMeans(sms.DescrStatsW(H1), sms.DescrStatsW(H3))
print(cm.tconfint_diff(usevar='unequal'))
cm = sms.CompareMeans(sms.DescrStatsW(H2), sms.DescrStatsW(H3))
print(cm.tconfint_diff(usevar='unequal'))

print('\n\ndifference across measure, same SATs')
cm = sms.CompareMeans(sms.DescrStatsW(S1), sms.DescrStatsW(H1))
print(cm.tconfint_diff(usevar='unequal'))
cm = sms.CompareMeans(sms.DescrStatsW(S2), sms.DescrStatsW(H2))
print(cm.tconfint_diff(usevar='unequal'))
cm = sms.CompareMeans(sms.DescrStatsW(S3), sms.DescrStatsW(H3))
print(cm.tconfint_diff(usevar='unequal'))


# In[ ]:


h = sns.boxplot(x="ID", y=df['delta'],
            hue="SAT",
            data=df, palette="Paired")
h.set(
    xlabel='Host', 
    ylabel='Range'
)
xticks = [.2, 2.2, 4.2, 6.2, 8.6,10.6,12.6,14.8,16.8,18.8]
h.set_xticks(xticks)
h.set_xticklabels(["ID 2", "ID 4", "ID 19", "ID 33", "ID 5", 'ID 9', 'ID 22','ID 12','ID 16','ID 17'])
sns.despine(offset=10, trim=True)
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.title('Killing Rate of Pathogen By Adaptive Immune Response')
#plt.rcParams["figure.figsize"] = (30,10)
plt.show()

h = sns.boxplot(x="SAT", y=df['delta'],
            hue="SAT",
            data=df, palette="Paired")
h.set(
    xlabel='Host', 
    ylabel='Range'
)
xticks = [.2, 2.5, 4.8]
h.set_xticks(xticks)
h.set_xticklabels(["SAT 1", "SAT 2", "SAT 3"])
sns.despine(offset=10, trim=True)
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.title('Killing Rate of Pathogen By Adaptive Immune Response')
#plt.rcParams["figure.figsize"] = (30,10)
#plt.rcParams["figure.figsize"] = (30,10)
plt.show()


# In[ ]:


S1, S2, S3 = dfSAT1SAA['delta'], dfSAT2SAA['delta'],dfSAT3SAA['delta']
H1, H2, H3 = dfSAT1Hapto['delta'], dfSAT2Hapto['delta'],dfSAT3Hapto['delta']


print('########################################################################')
print('\nKilling of Pathogen by Adaptive Response\n')

print('difference across SATs, same measure (SAA)')
cm = sms.CompareMeans(sms.DescrStatsW(S1), sms.DescrStatsW(S2))
print(cm.tconfint_diff(usevar='unequal'))
cm = sms.CompareMeans(sms.DescrStatsW(S1), sms.DescrStatsW(S3))
print(cm.tconfint_diff(usevar='unequal'))
cm = sms.CompareMeans(sms.DescrStatsW(S2), sms.DescrStatsW(S3))
print(cm.tconfint_diff(usevar='unequal'))
print('\ndifference across SATs, same measure (Hapto)')
cm = sms.CompareMeans(sms.DescrStatsW(H1), sms.DescrStatsW(H2))
print(cm.tconfint_diff(usevar='unequal'))
cm = sms.CompareMeans(sms.DescrStatsW(H1), sms.DescrStatsW(H3))
print(cm.tconfint_diff(usevar='unequal'))
cm = sms.CompareMeans(sms.DescrStatsW(H2), sms.DescrStatsW(H3))
print(cm.tconfint_diff(usevar='unequal'))

print('\n\ndifference across measure, same SATs')
cm = sms.CompareMeans(sms.DescrStatsW(S1), sms.DescrStatsW(H1))
print(cm.tconfint_diff(usevar='unequal'))
cm = sms.CompareMeans(sms.DescrStatsW(S2), sms.DescrStatsW(H2))
print(cm.tconfint_diff(usevar='unequal'))
cm = sms.CompareMeans(sms.DescrStatsW(S3), sms.DescrStatsW(H3))
print(cm.tconfint_diff(usevar='unequal'))


# In[ ]:



h = sns.boxplot(x="ID", y=df['b'],
            hue="SAT",
            data=df, palette="Paired")
h.set(
    xlabel='Host', 
    ylabel='Range'
)
xticks = [.2, 2.2, 4.2, 6.2, 8.6,10.6,12.6,14.8,16.8,18.8]
h.set_xticks(xticks)
h.set_xticklabels(["ID 2", "ID 4", "ID 19", "ID 33", "ID 5", 'ID 9', 'ID 22','ID 12','ID 16','ID 17'])
sns.despine(offset=10, trim=True)
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.title('Activation Rate of Adaptive Immune Response in the Presence of Pathogen')
#plt.rcParams["figure.figsize"] = (30,10)
plt.show()

h = sns.boxplot(x="SAT", y=df['b'],
            hue="SAT",
            data=df, palette="Paired")
h.set(
    xlabel='Host', 
    ylabel='Range'
)
xticks = [.2, 2.5, 4.8]
h.set_xticks(xticks)
h.set_xticklabels(["SAT 1", "SAT 2", "SAT 3"])
sns.despine(offset=10, trim=True)
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.title('Activation Rate of Adaptive Immune Response in the Presence of Pathogen')
#plt.rcParams["figure.figsize"] = (30,10)
plt.show()


# In[ ]:


S1, S2, S3 = dfSAT1SAA['b'], dfSAT2SAA['b'],dfSAT3SAA['b']
H1, H2, H3 = dfSAT1Hapto['b'], dfSAT2Hapto['b'],dfSAT3Hapto['b']


print('########################################################################')
print('\nKActivation of Adaptive Immune Response\n')


print('difference across SATs, same measure (SAA)')
cm = sms.CompareMeans(sms.DescrStatsW(S1), sms.DescrStatsW(S2))
print(cm.tconfint_diff(usevar='unequal'))
cm = sms.CompareMeans(sms.DescrStatsW(S1), sms.DescrStatsW(S3))
print(cm.tconfint_diff(usevar='unequal'))
cm = sms.CompareMeans(sms.DescrStatsW(S2), sms.DescrStatsW(S3))
print(cm.tconfint_diff(usevar='unequal'))
print('\ndifference across SATs, same measure (Hapto)')
cm = sms.CompareMeans(sms.DescrStatsW(H1), sms.DescrStatsW(H2))
print(cm.tconfint_diff(usevar='unequal'))
cm = sms.CompareMeans(sms.DescrStatsW(H1), sms.DescrStatsW(H3))
print(cm.tconfint_diff(usevar='unequal'))
cm = sms.CompareMeans(sms.DescrStatsW(H2), sms.DescrStatsW(H3))
print(cm.tconfint_diff(usevar='unequal'))

print('\n\ndifference across measure, same SATs')
cm = sms.CompareMeans(sms.DescrStatsW(S1), sms.DescrStatsW(H1))
print(cm.tconfint_diff(usevar='unequal'))
cm = sms.CompareMeans(sms.DescrStatsW(S2), sms.DescrStatsW(H2))
print(cm.tconfint_diff(usevar='unequal'))
cm = sms.CompareMeans(sms.DescrStatsW(S3), sms.DescrStatsW(H3))
print(cm.tconfint_diff(usevar='unequal'))


# In[ ]:


h = sns.violinplot(x="ID", y=df['TimePeak'],
            hue="SAT",
            data=df, palette="Paired")
h.set(
    xlabel='Host', 
    ylabel='Range'
)
xticks = [.2, 2.2, 4.2, 6.2, 8.6,10.6,12.6,14.8,16.8,18.8]
h.set_xticks(xticks)
h.set_xticklabels(["ID 2", "ID 4", "ID 19", "ID 33", "ID 5", 'ID 9', 'ID 22','ID 12','ID 16','ID 17'])
sns.despine(offset=10, trim=True)
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.title('Time to Peak Viral Load')
#plt.rcParams["figure.figsize"] = (30,10)
plt.show()

h = sns.boxplot(x="SAT", y=df['TimePeak'],
            hue="SAT",
            data=df, palette="Paired")
h.set(
    xlabel='Host', 
    ylabel='Range'
)
xticks = [.2, 2.5, 4.8]
h.set_xticks(xticks)
h.set_xticklabels(["SAT 1", "SAT 2", "SAT 3"])
sns.despine(offset=10, trim=True)
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.title('Time to Peak Viral Load (days)')
#plt.rcParams["figure.figsize"] = (30,10)
plt.show()


# In[ ]:


S1, S2, S3 = dfSAT1SAA['TimePeak'], dfSAT2SAA['TimePeak'],dfSAT3SAA['TimePeak']
H1, H2, H3 = dfSAT1Hapto['TimePeak'], dfSAT2Hapto['TimePeak'],dfSAT3Hapto['TimePeak']

print('########################################################################')
print('\nTme to Peak Viral Load\n')


print('difference across SATs, same measure (SAA)')
cm = sms.CompareMeans(sms.DescrStatsW(S1), sms.DescrStatsW(S2))
print(cm.tconfint_diff(usevar='unequal'))
cm = sms.CompareMeans(sms.DescrStatsW(S1), sms.DescrStatsW(S3))
print(cm.tconfint_diff(usevar='unequal'))
cm = sms.CompareMeans(sms.DescrStatsW(S2), sms.DescrStatsW(S3))
print(cm.tconfint_diff(usevar='unequal'))
print('\ndifference across SATs, same measure (Hapto)')
cm = sms.CompareMeans(sms.DescrStatsW(H1), sms.DescrStatsW(H2))
print(cm.tconfint_diff(usevar='unequal'))
cm = sms.CompareMeans(sms.DescrStatsW(H1), sms.DescrStatsW(H3))
print(cm.tconfint_diff(usevar='unequal'))
cm = sms.CompareMeans(sms.DescrStatsW(H2), sms.DescrStatsW(H3))
print(cm.tconfint_diff(usevar='unequal'))

print('\n\ndifference across measure, same SATs')
cm = sms.CompareMeans(sms.DescrStatsW(S1), sms.DescrStatsW(H1))
print(cm.tconfint_diff(usevar='unequal'))
cm = sms.CompareMeans(sms.DescrStatsW(S2), sms.DescrStatsW(H2))
print(cm.tconfint_diff(usevar='unequal'))
cm = sms.CompareMeans(sms.DescrStatsW(S3), sms.DescrStatsW(H3))
print(cm.tconfint_diff(usevar='unequal'))


# In[ ]:


h = sns.boxplot(x="ID", y=df['CumViral'],
            hue="SAT",
            data=df, palette="Paired")
h.set(
    xlabel='Host', 
    ylabel='Range'
)
xticks = [.2, 2.2, 4.2, 6.2, 8.6,10.6,12.6,14.8,16.8,18.8]
h.set_xticks(xticks)
h.set_xticklabels(["ID 2", "ID 4", "ID 19", "ID 33", "ID 5", 'ID 9', 'ID 22','ID 12','ID 16','ID 17'])
sns.despine(offset=10, trim=True)
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.title('Cum. Viral Load')
#plt.rcParams["figure.figsize"] = (30,10)
plt.show()

h = sns.boxplot(x="SAT", y=df['CumViral'],
            hue="SAT",
            data=df, palette="Paired")
h.set(
    xlabel='Host', 
    ylabel='Range'
)
xticks = [.2, 2.5, 4.8]
h.set_xticks(xticks)
h.set_xticklabels(["SAT 1", "SAT 2", "SAT 3"])
sns.despine(offset=10, trim=True)
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.title('Cum. Viral Load')
#plt.rcParams["figure.figsize"] = (30,10)
plt.show()


# In[ ]:


S1, S2, S3 = dfSAT1SAA['CumViral'], dfSAT2SAA['CumViral'],dfSAT3SAA['CumViral']
H1, H2, H3 = dfSAT1Hapto['CumViral'], dfSAT2Hapto['CumViral'],dfSAT3Hapto['CumViral']


print('########################################################################')
print('\nCum. Viral Load\n')


print('difference across SATs, same measure (SAA)')
cm = sms.CompareMeans(sms.DescrStatsW(S1), sms.DescrStatsW(S2))
print(cm.tconfint_diff(usevar='unequal'))
cm = sms.CompareMeans(sms.DescrStatsW(S1), sms.DescrStatsW(S3))
print(cm.tconfint_diff(usevar='unequal'))
cm = sms.CompareMeans(sms.DescrStatsW(S2), sms.DescrStatsW(S3))
print(cm.tconfint_diff(usevar='unequal'))
print('\ndifference across SATs, same measure (Hapto)')
cm = sms.CompareMeans(sms.DescrStatsW(H1), sms.DescrStatsW(H2))
print(cm.tconfint_diff(usevar='unequal'))
cm = sms.CompareMeans(sms.DescrStatsW(H1), sms.DescrStatsW(H3))
print(cm.tconfint_diff(usevar='unequal'))
cm = sms.CompareMeans(sms.DescrStatsW(H2), sms.DescrStatsW(H3))
print(cm.tconfint_diff(usevar='unequal'))

print('\n\ndifference across measure, same SATs')
cm = sms.CompareMeans(sms.DescrStatsW(S1), sms.DescrStatsW(H1))
print(cm.tconfint_diff(usevar='unequal'))
cm = sms.CompareMeans(sms.DescrStatsW(S2), sms.DescrStatsW(H2))
print(cm.tconfint_diff(usevar='unequal'))
cm = sms.CompareMeans(sms.DescrStatsW(S3), sms.DescrStatsW(H3))
print(cm.tconfint_diff(usevar='unequal'))


# In[ ]:
    
    
h = sns.boxplot(x="ID", y=df['I0'],
            hue="SAT",
            data=dfSAA,palette=['blue','green','red'])
h.set(
    xlabel='Host', 
    ylabel='Range'
)
xticks = [0, 1,2,3, 4,5,6,7,8,9]
h.set_xticks(xticks)
h.set_xticklabels(["ID 2", "ID 4", "ID 19", "ID 33", "ID 5", 'ID 9', 'ID 22','ID 12','ID 16','ID 17'])
sns.despine(offset=10, trim=True)
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.title('Initial SAA Level')
#plt.rcParams["figure.figsize"] = (30,10)
plt.show()

h = sns.boxplot(x="SAT", y=df['I0'],
            hue="SAT",
            data=dfSAA, palette=['blue','green','red'])
h.set(
    xlabel='Host', 
    ylabel='Range'
)
xticks = [-.2, 1,2.2]
h.set_xticks(xticks)
h.set_xticklabels(["SAT 1", "SAT 2", "SAT 3"])
sns.despine(offset=10, trim=True)
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.title('Initial SAA Level')
#plt.rcParams["figure.figsize"] = (30,10)
plt.show()

# In[ ]:
    
    
h = sns.boxplot(x="ID", y=df['I0'],
            hue="SAT",
            data=dfHapto,palette=['blue','green','red'])
h.set(
    xlabel='Host', 
    ylabel='Range'
)
xticks = [0, 1,2,3, 4,5,6,7,8,9]
h.set_xticks(xticks)
h.set_xticklabels(["ID 2", "ID 4", "ID 19", "ID 33", "ID 5", 'ID 9', 'ID 22','ID 12','ID 16','ID 17'])
sns.despine(offset=10, trim=True)
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.title('Initial Haptoglobin Level')
#plt.rcParams["figure.figsize"] = (30,10)
plt.show()

h = sns.boxplot(x="SAT", y=df['I0'],
            hue="SAT",
            data=dfHapto, palette=['blue','green','red'])
h.set(
    xlabel='Host', 
    ylabel='Range'
)
xticks = [-.2, 1,2.2]
h.set_xticks(xticks)
h.set_xticklabels(["SAT 1", "SAT 2", "SAT 3"])
sns.despine(offset=10, trim=True)
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.title('Initial Haptoglobin Level')
#plt.rcParams["figure.figsize"] = (30,10)
plt.show()


S1, S2, S3 = dfSAT1SAA['I0'], dfSAT2SAA['I0'],dfSAT3SAA['I0']
H1, H2, H3 = dfSAT1Hapto['I0'], dfSAT2Hapto['I0'],dfSAT3Hapto['I0']


print('########################################################################')
print('\nInitial Innate Measue Level\n')


print('difference across SATs, same measure (SAA)')
cm = sms.CompareMeans(sms.DescrStatsW(S1), sms.DescrStatsW(S2))
print(cm.tconfint_diff(usevar='unequal'))
cm = sms.CompareMeans(sms.DescrStatsW(S1), sms.DescrStatsW(S3))
print(cm.tconfint_diff(usevar='unequal'))
cm = sms.CompareMeans(sms.DescrStatsW(S2), sms.DescrStatsW(S3))
print(cm.tconfint_diff(usevar='unequal'))
print('\ndifference across SATs, same measure (Hapto)')
cm = sms.CompareMeans(sms.DescrStatsW(H1), sms.DescrStatsW(H2))
print(cm.tconfint_diff(usevar='unequal'))
cm = sms.CompareMeans(sms.DescrStatsW(H1), sms.DescrStatsW(H3))
print(cm.tconfint_diff(usevar='unequal'))
cm = sms.CompareMeans(sms.DescrStatsW(H2), sms.DescrStatsW(H3))
print(cm.tconfint_diff(usevar='unequal'))

