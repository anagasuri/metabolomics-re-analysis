#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 20 19:50:34 2024

@author: dabrahamsson
"""

# edited and rerun by AMRITA NAGASURI 
# ran once for pos and once for neg to create PosPVFStats2c_NEW.csv and NegPVFStats2c_NEW.csv

import numpy as np
import pandas as pd
import scipy
from scipy import stats
from scipy.stats import linregress
import seaborn as sns
import matplotlib as mlp
import matplotlib.pyplot as plt
import os

# output directory
output_dir = '/Users/amritanagasuri/Desktop/new Q-TOF Analysis/4. Cluster/1. All-sig/regress output'
os.makedirs(output_dir, exist_ok=True)

df = pd.read_csv('/Users/amritanagasuri/Desktop/new Q-TOF Analysis/pos, neg, corrected data/Neg_corrected_data_WITHREPS.csv')
ga = pd.read_csv('/Users/amritanagasuri/Desktop/new Q-TOF Analysis/final pos, neg, pheno/pheno-clean-NEWWITHREPS.csv')

dfa = df.loc[:, 'chem_id':'MS/MS spectrum']

base_name = sorted(set(col[:-3] for col in df.columns if col.endswith('_r1')))
dfb = pd.DataFrame()
for name in base_name:
    dfb[name] = df[f'{name}_r1'].add(df[f'{name}_r2'], fill_value=0) / 2
df = pd.concat([dfa, dfb], axis = 1)


ga['Sample Name'] = ga['Sample Name'].str.replace('_r1', '').str.replace('_r2', '')
ga = ga.drop_duplicates(subset=['Sample Name'])

df.columns.values
ga.columns.values

loc1 = df.loc[:, 'chem_id':'MS/MS spectrum']
loc3 = df.loc[:, df.columns.str.contains('BI')]

def stat(samples):
    df = pd.concat([loc1, samples], axis=1)
    
    print(df.shape)
    
    ids = samples.columns.values.tolist()
    mt = pd.melt(df, id_vars=['Alignment ID', 'Average Mz', 'Average Rt(min)'], value_vars=ids, value_name='Abundance')
    mt.columns = mt.columns.str.replace('variable', 'Sample Name')
    
    mt = pd.merge(mt, ga, on='Sample Name', how='left')
    # replaced mt['Type']
    mt['PVF'] = np.where(mt['age'] == 'C', 0, 1)
    mt['logA'] = np.log10(mt['Abundance'])
    
    def get_element(my_list, Position):
        return my_list[Position]
    
    mt['PVF'] = mt['PVF'].astype(str)
    mt['Sample Name'] = mt['Sample Name'].astype(str)
    mt['identifier'] = mt['PVF'] + '_' + mt['Sample Name']
    mt = mt.drop(['PVF','Sample Name'], axis=1)
    mt = mt.pivot_table('logA','Alignment ID','identifier')
    
    axisvalues = mt.columns.values
    axisvalues_mt = pd.DataFrame({'identifier': axisvalues})
    axisvalues_ = axisvalues_mt['identifier'].str.split('_').apply(lambda x: x[0])
    axisvalues_ = axisvalues_.astype(float)
    
    mt = mt.astype(float)
    def calc_slope(row):
        a = scipy.stats.linregress(axisvalues_ , row)
        return pd.Series(a._asdict())
    
    print (mt.apply(calc_slope,axis=1))
    
    mt = mt.join(mt.apply(calc_slope,axis=1))
    return(mt)


mts = stat(loc3)
mts = mts.reset_index()
mts = pd.concat([loc1, mts], axis=1)

#Fold-change
#mts.columns.values
#controls = mts.loc[:, mts.columns.str.startswith('0_')]
#cases = mts.loc[:, mts.columns.str.startswith('1_')]
#controls = 10**controls
#cases = 10**cases
#conav = controls.mean(axis=1)
#casav = cases.mean(axis=1)
#old = casav/conav
#mts['log2fold'] = np.log2(fold)

#Benjamini-Hochberg filtering p-values
mts = mts.sort_values(by='pvalue')

mts['R2'] = mts['rvalue']**2
#mts = mts[mts['R2'] > 0.1]
mts['rank'] = mts.reset_index().index + 1
mts['(I/m)Q'] = (mts['rank']/len(mts))*0.05
mts['(I/m)Q - p'] = mts['(I/m)Q'] - mts['pvalue']
mts = mts.sort_values(by='(I/m)Q - p', ascending =False)
mts['BH_sig'] = np.where(mts['(I/m)Q - p'] < 0, 0, 1)

mts = mts[mts['pvalue'] < 0.05]
mts = mts[mts['BH_sig'] == 1]
mts.to_csv(os.path.join(output_dir, 'NegPVFStats2c_NEW.csv'))
print('Significant after BH:', mts['BH_sig'].sum())
