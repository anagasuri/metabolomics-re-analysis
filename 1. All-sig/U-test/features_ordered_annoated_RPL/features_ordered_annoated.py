# -*- coding: utf-8 -*-
"""
Created on Tue Jul 22 16:51:36 2025

@author: Ji Xiaowen
"""
# is this the features ordered file generated from Clusterheatmap_sig_RPL_control_2.py

import pandas as pd

df1 = pd.read_csv('features_ordered_RPL.csv')

df2a = pd.read_csv('PosPVFStats2d.csv')
df2a.columns.values
df2a = df2a[['Alignment ID', 'U_stat', 'p_value', 'log2FC', 'cohens_d']]

df2b = pd.read_csv('NegPVFStats2d.csv')
df2b.columns.values
df2b = df2b[['Alignment ID', 'U_stat', 'p_value', 'log2FC', 'cohens_d']]

df2 = pd.concat([df2a, df2b], axis=0)


df = pd.merge(df1, df2, on = 'Alignment ID', how = 'left')

df3a = pd.read_csv('Pos_annoated_exogenous.csv')
df3a.columns.values
df3a = df3a[['Alignment ID','Adduct type', 'chemical_name']]

df3b = pd.read_csv('Neg_annoated_exogenous.csv')
df3b.columns.values
df3b = df3b[['Alignment ID','Adduct type', 'chemical_name']]

df3 = pd.concat([df3a, df3b], axis=0)


dfa = pd.merge(df, df3, on = 'Alignment ID', how = 'left')



df4 = pd.read_excel('pos_and_neg_endogenous_2.xlsx')
df4.columns.values
df4 = df4[['Alignment ID', 'chemical_name', 'mode']]

dfb = pd.merge(dfa, df4, on='Alignment ID', how='left', suffixes=('', '_df4'))
dfb['chemical_name'] = dfb['chemical_name'].combine_first(dfb['chemical_name_df4'])
dfb = dfb.drop(columns=['chemical_name_df4'])

dfb.to_csv('features_ordered_annotated.csv')





