# -*- coding: utf-8 -*-
"""
Created on Thu Jul 17 09:59:55 2025

@author: Ji Xiaowen
"""
# Amrita basically made these changes within the sig_loss_control.py file when trying to make it work so she did not use this in re-analysis 

import pandas as pd

df1 = pd.read_csv('pheno-clean.csv')
df1['Sample Name'] = df1['Sample Name'].str.replace(r'_r1', '').str.replace(r'_r2', '')
df1 = df1[['Sample Name', 'batch']]
df1 = df1.drop_duplicates(subset = 'Sample Name')


df2 = pd.read_excel('Name Match.xlsx')
df2['Sample Name'] = df2['Sample Name'].str.replace(r'_r1.d', '').str.replace(r'_r2.d', '')
df2 = df2[df2['Sample Name'].str.contains('BI')]
df2 = df2.drop_duplicates(subset=['Sample Name'])

df = pd.merge(df1, df2, on = 'Sample Name', how = 'left')

df3 = pd.read_excel('Demographics_EPL_final.xlsx')

df4 = pd.merge(df, df3, left_on='Sample Name_2', right_on='StudyID', how='left')

df4.to_csv('pheno-clean_2.csv')