import numpy as np
import pandas as pd
import scipy
from scipy import stats
import re

# Load data
df = pd.read_csv('Pos_removebatcheffect_corrected.csv', index_col=0)  # index is Alignment ID
ga = pd.read_csv('pheno-clean.csv')

base_name = sorted(set(col[:-3] for col in df.columns if col.endswith('_r1')))
dfb = pd.DataFrame()
for name in base_name:
    dfb[name] = df[f'{name}_r1'].add(df[f'{name}_r2'], fill_value=0) / 2
df = dfb

# Keep only sample columns
sample_data = df.copy()
sample_data.columns = sample_data.columns.astype(str)  # Ensure matching to ga

# Filter phenotype data to match sample columns
ga['Sample Name'] = ga['Sample Name'].str.replace('_r1', '').str.replace('_r2', '')
ga = ga.drop_duplicates(subset=['Sample Name'])
ga = ga[ga['Sample Name'].isin(sample_data.columns)]


df3 = pd.read_excel('Name Match.xlsx')
df3['Sample Name'] = df3['Sample Name'].str.replace(r'_r1.d', '').str.replace(r'_r2.d', '')
df3 = df3.drop_duplicates(subset=['Sample Name'])
df3 = df3[df3['Sample Name'].str.contains('BI', na=False)]


#change df1 name as well based on df3
rename_dict = dict(zip(df3['Sample Name'], df3['Sample Name_2']))
ga.loc[ga['Sample Name'].str.contains('BI', na=False), 'Sample Name'] = \
    ga['Sample Name'].replace(rename_dict)
ga = ga.drop_duplicates(subset=['Sample Name'])
ga['Sample Name'] = ga['Sample Name'].str.replace(r'(UC-CTRL-)(0+)(\d+)', r'\1\3', regex=True)
ga['Sample Name'] = ga['Sample Name'].str.replace(r'^.*?(UC-CTRL)', r'\1', regex=True)
ga['Sample Name'] = ga['Sample Name'].str.replace(r'CTRL-', 'CTRL', regex=True)
ga['Sample Name'] = ga['Sample Name'].str.replace(r'(?<=UC-)(0+)', '', regex=True)


# rename df columns name
rename_dict = dict(zip(df3['Sample Name'], df3['Sample Name_2']))
columns_to_rename = [col for col in df.columns if 'BI' in col]
filtered_rename_dict = {col: rename_dict[col] for col in columns_to_rename if col in rename_dict}
df = df.rename(columns=filtered_rename_dict)

# further change UC-Ctrl name 
df.columns = df.columns.str.replace(r'^.*?(UC-CTRL)', r'\1', regex=True)
df.columns = [re.sub(r'UC-(0+)(\d+)A', r'UC-\2A', col) if col.startswith("UC-") else col for col in df.columns]



# bioinfo match
df2 = pd.read_excel('Demographics_EPL_final.xlsx')
df2['StudyID'] = df2['StudyID'].str.replace(r'UC-0', 'UC-', regex=True)
df2['StudyID'] = df2['StudyID'].str.replace(r'(UC-CTRL-)(0+)(\d+)', r'\1\3', regex=True)
df2['StudyID'] = df2['StudyID'].str.replace(r'UC-CTRL-', 'UC-CTRL', regex=True)
df2 = df2.rename(columns = {'StudyID':'Sample Name'})
df2 = df2.drop_duplicates(subset=['Sample Name'])
df2 = df2[df2['Case'] == 1]
df2 = df2[df2['POC_karyotype'] != 99]

# Keep only samples with metadata
sample_data = df.loc[:, df.columns.isin(df2['Sample Name'])]



# Statistical test function
def stat(samples):
    df = samples.copy()
    df['Alignment ID'] = df.index
    mt = pd.melt(df, id_vars=['Alignment ID'], var_name='Sample Name', value_name='Abundance')
    
    mt = pd.merge(mt, df2, on='Sample Name', how='left')
    mt['PVF'] = np.where(mt['POC_karyotype'] == 0, 0, 1)
    mt['logA'] = np.log10(mt['Abundance'])

    mt['PVF'] = mt['PVF'].astype(str)
    mt['identifier'] = mt['PVF'] + '_' + mt['Sample Name']
    mt = mt.pivot_table('logA', 'Alignment ID', 'identifier')

    group_labels = pd.Series(mt.columns.str.split('_').str[0].astype(float), index=mt.columns)

    def calc_stats(row):
        group0 = row[group_labels == 0].dropna()
        group1 = row[group_labels == 1].dropna()

        if len(group0) < 2 or len(group1) < 2:
            return pd.Series({'U_stat': np.nan, 'p_value': np.nan, 'log2FC': np.nan, 'cohens_d': np.nan})

        U, p = scipy.stats.mannwhitneyu(group0, group1, alternative='two-sided')

        # log2FC (from log10)
        log2fc = (group1.mean() - group0.mean()) * np.log10(np.e) / np.log2(np.e)

        # Cohen's d
        mean_diff = group1.mean() - group0.mean()
        pooled_std = np.sqrt(((group0.std() ** 2 + group1.std() ** 2) / 2))
        cohens_d = mean_diff / pooled_std if pooled_std != 0 else np.nan

        return pd.Series({
            'U_stat': U,
            'p_value': p,
            'log2FC': log2fc,
            'cohens_d': cohens_d
        })

    results = mt.apply(calc_stats, axis=1)
    return results

# Run stats
results = stat(sample_data)
results = results.reset_index()

# Combine with original abundances
mts = pd.concat([df.reset_index(), results], axis=1)

from statsmodels.stats.multitest import multipletests
# Benjamini-Hochberg FDR
mts['p_adj'] = multipletests(mts['p_value'], method='fdr_bh')[1]
mts['BH_sig'] = np.where(mts['p_adj'] < 0.05, 1, 0)

# subset significant results
mts_2 = mts[mts['BH_sig'] == 1]

# Save
mts_2.to_csv('Aneuploid_PosPVFStats2d.csv', index=False)
