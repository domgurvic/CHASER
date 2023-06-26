import pandas as pd
import numpy as np
from scipy import stats

import matplotlib.pyplot as plt

import seaborn as sns

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
# from sklearn.metrics import silhouette_samples, silhouette_score
# from sklearn.cluster import KMeans
# from sklearn import datasets, decomposition
from sklearn.manifold import TSNE

#chem

from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, MACCSkeys, Descriptors, Descriptors3D, Draw, rdMolDescriptors, Draw, PandasTools, rdFingerprintGenerator
from rdkit.DataManip.Metric.rdMetricMatrixCalc import GetTanimotoSimMat, GetTanimotoDistMat


def mols_to_NHA(x): # input mmpa dataframe, grab smarts/smiles output number of heavy atoms
    return Chem.MolFromSmarts(x).GetNumHeavyAtoms()

def clean_mmpa_pairs_len(mmpa_df):
    temp=pd.DataFrame() # temp dataframe
    if 'LHS' not in mmpa_df.columns: # add LHS and RHS if not present
        mmpa_df = split_transition(mmpa_df, 'smirks')     # produce LHS and RHS
    else:
        temp['common_core_HA'] = mmpa_df['common_core'].apply(mols_to_NHA) # produce number of heavy atoms
        temp['LHS_HA'] = mmpa_df['LHS'].apply(mols_to_NHA)
        temp['RHS_HA'] = mmpa_df['LHS'].apply(mols_to_NHA)
        
        temp['len_check'] = np.where((temp['LHS_HA'] >= temp['common_core_HA']) & (temp['RHS_HA'] >= temp['common_core_HA'])
                     , 'fail', 'pass') # compare lengths of heavy atoms
        
        mmpa_df = mmpa_df.drop(temp[temp['len_check']=='fail'].index) # drop index that failed length check
        
        print('Initial number of transofrms: {} \nNumber fo transforms disqualified based on length discrepancy: {} \nRemaining number of transforms: {}'.format(len(temp[temp['len_check']=='fail']) +  len(mmpa_df) , len(temp[temp['len_check']=='fail']), len(mmpa_df)))
        # return temp to debug
    return mmpa_df

def stat_it_2(it):
    
    hts_pairs = it
    hts_smirk_repearts = hts_pairs.smirks.value_counts().index.tolist()
    hts_smirk_repearts_values = hts_pairs.smirks.value_counts().tolist()

    repeats_filtered = hts_smirk_repearts
    repeats_filtered=[y for x, y in zip(hts_smirk_repearts_values, hts_smirk_repearts) if x >= 3] # at least 3 pairs
    
    res_arr=[]
    
    print('Number of unique transforms: {} \nProcessing transforms:...\n '.format(len(repeats_filtered)))

    for i in range(len(repeats_filtered)):

        pair_eg = hts_pairs[hts_pairs['smirks'] == repeats_filtered[i]].reset_index(drop=True)

        t, p_2 = stats.ttest_rel(pair_eg.measurement_B,pair_eg.measurement_A) #Two-sided p-value / need only one side?

        new_row = {'smirks':hts_smirk_repearts[i], 'dof':len(pair_eg)-1,'t-stat':t, \
                   'p-val (t-test)':p_2, 'measurement_delta':pair_eg.measurement_B.mean() - pair_eg.measurement_A.mean(), 'std':pair_eg.measurement_delta.std(),  'sem':pair_eg.measurement_delta.sem(), }

        res_arr.append(new_row)
        
        if i % 1000 == 0 and i != 0:
            print(i)
    print('done!')
    return pd.DataFrame(res_arr, columns = ['smirks', 'dof' ,'t-stat', 'p-val (t-test)', 'measurement_delta', 'std',  'sem'])


def zero_in(dataframe,cutoff, pos_only=True):
    
    res = dataframe
    
    if pos_only==True:
        res_t_pos = res[res['t-stat']>0].reset_index(drop=True) # positive change
    else:
        res_t_pos=res
        
    res_t_pos_p_pos = res_t_pos[res_t_pos['p-val (t-test)']<cutoff] # significant change
    
    res_pos_t_sorted = res_t_pos_p_pos.sort_values(by='measurement_delta', ascending=False).reset_index(drop=True)
    
    print('Number of unique transforms where p-val < {} is {}'.format(cutoff, len(res_pos_t_sorted)))
    print('Split between {} positive transforms and {} negative transforms'.format(len(res_t_pos_p_pos[res_t_pos_p_pos['t-stat'] > 0]), len(res_t_pos_p_pos[res_t_pos_p_pos['t-stat'] < 0])))
    
    return res_pos_t_sorted