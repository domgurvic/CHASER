import pandas as pd
import numpy as np
from scipy import stats
import re
import os
import matplotlib.pyplot as plt

import seaborn as sns

#chem

from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, MACCSkeys, Descriptors, Descriptors3D, Draw, rdMolDescriptors, Draw, PandasTools, rdFingerprintGenerator
from rdkit.DataManip.Metric.rdMetricMatrixCalc import GetTanimotoSimMat, GetTanimotoDistMat

def split_transition(df, col):
    df['LHS'] = [re.split('>>',df[col].values[i])[0] for i in range(len(df)) ]
    df['RHS'] = [re.split('>>',df[col].values[i])[1] for i in range(len(df)) ]
    return df

def new_smarts():
#     print(os.getcwd())
    func_groups=pd.read_csv('smarts_fg.csv')
    
        #fetch all substructure definitions and calculate mosl for them
    print('Generating molecular objects from pre-defined substructures')
    mol_substructures=[]
    for substructure in func_groups.SMARTS:
        mol_substructures.append(Chem.MolFromSmarts(substructure))

    return mol_substructures,  func_groups.name.to_list()

def calculate_fractions_mk4(df):

#     mol_substructures, name_substructure = generate_mols_for_substructures()

#     name_substructure = name_substructure + ['smirks', 'target']
    
    mol_substructures, name_substructure = new_smarts()

    name_substructure = name_substructure + ['smirks', 'measurement' ,'target']


    frame_left=[]
    frame_right=[]

    print('Calcualting LHS+RHS matches')

    #for index in enumerate(df.LHS.values)):

    for index in range(len(df)):  

        #grab structure
        frame_temp_left=pd.DataFrame(0, index=range(1), columns=name_substructure)
        frame_temp_right=pd.DataFrame(0, index=range(1), columns=name_substructure)

        #turn it into mol 
        try:
            mol_target_left=Chem.MolFromSmarts(df.LHS.values[index])
            mol_target_left.UpdatePropertyCache()
            mol_target_left = Chem.AddHs(mol_target_left)
            Chem.SanitizeMol(mol_target_left)
        except TypeError:
            print('Error: ', index, target)

        try:
            mol_target_right=Chem.MolFromSmarts(df.RHS.values[index])
            mol_target_right.UpdatePropertyCache()
            mol_target_right = Chem.AddHs(mol_target_right)
            Chem.SanitizeMol(mol_target_right)
        except TypeError:
            print('Error: ', index, target)    

        if type(mol_target_right) != Chem.rdchem.Mol or type(mol_target_left) != Chem.rdchem.Mol:

            print('failed to MolObject: ', index)


        frame_temp_left['smirks'] = df.smirks.values[index]
        frame_temp_left['target'] = df.measurement_delta.values[index]    

        for sub_nr, sub in enumerate(mol_substructures):
            if mol_target_left.HasSubstructMatch(sub):
                frame_temp_left[name_substructure[sub_nr]] = [1]

        frame_temp_right['smirks'] = df.smirks.values[index]
        frame_temp_right['target'] = df.measurement_delta.values[index]    

        for sub_nr, sub in enumerate(mol_substructures):
            if mol_target_right.HasSubstructMatch(sub):
                frame_temp_right[name_substructure[sub_nr]] = [1]

        frame_left.append(frame_temp_left.values)
        frame_right.append(frame_temp_right.values)

    frame_left_df = pd.DataFrame(np.concatenate(frame_left), columns = name_substructure)
    # compare right hand side
    frame_right_df = pd.DataFrame(np.concatenate(frame_right), columns = name_substructure)

    diff = frame_right_df.iloc[:,:-3] - frame_left_df.iloc[:,:-3] 

    diff['smirks'] = frame_right_df['smirks']
    diff['target'] = frame_right_df['target']

    return diff.reset_index(drop=True), frame_left_df.reset_index(drop=True), frame_right_df.reset_index(drop=True)

def find_sig_feats_mk2(l_feats, r_feats, p_val):
    
    df_sig_feats=pd.DataFrame(columns=l_feats.columns)
    
    #df_feats_fr_only=df_feats.iloc[:,114:199]

    for name in df_sig_feats.columns[:-3]:
        t, p = stats.ttest_rel(r_feats[name],l_feats[name])
        if p < p_val:
            df_sig_feats[name]=pd.Series([p, t])

    df_sig_feats = df_sig_feats.dropna(axis=1)
    
    print('Found significant fractions: ', len(df_sig_feats.columns))

    return df_sig_feats

def results_arr(features_all, fr_sig_descriptors, r_feats, l_feats , fractions_to_drop):
    
    arr=[]

    #features_all_fr_only = features_all.iloc[:,114:199]
    #features_all_fr_only_target = features_all.iloc[:,114:200]


    #fr_descriptors=[]
    #for descriptor in features_all.columns:
    #    if re.search('fr_(.*?)', descriptor) is not None:
    #        fr_descriptors.append(descriptor)
    #

    #features_all_fr_only=features_all[fr_descriptors]
    
    features_all_fr_only = features_all.iloc[:,:-4]

    for name in fr_sig_descriptors: 
        if r_feats[name].mean() - l_feats[name].mean() == 0:
            print('{} has neutral correlation '.format(name))


        elif r_feats[name].mean() - l_feats[name].mean() < 0:

            print('{} has negative correlation '.format(name))


            name_loss= features_all_fr_only[features_all_fr_only[name]<0] # find all negative occurances of this named fraction

            name_loss = name_loss.drop(fractions_to_drop, axis = 1)

            # in all fractions that name is lost who are biggest gains:
            big_gain = pd.DataFrame(name_loss.mean().sort_values(ascending=False)).T.columns

            delta_p_p = round(features_all.iloc[name_loss.index].target.mean(), 2)
            
            delta_p_p_sem = round(features_all.iloc[name_loss.index].target.sem(), 2)
            
            delta_p_p_std = round(features_all.iloc[name_loss.index].target.std(), 2)
            
            #percentage_loss=[round(len(name_loss[name_loss[x]>0])/ len(name_loss) *100, 2 ) for x in big_gain[:10]]

            percentage_loss = [round(name_loss[x].sum() / len(name_loss) *100, 2 ) for x in big_gain[:10]]



            # Look at overlaping percentages

            if percentage_loss[1] == percentage_loss[2] == percentage_loss[3]:
                    big_gain=[big_gain[0], (big_gain[1] , big_gain[2], big_gain[3]), big_gain[4] ]
                    percentage_loss = [percentage_loss[0] , percentage_loss[2], percentage_loss[4]]

            elif percentage_loss[0] == percentage_loss[1] == percentage_loss[2]:
                big_gain=[( big_gain[0], big_gain[1] , big_gain[2]) , big_gain[3], big_gain[4]]
                percentage_loss = [percentage_loss[0] , percentage_loss[3], percentage_loss[4]]


                print('all gain')
                print(big_gain)
                print(percentage_loss)

            elif int(percentage_loss[0]) == int(percentage_loss[1]):
                big_gain=[(big_gain[0], big_gain[1]), big_gain[2], big_gain[3]]
                percentage_loss = [percentage_loss[0] , percentage_loss[2], percentage_loss[3]]

                print('first_gain')
                print(big_gain)
                print(percentage_loss)

            elif int(percentage_loss[1]) == int(percentage_loss[2]):
                big_gain=[big_gain[0], (big_gain[1] , big_gain[2]) , big_gain[3]]
                percentage_loss = [percentage_loss[0] , percentage_loss[2], percentage_loss[3]]


                print('second gain')
                print(big_gain)
                print(percentage_loss)



            if np.array(percentage_loss[:3]).sum() > 100:
                print('percentage_loss 100')


            arr.append([name, 'Negative', -delta_p_p, delta_p_p_sem, delta_p_p_std,  len(name_loss) , big_gain[0], percentage_loss[0], big_gain[1], percentage_loss[1], big_gain[2], percentage_loss[2]])




        elif r_feats[name].mean() - l_feats[name].mean() > 0:
            print('{} has positive correlation '.format(name))    

            name_gain = features_all_fr_only[features_all_fr_only[name]>0]


            name_gain = name_gain.drop(fractions_to_drop, axis = 1)
            

            big_loss = pd.DataFrame(name_gain.mean().sort_values(ascending=True)).T.columns

            delta_p_g = round(features_all.iloc[name_gain.index].target.mean(), 2)
            
            delta_p_g_sem = round(features_all.iloc[name_gain.index].target.sem(), 2)
            
            delta_p_g_std = round(features_all.iloc[name_gain.index].target.std(), 2)

            #percentage_gain=[round(len(name_gain[name_gain[x]<0])/ len(name_gain) *100, 2 ) for x in big_loss[:10]]

            percentage_gain = [round(name_gain[x].sum() / len(name_gain) *100, 2 ) for x in big_loss[:10]]


            # Look at overlaping percentages


            if percentage_gain[1] == percentage_gain[2] == percentage_gain[3]:
                    big_loss=[big_loss[0], (big_loss[1] , big_loss[2], big_loss[3]), big_loss[4] ]
                    percentage_gain = [percentage_gain[0] , percentage_gain[2], percentage_gain[4]]


                    print('1/2/3 loss')
                    print(big_loss)
                    print(percentage_gain)


            elif percentage_gain[0] == percentage_gain[1] == percentage_gain[2]:
                big_loss=[( big_loss[0], big_loss[1] , big_loss[2]) , big_loss[3], big_loss[4]]
                percentage_gain = [percentage_gain[0] , percentage_gain[3], percentage_gain[4]]

                print('0/1/2 loss')
                print(big_loss)
                print(percentage_gain)

            elif int(percentage_gain[0]) == int(percentage_gain[1]):
                big_loss=[(big_loss[0], big_loss[1]), big_loss[2], big_loss[3]]
                percentage_gain = [percentage_gain[0] , percentage_gain[2], percentage_gain[3]]

                print('first double loss')
                print(big_loss)
                print(percentage_gain)

            elif int(percentage_gain[1]) == int(percentage_gain[2]):
                big_loss=[big_loss[0], (big_loss[1] , big_loss[2]) , big_loss[3]]
                percentage_gain = [percentage_gain[0] , percentage_gain[2], percentage_gain[3]]

                print('second double loss')
                print(big_loss)
                print(percentage_gain)




            if np.array(percentage_gain[:3]).sum() < -100:
                print('percentage gain under -100')

            arr.append([name, 'Positive', delta_p_g,  delta_p_g_sem,delta_p_g_std , len(name_gain) , big_loss[0], percentage_gain[0], big_loss[1], percentage_gain[1], big_loss[2], percentage_gain[2]])



    res_neg = pd.DataFrame(arr, columns=[ 'Main fraction', 'Correlation', r'$\overline{\Delta P}$', 'sem','std', 'dof' ,'Opposite fraction 1', '% of opposite 1' ,'Opposite fraction 2', '% of opposite 2','Opposite fraction 3', '% of opposite 3'])

    res_neg = res_neg.sort_values(by='dof', ascending=False)

    res_neg = res_neg.iloc[res_neg[r'$\overline{\Delta P}$'].abs().argsort()[::-1]]


    return res_neg
