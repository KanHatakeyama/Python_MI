"""
this is a wrapper class to easily calculate fingerpritns or descriptors by RDKit
20200905 modif lib to show errors with broken molecules
"""

# import rdkit library
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors
import pandas as pd
import numpy as np
from numpy import inf

#make mol object from smiles
def mol_from_smiles(smiles):
    m=Chem.MolFromSmiles(smiles)
    if m is None:
        print("failed to purse: ", smiles )
        print("please recheck SMILES")
    return m

#draw chemical structure from smiles
def draw_SMILES(smiles):
    """
    smiles: SMILES (string)
    return: image
    """
    m=mol_from_smiles(smiles)
    return Draw.MolToImage(m) 


#calculate fingerprnt

#default FP function
def default_FP_Func(m):
    return AllChem.GetMorganFingerprintAsBitVect(m,2,nBits=512) 

#fingerprint class
class Fingerprint:
    def __init__(self,fp_func=default_FP_Func, str_mode=False):
        """
        fp_func: fingerprint function. you can set your favorite one
        str_mode: fingerprints will be retured as string if True
        """
        self.fp_func=fp_func
        self.str_mode=str_mode

    #calculate FP
    def calc(self,smiles):
        """
        smiles: SMILES
        return: fingerprint
        """
        try:
            m=mol_from_smiles(smiles)
            
            if m is None:
                print("invalid smiles!",smiles)
                return -1
            fp= self.fp_func(m)
            fp=fp.ToBitString()
        
        except:
            print("invalid smiles!", smiles)
            return -1
        

        if self.str_mode == False:
            fp=[int(i) for i in fp]
        return fp
    
    #calculate FP (for list-type input)
    def calc_list(self,ls,pandas_mode=True):
        """
        ls: list of smiles
        pandas_mode: if true, return dataframe otherwise return list 
        """
        res_list=[self.calc(i) for i in ls]
        
        if pandas_mode:
            df=pd.DataFrame(res_list)
            return df
        return [self.calc(i) for i in d_list]
        
        


#descriptor class
class RDKitDescriptors:
    def __init__(self):
        self.desc_list = [desc_name[0] for desc_name in Descriptors.descList]
        self.calculator = MoleculeDescriptors.MolecularDescriptorCalculator(self.desc_list)
        
    #calc descriptors    
    def calc(self,smiles,dict_mode=True):
        """
        smiles: smiles
        dict_mode: if true, return type will be a dict, otherwise, list
        """
        try:
            m=mol_from_smiles(smiles)
            
            if m is None:
                print("invalid smiles!",smiles)
                return -1
            
            descs=self.calculator.CalcDescriptors(m)
            descs=np.nan_to_num(descs)
            descs[descs > 10**5] = 0
            descs[descs < -10**5] = 0            
            
        except:
            print("invalid smiles!",smiles)
            return -1
        
        if dict_mode:
            desc_dict={k:v for k,v in zip(self.desc_list,descs)}
            
            return desc_dict
        else:
            return descs

    #calc descriptors from smiles list
    def calc_list(self,ls,pandas_mode=True):
        """
        ls: list of smiles
        pandas_mode: if true, return type will be dataframe, otherwise list
        """
        res_list=[self.calc(i,dict_mode=False) for i in ls]
        
        if pandas_mode:
            df=pd.DataFrame(res_list)
            df.columns=self.desc_list
            return df
        
        return res_list