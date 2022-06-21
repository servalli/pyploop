"""
Module for evaluation of binding sites of P-loop NTPases, see references [XXX] and [YYY] .
This file coontains functions for PDB processing"""


import os
import operator
import pandas as pd
from Bio.PDB import PDBList
import pandas as pd
import time
import re
import pandas as pd
from itertools import product
import os
import numpy as np
import Bio.PDB as PDB
from Bio.PDB import PDBList, MMCIFParser, Select, Selection, PDBIO, NeighborSearch
from Bio.PDB.MMCIF2Dict import MMCIF2Dict


def assign_arg_type(dists, d_sites_checked=""):
    distance_columns=[col for col in dists.columns if "dist" in col and col!="WBD-SerK+1_dist"]
    for col in distance_columns:
        dists[col]=dists[col].astype(float)
    #notype=(dists["TYPE_ARG"].isna())&(dists["include"].isna())

    dists.loc[(((dists["nh2-gamma-dist"]<4)&(dists["nh1-alpha-dist"]<4))|((dists["nh1-gamma-dist"]<4)&(dists["nh2-alpha-dist"]<4))),"TYPE_ARG"]="FORK."
    dists.loc[(dists["nh1-gamma-dist"]<4)&(dists["nh1-alpha-dist"]<4),"TYPE_ARG"]="NH1 weak."
    dists.loc[(dists["nh1-gamma-dist"]<3.2)&(dists["nh1-alpha-dist"]<3.2),"TYPE_ARG"]="NH1."
    dists.loc[(dists["nh2-gamma-dist"]<4)&(dists["nh2-alpha-dist"]<4),"TYPE_ARG"]="NH2 weak."
    dists.loc[(dists["nh2-gamma-dist"]<3.2)&(dists["nh2-alpha-dist"]<3.2),"TYPE_ARG"]="NH2."
    dists.loc[(((dists["nh2-gamma-dist"]<4)|(dists["nh1-gamma-dist"]<4))&((dists["nh1-alpha-dist"]>4)&(dists["nh2-alpha-dist"]>4))),"TYPE_ARG"]="ONLY GAMMA weak."
    dists.loc[(((dists["nh2-gamma-dist"]<3.2)|(dists["nh1-gamma-dist"]<3.2))&((dists["nh1-alpha-dist"]>4)&(dists["nh2-alpha-dist"]>4))),"TYPE_ARG"]="ONLY GAMMA."
    dists.loc[(((dists["nh2-gamma-dist"]>4)&(dists["nh1-gamma-dist"]>4))),"TYPE_ARG"]="NONE."
    if type(d_sites_checked)!=str:
        for ind,row in d_sites_checked.iterrows():
            dists.loc[(dists.PDBID==row.PDB)&((dists.nuc_chain+dists.nuc_id.astype(str))==row.nuc),"TYPE_ARG"]=row.argtype




def assign_mono_type(dists):
    dists.loc[(dists["nz-gamma-dist"]<4),"LYS-site"]="G"
    dists.loc[(dists["nz-gamma-dist"]<4)&(dists["nz-alpha-dist"]<4),"LYS-site"]="AG"

    dists.loc[(dists["LYS-site"].isna()),"LYS-site"]="NONE"
    
    
    dists.loc[(dists["asn_ne2-gamma-dist"]<4),"asn_TYPE"]="(only G)"
    dists.loc[(dists["asn_ne2-gamma-dist"]<4)&(dists["asn_ne2-alpha-dist"]<4),"asn_TYPE"]="AG_weak"
    dists.loc[(dists["asn_ne2-gamma-dist"]<3.2)&(dists["asn_ne2-alpha-dist"]<3.2),"asn_TYPE"]="AG"
    dists.loc[(dists["asn_ne2-gamma-dist"]>4),"asn_TYPE"]="NONE"
    dists.loc[(dists["asn_ne2-gamma-dist"].isna()),"asn_TYPE"]="NO_ASN"


def assign_type(dists):
    

    
    dists.loc[dists["LYS-site"]=="AG","AG-site"]="LYS"
    dists.loc[dists.TYPE_ARG.str.contains("NH"),"AG-site"]="ARG"
    dists.loc[dists.TYPE_ARG.str.contains("FORK"),"AG-site"]="ARG"
    dists.loc[dists.asn_TYPE.str.contains("AG"),"AG-site"]="ASN"
    dists.loc[dists["AG-site"].isna(),"AG-site"]="NONE"
    
    #GAMMA
    g_main_k_g=(dists["LYS-site"]=="G")
    
    g_main_r_g=dists.TYPE_ARG.str.contains("GAMMA")
    

    
    n_lys_g=dists["Surr_N_SCh"].str.count("LYS")+g_main_k_g.astype(int)   
    
    n_arg_g=dists["Surr_N_SCh"].str.count("ARG")+g_main_r_g.astype(int)
    
    comma=(np.sign(n_arg_g)*np.sign(n_lys_g))*(n_lys_g.astype(str).str.slice(0,0)+", ")
    y=(n_arg_g.astype(str)+"ARG")*np.sign(n_arg_g)+comma+(n_lys_g.astype(str)+"LYS")*np.sign(n_lys_g)
    
    
    dists.loc[:,"G-site"]=y.copy()
    dists.loc[dists["G-site"]=="","G-site"]="NONE"



