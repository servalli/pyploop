"""
Module for evaluation of binding sites of P-loop NTPases, see references [XXX] and [YYY] .
This file coontains functions for type assignment for finger (stimulator) residues"""

import numpy as np

def assign_arg_type(dists, d_sites_checked="",strong=3.2,weak=4):
    """Assign interaction type to main Arg (closest to beta-phosphate) in each site

    Args:
        dists (Pandas DataFrame): distances dataframe
        d_sites_checked (optional): Table with binding sites previously checked.
        strong (float): strong H-bond distance treashold. Defaults to 3.2.
        weak (float): Weak H-bond distance treashold. Defaults to 4.
    """    
    distance_columns=[col for col in dists.columns if "dist" in col and col!="WBD-SerK+1_dist"]
    for col in distance_columns:
        dists[col]=dists[col].astype(float)
    

    dists.loc[(((dists["nh2-gamma-dist"]<weak)&(dists["nh1-alpha-dist"]<weak))|((dists["nh1-gamma-dist"]<weak)&(dists["nh2-alpha-dist"]<weak))),"TYPE_ARG"]="FORK."
    dists.loc[(dists["nh1-gamma-dist"]<weak)&(dists["nh1-alpha-dist"]<weak),"TYPE_ARG"]="NH1 weak."
    dists.loc[(dists["nh1-gamma-dist"]<strong)&(dists["nh1-alpha-dist"]<strong),"TYPE_ARG"]="NH1."
    dists.loc[(dists["nh2-gamma-dist"]<weak)&(dists["nh2-alpha-dist"]<weak),"TYPE_ARG"]="NH2 weak."
    dists.loc[(dists["nh2-gamma-dist"]<strong)&(dists["nh2-alpha-dist"]<strong),"TYPE_ARG"]="NH2."
    dists.loc[(((dists["nh2-gamma-dist"]<weak)|(dists["nh1-gamma-dist"]<weak))&((dists["nh1-alpha-dist"]>weak)&(dists["nh2-alpha-dist"]>weak))),"TYPE_ARG"]="ONLY GAMMA weak."
    dists.loc[(((dists["nh2-gamma-dist"]<strong)|(dists["nh1-gamma-dist"]<strong))&((dists["nh1-alpha-dist"]>weak)&(dists["nh2-alpha-dist"]>weak))),"TYPE_ARG"]="ONLY GAMMA."
    dists.loc[(((dists["nh2-gamma-dist"]>weak)&(dists["nh1-gamma-dist"]>weak))),"TYPE_ARG"]="NONE."
    if type(d_sites_checked)!=str:
        for ind,row in d_sites_checked.iterrows():
            dists.loc[(dists.PDBID==row.PDB)&((dists.nuc_chain+dists.nuc_id.astype(str))==row.nuc),"TYPE_ARG"]=row.argtype




def assign_mono_type(dists, strong=3.2,weak=4):
    """Assign interaction type to main Lys and Asn (closest to beta-phosphate) in each site

    Args:
        dists (Pandas DataFrame): distances dataframe
        strong (float): strong H-bond distance treashold. Defaults to 3.2.
        weak (float): Weak H-bond distance treashold. Defaults to 4.
    """
    dists.loc[(dists["nz-gamma-dist"]<weak),"LYS-site"]="G"
    dists.loc[(dists["nz-gamma-dist"]<weak)&(dists["nz-alpha-dist"]<weak),"LYS-site"]="AG"

    dists.loc[(dists["LYS-site"].isna()),"LYS-site"]="NONE"
    
    
    dists.loc[(dists["asn_ne2-gamma-dist"]<weak),"asn_TYPE"]="(only G)"
    dists.loc[(dists["asn_ne2-gamma-dist"]<weak)&(dists["asn_ne2-alpha-dist"]<weak),"asn_TYPE"]="AG_weak"
    dists.loc[(dists["asn_ne2-gamma-dist"]<strong)&(dists["asn_ne2-alpha-dist"]<strong),"asn_TYPE"]="AG"
    dists.loc[(dists["asn_ne2-gamma-dist"]>weak),"asn_TYPE"]="NONE"
    dists.loc[(dists["asn_ne2-gamma-dist"].isna()),"asn_TYPE"]="NO_ASN"


def assign_type(dists):   
    """Assign general type of stimulator interaction based on fields for main Arg, Lys, Asn
    and additional interacting positively charged residues

    Args:
        dists (Pandas DataFrame): distance dataframe
    """
    
    dists.loc[dists["LYS-site"]=="AG","AG-site"]="LYS"
    dists.loc[dists.TYPE_ARG.str.contains("NH"),"AG-site"]="ARG"
    dists.loc[dists.TYPE_ARG.str.contains("FORK"),"AG-site"]="ARG"
    dists.loc[dists.asn_TYPE.str.contains("AG"),"AG-site"]="ASN"
    dists.loc[dists["AG-site"].isna(),"AG-site"]="NONE"
    
    #GAMMA
    g_main_k_g=(dists["LYS-site"]=="G")
    
    g_main_r_g=dists.TYPE_ARG.str.contains("GAMMA")
    

    n_arg_g=dists["Surr_N_SCh"].str.count("ARG")+g_main_r_g.astype(int)
    
    n_lys_g=dists["Surr_N_SCh"].str.count("LYS")+g_main_k_g.astype(int)      
     
    comma=(np.sign(n_arg_g)*np.sign(n_lys_g))*(n_lys_g.astype(str).str.slice(0,0)+", ")
    y=(n_arg_g.astype(str)+"ARG")*np.sign(n_arg_g)+comma+(n_lys_g.astype(str)+"LYS")*np.sign(n_lys_g)
    
    
    dists.loc[:,"G-site"]=y.copy()
    dists.loc[dists["G-site"]=="","G-site"]="NONE"



