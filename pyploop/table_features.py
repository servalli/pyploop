#!/usr/bin/env python
__doc__ = """
Module for evaluation of binding sites of P-loop NTPases, see references [XXX] and [YYY] .
This file coontains functions for annotation of a table with calculated distances 
(Protein family assignment, site quality check).
More at https://github.com/servalli/pyploop."""

import re

import pandas as pd
import numpy as np

def mark_repeated_nmr(dists):
    """Mark repeated NMR models .
    In-place modification of the DataFrame
    Args:
        dists : Pandas DataFrame
    """    #inplace
    for model in range(1,max(dists.model)):
        nmrs=dists[~dists["include"].isin(["NO_PLOOP","KPLDIST>5","NOMG","NMR_RED"])&(dists["method"]=="solution nmr")&(dists["model"]==model)]
        #print(model)
        #display(nmrs)
        if nmrs.shape[0]>0:
            for indr, row in nmrs.iterrows():
                p,n,i=row[["PDBID","nuc_chain","nuc_id"]]
                #print(p,n,i)
                models_bef=dists[(dists["model"]<model)&(dists["PDBID"]==p)&(dists["nuc_chain"]==n)&(dists["nuc_id"]==i)]
                if models_bef.shape[0]>0:
                    dists.loc[(dists["PDBID"]==p)&(dists["nuc_chain"]==n)&(dists["nuc_id"]==i)&(~dists["include"].isin(["NO_PLOOP","KPLDIST>5","NOMG"]))&(dists["method"]=="solution nmr")&(dists["model"]>=model),"include"]="NMR_RED"
                    
                    
def mark_lowqual_sites(dists):
    """Mark sites inappropriate for further processing in table with P-loop NTPase 
    binding sites information. In-place modification of the DataFrame
    Args:
        dists : Pandas DataFrame
    """    
    dists.loc[(dists["ploopk-dist"].isna()),"include"]="NO_PLOOP_R" #P-loop residue was not found
    dists.loc[(dists["ploopk-dist"]>5)&(dists["include"].isna()),"include"]="KPLDIST>5" #It was found but is located at 5A or more from beta-phosphate
    dists.loc[(dists["mg_in_site_4a_from_beta"]=="NONE")&(dists["include"].isna()),"include"]="NOMG" #No divalent cation in 4A from beta-phosphate
    dists.loc[dists.resolution>5,"include"]="BAD_RES" #Bad resolution
    dists.loc[((dists["nh1-alpha-dist"]=="")|(dists["nh2-alpha-dist"]==""))&dists["include"].isna(),"include"]="FAULTY_ARG"
    ###Closest Arg residue lacks one of the terminal NH2s, this will produce errors later
    mark_repeated_nmr(dists) #exclude NMR models apart from the first good one
    
                    
def assign_uniprot(dists,pdb,chain,row, d_up_mapping,d_up_info, indexr):
    """Assign Uniprot Ids to each binding site, if available

    Args:
        dists : distance Pandas DataFrame
        pdb (str): PDB ID
        chain (str): PDB chain
        row : current row in dataframe
        d_up_mapping (df): PDB - Uniprot mapping
        d_up_info (df): Uniprot additional information
        indexr: current row index

    Returns:
        uniprot_id, Walker A Lys id
    """    
    ids_m=d_up_mapping.loc[(d_up_mapping['PDB']==pdb)&(d_up_mapping['CHAIN']==chain)]
    if not pd.isnull(chain) and not type(row["gly13-id"]==str):
        ploop_n=row["gly13-id"]+3
    else:
        ploop_n=np.nan
    upid="None"    
    if ids_m.shape[0]>0: 
        if ids_m.shape[0]>1 : #If several ids are mapped to this chain and ploop  region is known, get the right id
            if ploop_n!="":
                if ids_m.loc[(ids_m.RES_BEG<=ploop_n)&(ids_m.RES_END>=ploop_n),'SP_PRIMARY'].shape[0]>0:
                    upid=ids_m.loc[(ids_m.RES_BEG<=ploop_n)&(ids_m.RES_END>=ploop_n),'SP_PRIMARY'].values[0]
        if upid=="None":    
            upid=ids_m.SP_PRIMARY.values[0]
            #print( upid)
  
        
    if upid!="None":
        upinfo=d_up_info.loc[d_up_info['Entry'].values==upid,:]
        if upinfo.shape[0]>0:
            uname=upinfo['Entry name'].values[0]
            p_longname=upinfo['Protein names'].values[0]
            pname=re.split('[\(\)\[\]]',p_longname)[0]
            gname=upinfo['Gene names'].values[0]
            #pfams=upinfo["Cross-reference (Pfam)"].values[0]
            dists.loc[indexr,['Uniprot_Ac','Uniprot_Id','Protein_name_up','Gene_name_up']]=upid,uname,pname,gname
    else:
        print("No Uniprot ID, ", pdb, chain)
    return upid,ploop_n

def get_pfamid(pdb, chain, ploop_n, d_pfam_hmm, pfams_listed_as_pl):
    """Get Pfam id with sequence coordinate check

    Args:
        pdb (str): PDB Id
        chain (str): PDB chain id
        ploop_n (int): Walker A Lys id
        d_pfam_hmm (Pandas DataFrame): Pfam profiles HMMer to PDB
        pfams_listed_as_pl (set): Pfam domains belonging to CL0023

    Returns:
        pfam id, pfam domain name , pfam domain description, comment on how the ID was mapped
    """    
    pf_i=d_pfam_hmm[(d_pfam_hmm.PDB_ID==pdb.upper())&(d_pfam_hmm.CHAIN_ID==chain)]
    pfam=""

    if pf_i.shape[0]>0:
        pfams=set(pf_i.PFAM.to_list()).intersection(pfams_listed_as_pl)
        if not pd.isnull(ploop_n):

            pf_bycoord=pf_i[(pf_i.start<=ploop_n)&(pf_i.end>=ploop_n)&pf_i.PFAM.isin(pfams)]

            #Pfam df was already sorted by value, so then
            if pf_bycoord.shape[0]>0:
                pfam, pfam_n, pfam_descr=pf_bycoord[["PFAM","PFAM_Name","PFAM_desc"]].values[0]
                comment="Ok"
        else:
            
            if len(pfams)>0:
                pfam, pfam_n, pfam_descr=pf_i.loc[pf_i.PFAM.isin(pfams),["PFAM","PFAM_Name","PFAM_desc"]].values[0]               

                comment="coordinate missmatch,picking best P-loop dom"
            
    if pfam=="": 
        pfam, pfam_n, pfam_descr, comment="None","None","None","not found"
            
    return pfam, pfam_n, pfam_descr, comment    

def add_identifiers(dists, d_up_mapping, d_up_info, d_pfam_hmm, d_pfam_mapping, pf_to_suprf,ploop_all_chains,pdbs_ploop, *args):
    """Add Uniprot, Pfam and superfamily Information. 
    Using several mappings to make sure we can annotate as many sites as we can

    Args:
        dists (pandas DataFrame): distance dataframe
        d_up_mapping (pandas DataFrame): Uniprot-PDB mapping
        d_up_info (pandas DataFrame): Uniprot info for each accession
        d_pfam_hmm (pandas DataFrame): HMMer output from Pfam against PDB
        d_pfam_mapping (pandas DataFrame): _description_
        pf_to_suprf (pandas DataFrame): Pfam to P-loop classes
        ploop_all_chains (pandas DataFrame): P-loop list from IPR027417
        pdbs_ploop (pandas DataFrame): Another mapping
    """    
    all_pdbs_ploop=set(ploop_all_chains.ind_ch.to_list())
    pfams_listed_as_pl=set(pf_to_suprf.pf_id.to_list())
    
    for indexr, row in dists.iterrows():

        pdb=row['PDBID']
        chain=row['gly13-chain']
        if pd.isnull(chain):
            chain=row['nuc_chain']
        #print pdb
        #if pdb in all_pdbs_mapped:
        upid,ploop_n=assign_uniprot(dists,pdb,chain,row, d_up_mapping,d_up_info, indexr)
        pfam=""
        pfam, pfam_n, pfam_descr, pfc=get_pfamid(pdb, chain, ploop_n, d_pfam_hmm, pfams_listed_as_pl)


        if pfam=="None":
            #Check if other chains with similar sequence were mapped

            pf_alt_chain=d_pfam_hmm[(d_pfam_hmm.PDB_ID==pdb.upper())].CHAIN_ID.to_list()
            alt_chains=d_pfam_mapping[(d_pfam_mapping.PDB==pdb)&(d_pfam_mapping.SP_PRIMARY==upid)].CHAIN.to_list()
            alt_chains.extend(pf_alt_chain)
            if len(alt_chains)>0:
                sim_site=dists[(dists.PDBID==pdb)&(dists["gly13-chain"].isin(alt_chains))&(dists["gly13-id"]==row["gly13-id"])]
                if sim_site.shape[0]>0:
                    chain_similar=sim_site["gly13-chain"].values[0]
                    pfam, pfam_n, pfam_descr, pfc=get_pfamid(pdb, chain_similar, ploop_n, d_pfam_hmm, pfams_listed_as_pl)
        #If this doesn't help, use mapping by Id     
        if pfam=="None":

            pfams=d_pfam_mapping[(d_pfam_mapping.SP_PRIMARY==upid)].PFAM_ID.values
            if len(pfams_listed_as_pl.intersection(pfams))>0:
                pfam=list(pfams_listed_as_pl.intersection(pfams))[0]
                pfam_n,pfam_descr=d_pfam_hmm.loc[d_pfam_hmm.PFAM==pfam,["PFAM_Name","PFAM_desc"]].values[0]
                pfc="assigned by UP id"
        #If this doesn't help, try to use another mapping by Id 
        if pfam=="None":
            if pdbs_ploop.loc[(pdbs_ploop['pdb_id']==pdb)&(pdbs_ploop['CHAIN_ID']==chain)].shape[0]>0:
                pfam=pdbs_ploop.loc[(pdbs_ploop['pdb_id']==pdb)&(pdbs_ploop['CHAIN_ID']==chain)]["PFAM_ACC"].values[0]
                pfam_n,pfam_descr=d_pfam_hmm.loc[d_pfam_hmm.PFAM==pfam,["PFAM_Name","PFAM_desc"]].values[0]

                pfc="assigned from pdbmap"
        dists.loc[indexr,["pfam_acc",'pfam_domain','domain_name',"pfam_comm"]]=pfam, pfam_n, pfam_descr,pfc        
        if pfam in pfams_listed_as_pl:
            sfam=pf_to_suprf.loc[pf_to_suprf['pf_id']==pfam,"superfamily"].values[0]

            dists.loc[indexr,["superfamily"]]=sfam
        else:
            #Give up
            if f"{pdb}_{chain}" in all_pdbs_ploop:
                dists.loc[indexr,["superfamily"]]="Unassigned P-loop (PFAM not mapped)"   
    return 1              

reduce_cols=[
  #'Unnamed: 0',
  'superfamily',
 'pfam_acc',
 'pfam_domain',

 'Uniprot_Id',
 'Protein_name_up',
 'PDBID',
 'resolution',
 'method',
 "nucleotide_type",
 "nucleotide_id",
  'model',
'mg_in_site_4a_from_beta',
 'AG-site',
    'G-site',
 
 "lys_id",
 'LYS-site',


 'nz-alpha-atom',
 'nz-alpha-dist',
 'nz-gamma-atom',
 'nz-gamma-dist',
 
  'arg_id',
'TYPE_ARG',
 'nh1-alpha-atom',
 'nh1-alpha-dist',
 'nh1-gamma-atom',
 'nh1-gamma-dist',
 'nh2-alpha-atom',
 'nh2-alpha-dist',
 'nh2-gamma-atom',
 'nh2-gamma-dist',
 'ne-alpha-atom',
 'ne-alpha-dist',
 'ne-gamma-atom',
 'ne-gamma-dist',
 
 'asn_ID',
 'asn_TYPE',
 
 'asn_ne2-alpha-atom',
 'asn_ne2-alpha-dist',
 'asn_ne2-gamma-atom',
 'asn_ne2-gamma-dist',
 'Surr_N_BB',
 'Surr_N_SCh',
 'gly13',

 'nuc-to-g13-atom',
 'dist-gly13',
 'lys-ploop-info',
 'ploopk-dist',
 'SerK+1-Mg',
'WB-Asp/Glu',
"WBD-SerK+1_dist",
             
"WBD-Mg",
'is_hydro',
'preceding_res',             
"water_present"            
             
 ]

rename_cols={
  'asn_ne2-alpha-atom':'asn_nd2-alpha-atom',
 'asn_ne2-alpha-dist':'asn_nd2-alpha-dist',
 'asn_ne2-gamma-atom':'asn_nd2-gamma-atom',
 'asn_ne2-gamma-dist':'asn_nd2-gamma-dist',}
                                         
                
def format_for_excel(dists_selected, reduce_cols=reduce_cols, rename_cols=rename_cols):
    """Quick formatting of the dataframe to more easily readable"""

    r=dists_selected.copy()
    r["nucleotide_type"]=r.nuc_type

    r.loc[r.nucleotide_type.isin(["ADP","GDP"]),"nucleotide_type"]=r.loc[r.nucleotide_type.isin(["ADP","GDP"]),"nucleotide_type"]+"-"+r.loc[r.nucleotide_type.isin(["ADP","GDP"]),"nuc_gamma_moiety_type"]

    r["nucleotide_id"]=r.nuc_chain+r.nuc_id.astype(str)


    r["lys_id"]=r.Lys_res_ch+r.Lys_res_id.astype(str)

    r["arg_id"]=r.Arg_res_ch+r.Arg_res_id.astype(str)

    r["gly13"]=r['gly13-type']+" "+r['gly13-chain']+r['gly13-id'].astype(str)

    r=r[reduce_cols]
    r.rename(columns=rename_cols)
    return r

