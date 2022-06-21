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





def mark_repeated_nmr(dists):
    #inplace
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
    #inplace
    dists.loc[(dists["ploopk-dist"].isna()),"include"]="NO_PLOOP_R"
    dists.loc[(dists["ploopk-dist"]>5)&(dists["include"].isna()),"include"]="KPLDIST>5"
    dists.loc[(dists["mg_in_site_4a_from_beta"]=="NONE")&(dists["include"].isna()),"include"]="NOMG"
    dists.loc[dists.resolution>5,"include"]="BAD_RES"
    dists.loc[((dists["nh1-alpha-dist"]=="")|(dists["nh2-alpha-dist"]==""))&dists["include"].isna(),"include"]="FAULTY_ARG"
    ###Closest Arg residue lacks one of the terminal NH2s, this will produce errors later
    mark_repeated_nmr(dists)
                    
def assign_uniprot(dists,pdb,chain,row, d_up_mapping,d_up_info, indexr):
        ids_m=d_up_mapping.loc[(d_up_mapping['PDB']==pdb)&(d_up_mapping['CHAIN']==chain)]
        if not pd.isnull(chain):
            ploop_n=row["gly13-id"]+3
        else:
            ploop_n=pd.nan
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
                pfams=upinfo["Cross-reference (Pfam)"].values[0]


                dists.loc[indexr,['Uniprot_Ac','Uniprot_Id','Protein_name_up','Gene_name_up']]=upid,uname,pname,gname
        else:
            print("No Uniprot ID, ", pdb, chain)
        return upid,ploop_n

def get_pfamid(pdb, chain, ploop_n, d_pfam_hmm, pfams_listed_as_pl):
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

def add_identifiers(dists, d_ip_mapping, d_up_info, d_pfam_hmm, d_pfam_mapping, pf_to_suprf):
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
            """pf_alt_chain=d_pfam_hmm[(d_pfam_hmm.PDB_ID==pdb.upper())].CHAIN_ID.to_list()
            if len(pf_alt_chain)>0:
                sim_site=dists[(dists.PDBID==pdb)&(dists["gly13-chain"].isin(pf_alt_chain))&(dists["gly13-id"]==row["gly13-id"])]
                if sim_site.shape[0]>0:
                    chain_similar=sim_site["gly13-chain"].values[0]
                    pfam, pfam_n, pfam_descr, pfc=get_pfamid(pdb, chain_similar, ploop_n, d_pfam_hmm, pfams_listed_as_pl)

                """
            pf_alt_chain=d_pfam_hmm[(d_pfam_hmm.PDB_ID==pdb.upper())].CHAIN_ID.to_list()
            alt_chains=d_pfam_mapping[(d_pfam_mapping.PDB==pdb)&(d_pfam_mapping.SP_PRIMARY==upid)].CHAIN.to_list()
            alt_chains.extend(pf_alt_chain)
            if len(alt_chains)>0:
                sim_site=dists[(dists.PDBID==pdb)&(dists["gly13-chain"].isin(alt_chains))&(dists["gly13-id"]==row["gly13-id"])]
                if sim_site.shape[0]>0:
                    chain_similar=sim_site["gly13-chain"].values[0]
                    pfam, pfam_n, pfam_descr, pfc=get_pfamid(pdb, chain_similar, ploop_n, d_pfam_hmm, pfams_listed_as_pl)
                """display(chain_similar, sim_site, chain)

                alt_chains=set(d_pfam_mapping[(d_pfam_mapping.PDB==pdb)&(d_pfam_mapping.SP_PRIMARY==upid)].CHAIN.to_list())
                chains_same=list((alt_chains-set(chain)).intersection(set(pf_alt_chain)))
                if len(chains_same)!=0:
                    pfam, pfam_n, pfam_descr, pfc=get_pfamid(pdb, chains_same[0], ploop_n, d_pfam_hmm, pfams_listed_as_pl)

                else:
                    if len(alt_chains)>0:
                        pfams=d_pfam_mapping[(d_pfam_mapping.PDB==pdb)&(d_pfam_mapping.SP_PRIMARY==upid)].PFAM_ID.values
                        if len(pfams_listed_as_pl.intersection(pfams))>0:
                            pfam=list(pfams_listed_as_pl.intersection(pfams))[0]
                            pfam_n,pfam_descr=d_pfam_hmm.loc[d_pfam_hmm.PFAM==pfam,["PFAM_Name","PFAM_desc"]].values[0]
                            pfc="assigned without coord check"""



        if pfam=="None":

            pfams=d_pfam_mapping[(d_pfam_mapping.SP_PRIMARY==upid)].PFAM_ID.values
            if len(pfams_listed_as_pl.intersection(pfams))>0:
                pfam=list(pfams_listed_as_pl.intersection(pfams))[0]
                pfam_n,pfam_descr=d_pfam_hmm.loc[d_pfam_hmm.PFAM==pfam,["PFAM_Name","PFAM_desc"]].values[0]
                pfc="assigned by UP id"


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
            if f"{pdb}_{chain}" in all_pdbs_ploop:
                dists.loc[indexr,["superfamily"]]="Unassigned P-loop (PFAM not mapped)"