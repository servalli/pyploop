#!/usr/bin/env python
__doc__ ="""
Module for evaluation of binding sites of P-loop NTPases, see references [XXX] and [YYY] .
This file coontains functions for PDB processing.
More at https://github.com/servalli/pyploop."""



import os
import pandas as pd
from Bio.PDB import PDBList
import pandas as pd
import time
import pandas as pd
import os
from Bio.PDB import PDBList, MMCIFParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

from pyploop.process_nucleotide import process_nucleotides, process_one_nucleotide



def download_pdbs(compound_pdb, ploop_wdir, pdb_dir, log_dir):
    
    
    """
    Download PDBs for P-loop NTPases containing each type of NTP-analog 
    and save those in separate folders inside the working directory.

    Args:
        compound_pdb (dict): Dictionary with structures assigned to compounds, ex: {'ATP': ['1A82']}
        ploop_wdir (str): Work directory location
        pdb_dir (str): Subdirectory where to save PDBs
        log_dir (str): Subdirectory where to save log

    Returns:
        Pandas dataframe with structures downloaded
    """    
    pdbl = PDBList(verbose=False)

    log_p=os.path.join(ploop_wdir,log_dir,"PDB_retrieval.txt")
    curr_time=time.strftime("%Y-%m-%d %H:%M")
    with open(log_p,"w") as log_f:
        log_f.write(f"# 3D structures (mmCif) of P-loop NTPases with NTPs/analogs downloaded, job started at {curr_time} \n")
    colnames=["PDBID","nucl"]
    downloaded=pd.DataFrame(columns=colnames)

    for c_id in compound_pdb.keys():
        success=[]
        save_dir=os.path.join(ploop_wdir,pdb_dir, f"PDB{c_id}","full")
        if not os.path.isdir(save_dir): os.makedirs(save_dir) 
        for i in compound_pdb[c_id]:

            fname=pdbl.retrieve_pdb_file(i,pdir=save_dir, file_format="mmCif")
            fmt=os.path.join(save_dir,f"{i.lower()}.cif")
            if fname==fmt:
                success.append(i)
            else:
                print(fname,i)
        downloaded = pd.concat([downloaded, pd.DataFrame([(i,c_id) for i in success],columns=colnames)])
    with open(log_p,"a") as log_f:
        downloaded.to_csv(log_f,sep="\t",index=False)
    return downloaded




def get_terminals(s,atom,rtype):
    """Retrieve terminals of residues of a specific type

    Args:
        s (_type_): _description_
        atom (str): atom name, i.e. "NZ" for Lys
        rtype (_type_): _description_

    Returns:
        _type_: _description_
    """    
    
    all_term_Ns=[re[atom] for re in s.get_residues() if re.resname==rtype and re.has_id(atom) ]
    return all_term_Ns

def get_structinfo(path):
    mmcif_dict = MMCIF2Dict(path)
    
    method=mmcif_dict["_exptl.method"][0].lower()
    pdb_date=mmcif_dict["_pdbx_database_status.recvd_initial_deposition_date"][0] 
    
    if method!="electron microscopy":
        if "_reflns.d_resolution_high" in mmcif_dict:
            resolution=mmcif_dict["_reflns.d_resolution_high"][0]
            if resolution=="?" and "_refine.ls_d_res_high" in mmcif_dict:
                resolution=mmcif_dict["_refine.ls_d_res_high"][0]                
        elif "_refine.ls_d_res_high" in mmcif_dict:
            resolution=mmcif_dict["_refine.ls_d_res_high"][0]
            
            
        else:
            
            resolution=pd.NA
            
       
        
    else:
        resolution=mmcif_dict["_em_3d_reconstruction.resolution"][0]
    if type(resolution)==str:
        if resolution=="?":
            resolution=pd.NA
        else:
            resolution=float(resolution)
    
    p_name=mmcif_dict["_struct.title"][0] 
    
    has_water='HOH' in mmcif_dict["_pdbx_entity_nonpoly.comp_id"]
    
    return method, resolution, p_name, pdb_date, has_water

def process_pdb_dir_classic(p_dir,comps, comps_F, alpha, beta, gamma_F, gamma_N, gamma_types,ploop_all_chains):
    """_summary_

    Args:
        p_dir (str): Folder with PDBs
        comps (list): compound types list
        comps_F (_type_): compound types for which to look for gamma-phosphate mimic
        alpha (list): _description_
        beta (list): _description_
        gamma_F (list): _description_
        gamma_N (list): _description_
        gamma_types (list): gamma mimic compounds
        ploop_all_chains (DataFrame): Interpro mapping (IPR027417)

    Returns:
        _type_: _description_
    """    
    gia_pdb={comp:(sorted([w for w in os.listdir(os.path.join(p_dir,"PDB"+comp,"full")) if w.endswith(".cif")]) if os.path.isdir(os.path.join(p_dir,"PDB"+comp,"full")) else []) for comp in comps}       
    #Dict with structure filenames
    
    
    asn_cols=[
        "asn_ID","asn_TYPE",
        "asn_ne2-alpha-atom","asn_ne2-alpha-dist","asn_ne2-gamma-atom","asn_ne2-gamma-dist"]

    cols=[
        #protein description
          'superfamily','pfam_acc','pfam_domain','domain_name',
          'Uniprot_Ac','Uniprot_Id',
          'Protein_name_up','Gene_name_up',
        #PDB info
          "PDBID", "model",
        #Nucleotide info
          "nuc_type","nuc_chain","nuc_id",
          "nuc_gamma_moiety_type","nuc_gamma_moiety_chain","nuc_gamma_moiety_id",
        #info on site
          "mg_in_site_4a_from_beta",
          "AG-site","G-site",
          "TYPE_ARG","LYS-site","asn_TYPE",
          
        #Closest non-P-loop Lys
          "Lys_res_ch","Lys_res_id",
          "nz-alpha-atom","nz-alpha-dist","nz-gamma-atom","nz-gamma-dist",
        #Closest Arg 
          
          "Arg_res_ch","Arg_res_id",
          "nh1-alpha-atom","nh1-alpha-dist","nh1-gamma-atom","nh1-gamma-dist",
          "nh2-alpha-atom","nh2-alpha-dist","nh2-gamma-atom","nh2-gamma-dist",
          "ne-alpha-atom","ne-alpha-dist","ne-gamma-atom","ne-gamma-dist",
        #Closest Asn  
        "asn_ID",
        "asn_ne2-alpha-atom","asn_ne2-alpha-dist","asn_ne2-gamma-atom","asn_ne2-gamma-dist",
               
        #Other interactions of gamma-phosphate  with surrounding residues
          
          "Surr_N_BB","Surr_N_SCh",
        #P-loop K-3 backbone to gamma
          "gly13-chain","gly13-type","gly13-id","nuc-to-g13-atom", "dist-gly13",
        #P-loop residue availability
            "lys-ploop-info","ploopk-dist",
        #PDB info
        'resolution',"method","pdbname","pdbdate","include",
        #WB-WA interaction and Mg-coordination
        'SerK+1-Mg','WB-Asp/Glu',"WBD-SerK+1_dist","WBD-Mg",
        
        "D+3_CA..K-3_N","D+3_N..K-3_N", "D+4_CA..K-3_N", "D+4_N..K-3_N","rASA_WBD,%","SASA_WBD",
        "water_present"]
    
    
    dists=pd.DataFrame(columns=cols,index=[0])
    
    df_ind=-1
    ik=0
    epa=MMCIFParser(QUIET=True)
    for ind,comp in enumerate(comps):
        print( comp, len(gia_pdb[comp]), "structures of P-loop NTPases available")
        for T,pdb in enumerate(gia_pdb[comp]):
            if T%50==0:
                print( T, pdb[:4],comp)
            pid=pdb[:4]
            path=os.path.join(p_dir,"PDB"+comp,"full", pdb)
            #print(pdb)
            stru=epa.get_structure(pdb, path)
            ###description

            method, resolution, p_name, pdb_date, has_water=get_structinfo(path)
            if type(resolution)==float:
                if resolution>5:   #Skip low-res structures alltogether
                    continue
            ###
            
            ###terminals of relevant atoms
            for s in stru.get_models():
                all_term_ND2s_ASN=get_terminals(s,atom="ND2",rtype="ASN")
                all_term_NZs_LYS=get_terminals(s,atom="NZ",rtype="LYS")

                ###ARGS
                all_term_NH1s_ARG=get_terminals(s,atom="NH1",rtype="ARG")
                
                all_term_NH2s_ARG=get_terminals(s,atom="NH2",rtype="ARG")
                
                #negative
                all_term_acid=get_terminals(s,atom="OE1",rtype="GLU")+get_terminals(s,atom="OE2",rtype="GLU")+get_terminals(s,atom="OD1",rtype="ASP")+get_terminals(s,atom="OD2",rtype="ASP")

                ions=[atm for atm in s.get_atoms() if atm.element in ["MG","MN","CA","SR"]]
                
                
                res_all=[res for res in s.get_residues()]
                for ch in s.get_chains():
                    
                    #### Only use chains listed as P-loop NTPases in INTERPRO
                    #### Nucleotides erroneously assigned to activator and not P-loop chain will thus be skipped
                    if pdb[:4]+"_"+ch.get_id() in ploop_all_chains["ind_ch"].values:
                        ##List nucleotides
                        nuc=[res for res in ch.get_residues() if res.get_resname()==comp]
                        ####CHECK FOR INCOMPLETE NUCLEOTIDES AND SKIP THEM
                        nuc_g=[nc for nc in nuc if "O1B" in nc]
                        bad=set(nuc)-set(nuc_g)
                        if len(bad)>0:
                            print("Warning! Truncated nucleotides (no beta-phosphate)", pdb[:4], ch.get_id(), [nu.get_resname()+str(nu.get_id()[1]) for nu in bad])
                        nuc=nuc_g
                        if len(nuc)>0:
                            df_ind=process_nucleotides(dists,
                        s,nuc, comp,df_ind,
                        resolution,method,p_name,pdb_date,has_water, pdb,
                        alpha,beta,gamma_F,gamma_N,comps_F,
                       all_term_NZs_LYS, all_term_NH1s_ARG, all_term_NH2s_ARG, all_term_ND2s_ASN, all_term_acid, ions,
                       asn_cols)
                       
                       


    return dists




def process_pdb_dir(p_dir,comps, comps_F, alpha, beta, gamma_F, gamma_N, gamma_types,ploop_all_chains):
    """_summary_

    Args:
        p_dir (str): Folder with PDBs
        comps (list): compound types list
        comps_F (_type_): compound types for which to look for gamma-phosphate mimic
        alpha (list): _description_
        beta (list): _description_
        gamma_F (list): _description_
        gamma_N (list): _description_
        gamma_types (list): gamma mimic compounds
        ploop_all_chains (DataFrame): Interpro mapping (IPR027417)

    Returns:
        _type_: _description_
    """    
    gia_pdb={comp:(sorted([w for w in os.listdir(os.path.join(p_dir,"PDB"+comp,"full")) if w.endswith(".cif")]) if os.path.isdir(os.path.join(p_dir,"PDB"+comp,"full")) else []) for comp in comps}       
    #Dict with structure filenames
    
    
    asn_cols=[
        "asn_ID","asn_TYPE",
        "asn_ne2-alpha-atom","asn_ne2-alpha-dist","asn_ne2-gamma-atom","asn_ne2-gamma-dist"]

    cols=[
        #protein description
          'superfamily','pfam_acc','pfam_domain','domain_name',
          'Uniprot_Ac','Uniprot_Id',
          'Protein_name_up','Gene_name_up',
        #PDB info
          "PDBID", "model",
        #Nucleotide info
          "nuc_type","nuc_chain","nuc_id",
          "nuc_gamma_moiety_type","nuc_gamma_moiety_chain","nuc_gamma_moiety_id",
        #info on site
          "mg_in_site_4a_from_beta",
          "AG-site","G-site",
          "TYPE_ARG","LYS-site","asn_TYPE",
          
        #Closest non-P-loop Lys
          "Lys_res_ch","Lys_res_id",
          "nz-alpha-atom","nz-alpha-dist","nz-gamma-atom","nz-gamma-dist",
        #Closest Arg 
          
          "Arg_res_ch","Arg_res_id",
          "nh1-alpha-atom","nh1-alpha-dist","nh1-gamma-atom","nh1-gamma-dist",
          "nh2-alpha-atom","nh2-alpha-dist","nh2-gamma-atom","nh2-gamma-dist",
          "ne-alpha-atom","ne-alpha-dist","ne-gamma-atom","ne-gamma-dist",
        #Closest Asn  
        "asn_ID",
        "asn_ne2-alpha-atom","asn_ne2-alpha-dist","asn_ne2-gamma-atom","asn_ne2-gamma-dist",
               
        #Other interactions of gamma-phosphate  with surrounding residues
          
          "Surr_N_BB","Surr_N_SCh",
        #P-loop K-3 backbone to gamma
          "gly13-chain","gly13-type","gly13-id","nuc-to-g13-atom", "dist-gly13",
        #P-loop residue availability
            "lys-ploop-info","ploopk-dist",
        #PDB info
        'resolution',"method","pdbname","pdbdate","include",
        #WB-WA interaction and Mg-coordination
        'SerK+1-Mg','WB-Asp/Glu',"WBD-SerK+1_dist","WBD-Mg",
        
        "D+3_CA..K-3_N","D+3_N..K-3_N", "D+4_CA..K-3_N", "D+4_N..K-3_N","rASA_WBD,%","SASA_WBD",
        "water_present"]
    
    
    dists=pd.DataFrame(columns=cols,index=[0])
    
    df_ind=-1
    ik=0
    epa=MMCIFParser(QUIET=True)
    for ind,comp in enumerate(comps):
        print( comp, len(gia_pdb[comp]), "structures of P-loop NTPases available")
        for T,pdb in enumerate(gia_pdb[comp]):
            if T%50==0:
                print( T, pdb[:4],comp)
            pid=pdb[:4]
            path=os.path.join(p_dir,"PDB"+comp,"full", pdb)
            #print(pdb)
            stru=epa.get_structure(pdb, path)
            ###description

            method, resolution, p_name, pdb_date, has_water=get_structinfo(path)
            if type(resolution)==float:
                if resolution>5:   #Skip low-res structures alltogether
                    continue
            ###
            
            
            for s in stru.get_models():               
                
                
                for ch in s.get_chains():
                    
                    #### Only use chains listed as P-loop NTPases in INTERPRO
                    #### Nucleotides erroneously assigned to activator and not P-loop chain will thus be skipped
                    if pdb[:4]+"_"+ch.get_id() in ploop_all_chains["ind_ch"].values:
                        ##List nucleotides
                        nuc=[res for res in ch.get_residues() if res.get_resname()==comp]
                        ####CHECK FOR INCOMPLETE NUCLEOTIDES AND SKIP THEM
                        nuc_g=[nc for nc in nuc if "O1B" in nc]
                        bad=set(nuc)-set(nuc_g)
                        if len(bad)>0:
                            print("Warning! Truncated nucleotides (no beta-phosphate)", pdb[:4], ch.get_id(), [nu.get_resname()+str(nu.get_id()[1]) for nu in bad])
                        nuc=nuc_g
                        if len(nuc)>0:
                            for nucleotide in nuc:
                                df_ind_stored=df_ind                                
                                df_ind=process_one_nucleotide(dists,
                            s,nucleotide, comp,df_ind,
                            resolution,method,p_name,pdb_date,has_water, pdb,
                            alpha,beta,gamma_F,gamma_N,comps_F,
                        asn_cols)
                                if df_ind_stored>df_ind:
                                    #Skipped Nucleotides
                                    df_ind=df_ind_stored
                       
                       


    return dists