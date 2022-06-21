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



def get_oneletter(three):
    aa_dict={
        "ALA":"A",
        "ASP":"D",
        "ASN":"N",
        "ARG":"R",
        "CYS":"C",
        "GLY":"G",
        "GLU":"E",
        "GLN":"Q",
        "HIS":"H",
        "ILE":"I",
        "LEU":"L",
        "LYS":"K",
        "MET":"M",
        "PRO":"P",
        "PHE":"F",
        "SER":"S",
        "THR":"T",
        "TYR":"Y",
        "TRP":"W",
        "VAL":"V",
        
        
    }
    return aa_dict[three]





def get_gamma(s, nc, gamma_types=""):
    beta=["O1B","O2B","O3B","N3B","C3B"]
    alpha=["O1A","O2A","O3A"] 
    gamma=["O1G","O2G","O3G","S1G"]
    gamma_F=["F1","F2","F3","F4"]
    
    if gamma_types=="":
        gamma_types=["ALF","AF3","BEF","MGF","VO4"]
    
    gamma_p=[res for res in s.get_residues() if res.get_resname() in gamma_types]
    
    
    gamma_dict={}
    #ik=ik+1
    #print ik
    dist_g_b=[list((gm, min([atg-nc[atbet] for atbet in beta if atbet in nc]))) for gm in gamma_p for atg in gm.get_atoms()]
    if len(dist_g_b)>0:
        g_cl=min(dist_g_b, key= lambda elem: elem[1])
        #print "Gamma", g_cl, [a for a in g_cl[0].get_atoms()]
        if g_cl[1]<10:
            gamma_dict[nc]=g_cl[0]
    if len(gamma_dict)>0:
        return gamma_dict[nc]
    else:
        return 1


def get_WB_Asp(ch,lys_id,nxt, mg_at,all_term_acid, check_hydro=False,distance_thres=5): 
    #Find Walker B Asp and evaluate Mg binding
    not_allowed=['E', 'D', 'S', 'T', 'Y', 'K', 'R', 'H']
    
    WB_found=False
    hydro=False
    preceding=""
    ser=ch[lys_id+1]
    if "OG1" in ser:
        ser_atom=ser["OG1"]
    else:
        if "OG" in ser:
            ser_atom=ser["OG"]
        else:
            return "SER_NO_TERMINAL_HYDROXYL", "", "","","",""
    
    dst_to_a={ox_at:ser_atom-ox_at for ox_at in all_term_acid}
    
    dst_to_a=dict(sorted(dst_to_a.items(),key=operator.itemgetter(1)))
    
    for atom, dist in dst_to_a.items():
        if WB_found or dist>=distance_thres:
            break
        resn=atom.get_parent()
        rid=resn.get_id()[1]
        #WB_Asp_dist=dist
        aspch=resn.get_parent() #should be the same as ch but who knows as ploop domain is supposed to be one chain 
        #if not all([aspch[rid-x].get_resname() not in not_allowed for x in range (1,4) if rid-x in aspch]):
        preceding="".join([get_oneletter(aspch[rid-x].get_resname()) for x in range (1,4) if rid-x in aspch])
        hydro=not any([x in preceding for x in not_allowed])
        if  hydro or not check_hydro:
            
            WB_found=True
            WB_Asp_dist=dist
            WB_Asp_res=resn
            WB_Asp_info=f"{WB_Asp_res.get_resname()}_{WB_Asp_res.get_parent().get_id()}{WB_Asp_res.get_id()[1]}"
    
            two_terminal_ox=[atom for atom in WB_Asp_res.get_atoms() if atom.name in ["OD1","OD2","OE1","OE2"]]
            if mg_at!="NONE":        
                ser_to_mg=ser_atom-mg_at
        
                WB_to_mg=min([mg_at-ox for ox in two_terminal_ox])             
            else:
                ser_to_mg="NO_MG"
                WB_to_mg="NO_MG" 
        #print (preceding,dist, WB_found)       
    if not WB_found:
        WB_Asp_info, WB_to_mg, WB_Asp_dist = "NOT_FOUND", "", ""
        if mg_at!="NONE":        
            ser_to_mg=ser_atom-mg_at
        else:
            ser_to_mg="NO_MG"
    
               
                     
    return WB_Asp_info, WB_Asp_dist, ser_to_mg, WB_to_mg, hydro, preceding
                                
def find_Lys(beta_atoms, all_term_NZs_LYS,all_term_acid, mg_at):
        
        ploop_lys_found=False
        finger_lys_found=False
        ploop_lys=""
        finger_lys=""
        ser_props=["","","","","",""]
        lys_dsts={(NZ,atb):NZ-atb for atb in beta_atoms for NZ in all_term_NZs_LYS}
        if len(lys_dsts)>0:
            lys_dsts={k: v for k, v in sorted(lys_dsts.items(), key=lambda item: item[1])}
            for lysnum, ((NZ,atb), lysdist) in enumerate(lys_dsts.items()):

                lys_id=NZ.get_parent().get_id()[1]                    
                ch=NZ.get_parent().get_parent()
                if not ploop_lys_found:
                    
                    if ch.has_id(lys_id+1): #If next residue number exist in the structure
                        nxt=ch[lys_id+1].get_resname()
                        if nxt=="SER" or nxt=="THR":
                            ploop_lys=NZ.get_parent()
                            ser_props=get_WB_Asp(ch,lys_id,nxt, mg_at,all_term_acid)
                            
                            ploop_lys_found=True

                        else:
                            finger_lys=NZ.get_parent()
                            finger_lys_found=True

                    else:
                        finger_lys=NZ.get_parent()
                        finger_lys_found=True 

                else:
                    if not finger_lys_found and (lys_id!=ploop_lys.get_id()[1] or ch.get_id()!=ploop_lys.get_parent().get_id()) :
                        finger_lys=NZ.get_parent()
                        finger_lys_found=True

                if finger_lys_found and ploop_lys_found:
                    break
        
        return ploop_lys, finger_lys, ser_props
    
def calc_Asn(alpha_atoms, beta_atoms, gamma_atoms,all_term_ND2s_ASN):
    asn_dsts={(NE2,atb):NE2-atb for atb in beta_atoms for NE2 in all_term_ND2s_ASN}
    if len(asn_dsts)>0:
        
            asn_closest_ne,atb_cl=min(asn_dsts, key=asn_dsts.get)
            
            dst_to_a={at:asn_closest_ne-at for at in alpha_atoms}
            dst_to_g={at: asn_closest_ne-at for at in gamma_atoms}           
            
            a_closest=min(dst_to_a, key=dst_to_a.get)
            g_closest=min(dst_to_g, key=dst_to_g.get)
            
            #asn_cols=["asn_ID","asn_TYPE","asn_ne2-alpha-atom","asn_ne2-alpha-dist","asn_ne2-gamma-atom","asn_ne2-gamma-dist"]

            asn_line=f"{asn_closest_ne.get_parent().get_parent().get_id()}{asn_closest_ne.get_parent().get_id()[1]}","",a_closest.id,dst_to_a[a_closest],g_closest.id,dst_to_g[g_closest]
            return asn_line,asn_closest_ne.get_parent()
    else:
        return "",""


def calc_ag_distances(residue,atom_name, alpha_atoms,gamma_atoms):
    """
    Calculate shortest distances from <atom_name> atom of <residue> (Biopython PDB Residue object)
    to alpha-phosphate oxygens (supplied as Biopython PDB Atom object list <alpha_atoms>)
    and to gamma-phosphate oxygens (idem, as <gamma_atoms>)
    """
    #if len(gamma_atoms)>0:
    if residue.has_id(atom_name):
        alpha_atom,alpha_dist=min([(at,residue[atom_name]-at) for at in alpha_atoms], key=lambda x: x[1])     
        gamma_atom,gamma_dist=min([(at,residue[atom_name]-at) for at in gamma_atoms], key=lambda x: x[1])

        return alpha_atom.get_name(),alpha_dist,gamma_atom.get_name(),gamma_dist
    else: return "","","",""
    
    

def list_all_Ns(gamma_atoms, distlimit, excluded_sidechain, excluded_backbone):
    """
    Return two lists of tuples with information on backbone and sidechain nitrogens interacting
    (located closer than <distlimit> Angstrom) with gamma-phosphate oxygens <gamma_atoms>,
    excluding residues provided in lists <excluded_sidechain>, <excluded_backbone>.
    Output lists consist of tuples (interacting_nitrogen_atom,interacting_gamma_atom, distance)
    """
    nc=gamma_atoms[0].get_parent()
    model = nc.get_parent().get_parent()

    atom_list = [atom for atom in model.get_atoms() if atom.element == 'N' and atom.get_parent()!=nc]
    backbone_list= [atom for atom in atom_list if atom.get_name()=="N" and atom.get_parent() not in excluded_backbone]
    if len(backbone_list)==0:
        return "",""
    sidechain_list=[atom for atom in atom_list if atom.get_name()!="N" and atom.get_parent() not in excluded_sidechain]
    ns_bb = NeighborSearch(backbone_list)
    ns_sch= NeighborSearch(sidechain_list)

    nearby_nitro_bb = [(atm,ga,atm-ga) for ga in gamma_atoms for atm in ns_bb.search(ga.coord, distlimit, 'A')]
    nearby_nitro_sch = [(atm,ga,atm-ga) for ga in gamma_atoms for atm in ns_sch.search(ga.coord, distlimit, 'A')]
    sidechains={}
    for atm,ga,dist in nearby_nitro_sch:
        if atm.get_parent() in sidechains:
            sidechains[atm.get_parent()].append((atm,ga,dist))
        else:
            sidechains[atm.get_parent()]=[(atm,ga,dist)]
        

    return nearby_nitro_bb, sidechains

    
def process_nucleotides(dists,
                        s,nuc, comp,df_ind,
                        resolution,method,p_name,pdbdate,has_water, pdb,
                        alpha,beta,gamma_F,gamma_N,comps_F,
                       all_term_NZs_LYS, all_term_NH1s_ARG, all_term_NH2s_ARG, all_term_ND2s_ASN, all_term_acid, ions,
                       asn_cols
                       
                       ):

    ##Calculate the dists for each chain
    for nc in nuc:
        if comp in comps_F:
            gamma=get_gamma(s, nc)
            if type(gamma)==int:
                print(f"Gamma not found! {pdb[:4]} {nc.get_resname()} {nc.get_parent().get_id()}{nc.get_id()[1]}, skipping")
                continue
            gamma_atoms=[at for at in gamma.get_atoms() if at.name in gamma_F]
                
        else:
            gamma_atoms=[at for at in nc.get_atoms() if at.name in gamma_N]
            if len(gamma_atoms)==0:
                print(f"Gamma not found! {pdb[:4]} {nc.get_resname()} {nc.get_parent().get_id()}{nc.get_id()[1]}, skipping")
                continue


        beta_atoms=[at for at in nc.get_atoms() if at.name in beta] 
        alpha_atoms=[at for at in nc.get_atoms() if at.name in alpha]
        
        ##Check MG ion in 5A radius from beta phopsphates
        ion_dsts={(M,atb):M-atb for atb in beta_atoms for M in ions}
        if len(ion_dsts)>0:        
            mg_at,atb_cl=min(ion_dsts, key=ion_dsts.get)
            if mg_at-atb_cl>5: mg_at="NONE" 
        else:
            mg_at="NONE"
                     
        if mg_at!="NONE":
            mg_at_e=mg_at.element
        else:
            mg_at_e="NONE"
            
        #################
        df_ind=df_ind+1
        ################
        
        #Write Mg search results
        dists.loc[df_ind,"mg_in_site_4a_from_beta"]=mg_at_e #this should create the line
        dists.loc[df_ind,"water_present"]=has_water
        
        
        #Write General struct props
        dists.loc[df_ind,['resolution',"method","pdbname","pdbdate","PDBID"]]=resolution,method,p_name,pdbdate,pdb[:4]

        #write Nucleotide info
        
        dists.loc[df_ind,["nuc_type","nuc_chain","nuc_id","model"]]=comp,nc.get_parent().get_id(),nc.get_id()[1],nc.get_parent().get_parent().get_id()

        #write info on gamma-mimic, if any
                  
        if comp in comps_F:
            dists.loc[df_ind,["nuc_gamma_moiety_type","nuc_gamma_moiety_chain","nuc_gamma_moiety_id"]]=gamma.get_resname(),gamma.get_parent().get_id(),gamma.get_id()[1]        
  
                  
        #####LYS
        ploop_lys, finger_lys,ser_props=find_Lys(beta_atoms, all_term_NZs_LYS,all_term_acid, mg_at) 
        
        
        
        dists.loc[df_ind,['WB-Asp/Glu',"WBD-SerK+1_dist",'SerK+1-Mg',"WBD-Mg","is_hydro","preceding_res"]]=ser_props

        #P-loop Lys
        if ploop_lys!="":          
            ploop_ch=ploop_lys.get_parent().get_id() #It should be the same as nucleotide chain. But in real PDBs...
            ploop_id=ploop_lys.get_id()[1]   
            ploop_dist=min([ploop_lys["NZ"]-atb for atb in beta_atoms])
            dists.loc[df_ind,["lys-ploop-info","ploopk-dist"]]= f"{ploop_ch}{ploop_id}",ploop_dist 
        #Finger Lys
        if finger_lys!="":
            Kdist=calc_ag_distances(finger_lys,"NZ", alpha_atoms,gamma_atoms)         
            dists.loc[df_ind,["Lys_res_ch","Lys_res_id",
                              "nz-alpha-atom","nz-alpha-dist",
                              "nz-gamma-atom","nz-gamma-dist"
                             ]]= finger_lys.get_parent().get_id(),finger_lys.get_id()[1],*Kdist        
                 

        
        #####ASN
        asn_line,finger_asn=calc_Asn(alpha_atoms, beta_atoms, gamma_atoms,all_term_ND2s_ASN)
        if asn_line!="":
            dists.loc[df_ind,asn_cols]=asn_line
        
        
        ####ARG
        if len(all_term_NH1s_ARG)>0:


            arg_dsts={(NH,atb):NH-atb for atb in beta_atoms for NH in all_term_NH1s_ARG+all_term_NH2s_ARG}
            

            arg_closest_nhx,atb_cl=min(arg_dsts, key=arg_dsts.get)
            finger_arg=arg_closest_nhx.get_parent()
            dists.loc[df_ind,["Arg_res_ch","Arg_res_id"]]=finger_arg.get_parent().get_id(),finger_arg.get_id()[1]

            #NH1
            Rdist1=calc_ag_distances(finger_arg,"NH1", alpha_atoms,gamma_atoms)     
            dists.loc[df_ind,["nh1-alpha-atom","nh1-alpha-dist","nh1-gamma-atom","nh1-gamma-dist"]]=Rdist1          

            #NH2
            Rdist2=calc_ag_distances(finger_arg,"NH2", alpha_atoms,gamma_atoms)
            dists.loc[df_ind,["nh2-alpha-atom","nh2-alpha-dist","nh2-gamma-atom","nh2-gamma-dist"]]=Rdist2 

            #NE          
            Rdiste=calc_ag_distances(finger_arg,"NE", alpha_atoms,gamma_atoms)
            dists.loc[df_ind,["ne-alpha-atom","ne-alpha-dist","ne-gamma-atom","ne-gamma-dist"]]=Rdiste           
        
        else: finger_arg=""
        #####K-3 residue
        if ploop_lys!="":
            gly13=ploop_lys.get_parent()[ploop_id-3]
            nuctog13a,g13_dist=min([(at,gly13["N"]-at) for at in gamma_atoms], key=lambda x: x[1])          
            dists.loc[df_ind,["gly13-chain","gly13-type","gly13-id",
                              "nuc-to-g13-atom", "dist-gly13"]]=ploop_ch,gly13.get_resname(),ploop_id-3,nuctog13a.get_name(),g13_dist
            excluded_bb=[gly13]
        else:
            excluded_bb=[]
       
 
        ####Additional interactions with other residues

        excluded_sch=[res for res in [ploop_lys,finger_lys,finger_arg,finger_asn] if type(res)!=str]
               
        bb,sch=list_all_Ns(gamma_atoms, 4, excluded_sch,excluded_bb)

        if len(bb)!=0:
            ln_bb="; ".join([f"{N.get_parent().get_resname()}_{N.get_parent().get_parent().get_id()}{N.get_parent().get_id()[1]}..{Og.get_name()}({dist:.2f})" for N,Og,dist in bb])
        else:
            if bb=="":
                print(f"Warning! cannot find backbone nitrogens. {pdb[:4]} {nc.get_resname()} {nc.get_parent().get_id()}{nc.get_id()[1]}")
            ln_bb=""
        
        if len(sch)>0:            
            #ln_sch="; ".join([f"{N.get_parent().get_resname()}_{N.get_parent().get_parent().get_id()}{N.get_parent().get_id()[1]}_{N.get_name()}..{Og.get_name()}({dist:.2f})" for N,Og,dist in sch])          
            sep2=","
            cont_temp='{}..{}({:.2f})'
            ln_sch="; ".join([f"{res.get_resname()}_{res.get_parent().get_id()}{res.get_id()[1]}_{sep2.join([cont_temp.format(N.get_name(),Og.get_name(),dist) for (N,Og,dist) in valuelist if (float(dist)<=3.2 or dist==min(valuelist, key=lambda x:float(x[2]))[2])])}" for res,valuelist in sch.items()])          
            
        else:
            ln_sch=""
        dists.loc[df_ind,["Surr_N_BB","Surr_N_SCh"]]=ln_bb,ln_sch
    return df_ind 
