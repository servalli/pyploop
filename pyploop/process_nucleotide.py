#!/usr/bin/env python
__doc__ = """
Module for evaluation of binding sites of P-loop NTPases, see references [XXX] and [YYY] .
This file coontains functions for processing of an individual binding site.
More at https://github.com/servalli/pyploop."""

from copy import copy
import operator
from Bio.PDB import NeighborSearch, Selection
from Bio.PDB.Residue import Residue 
from Bio.PDB.Chain import Chain 
from Bio.PDB.Model import Model
from Bio.PDB.SASA import ShrakeRupley
import Bio.PDB

### DEFAULTS ###
DWB_THRES=6

SR_POINTS=200
SV_RADIUS=1.4
SASA_ONLY_PROTEIN=True
SASA_DISTLIMIT=8
AAS=['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
MAX_EXP={
"A" : 106,
"C" : 135,
"D" : 163,
"E" : 194,
"F" : 197,
"G" : 84,
"H" : 184,
"I" : 169,
"K" : 205,
"L" : 164,
"M" : 188,
"N" : 157,
"P" : 136,
"Q" : 198,
"R" : 248,
"S" : 130,
"T" : 142,
"V" : 142,
"W" : 227,
"Y" : 222,
}
RADII_DICT={
        "H": 1.200,
        "D": 1.20, #This might be why the calculation was failing. The radius might be off.
        "HE": 1.400,
        "LI": 1.81,  
        "BE" :1.53,
        "C": 1.700,
        "N": 1.550,
        "O": 1.520,
        "F": 1.470,
        "NA": 2.270,
        "MG": 1.730,
        "AL": 1.84,
        "SI": 2.10,
        "P": 1.800,
        "S": 1.800,
        "CL": 1.750,
        "K": 2.750,
        "CA": 2.310,
        "NI": 1.630,
        "CU": 1.400,
        "ZN": 1.390,
        "SE": 1.900,
        "BR": 1.850,
        "CD": 1.580,
        "I": 1.980,
        "HG": 1.550,
        #incertain about following  
        "V": 2.05,
        "MN": 2.05,
        "FE": 2.05,
        "CO": 2.0,
    }

### Biopython disorder fix ###

def get_unpacked_list(self):
     """
     from https://github.com/biopython/biopython/issues/455
     Returns all atoms from the residue,
     in case of disordered, keep only first alt loc and remove the alt-loc tag
     """
     atom_list = self.get_list()
     undisordered_atom_list = []
     for atom in atom_list:
         if atom.is_disordered():
             atom.altloc=" "
             undisordered_atom_list.append(atom)
         else:
             undisordered_atom_list.append(atom)
     return undisordered_atom_list
 
Bio.PDB.Residue.Residue.get_unpacked_list = get_unpacked_list

### FUNCTIONS ###

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

    dist_g_b=[list((gm, min([atg-nc[atbet] for atbet in beta if atbet in nc]))) for gm in gamma_p for atg in gm.get_atoms()]
    if len(dist_g_b)>0:
        g_cl=min(dist_g_b, key= lambda elem: elem[1])
        
        if g_cl[1]<10:
            gamma_dict[nc]=g_cl[0]
    if len(gamma_dict)>0:
        return gamma_dict[nc]
    else:
        return 1


def exclude_residues(entity,resnames=[], onlyprotein=True):
    
    selection=[r for r in Selection.unfold_entities(entity, 'R') if r.resname not in resnames]        
    if onlyprotein:
        selection=[r for r in selection if r.resname in AAS]        
    return selection


def get_res_around(res, entity, distlimit=7, atname="CA", ban=["HOH"],onlyprotein=True):
    """
    Wrapper to get residues within {distlimit} from atom {atname} of residue {res}
    """
    
    if len(ban)>=0 or onlyprotein:
        selection=exclude_residues(entity,ban,onlyprotein)        
        selection=Selection.unfold_entities(selection, 'A')
    else:
        selection=Selection.unfold_entities(entity, 'A')    
    ns = NeighborSearch(selection)
   
    nearby_res = ns.search(res[atname].coord, distlimit, 'R')
    return nearby_res



def get_pocket(nucleotide, model, radius=15, nucatom="O1B"):
    #nucleotide in nuc
    #model is s
    #(res, entity, distlimit=7, atname="CA", ban=["HOH"],onlyprotein=True):
    """
    Only look at residues relatively around the nucleotide """
    
    selection=Selection.unfold_entities(model, 'A')    
    ns = NeighborSearch(selection)
   
    pocket_res = ns.search(nucleotide[nucatom].coord, radius, 'R')
    pocket_res=sorted(pocket_res,key=lambda x:int(x.id[1]))
    #pocket_res=undisorder(pocket_res)        
    pocket=copy_selected_residues(pocket_res)
    
    return pocket
    
    

def copy_selected_residues(residue_list):
    #Preserve the object structure, expecting a reasonably short residue list
    unique_chains={r.get_parent().id for r in residue_list}    
    
    copy_model=Model("tmp")
    for chain in unique_chains:
        copy_model.add(Chain(chain))
    for res in residue_list:
        copy_model[res.get_parent().id].add(res.copy())    
    return copy_model

def quick_sasa(res, entity,distlimit=8,radius=SV_RADIUS,solvent=["HOH"],atname="CA",onlyprotein=SASA_ONLY_PROTEIN):
    """SASA for a residue via Shrake-Rupley, using immediate neigbouring residues only"""
    
    if distlimit==0:
        distlimit=0.001
    residues=get_res_around(res, entity,distlimit,atname,solvent, onlyprotein)   
    
    tmp_model=copy_selected_residues(residues)
    
    sr = ShrakeRupley(probe_radius=radius, radii_dict=RADII_DICT, n_points=SR_POINTS)
    
    sr.compute(tmp_model, level="R")
    sasa=tmp_model[res.get_parent().id][res.id].sasa
    return sasa


def rel_sasa(res, entity,  distlimit=SASA_DISTLIMIT, radius=SV_RADIUS, solvent=["HOH"], atname="CB", onlyprotein=SASA_ONLY_PROTEIN):
    
    sasa=quick_sasa(res, entity, distlimit, radius, solvent, atname, onlyprotein)
    return sasa*100/MAX_EXP[get_oneletter(res.resname)]

def WBD_sasa(WB_Asp_res, ch):
    if type(WB_Asp_res)!=str:
        if WB_Asp_res.resname=="ASP": #Ignore GLU
            sasa=quick_sasa(WB_Asp_res, ch.get_parent(), distlimit=SASA_DISTLIMIT,onlyprotein=SASA_ONLY_PROTEIN)
            relsasa=sasa*100/MAX_EXP["D"]
            return sasa,relsasa
    return "N/A","N/A"

def WBD_sasa2(WB_Asp_res, ch):
    if type(WB_Asp_res)!=str:
        if WB_Asp_res.resname=="ASP": #Ignore GLU
            sasa=quick_sasa(WB_Asp_res, ch.get_parent(), distlimit=SASA_DISTLIMIT,onlyprotein=SASA_ONLY_PROTEIN)
            relsasa=sasa*100/MAX_EXP["D"]
            return sasa,relsasa
    return "N/A","N/A"




def Ser_Mg(mg_at,ser_atom, WB_found,WB_Asp_res=""):
    if mg_at!="NONE":        
        ser_to_mg=ser_atom-mg_at
        if  WB_found:
            two_terminal_ox=[atom for atom in WB_Asp_res.get_atoms() if atom.name in ["OD1","OD2","OE1","OE2"]]
            WB_to_mg=min([mg_at-ox for ox in two_terminal_ox])             
    else:
        ser_to_mg="NO_MG"
        WB_to_mg="NO_MG" 
    return ser_to_mg, WB_to_mg    

def measure_WA_WB_crest(ch,lys_id, WB_Asp_res):
    if lys_id-3 in ch:
        k_3=ch[lys_id-3]
    else:
        print("Problem with Walker A or disorder")
        return ["WA_CONTINUITY_ERROR"]*4
        

    
    ch_WBD=WB_Asp_res.get_parent() #For real, should be same as ch, kept separate in case of buggy PDB files
    n_WBD=WB_Asp_res.get_id()[1]
    distances=[]
    for follows in [3,4]:
        if ch_WBD.has_id(n_WBD+follows):
            r=ch_WBD[n_WBD+follows]
            dist_1res=[r[atom]-k_3["N"] for atom in ("CA","N")]
        else:
            dist_1res=["WB_CONTINUITY_ERROR"]*2    
        distances.extend(dist_1res)
    
    #D+3_CA, D+3_N, #D+4_CA, D+4_N 
    return distances           
            
        
    

def get_WB_Asp(ch,lys_id,nxt, mg_at,all_term_acid, check_hydro=True,distance_thres=5, recursed=False): 
    #Find Walker B Asp and evaluate Mg binding
    #This is suboptimal and needs rewriting via NeighborSearch
    
    #Output fields:
    #return WB_Asp_info, WB_Asp_dist, ser_to_mg, WB_to_mg, hydro, preceding, *crest_distances {four of them}, sasa_rel, sasa =
    #WB_Asp_info, *[""]*11, 
    
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
            WB_Asp_info = "SER_NO_TERMINAL_HYDROXYL"
            return WB_Asp_info, *[""]*11       
            
    
    dst_to_a={ox_at:ser_atom-ox_at for ox_at in all_term_acid}
    
    dst_to_a=dict(sorted(dst_to_a.items(),key=operator.itemgetter(1)))
    
    for atom, dist in dst_to_a.items():
        if WB_found or dist>=distance_thres:
            break
        resn=atom.get_parent()
        rid=resn.get_id()[1]        
        aspch=resn.get_parent() #should be the same as ch as ploop domain is supposed to be one chain but who knows, there are different kinds of bs in structures
        #if not all([aspch[rid-x].get_resname() not in not_allowed for x in range (1,4) if rid-x in aspch]):
        #print(aspch.id, check_hydro, [f"{r.id[1]} {r.resname}" for r in aspch.get_residues()])
        
        preceding="".join([get_oneletter(aspch[rid-x].get_resname()) for x in range (1,4) if rid-x in aspch])
        hydro=not any([x in preceding for x in not_allowed])
        if  hydro or not check_hydro:
            
            WB_found=True
            WB_Asp_dist=dist
            WB_Asp_res=resn
            rtype=WB_Asp_res.get_resname()
            rchain=WB_Asp_res.get_parent().get_id()
            rnum=WB_Asp_res.get_id()[1]
            WB_Asp_info=f"{rtype}_{rchain}{rnum}"            
            crest_distances=measure_WA_WB_crest(ch,lys_id, WB_Asp_res) 
            sasa, sasa_rel=WBD_sasa(WB_Asp_res, ch)
            ser_to_mg,  WB_to_mg=Ser_Mg(mg_at,ser_atom, WB_found,WB_Asp_res)
    if not WB_found:
        #If we ran out of Asp residues in 5A radius, but could not find WB-Asp, we might be dealing with a protein
        # with a polar substitution in hhhhD/E, such as many Ras-like proteins. We try to rerun search without hydro check:
        if not recursed:
            #Trying again with no hy-check
            output=get_WB_Asp(ch,lys_id,nxt, mg_at,all_term_acid, 
                              check_hydro=False, recursed=True)
            WB_Asp_info, WB_Asp_dist, *misc = output
            #If we still couldnt find anything, give up and announce there is no acidic residue around. Probably, open catalytic site.
            if type(WB_Asp_dist)==str:
                WB_Asp_info = "NOT_FOUND",
                return WB_Asp_info, *[""]*11 
            else:
                ser_to_mg, WB_to_mg, hydro, preceding,*crest_distances, sasa_rel, sasa = misc
        else:
            WB_Asp_info = "NOT_FOUND",
            return WB_Asp_info, *[""]*11        
    #ser_to_mg,  WB_to_mg=Ser_Mg(mg_at,ser_atom, WB_found,WB_Asp_res)    
              
    output= WB_Asp_info, WB_Asp_dist, ser_to_mg, WB_to_mg, hydro, preceding,*crest_distances, sasa_rel, sasa                
    return output
                                
def find_Lys(beta_atoms, all_term_NZs_LYS,all_term_acid, mg_at):
        
        ploop_lys_found=False
        finger_lys_found=False
        ploop_lys=""
        finger_lys=""
        ser_props=[""]*12
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
                            
                            ser_props=get_WB_Asp(ch,lys_id,nxt, mg_at,all_term_acid, distance_thres=DWB_THRES)
                            
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
    #else:
    #    print(residue, gamma_atoms)
    

def list_all_Ns(pocket, nc,gamma_atoms, distlimit, excluded_sidechain, excluded_backbone):
    """
    Return two lists of tuples with information on backbone and sidechain nitrogens interacting
    (located closer than <distlimit> Angstrom) with gamma-phosphate oxygens <gamma_atoms>,
    excluding residues provided in lists <excluded_sidechain>, <excluded_backbone>.
    Output lists consist of tuples (interacting_nitrogen_atom,interacting_gamma_atom, distance)
    """
    #nc=gamma_atoms[0].get_parent()
    #pocket = nc.get_parent().get_parent()

    atom_list = [atom for atom in pocket.get_atoms() if atom.element == 'N' and atom.get_parent()!=nc]
    
    backbone_list= [atom for atom in atom_list if atom.get_name()=="N" and atom.get_parent() not in excluded_backbone]
    if len(backbone_list)==0:
        return "",""
    sidechain_list=[atom for atom in atom_list if atom.get_name()!="N" and atom.get_parent() not in excluded_sidechain]
    
    ns_bb = NeighborSearch(backbone_list)    
    
    nearby_nitro_bb = [(atm,ga,atm-ga) for ga in gamma_atoms for atm in ns_bb.search(ga.coord, distlimit, 'A')]
    
    if len(sidechain_list)>0:
        #print(excluded_sidechain)
        ns_sch= NeighborSearch(sidechain_list)
        #print(ns_sch.search(gamma_atoms[0].coord, distlimit, 'A'))
        nearby_nitro_sch = [(atm,ga,atm-ga) for ga in gamma_atoms for atm in ns_sch.search(ga.coord, distlimit, 'A')]
        #print(nearby_nitro_sch)
        sidechains={}
        for atm,ga,dist in nearby_nitro_sch:
            if atm.get_parent() in sidechains:
                sidechains[atm.get_parent()].append((atm,ga,dist))
            else:
                sidechains[atm.get_parent()]=[(atm,ga,dist)]
    else:
        sidechains=dict()        
    #print(sidechains)
    return nearby_nitro_bb, sidechains

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

def process_one_nucleotide(dists,
                        model,nucleotide, comp,df_ind,
                        resolution,method,p_name,pdbdate,has_water, pdb,
                        alpha,beta,gamma_F,gamma_N,comps_F,asn_cols
                       
                       ):

    nc=nucleotide
    nuc_chain=nc.get_parent().get_id()
    nuc_fullid=nc.get_id()    
    model_id=nc.get_parent().get_parent().get_id() #Only used for our table
    
    
    s=get_pocket(nucleotide, model, radius=24, nucatom="O1B")
    nc=s[nuc_chain][nuc_fullid]

    
    
    all_term_ND2s_ASN=get_terminals(s,atom="ND2",rtype="ASN")
    all_term_NZs_LYS=get_terminals(s,atom="NZ",rtype="LYS")

    ###ARGS
    all_term_NH1s_ARG=get_terminals(s,atom="NH1",rtype="ARG")
                
    all_term_NH2s_ARG=get_terminals(s,atom="NH2",rtype="ARG")
                
    #negative
    all_term_acid=get_terminals(s,atom="OE1",rtype="GLU")+get_terminals(s,atom="OE2",rtype="GLU")+get_terminals(s,atom="OD1",rtype="ASP")+get_terminals(s,atom="OD2",rtype="ASP")

    ions=[atm for atm in s.get_atoms() if atm.element in ["MG","MN","CA","SR"]]
                
    
    
    if comp in comps_F:
            gamma=get_gamma(s, nc)
            if type(gamma)==int:
                print(f"Gamma not found! {pdb[:4]} {nc.get_resname()} {nc.get_parent().get_id()}{nc.get_id()[1]}, skipping")
                return 1
            gamma_atoms=[at for at in gamma.get_atoms() if at.name in gamma_F]
             
    else:
        gamma_atoms=[at for at in nc.get_atoms() if at.name in gamma_N]
        if len(gamma_atoms)==0:
            print(f"Gamma not found! {pdb[:4]} {nc.get_resname()} {nc.get_parent().get_id()}{nc.get_id()[1]}, skipping")
            return 1


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
        
    dists.loc[df_ind,["nuc_type","nuc_chain","nuc_id","model"]]=comp,nc.get_parent().get_id(),nc.get_id()[1],model_id

    #write info on gamma-mimic, if any
                  
    if comp in comps_F:
        dists.loc[df_ind,["nuc_gamma_moiety_type","nuc_gamma_moiety_chain","nuc_gamma_moiety_id"]]=gamma.get_resname(),gamma.get_parent().get_id(),gamma.get_id()[1]        
  
                  
    #####LYS
    ploop_lys, finger_lys,ser_props=find_Lys(beta_atoms, all_term_NZs_LYS,all_term_acid, mg_at) 
        
    ser_prop_list=['WB-Asp/Glu',"WBD-SerK+1_dist",'SerK+1-Mg',"WBD-Mg","is_hydro","preceding_res","D+3_CA..K-3_N","D+3_N..K-3_N", "D+4_CA..K-3_N", "D+4_N..K-3_N","rASA_WBD,%","SASA_WBD"]
    #print( len(ser_prop_list),len(ser_props),ser_props)
    dists.loc[df_ind,ser_prop_list]=ser_props
    #WB_Asp_info, WB_Asp_dist, ser_to_mg, WB_to_mg, hydro, preceding,*crest_distances, sasa_rel, sasa   
    #P-loop Lys
    if ploop_lys!="":          
        ploop_ch=ploop_lys.get_parent().get_id() #It should be the same as nucleotide chain. But in real PDBs...
        ploop_id=ploop_lys.get_id()[1] 
        #print(ploop_id,pdb,ploop_ch)  
        ploop_dist=min([ploop_lys["NZ"]-atb for atb in beta_atoms])
        ploop_gamma=min([ploop_lys["NZ"]-atb for atb in gamma_atoms])
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
    
    excluded_bb=[]
    if ploop_lys!="":
        if ploop_id-3 in ploop_lys.get_parent():
            gly13=ploop_lys.get_parent()[ploop_id-3]
            nuctog13a,g13_dist=min([(at,gly13["N"]-at) for at in gamma_atoms], key=lambda x: x[1])
            nuctog13a_beta,g13_dist_bridg=min([(at,gly13["N"]-at) for at in beta_atoms], key=lambda x: x[1])          
            dists.loc[df_ind,["gly13-chain","gly13-type","gly13-id",
                            "nuc-to-g13-atom", "dist-gly13","dist-gly13-bridging"]]=ploop_ch,gly13.get_resname(),ploop_id-3,nuctog13a.get_name(),g13_dist,g13_dist_bridg
            excluded_bb=[gly13]
        else:
            
            dists.loc[df_ind,"gly13-id"]="ERROR_RETRIEVING_K-3"
              
 
    ####Additional interactions with other residues

    excluded_sch=[res for res in [ploop_lys,finger_lys,finger_arg,finger_asn] if type(res)!=str]
               
    bb,sch=list_all_Ns(s, nc, gamma_atoms, 4, excluded_sch,excluded_bb)
    
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
        #print(sch, len(sch))
        residues=[f"{res.get_resname()}_{res.get_parent().get_id()}{res.get_id()[1]}_{sep2.join([cont_temp.format(N.get_name(),Og.get_name(),dist) for (N,Og,dist) in valuelist if (float(dist)<=3.2 or dist==min(valuelist, key=lambda x:float(x[2]))[2])])}" for res,valuelist in sch.items()]
        ln_sch="; ".join(residues)          
        #print(residues,ln_sch)    
    else:
        ln_sch=""
    dists.loc[df_ind,["Surr_N_BB","Surr_N_SCh"]]=ln_bb,ln_sch
    return df_ind 
    

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
        
        ser_prop_list=['WB-Asp/Glu',"WBD-SerK+1_dist",'SerK+1-Mg',"WBD-Mg","is_hydro","preceding_res","D+3_CA..K-3_N","D+3_N..K-3_N", "D+4_CA..K-3_N", "D+4_N..K-3_N","rASA_WBD,%","SASA_WBD"]
        #print( len(ser_prop_list),len(ser_props),ser_props)
        dists.loc[df_ind,ser_prop_list]=ser_props
    #WB_Asp_info, WB_Asp_dist, ser_to_mg, WB_to_mg, hydro, preceding,*crest_distances, sasa_rel, sasa   
        #P-loop Lys
        if ploop_lys!="":          
            ploop_ch=ploop_lys.get_parent().get_id() #It should be the same as nucleotide chain. But in real PDBs...
            ploop_id=ploop_lys.get_id()[1] 
            #print(ploop_id,pdb,ploop_ch)  
            ploop_dist=min([ploop_lys["NZ"]-atb for atb in beta_atoms])
            ploop_gamma=min([ploop_lys["NZ"]-atb for atb in gamma_atoms])
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
            pch=ploop_lys.get_parent()
            if pch.has_id(ploop_id-3):
                gly13=pch[ploop_id-3]
                nuctog13a,g13_dist=min([(at,gly13["N"]-at) for at in gamma_atoms], key=lambda x: x[1])
                nuctog13a_beta,g13_dist_bridg=min([(at,gly13["N"]-at) for at in beta_atoms], key=lambda x: x[1])          
                dists.loc[df_ind,["gly13-chain","gly13-type","gly13-id",
                                "nuc-to-g13-atom", "dist-gly13","dist-gly13-bridging"]]=ploop_ch,gly13.get_resname(),ploop_id-3,nuctog13a.get_name(),g13_dist,g13_dist_bridg
                excluded_bb=[gly13]
            else:
                dists.loc[df_ind,"gly13-type"]="WA_CONTINUITY_ERROR"
                excluded_bb=[]
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