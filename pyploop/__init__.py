import os,time
import pandas as pd
from Bio.PDB import PDBList
__all__ =["download_pdbs"]
from . import describe_finger,process_nucleotide, process_pdb, table_features
def download_pdbs(compound_pdb, ploop_wdir, pdb_dir, log_dir):

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