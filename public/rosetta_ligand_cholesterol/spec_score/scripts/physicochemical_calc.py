#!/usr/bin/env python3
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.
'''
Functions and executable to calculate the structural similarity between the reference set and query cholesterol binding site.

Author: Brennica Marlow
'''


from re import S
import numpy as np
import pandas as pd
import jellyfish as jf
from scipy.spatial.distance import cosine
import MDAnalysis as mda
import MDAnalysis.lib as mdal
from sklearn.preprocessing import MinMaxScaler
scaler1 = MinMaxScaler(feature_range=(0,1))


import collections
import glob2


pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)

## import hydrophobicity and bulkiness scales
hbs = pd.read_csv("../reference_files/hydro_bulk_scale.csv")
## import reference files
ref_res_int = pd.read_csv("../reference_files/res_interface.csv",sep=",",index_col="clr")
ref_hydro = pd.read_csv("../reference_files/hydropathy.csv",sep=",",index_col="clr")
ref_bulk = pd.read_csv("../reference_files/bulkiness.csv",sep=",",index_col="clr")


def hydro_bulk_type(df):
    df.loc[df["resn"].str.contains("TYR|TRP|PHE"),"type"] = "A"
    df.loc[df["resn"].str.contains("LEU|ALA|MET|VAL|ILE|GLY"),"type"] = "N"
    df.loc[df["resn"].str.contains("ASN|GLN|SER|THR|PRO|CYS"),"type"] = "U"
    df.loc[df["resn"].str.contains("LYS|ARG|GLU|ASP|HIS"),"type"] = "C"
    df = df.set_index("type")
    return df

def calc_res(r,q):
    ## ref
    rtot = r.sum(axis=0)
    rh,rb,rt =  rtot[:2],rtot[1:20],rtot[-8]
    ## query
    qtot = q.sum(axis=0)
    qh,qb,qt =  qtot[:2],qtot[1:20],qtot[-8]
    
    jh = jf.jaro_winkler_similarity(rh,qh)
    jb = jf.jaro_winkler_similarity(rb,qb)
    jt = jf.jaro_winkler_similarity(rt,qt)
    x,y,z = 0.17,0.66,0.17
    return jh*x+jb*y+jt*z

def calc_hydro(r,q):
    ## ref
    rtot = r.to_numpy()
    rh,rb,rt =  rtot[:2],rtot[1:20],rtot[-8]
    ## query
    qtot = q.to_numpy()
    qh,qb,qt =  qtot[:2],qtot[1:20],qtot[-8]
    
    ch = 1-cosine(rh,qh)
    cb = 1-cosine(rb,qb)
    ct = 1-cosine(rt,qt)
    x,y,z = 0.2,0.67,0.13
    return ch*x+cb*y+ct*z
    
def calc_bulk(r,q):
    ## ref
    rtot = r.to_numpy()
    rh,rb,rt =  rtot[:2],rtot[1:20],rtot[-8]
    ## query
    qtot = q.to_numpy()
    qh,qb,qt =  qtot[:2],qtot[1:20],qtot[-8] 
    ch = 1-cosine(rh,qh)
    cb = 1-cosine(rb,qb)
    ct = 1-cosine(rt,qt)
    x,y,z = 0.06,0.94,0.01
    return ch*x+cb*y+ct*z




def physico_chemical_calc(pdbss,querys):
    # Determine interaction site
    pdbs = sorted(glob2.glob(pdbss+"/*"))
    querynames = [x.split(',')[0] for x in open(querys).readlines()]
    querynames = ' '.join(querynames).split()


    ## Count interactions around CLR atom types
    isite_dict = {}
    for p in pdbs:
        
        query_name = str(p.split("/")[-1])

        if query_name in querynames:
            resid1 = []
            resn1 = []
            dist1 = []
            clr1 = []
            infile = mda.Universe(p)
            prot = infile.select_atoms("protein and not backbone")
            clr = infile.select_atoms("resname CLR")
            prt = prot.positions
            clr_atom = clr.atoms.positions
            pairs,dists = mdal.distances.capped_distance(prt,clr_atom,max_cutoff=6,return_distances=True)
            dDist = pd.DataFrame()

            for k,[i,j]  in enumerate(pairs):
                coord1=prt[i]
                coord2=clr_atom[j]
                dist=dists[k]
                resid = mda.core.groups.AtomGroup([prot.atoms[i]]).residues.resids
                resn = mda.core.groups.AtomGroup([prot.atoms[i]]).residues.resnames
                atoms_name = mda.core.groups.AtomGroup([clr.atoms[j]]).names
                resid1.append(resid)
                resn1.append(resn)
                clr1.append(atoms_name)
                dist1.append(dist)
            dDist = pd.DataFrame(np.column_stack([resid1, resn1, dist1,clr1]),columns=["resid","resn","dist","clr"])
            dDist["pdb"] = query_name
            isite_dict[query_name] = dDist

    isite_dicts = collections.OrderedDict(sorted(isite_dict.items())) 




    ## Reformat interaction site dictionary
    plist = list(isite_dicts.keys())
    hydro_dict = {}
    res_int_dict = {}
    bulk_dict = {}


    atoms = ['O1', 'C3', 'C1', 'C2', 'C4', 'C5', 'C19','C10','C6', 'C7', 'C8', 'C9', 
    'C11', 'C12', 'C13','C14', 'C15', 'C16', 'C17' ,'C18' ,'C20', 'C21', 'C22','C23', 'C24', 'C25', 'C26', 'C27']



    for p in plist:
        temp = isite_dicts.get(p)
        cdata = pd.DataFrame.from_dict(temp)
        gdata = cdata.groupby(["resid","resn","clr","pdb"],as_index=False)
        cdata.update(gdata.transform(np.min))
        gdata1 = gdata.nth(0).reset_index(drop=True)
        gdata1["resid_resn"] = gdata1["resid"].astype(str)+"_"+gdata1["resn"]
        gdata2 = gdata1.drop(columns=["resid","resn","pdb"])
        gdata3 = gdata2.pivot(index="resid_resn",columns="clr",values="dist")
        
        for k in atoms:
            if k not in gdata3.columns:
                gdata3[k] = np.nan
        
        gdata3 = gdata3[atoms]
        gdata4 = gdata3.copy()

        gdata4["type"] = gdata4.index
        gdata4.loc[gdata4["type"].str.contains("LEU|ALA|MET|VAL|ILE|GLY"),"type"] = "N"
        gdata4.loc[gdata4["type"].str.contains("TYR|TRP|PHE"),"type"] = "A"
        gdata4.loc[gdata4["type"].str.contains("ASN|GLN|SER|THR|PRO|CYS"),"type"] = "U"
        gdata4.loc[gdata4["type"].str.contains("LYS|ARG|GLU|ASP|HIS"),"type"] = "C"
        gdata4 = gdata4.set_index("type")
        
        gdata44 = gdata4.where(gdata4 == gdata4.min(), other=0).astype(float)
        gdata5 = gdata44.groupby(level=0).sum().T
        for r in ["N","A","U","C"]:#"S",
            if r not in gdata5.columns:
                gdata5[r] = 0
        

        gdata5.loc[gdata5["C"]  > 1,"C"]= 0.25
        gdata5.loc[gdata5["U"]  > 1,"U"]= 0.5
        gdata5.loc[gdata5["A"]  > 1,"A"]= 0.75
        gdata5.loc[gdata5["N"] > 1,"N"] = 1
        gdata5 = gdata5[["N","A","U","C"]]
        
        
        ## residue interface
        gdata6 = gdata5.copy()
        gdata6["res"] = gdata6.sum(axis=1)
        gdata6 = gdata6[['res',"N","A","U","C"]]
        gdata66 = gdata6.drop(columns=["N","A","U","C"]) 
        gdata66[gdata66 == 1.0] = "N" 
        gdata66[gdata66 == 0.75] = "A" 
        gdata66[gdata66 == 0.5] = "U" 
        gdata66[gdata66 == 0.25] = "C" 
        gdata66[gdata66 == 0.0] = "0"
        gdata666 = gdata66.astype(str)
        res_int_dict[p] = gdata66.T
        
        
        ## hydrophobicity
        gdata7 = gdata3.copy()
        gdata77 = gdata7.where(gdata7 == gdata7.min(), other=0).astype(float)
        gdata77["type"] = gdata77.index
        gdata77["resn"] = gdata77["type"].str.split("_",expand=True)[1]#.head(4))
        
    
        gdata77 = gdata77.drop(columns=["type"])
        gdata777 = gdata77.merge(hbs[["resn","scaled UHS"]],on="resn")
        for c in gdata777.columns:
            if c not in ["resn","scaled UHS"]:
                gdata777[c] = np.where(gdata777[c] > 0, gdata777["scaled UHS"],0)
        gdata777 = hydro_bulk_type(gdata777)
        
        gdata777 = gdata777.drop(columns=["scaled UHS","resn"])
        gdata8 = gdata777.groupby(level=0).sum().T
        for r in ["N","A","U","C"]:#"S",
            if r not in gdata8.columns:
                gdata8[r] = 0
    
        gdata9 = gdata8.copy()
        gdata9["res"] = gdata9.sum(axis=1)
        gdata9 = gdata9[['res',"N","A","U","C"]]#'S',
        gdata9 = gdata9.drop(columns=["N","A","U","C"])#'S', 
        hydro_dict[p] = gdata9.T
        
        ## bulkiness
        gdata10 = gdata3.copy()
        gdata111 = gdata10.where(gdata10 == gdata10.min(), other=0).astype(float)
        gdata111["type"] = gdata111.index
        gdata111["resn"] = gdata111["type"].str.split("_",expand=True)[1]#.head(4))
        
    
        gdata111 = gdata111.drop(columns=["type"])
        gdata1111 = gdata111.merge(hbs[["resn","scaled volume (A^3)"]],on="resn")

        for c in gdata1111.columns:
            if c not in ["resn","scaled volume (A^3))"]:
                gdata1111[c] = np.where(gdata1111[c] > 0, gdata1111["scaled volume (A^3)"],0)
        gdata1111 = hydro_bulk_type(gdata1111)
        gdata1111 = gdata1111.drop(columns=["scaled volume (A^3)","resn"])
        gdata12 = gdata1111.groupby(level=0).sum().T
        for r in ["N","A","U","C"]:#"S",
            if r not in gdata12.columns:
                gdata12[r] = 0
    
        gdata13 = gdata12.copy()
        gdata13["res"] = gdata13.sum(axis=1)
        gdata13 = gdata13[['res',"N","A","U","C"]]#'S',
        gdata13 = gdata13.drop(columns=["N","A","U","C"])#'S', 
        bulk_dict[p] = gdata13.T


    ##### create total dataframe from new dictionary ####
    ## residue interface
    slist = [v for k,v in res_int_dict.items()]
    sdata = pd.concat(slist ,axis=0,keys=res_int_dict.keys())
    sdata = sdata.droplevel(1).T

    ## hydrophobicity
    ulist = [ v for k,v in hydro_dict.items()] 
    udata = pd.concat(ulist ,axis=0,keys=hydro_dict.keys())
    udata = udata.droplevel(1).T

    ## bulkiness
    blist = [ v for k,v in bulk_dict.items()] 
    bdata = pd.concat(blist ,axis=0,keys=bulk_dict.keys())
    bdata = bdata.droplevel(1).T



    ###### CALCULATE RESIDUE TYPE SIMILARITY MATCHES #####
    subs = sdata.columns.tolist()
    res = pd.DataFrame()
    ress = []
    for s in subs:
        ds = sdata[s]
        ref_res_int1 = pd.concat([ds,ref_res_int],axis=1)
        for column_name in ref_res_int1.columns:
            ress.append((s,column_name,calc_res(ds,ref_res_int1[column_name])))
    res = pd.DataFrame(ress,columns=["query","ref","jws1"])
    res["res_int"] = scaler1.fit_transform(res["jws1"].values.reshape(-1,1))
    res.drop(columns=["jws1"],inplace=True)


    ###### CALCULATE HYDROPHOBICITY ######
    hydros = []
    hydro = pd.DataFrame()
    for s in subs:
        ds = udata[s]
        ref_hydro1 = pd.concat([ds,ref_hydro],axis=1)
        for column_name in ref_hydro1.columns:
            hydros.append((s,column_name,calc_hydro(ds,ref_hydro1[column_name])))
    hydro = pd.DataFrame(hydros,columns=["query","ref","cos1"])
    hydro["hydro"] = scaler1.fit_transform(hydro["cos1"].values.reshape(-1,1))
    hydro.drop(columns=["cos1"],inplace=True)

    ###### CALCULATE BULKINESS MATCHES ######
    bulk = pd.DataFrame()
    bulks = []
    for s in subs:
        ds = bdata[s]
        ref_bulk1 = pd.concat([ds,ref_bulk],axis=1)
        for column_name in ref_bulk1.columns:
            bulks.append((s,column_name,calc_bulk(ds,ref_bulk1[column_name])))
    bulk = pd.DataFrame(bulks,columns=["query","ref","cos1"])
    bulk["bulk"] = scaler1.fit_transform(bulk["cos1"].values.reshape(-1,1))
    bulk.drop(columns=["cos1"],inplace=True)


    stats_full = res.copy()
    stats_full["hydro"] = hydro["hydro"]
    stats_full["bulk"] = bulk["bulk"]

    return stats_full


