#!/usr/bin/env python3
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.
'''
Functions and executable to calculate the rate of evolution of a query cholesterol binding site.

Author: Brennica Marlow
'''

import numpy as np
import pandas as pd

from sklearn.preprocessing import MinMaxScaler
scaler1 = MinMaxScaler(feature_range=(0,1))
from scipy.stats import ranksums
import MDAnalysis as mda
import MDAnalysis.lib as mdal
import collections
import glob2
import argparse 

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)


def smd(x,y):
    ## takes two arrays
    out = (np.median(x)-np.median(y))/np.sqrt((np.std(x)**2 + np.std(y)**2)/2)
    return out



def conserv_calc(pdbss,querys,lipss,nets,consvs=None,evors=None):
    # Determine interaction site
    pdbs = sorted(glob2.glob(pdbss+"/*"))
    querynames = [x.split(',')[0] for x in open(querys).readlines()]
    querynames = ' '.join(querynames).split()


    isite_dict = {}
    for p in pdbs:
        resn = []
        resid = []
        query_name = str(p.split("/")[-1]) 
        if query_name in querynames:
            infile = mda.Universe(p)
            prot = infile.select_atoms("protein and not backbone")
            clr = infile.select_atoms("resname CLR")
            try:
                nearClr = mdal.NeighborSearch.AtomNeighborSearch(clr)
                nearProt = mdal.NeighborSearch.AtomNeighborSearch(prot)

                closeClr = nearClr.search(prot,6)
                res = nearProt.search(closeClr,6)

                for r in res.residues:
                    resn.append(r.resname)
                    resid.append(r.resid)

            except AttributeError:
                pass

            rf = pd.DataFrame(np.column_stack([resid,resn]),columns=["resid","resn"])
            rf["pdb"] = query_name
           
            isite_dict[query_name] = rf

    isite_dicts = collections.OrderedDict(sorted(isite_dict.items()))   
    isite_df = pd.concat(isite_dicts,axis=0).reset_index(level=0,drop=True)



    ## conservation calculation
    cf = pd.DataFrame()
    if consvs:
        ## consurf
        c_names = ["resid","SEQ","3LATOM","SCORE","COLOR","CONFIDENCE INTERVAL","CONFIDENCE INTERVAL COLORS","MSA DATA","RESIDUE VARIETY"]
        cf = pd.read_csv(consvs,header=None,skiprows=15,sep="\s+",skipfooter=4,engine="python",index_col=False,names=c_names)
    elif evors:
        ## evorator
        cf = pd.read_csv(evors,sep=",")
        cf["pdb"] = cf["pdb_position"].str.split("_",expand=True)[1]
        cf["Unnamed: 0"] = cf["Unnamed: 0"] + 1
        cf = cf.rename(columns={"Unnamed: 0":"resid","Score":"SCORE"})
    else:
        print ("ERROR: No conservation file provided")


    ## mp_lipid_acc pdb
    resid_lip = []
    x = []
    y = []
    z = []
    lips = []
    with open(lipss,"r") as f:
        for lines in f:
            if "ATOM" in lines:
                resid_lip.append(lines[23:28].strip())
                x.append(lines[31:38].strip())
                y.append(lines[39:46].strip())
                z.append(lines[47:54].strip())
                lips.append(lines[60:67].strip())
    lips_data = list(zip(resid_lip,x,y,z,lips))
    lf = pd.DataFrame(columns=["resid","x","y","z","lips"],data=lips_data).astype(float)
    lf1 = lf.groupby("resid",as_index=True)[["x","y","z","lips"]].apply(np.mean).round(3).reset_index()
  

    ## netsurfp-3.0 
    n_names = ["class","name","seq","resid","RSA","ASA","null","palpha","pbeta","pcoil"]
    nf = pd.read_csv(nets,skiprows=21,sep="\s+",header=None,names=n_names)

    ## merge
    dcl = cf.merge(lf1, on="resid")    
    dcln = dcl.merge(nf,on="resid")
    dcln["resid"] = dcln["resid"].astype(str)
    dcln["norm"] = scaler1.fit_transform(dcln.SCORE.values.reshape(-1,1))
    


    ## calculate rate of evolution
    roes = []
    for key in querynames:
        isite_df1 = isite_df[isite_df["pdb"] == key]

        inter = isite_df1.merge(dcln,on="resid",how="left")
        ilist = inter.resid.tolist()
        interL = inter.norm.to_numpy()

        surf = dcln[~dcln["resid"].isin(ilist)]
        surfL = surf[(surf["class"] == "E") & (surf["lips"] == 50.00)].norm.to_numpy() ## NOTE: sometimes mp_lipid_acc gives you 50.0 instead of 50.00 

        w,p = ranksums(surfL,interL,equal_var="False",alternative="greater")
        roes.append((key,smd(interL,surfL),w,p))

    roe = pd.DataFrame(roes,columns=["query","roe","statistic","pvalue"])  

    return roe

