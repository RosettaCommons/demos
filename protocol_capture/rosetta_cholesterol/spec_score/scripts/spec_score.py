#!/usr/bin/env python3
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.
'''
Functions and executable to calculate the likelihood a query cholesterol binding site is specific.

Author: Brennica Marlow
'''


import numpy as np
import pandas as pd
import argparse
from sklearn.preprocessing import MinMaxScaler
scaler1 = MinMaxScaler(feature_range=(0,1))

import conservation_calc
import physicochemical_calc

parser = argparse.ArgumentParser()
parser.add_argument('-q', '--query', required=True, help="query pdbs with their respective buried areas")
parser.add_argument('-p', '--pdbs', required=True, help="path to docked pdbs")
parser.add_argument('-c', '--consv', required=True, help="path to consurf file")
parser.add_argument('-l', '--lips', required=True, help="lipid accessible pdb file")
parser.add_argument('-n', '--net', required=True, help="path netsurfp file")
args = parser.parse_args()

## import reference volume scale
sv = pd.read_csv("../reference_files/volume_scale.csv")

pc = physicochemical_calc.physico_chemical_calc(args.pdbs,args.query)
re = conservation_calc.conserv_calc(args.pdbs,args.query,args.lips,args.net,args.consv)
pv = pd.read_csv(args.query,header=None,names=["pdb","buried_area"])

psv = pd.concat([pv,sv[["pdb","buried_area"]]],axis=0,ignore_index=True)
psv["scaled_volume"] = scaler1.fit_transform(psv["buried_area"].values.reshape(-1,1))

a,b,c,d,e,f = 0.55,0.4,0.2,0.2,0.2,0.45
querys = pv["pdb"].unique().tolist()
spec_score = pd.DataFrame()


for q in querys:
    pc1 = pc.loc[pc["query"] == q].copy()
    pc1["Rank"] = pc1[["res_int","hydro","bulk"]].apply(tuple,axis=1)\
             .rank(method='average',ascending=False).astype(int)
    pc2 = pc1.loc[pc1["Rank"] == 2].copy()
    re1 = re.loc[re["query"] == q]
    psv1 = psv.loc[psv["pdb"] == q]
    
    pc2["Escore"] = re1["roe"].values[0]
    pc2["pvalue"] = re1["pvalue"].values[0]
    pc2["volume"] = psv1["scaled_volume"].values[0]
    pc2["Pscore"] = (pc2["res_int"]*b)+(pc2["hydro"]*c)+(pc2["bulk"]*d)+(pc2["volume"]*e)
    pc2["specifcity_score"] = (1/(1+np.exp(pc2["Escore"])))*a +(pc2["Pscore"]*f)    

    spec_score = spec_score.append(pc2,ignore_index=True)


spec_score2 = spec_score[['query','Escore','Pscore', 'specifcity_score']]
spec_score2 = spec_score2.round(2)
spec_score2.to_csv("specificity.csv",sep=",",index=False)

