~/Rosetta/main/source/bin/mp_transform.macosclangrelease \
-database ~/Rosetta/main/database \
-in:file:s input/4P79_A.pdb \
-in:membrane \
-mp:setup:spanfiles input/4P79_tr_A.span \
-score:weights mpframework_smooth_fa_2012.wts \
-mp:setup:transform_into_membrane true \
-mp:transform:optimize_embedding true \
-out:pdb true \
-out:file:scorefile score_transform_optimize.sc \
