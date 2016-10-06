~/Rosetta/main/source/bin/score_jd2.macosclangrelease \
-database ~/Rosetta/main/database \
-in:file:s input/4P79_A.pdb \
-in:membrane \
-mp:setup:spanfiles input/4P79_tr_A.span \
-score:weights mpframework_smooth_fa_2012.wts \
-mp:setup:transform_into_membrane true \
-out:pdb true \
-out:file:scorefile score_transform.sc \
