~/Rosetta/main/source/bin/score_jd2.macosclangrelease \
-database ~/Rosetta/main/database \
-in:file:s input/4P79_tr_A.pdb \
-in:membrane \
-mp:setup:spanfiles input/4P79_tr_A.span \
-score:weights mpframework_smooth_fa_2012.wts \
-out:file:scorefile score_as-is.sc \
