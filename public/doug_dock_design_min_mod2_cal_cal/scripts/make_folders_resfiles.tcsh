#!/bin/tcsh

#foreach i ( 601 602 603 604 605 606 607 608 609 610 611 )
foreach i ( 610 )
    #foreach j ( ABA APA HLU HPR MAL MPA NLU NVL A04 A05 A06 A07 A12 A20 A30 A31 A33 A34 A43 A44 A45 A48 A68 A69 A78 A80 A82 A83 A84 A91 A92 A94 B00 B01 B02 B03 B04 B05 B06 B07 B12 B19 B21 B27 B28 B30 B31 B35 B36 B38 B40 B44 B47 B48 B49 B50 B51 B53 B54 B56 B57 B58 B59 B60 B61 B62 B63 B67 B74 B92 B93 B94 B95 B96 B97 B99 C00 C01 C02 C03 C04 C05 C11 C12 C15 C16 C20 C26 C27 C30 C36 C40 C41 C42 C43 C53 C54 C55 C60 C61 C80 C81 C83 C84 C85 C86 C87 C88 C89 C90 C91 C92 C93 C94 )
    foreach j ( MPA )
	# make some unique id
	set label = ${i}_${j}

	# make a directory
	mkdir pos_$label

	# make a unique resfile
	sed "s/XXX/$j/g" resfile_pos_${i}.template > resfile_pos_$label
	
    end
end
