# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
#
# @brief Script to visualize membrane planes in pymol
# @author JKLeman (julia.koehler1982@gmail.com)

from pymol import cmd

# run scripts
run membrane_planes.py

# basic settings and colorings
util.cbc
bg white
show cartoon
select mem, resn MEM
select emb, resn EMB
show sphere, mem
show sphere, emb
alter emb, vdw=1.5
rebuild
set_view (\
     0.958278537,   -0.025210002,    0.284718215,\
     0.285831660,    0.085169815,   -0.954486012,\
    -0.000186668,    0.996048152,    0.088822074,\
    -0.000005720,    0.000001855, -301.575897217,\
     3.013955355,   -3.607778311,   -1.902338266,\
   241.575897217,  361.575897217,  -20.000000000 )


# for proteins without membrane residue
with_mem = []
iterate (pdb and resn MEM), with_mem.append("MEM")
if "MEM" not in with_mem:

	#create the pseudoatoms
	pseudoatom a1, pos=[-100, -100, -15]
	pseudoatom a2, pos=[-100,  100, -15]
	pseudoatom a3, pos=[ 100, -100, -15]
	pseudoatom a4, pos=[ 100,  100, -15]

	pseudoatom a5, pos=[-100, -100,  15]
	pseudoatom a6, pos=[-100,  100,  15]
	pseudoatom a7, pos=[ 100, -100,  15]
	pseudoatom a8, pos=[ 100,  100,  15]

	draw_plane origin_lo, a1, a2, a3, a4
	draw_plane origin_up, a5, a6, a7, a8

	#set plane transparency, only shows up after 'ray'
	set cgo_transparency, 0.5, origin_lo
	set cgo_transparency, 0.5, origin_up

	# remove pseudoatoms
	delete a1
	delete a2
	delete a3
	delete a4
	delete a5
	delete a6
	delete a7
	delete a8

else:
	# set membrane plane length
	width = 100

	# Read in center
	center_list = []
	cmd.iterate_state(1, "pdb and resn MEM and name CNTR", "center_list.append((x,y,z))" )
	center = XYZCoord(center_list[0][0], center_list[0][1], center_list[0][2])

	# Read in normal position
	normal_list = []
	cmd.iterate_state(1, "pdb and resn MEM and name NORM", "normal_list.append((x, y, z))" )
	normalp = XYZCoord(normal_list[0][0], normal_list[0][1], normal_list[0][2])

	# compute normal vector, leaflet thickness is 15A
	normal = subtract( normalp, center )
	normal = normalize( normal, 15 )

	# get upper and lower center point along normal
	upper_centerp = add( center, normal )
	lower_centerp = subtract( center, normal )
	print upper_centerp
	print lower_centerp

	# get a vector perpendicular (in membrane plane) to normal
	v1 = XYZCoord()
	v1.x = normal.z
	v1.y = normal.z
	v1.z = -normal.x - normal.y
	v1n = normalize( v1, width )

	# get vector perpendicular (in membrane plane) to v1 and normal
	v2 = XYZCoord()
	v2 = cross( normal, v1 )
	v2n = normalize( v2, width )

	# get 4 points defining upper plane
	p1 = add( upper_centerp, v1n )
	p2 = add( upper_centerp, v2n )
	p3 = subtract( upper_centerp, v2n )
	p4 = subtract( upper_centerp, v1n )

	# get 4 points defining the lower plane 
	q1 = add( lower_centerp, v1n )
	q2 = add( lower_centerp, v2n )
	q3 = subtract( lower_centerp, v2n )
	q4 = subtract( lower_centerp, v1n )

	# create pseudoatoms for lower and upper plane
	f=open('tmp.pml', 'w')
	f.write("pseudoatom p1, pos=[" + str(p1.x) + ", " + str(p1.y) + ", " + str(p1.z) + "]\n")
	f.write("pseudoatom p2, pos=[" + str(p2.x) + ", " + str(p2.y) + ", " + str(p2.z) + "]\n")
	f.write("pseudoatom p3, pos=[" + str(p3.x) + ", " + str(p3.y) + ", " + str(p3.z) + "]\n")
	f.write("pseudoatom p4, pos=[" + str(p4.x) + ", " + str(p4.y) + ", " + str(p4.z) + "]\n")

	f.write("pseudoatom q1, pos=[" + str(q1.x) + ", " + str(q1.y) + ", " + str(q1.z) + "]\n")
	f.write("pseudoatom q2, pos=[" + str(q2.x) + ", " + str(q2.y) + ", " + str(q2.z) + "]\n")
	f.write("pseudoatom q3, pos=[" + str(q3.x) + ", " + str(q3.y) + ", " + str(q3.z) + "]\n")
	f.write("pseudoatom q4, pos=[" + str(q4.x) + ", " + str(q4.y) + ", " + str(q4.z) + "]\n\n")
	f.close()

	# create pseudoatoms
	run tmp.pml

	# draw the plane
	draw_plane arbtr_up, p1, p2, p3, p4
	draw_plane arbtr_lo, q1, q2, q3, q4

	# set transparency
	set cgo_transparency, 0.5, "arbtr_lo"
	set cgo_transparency, 0.5, "arbtr_up"

	# remove pseudoatoms
	delete p1
	delete p2
	delete p3
	delete p4
	delete q1
	delete q2
	delete q3
	delete q4

