PymolMover Protocol Capture
===========================
KEYWORDS: MEMBRANES, ANALYSIS

Visualization with the PyMOLMover in any type of protocol. In this example, we use the mutate_relax 
application. To run this, adjust the paths in 

	cmd_mutate_relax.sh

Then open PyMOL and run the cmd script in a regular terminal outside of PyMOL. The PyMOL 
window will automatically be updated with the current structure at each frame. If you want to run this
for your protocol of choice, just add the option

	-show_simulation_in_pymol 0 
	
to your options file. If you want to make a movie and want PyMOL to save the individual frames, also add

	-keep_pymol_simulation_history true
	
to your options file.

