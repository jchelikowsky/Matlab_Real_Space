# Matlab_Real_Space
Electronic structure code for atoms, molecules, clusters and nanocrystals.  

# Author's Note
For background material see:  Introductory Quantum Mechanics with MATLAB: For Atoms, Molecules, Clusters, and Nanocrystals by James R. Chelikowsky.  The book is now available on Amazon.  Some details are given below, but specifics of how the code works are outlined in the book.   

We will be updating this site on a regular basis to illustrate how the MatLab code be used.   Since this is the initial offering, there will likely be "bugs."  In particular, the code is set up to sacrifice "accuracy" for "speed."  Still, the results should be reasonable, albeit not at the accuracy for a "research code.  

# Starting the Program

Run MATLAB and then `RSDFT` inside the RSDFT/ directory for a GUI. Running`main` will instead run without any GUI. The included help assumes use of the GUI.

## Entering coordinates of the system of interest.

Use the drop-down menu at the top left to select the element of the next atom.  Input the x,y,z coordinates of the atom in the three  text boxes to the right.  The coordinates need to be in atomic units!  (1 a.u. is 0.5292 A).   Click 'Add Atom'.  The atom will be
added to the list of atoms and will also be shown in the molecule visualization.  The molecule visualization can be rotated for a better view by left clicking and dragging. If a mistake is made, select the offending atoms from the list and click the 'Delete Selected Atom(s)' button.  Individual atoms can be done too, e.g., putting the atom at the origin.  Only elements in the first two rows of the periodic table are available. 

## Saving your entered coordinates. 
Once the coordinates for an atom, cluster, molecule or nanocrystal are enetered is done, the coordinates can be saved by clicking the
'Save' button or selecting 'Save' from the file menu.   Give the system a name and choose the file extension to use.  '.mat' is binary and can only be read by MATLAB, '.dat' is a text format and can be viewed by any text editor.

## Loading a previously saved system. 

Click the button marked 'Load' or select 'Load' from the file menu.  Choose either '.mat' or '.dat' extension and select a previously created system. The atoms that make up the system will be added to the list and visualization. 

## Preloaded molecules.   

Pre-loaded molecular coordinates are available for some representaive molecules.  Solutions for these molecules have been checked and should work.  Other molecular systems will likely work, but there can be come issues.  For example, the silicon or carbon dimer may not converge for some bond lengths.   


## Solving the problem

Once the 'Start' button has been clicked, the Progress Window will open.  This window will show the current progress being made on the
calculations in two ways.  At the bottom of the  window will be a progress bar.  NOTE: the progress bar determines the percentage completed by using a specified error tolerance in the self-consistent field (SCF).  `Because the SCF error can increase (or go the wrong way), the progress bar can go backwards!

In addition to the "progress bar," text will be displayed to show the progress in the text window.  The text shows the SCF iteration and what the current error is.   Once the SCF error decreases below a specified value, the code will stop and print out the occupied eigenvalues and some empty states.  Various contributions to the total energy of the system are also listed.  More details, such as timings and the eignevalues for each SCF iteration are given in the file rsdft.out.  


The RSDFT code is designed to run in a "simple mode."  Some of the optimizations used may not be compatible with all 
computers.  This is why  a flag is made available to the user, allowing the user to decide which optimizations to use.  Level 0 
optimization only uses MATLAB code.  This level should  be compatible with recent versions of MATLAB, at least after 2014.  The next step up is level 1 and this level allows the use of mex files.  Mex files are binary files that have been compiled from C code.  These mex files will  execute faster than their m file counterparts.  The downside is that compatibility is not guaranteed.

# Settings

The values set in RSDFTsettings.m are the values that these settings will use when RSDFT is first started up.  Many of these values are set up to produce "reasonable results."  There can be issues where the preloaded values may not be optimal for the computation. For example, the code will stop when a tolerance value for the SCF iteration is met.  If one is looking to compute small energy differences, e.g., in bond bending, the accuracy many need to be improved.  This can be done by decreasing the "tol" parameter in this file.   Many of these settings are somewhat technical.  Other settings control the display of charge densities:

* Enable Charge Density Visualization: Controls whether a visualization of the molecule's charge density is shown after execution is completed.

* useIsosurface: {1 | [0]} - determines whether visualization uses isosurface rendering methods or contour rendering.

# Visualization of Charge Density
Explaination of Charge Density Visualization.  If the check box enabling the visualization is checked, once the calculations are finished, a new window showing a visualization of the charge density of the current molecule will open.  The way to interpret
the graph is that the graph is showing a slice through the volume that is parallel to the xy plane. The z axis and colors 
of the graph represent charge density values.  By using the  scroll bar on the right, different slices can be viewed. The scroll bar on the bottom changes the y axis and color scale.

# Advanced Settings (Unsupported)
In constructing this code, several sophisticated algorithms for the eigensolver operation were examined.  These algorithms may or may not work. 

* Precondition CG: Control whether or not preconditioning is used.

* Use Adpative Scheme: Controls whether or not certain perameters are changed during execution.  This is mainly used to speedup
execution, but this speed up has two potential downsides.  The first is that the final eigenvalues are slightly different and that under the wrong conditions, execution might take longer so use this option with care.

* Diagonalization Method: Use drop-down menu to select which diagonalization method(s) are used.
0 --> Lanczos 1st step and chebyshev filtering thereafter
1 --> Lanczos all the time
2 --> Full-Chebyshev subspace iteration first step
      chebyshev filtering thereafter.

* Polynomial Degree for Chebyshev filtering: Used to specify the degree of the Chebyshev polynomial.  10 is a good mix of speed and
accuracy.

## Potential problems with MEX
One potential problem is the lack of a MATLAB compatible compiler.  The first time the user attempts to compile the mex files, if MATLAB has not been setup, the user will be prompted to choose a compiler that is installed on the current computer.  If a compatible 
compiler cannot be found, the compiling cannot continue. RSDFT is shipped with precompiled mex files, but if these do not execute properly, then the optimization level must be set to 0.

There are two main reasons that a mex file might have compatibility problems.The first is that the mex file
was compiled for a different operating systemor architecture. The second problem is that mex files are not guaranteed to
be compatible on a different version of MATLAB than the  version it was compiled on.  The only solution to these problems is to recompile the mex files.  This can be done two ways.  Either run the script, compileMexFiles.m, or from the first graphical user interface window, click  the 'compile C code' button. 

The user has two ways to control the optimization level. If the graphical user interface is being used, on the 
first window, use the drop down menu near the bottom center of the window.  If the text based interface is in use, go
into the file, settings.m, and near the bottom should be a statement, OPTIMIZATIONLEVEL= #.  Replace the number with the desired optimization level.  Note that the value used set to the variable OPTIMIZATIONLEVEL in settings.m is used as the default value in the graphical user interface.


