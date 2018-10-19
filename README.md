# Matlab_Real_Space
Electronic structure code for molecules and clusters

# Author's Note
TO-DO - Associated with <Book>

# Starting the Program

Run MATLAB and then `RSDFT` inside the RSDFT/ directory for a GUI. Running
`main` will instead run without any GUI, and should be 
considered a fall back method. The included help assumes use of the GUI.

## Manually Creating a Molecule

Use the drop-down menu at the top left
to select the element of the next atom.  Input
the xyz coordinates of the atom in the 3 text boxes
to the right.  Click 'Add Atom'.  The atom will be
added to the list of atoms and will also be shown
in the molecule visualization.  The molecule visualization
can be rotated for a better view by left clicking and dragging.
If a mistake is made, select the offending atoms from the list
and click the 'Delete Selected Atom(s)' button.

## Saving a Molecule
Once the molecule is done, it can be saved by clicking the
'Save' button or selecting 'Save' from the file menu.
Give the molecule a name and choose the file extension to
use.  '.mat' is binary and can only be read by MATLAB, '.dat'
is a text format and can be viewed by any text editor.

## Loading a Previously Saved Molecule

Click the button marked 'Load' or select 'Load'
from the file menu.  Choose either '.mat' or '.dat'
extension and select a previously created molecule.
The atoms that make up the molecule will be added to
the list and visualization. 

## Solving the problem

Once the 'Start' button has been clicked,
the Progress Window will open.  This window
will show the current progress being made on the
calculations in two ways.  At the bottom of the 
window will be a progress bar.  NOTE: the progress
bar determines the percentage completed by using 
the SCF error.  Because the SCF error can increase, 
the progress bar can go backwards.

The second method provides more information.
The text shows at what step the calculations are at.
It also shows the current diagonalization method and
the subsequent error.  Once the calculations are finished,
the eigenvalues and energy statistics will be printed.

The most recent modifications to RSDFT
have the purpose of decreasing the 
time it takes to run the calculations. Some of the 
optimizations used may not be compatible with all 
computers, may not speed up the program much, or in 
the worst case, may slow down execution.  This is why 
a flag is made available to the user, allowing the
user to decide which optimizations to use.  Level 0 
optimization only uses MATLAB code.  This level should 
be compatible with all recent versions of MATLAB.  
The next step up is level 1 and this level allows the 
use of mex files.  Mex files are binary files that have 
been compiled from C code.  These mex files will 
execute faster than their m file counterparts.  The 
downside is that compatibility is not guaranteed. 
*MEX Files are currently unsupported* 

# Settings

The values set in RSDFTsettings.m are the
values that these settings will use when
RSDFT is first started up.

* Enable Charge Density Visualization: Controls
whether a visualization of the molecule's charge
density is shown after execution is completed.

* useIsosurface: {1 | [0]} - determines whether visualization uses isosurface
rendering methods or contour rendering.

# Visualization of Charge Density
Explaination of Charge Density Visualization

If the check box enabling the visualization is
checked, once the calculations are finished, a
new window showing a visualization of the charge density
of the current molecule will open.  The way to interpret
the graph is that the graph is showing a slice through the
volume that is parallel to the xy plane. The z axis and colors 
of the graph represent charge density values.  By using the 
scroll bar on the right, different slices can be viewed.  
The scroll bar on the bottom changes the y axis and color scale.

# Advanced Settings (Unsupported)
* Precondition CG: Control whether or not
preconditioning is used.

* Use Adpative Scheme: Controls whether
or not certain perameters are changed during
execution.  This is mainly used to speedup
execution, but this speed up has two potential
downsides.  The first is that the final eigenvalues 
are slightly different and that under the wrong
conditions, execution might take longer so use this
option with care.

* Diagonalization Method: Use drop-down menu
to select which diagonalization method(s) are
used.
0 --> Lanczos 1st step and chebyshev filtering thereafter
1 --> Lanczos all the time
2 --> Full-Chebyshev subspace iteration first step
      chebyshev filtering thereafter.

* Polynomial Degree for Chebyshev filtering:
Used to specify the degree of the Chebyshev
polynomial.  10 is a good mix of speed and
accuracy.

## Potential problems with MEX
The one potential problem is the lack of a MATLAB 
compatible compiler.  The first time the user attempts 
to compile the mex files, if MATLAB has not been setup,
the user will be prompted to choose a compiler that is
installed on the current computer.  If a compatible 
compiler cannot be found, the compiling cannot continue.
RSDFT is shipped with precompiled mex files, but if these
do not execute properly, then the optimization level must
be set to 0.

There are two main reasons that a mex file might have 
compatibility problems.The first is that the mex file
was compiled for a different operating systemor architecture.
The second problem is that mex files are not guaranteed to
be compatible on a different version of MATLAB than the 
version it was compiled on.  The only solution to these
two problems is to recompile the mex files.  This can
be done two ways.  Either run the script, compileMexFiles.m,
or from the first graphical user interface window, click 
the 'compile C code' button. 

The user has two ways to control the optimization level.
If the graphical user interface is being used, on the 
first window, use the drop down menu near the bottom center
of the window.  If the text based interface is in use, go
into the file, settings.m, and near the bottom should be
a statement, OPTIMIZATIONLEVEL= #.  Replace the number 
with the desired optimization level.  Note that the value
 used set to the variable OPTIMIZATIONLEVEL in settings.m 
is used as the default value in the graphical user interface.
Progress Window Information

