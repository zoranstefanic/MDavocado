# Methodology of MDavocado plots

MDavocado diagrams can be produced from single molecular dynamics simulation by
runinng a simple Python script with the name of the topology and trajectory
file as inputs: 

`bash 
./MDavocado.py topology.file trajectory.file 
`

This has to be run inside the Python virtual environment where all the
necessary Python packages are installed (more detailed instructions on how to
set up this virtual environment are given on the [install.md](install.md)  

The net result of the calculations are dynamic images for each of the amino
acids in the protein, conveniently grouped by chain and named A.gif, B.gif etc.
Internally this script is relying on MDanalysis [^1,2] to read the trajectory and
toplogy file and to extract the Ramachandran angles for all the amino acids.
The main component of MDavocado program is a Python class **RamachandranPlots**
which is initialized by MDanalysis **Universe**. This class contains a number of
properties (number of frames, number of residues, number of chains etc.) which
are all initialized automatically from the **Universe**.  The only parameter which
user may need to adjust is a number of parts *N* to devide the trajectory into,
and this is set to 10 by default. So, for example, if the trajectory has 400000
frames in total, it will be devided into 10 batches of 40000 frames, and each
dynamical figure will be made by concatenating 10 images into a movie (and
consequently each one of the images will have 40000 points). 

The major obstacle in making one of those images is drawing this many points on
a single plot.  This is accomplished by using another Python library called
[Datashader](https://datashader.org/):w[^3] which is able to draw scatterplots of
arbitrary number of points. It does so by dividing the area of the plot into a
number of bins, and then counting the number of points that fall into each bin
(this is called the aggregation step), and representing the resulting density
of points by a sutable color mapping (where for example red shades would mean a
high density and blue shades low density of points). All the parameters for an
individual plot are determined automatically by Datashader, giving an optimal
image for any number of points, thus avoiding common problems such as
overplotting, oversaturation, undersampling etc. This step is essential, as
otherwise it would be impossible to draw scatterplots with 40000 points without
point overlaps. For each amino acid there will be *N* images, where *N* is the
number of equal parts that the trajectory is devided into.  


In the final step, for each amino acid in the protein, N images are combined
together into a movie in form of the animated GIF, using a method *create_gifs*
of  **RamachandranPlots** class. After this all animated gifs for one protein chain
are combined to an rectangular panel using montage class and named A.gif, B.gif
and so on. This is the final result of the analysis giving dynamic MDavocado
diagrams as the visual represtation of the entire trajectory.  The MDavocado
program requires only two input files: the topology and the trajectory file of
the MD simulation. It uses MDAnalysis to extract the Ramachandran angles
from the trajectories. The Ramachandran angles are then sliced into *N* equal
segments (by default 10) which are stored into n Pandas[^4,5] dataframes which
contain amino-acids as columns and timesteps as rows. For each segment and each
amino acid a φ-ψ plot is drawn using aggregation in the program Datashader.
Then all the images are converted into gif images and grouped by chains using
convert and montage utilities from ImageMagick[^6] representing time evolution of
each particular chain.

[^1]: Gowers, R.; Linke, M.; Barnoud, J.; Reddy, T.; Melo, M.; Seyler, S.;
Domanski, J.; Dotson, D.; Buchoux, S.; Kenney, I.; Beckstein, O. MDAnalysis: A
python package for the rapid analysis of molecular dynamics simulations.
Proceedings of the 15th Pythonin Science Conference. 2016; 
[^2]: Michaud-Agrawal, N.; Denning, E. J.; Woolf, T. B.; Beckstein, O. MDAnalysis: a
toolkit for the analysis of molecular dynamics simulations. J. Comput. Chem.
2011, 32, 2319–2327.  
[^3]: Datashader. https://datashader.org/# (accessed 2024-05-31).
[^4]: pandas development team, T. Pandas-dev/pandas: Pandas. 2020;
https://doi.org/10.5281/zenodo.3509134 (accessed 2024-04-10); 
[^5]: Wes McKinney, Data Structures for Statistical Computing in Python. Proceedings of
the 9th Python in Science Conference. 2010; pp 56 – 61.  
[^6]: ImageMagick Studio LLC, ImageMagick. https://imagemagick.org.

