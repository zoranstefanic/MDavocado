![MDavocado_logo](./MDavocado.png)

MDavocado
---------

This is a repository of the MDavocado utility. It is a set of programes
written in Python which enable quick visualization of Molecular dynamics
trajectories. It can be considerred as an extension for the well known
[MDanalysis](https://www.mdanalysis.org/) suite of programs. 

What does the name stand for?
----------------------------

The name is the acronym for: MD Analysis and Visualization Of Correlated Angular Diagrams.  


#Installation instructions

Perhaps the easiest way is start using MDavocado is to clone this repository,
and then folow these simple steps.

First we make a clean new conda environment (here we have named this conda environment **MDavocado**, but
you can call it anything you like). All the dependences will be installed from
the packages file into that newly created conda environment. 

```bash
# To exactly replicate this conda environment use this command:
conda create --name MDavocado --file packages.txt
```

Once this is done we activate the new environment. 
```bash
# Activate new conda environment
conda activate MDavocado
```

Now that we are in the new environment we execute the command:
```bash
./MDavocado.py topology.file trajectory.file
```

This will make a number of gif images named A.gif, B.gif ... Each one is a visual representation of the whole trajectory of one chain in your protein.
It may take anywhere from few minutes for smaller trajectories, to a few hours for a long trajectories and bigger proteins, for this analysis to finish.
That is it!

#More detailed installation instructions

##Installing MDAnalysis

```bash
# Then we have to install [MDAnalysis](https://www.mdanalysis.org/).

conda config --add channels conda-forge
conda install mdanalysis

# List packages installed by MDAnalysis
conda list --export > packages.txt

# To exactly replicate this conda environment use this command:
conda create --name MDavocado --file packages.txt
```

##Installing Ruptures
If you would like to use change point detection then install [Ruptures](https://centre-borelli.github.io/ruptures-docs/) library.
Anyway the main command imports this so we need to install it too.

```bash
conda install ruptures
```

##Installing Datashader
[Datashader](https://datashader.org/) we deffinitely need to install as the ploting of avocado diagrams depends on it.

```bash
conda install datashader
```
Trying to import Datashader complains about Dask, so we install [Dask](https://www.dask.org/) too.

```bash
conda install dask
```
