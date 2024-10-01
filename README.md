![MDavocado_logo](./MDavocado.png)

# MDavocado

This is a repository of the MDavocado utility. It is a set of programes
written in Python which enable quick visualization of Molecular dynamics
trajectories. It can be considerred as an extension for the well known
[MDanalysis](https://www.mdanalysis.org/) suite of programs. 

## What does the name stand for?

The name is the acronym for: **MD** **A**nalysis and **V**isualization **O**f **C**orrelated **A**ngular **D**iagrams.  

# Methodology

The methodology behind MDavocado diagrams is explained in more detail in [methodology.md](doc/methodology.md).
You can watch the presentation of MDavocado program held at PyCon DE & PyData Berlin 2024 [here](https://www.youtube.com/watch?v=2z66DLkue9c).

[![Watch the video](https://img.youtube.com/vi/2z66DLkue9c/0.jpg)](https://www.youtube.com/watch?v=2z66DLkue9c)

# Installation instructions

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

# More detailed installation instructions

## Installing MDAnalysis


```bash
# Then we have to install [MDAnalysis](https://www.mdanalysis.org/).

conda config --add channels conda-forge
conda install mdanalysis

# List packages installed by MDAnalysis
conda list --export > packages.txt

# To exactly replicate this conda environment use this command:
conda create --name MDavocado --file packages.txt
```

## Installing Ruptures

If you would like to use change point detection then install [Ruptures](https://centre-borelli.github.io/ruptures-docs/) library.
Anyway the main command imports this so we need to install it too.

```bash
conda install ruptures
```

## Installing Datashader
[Datashader](https://datashader.org/) we deffinitely need to install as the ploting of avocado diagrams depends on it.

```bash
conda install datashader
```
Trying to import Datashader complains about Dask, so we install [Dask](https://www.dask.org/) too.

```bash
conda install dask
```
## Funding

The development of this project is financed by [Croatian Science Foundation](https://hrzz.hr/en/) under the project grant IP-2019-04-6764. It is a part of [ALOKOMP project](https://alokomp.irb.hr/).

## Citation
A detailed description of the methodology and usage of MDavocado has been published in the Journal of Chemical Information and Modeling in the paper titled *MDavocado: Analysis and Visualization of Protein Motion by Time-Dependent Angular Diagrams* (10.1021/acs.jcim.4c00650).
