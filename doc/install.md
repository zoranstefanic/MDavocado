#Installation

To start using MDavocado the follwing dependances need to be installed first:

It is easiest to make a new conda environment:

```bash

# Create new conda environment named 'MDavocado'
conda create -n MDavocado

# Activate new conda environment
conda activate MDavocado
```

##Installing MDAnalysis
Then we have to install [MDAnalysis](https://www.mdanalysis.org/).

```bash
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
