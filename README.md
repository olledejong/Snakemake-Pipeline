# Snakemake pipeline implementation of FAST-GBS

An implementation of the FAST-GBS pipeline using Snakemake.

## Installation

*Described only for unix systems*

To be able to use this pipeline, some things need to be set-up.

First, clone this repository. Then, within the cloned folder, perform this command to create a virtual environment:

```
python3 -m venv venv
```

Next, we need to activate the virtual environment.

```
source venv/bin/activate
```

Now the virtual environment is activated.

Now Snakemake needs to be installed:

```
pip install snakemake
```

If this is done, go back to the root of the cloned project, and run the snakemake command like so:

(You should determine yourself how many cores you want to use)

```
Snakemake --snakefile Snakefile --cores 6
```
