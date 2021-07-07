# LingPy Tutorial

## Status

The last time this tutorial has been tested was on 07/07/2021 with LingPy Version 2.6.8 and Python 3.6. All dependencies (`pip freeze`) are listed in the file `requirements.txt`.

## Preliminaries

The first version of this tutorial had certain shortcomings, which we have now tried to overcome by providing an updated version. 
If you want to follow the tutorial in its old version, you should make sure to have jupyter notebooks installed, as well as the following packages, which you best install via commandline with the help of `pip`:

```shell
$ pip install lingpy
$ pip install python-igraph
$ pip install segments
```

If you run into trouble installing the `igraph` package, we ask you kindly to turn to the official website of the package at [igraph.org](https://igraph.org) and follow installations issues over there. As long as LingPy works on your system, you can also run the code without the `infomap` clustering algorithm from the `igraaph` package. In this case, simply make sure that you replace the keyword `cluster_method="infomap"` by `cluster_method="upgma"`. 

## Run the Tutorial with Python

After having created a fresh virtual environment, all you need to do to run this tutorial is to install the packages (see *Preliminaries* above) and then run the Python script called `notebook.py`.

```shell
$ python notebook.py
```

This will run all the code listed in the tutorial.

## Run the Tutorial with Jupyter

If you have jupyter notebooks installed on your machine, just open a terminal in this folder and type:

```shell
$ jupyter notebook notebook.ipynb
```

## Inspecting the Tutorial in HTML

If you do not have jupyter notebooks installed, just double-click on the file `notebook.html`, where you can see the tutorial (but won't be able to run it interactively).

