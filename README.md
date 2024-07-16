# PyMol-QTAIM

Visualisation of Quantum Topology of Atoms in Molecules (QTAIM) atomic basins through PyMol GUI.
This repository contains the Plugin 'pymol_qtaim.py' for PyMol (tested on version 2.6.0 Open-Source).
This works with outputs from the QTAIM program AIMAll.

# Usage

To load the `pymol_qtaim.py` script as a PyMol plugin, please follow this [link](https://pymolwiki.org/index.php/Plugins).
Otherwise, simply type

- `run path/to/pymol_qtaim.py`
  in the PyMol console to load it in the current PyMol session.
  These methods will add a function named `qtaim_visualiser` to the standard functions of PyMol.

_Note_: To properly run the script you will need an initial object (e.g. xyz file of the geometry on which AIMAll was run) and the corresponding `_atomicfiles_` folder output from an AIMAll calculation with `-iaswrite=true`. The files output from AIMAll for visualisation are named `.iasviz`. The only files that can be actually visualised are the ones obtained with the `-bim=proaim` parameters (i.e. basin integration method : ProAIM).

You can install a python environment from the pymol python, and then add any relevant packages. So you can do

```
....path_to_pymol/bin/python -m venv ~/.venv/pymol --system-site-packages
```

That will give you access to pymol and its installed packages. After that you can just pip install whatever packages are needed (scipy, etc.). Then you can source the venv and run any of the scripts that import pymol without getting errors. Note that pymol comes with `python3.7`.

In addition, if you want to call the script you have written directly from Pymol (i.e. adding it as a plugin or running it directly), you should update the `PYTHONPATH` environment variable.

Doing 
```
export PYTHONPATH=$PYTHONPATH:~/.venv/pymol/lib/python3.7/site-packages
```
will make any packages installed in the venv environment available to Pymol. You have the option to either add this to your `.bashrc` file, or just add it directly to the current terminal from which you will launch Pymol.

## Variables

The variables to the main `qtaim_visualiser` functions are:

1. **selection** : object/geometry uploaded in PyMol (geometry for which the AIMAll calculation was done is required).
2. **file** : path to the .int file of a specific atom contained in an \_atomicfiles folder output from AIMAll. Note that if a path to the \_atomicfiles folder is passed directly then all the .iasviz files will be read automatically. Default is None.
3. **main_color** : color in either normal text (uses PIL library) or HEX(#) form. The default colors are obtained directly from the selection object.
4. **cp_color** : color in either normal text or HEX(#) form. The default color for critical points is green.
5. **iso_rho** : IsoDensity Surface electron density cut-off. Only three values are available from AIMAll (1e-3,2e-3,4e-4). Default is 1e-3.
6. **ias_rho** : Interatomic Surface electron density cut-off. The lower the value (i.e. closer to zero) the longer the gradient paths shown. Default will be set depending on iso_rho. This value can be different than iso_rho.
7. **transparency** : default is 0.0. Any integer between 0 and 1 will set some transparency on the QTAIM objects.

Note: in the current version some functionalities are missing and will be updated soon.

## Example images

Some example AIMAll output files for testing are put in the `tests` folder.
Here are some images generated with the PyMol-QTAIM visualiser and ray tracing functionality of PyMol.

| ![alt text](https://github.com/popelier-group/PyMol-QTAIM/blob/main/iasmesh_points_imidazole.png) |
| :-----------------------------------------------------------------------------------------------: |
|       <b>Figure showing an imidazole QTAIM interatomic surfaces and isodensity surfaces</b>       |

| ![alt text](https://github.com/popelier-group/PyMol-QTAIM/blob/main/iasmesh_points_water.png) |
| :-------------------------------------------------------------------------------------------: |
|    <b>Figure showing a water dimer QTAIM interatomic surfaces and isodensity surfaces</b>     |

# To do

- Add a QTWidget in PyMol GUI.

# License

The MIT License makes this plugin available for everyone. You are more than welcome to help with the development of this repository.
Please cite this github page if you use the plugin for your studies/research.
