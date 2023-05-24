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

_Note_: To properly run the script you will need an initial object (e.g. xyz file of the geometry on which AIMAll was run) and the corresponding `_atomicfiles_` folder output from an AIMAll calculation with `-iaswrite=true`.

## Variables

The variables to the main `qtaim_visualiser` functions are:

1. **selection** : object/geometry uploaded in PyMol (geometry for which the AIMAll calculation was done is required).
2. **file** : path to the .int file of a specific atom contained in an \_atomicfiles folder output from AIMAll. Note that if a path to the \_atomicfiles folder is passed directly then all the .int files will be read automatically.
3. **color** : color in either normal text or (r,g,b) form. The default colors are obtained directly from the selection object.
4. **transparency** : default is 0.0. Any integer between 0 and 1 will set some transparency on the QTAIM objects.

Note: in the current version some functionalities are missing and will be updated soon.

# Example images

Some example AIMAll output files for testing are put in the `tests` folder.
Here are some images generated with the PyMol-QTAIM visualiser and ray tracing functionality of PyMol.

| ![alt text](https://github.com/popelier-group/PyMol-QTAIM/blob/main/iasmesh_points_imidazole.png) |
| :-----------------------------------------------------------------------------------------------: |
|       <b>Figure showing an imidazole QTAIM interatomic surfaces and isodensity surfaces</b>       |

| ![alt text](https://github.com/popelier-group/PyMol-QTAIM/blob/main/iasmesh_points_water.png) |
| :-------------------------------------------------------------------------------------------: |
|     <b>Figure showing an imidazole QTAIM interatomic surfaces and isodensity surfaces</b>     |

# To do

- Add a QTWidget in PyMol GUI.
- Add auto-fill capabilities.

# License

The MIT License makes this plugin available for everyone. You are more than welcome to help with the development of this repository.
Please cite this github page if you use the plugin for your studies/research.
