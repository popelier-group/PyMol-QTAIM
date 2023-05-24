# PyMol-QTAIM

Visualisation of Quantum Topology of Atoms in Molecules (QTAIM) atomic basins through PyMol GUI.
This repository contains the Plugin 'pymol_qtaim.py' for PyMol (tested on version 2.5.0 Open-Source).
This works with outputs from the QTAIM program AIMAll.

# Usage

To load the `pymol_qtaim.py` script as a PyMol plugin, please follow this [link](https://pymolwiki.org/index.php/Plugins).
Otherwise, simply type

- `run path/to/pymol_qtaim.py`
  in the PyMol console to load it in the current PyMol session.
  These methods will add a function named "qtaim_visualiser" to the standard functions of PyMol.

_Note_: To properly run the script you will need an initial object (e.g. xyz file of the geometry on which AIMAll was run) and the corresponding _\_atomicfiles_ folder output from an AIMAll calculation with _-iaswrite=true_.

# To do

- Add a QTWidget in PyMol GUI

# License

The MIT License makes this plugin available for everyone. You are more than welcome to help with the development of this repository.
Please cite this github page if you use the plugin for your studies/research.
