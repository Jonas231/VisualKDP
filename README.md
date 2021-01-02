# VisualKDP

Visual KDP is intended to be a python backend for simulations done with Klein's optical design software (www.ecalculations.com, KDP-2).
As such it serves as a playground for KDP-2. Functions are provided which allow the user to set up optical system files and analyses as .dat files, which can than be input in KDP-2 using the command "input file filename.dat". This is similar to writing KDP-2-macros, however it allows to do repetitive tasks using for-loops in python (The macro language of KDP-2 does not provide loops!). The results are stored by KDP-2 in .dat-files and can be loaded into python with special "read"-functions for nicer visualizations with python (e.g. spot diagrams, aberration fans, footprints etc.). 3D visualizations of KDP optical setup files (parsed to pyrate) can be done with matplotlib and pyvista. 3D-visualization is a major deficiency of KDP-2 as it is.

At the same time this project is intended to learn to use KDP-2 and eventually to understand its source code (in the far future ;-) ), such that additional features such as additional operands can be added. For this purpose this projects aims at a nicer documentation (using sphinx) and automated ways to explore the structure of the Fortran-90 legacy code of KDP-2 provided by J. E. Klein.

Features will continuously be added by the author.

Structure:
- VisualKDP: Playground for KDP-2 with python ("macro" and "reading") functions.
- VisualOptalix: Comparison with Optalix. Parsers to pyrate and KDP, because Optalix provides a very extensive
  library of optical systems, see: http://optenso.com/hos/designfiles.html)
- VisualOpticsPy: Adds some functions to pyrate for display of optical setups with pyvista.
