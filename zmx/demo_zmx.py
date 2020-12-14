#!/usr/bin/env/python
"""
Pyrate - Optical raytracing based on Python

Copyright (C) 2014-2020
               by     Moritz Esslinger moritz.esslinger@web.de
               and    Johannes Hartung j.hartung@gmx.net
               and    Uwe Lippmann  uwe.lippmann@web.de
               and    Thomas Heinze t.heinze@uni-jena.de
               and    others

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
"""
#warnings.filterwarnings("ignore")
import sys
import logging
import argparse

import matplotlib.pyplot as plt


from pyrateoptics.raytracer.localcoordinates import LocalCoordinates
from pyrateoptics.raytracer.material.material_isotropic import\
    ConstantIndexGlass
from pyrateoptics.raytracer.globalconstants import standard_wavelength
from pyrateoptics.raytracer.io.zmx import ZMXParser

from pyrateoptics.sampling2d import raster
from pyrateoptics import draw

from pyrateoptics.raytracer.analysis.optical_system_analysis import\
    OpticalSystemAnalysis


logging.basicConfig(level=logging.DEBUG)

# download ZMX files from e.g.:
# http://astro.dur.ac.uk/~rsharp/opticaldesign/
# some good demonstration of coordinate breaks is: FIELDROTATOR-LECT5.ZMX

#file = "BAM25-150-D-A-532.zmx"

parser = argparse.ArgumentParser(description="Read ZMX files.")
parser.add_argument("file", help="File to be interpreted",
                    type=str)
parser.add_argument("--bundletype", nargs='?', help="Bundle type",
                    type=str, default='collimated')
parser.add_argument("--epd", nargs='?', help="Entrance pupil diameter",
                    type=float, default=1.0)
parser.add_argument("--numrays", nargs='?', help="Number of rays",
                    type=int, default=11)
parser.add_argument("--showspot", help="Show spot diagram?",
                    action="store_true")
parser.add_argument("--anglex", help="Angle", type=float, default=0.0)
parser.add_argument("--reverse", help="Send light in reverse direction?",
                    action="store_true")
parser.add_argument("--do_not_draw_surfaces",
                    help="List of surfaces not to be drawn",
                    type=str, default="")
parser.add_argument("--do_not_draw_raybundles",
                    help="List of raybundles not to be drawn",
                    type=str, default="")

parsed = parser.parse_args()

# TODO: add materials via command line

show_spot = parsed.showspot
file_to_read = parsed.file
#file_to_read = "BAM25-150-D-A-532.zmx"
enpd = parsed.epd
num_rays = parsed.numrays
bundletype = parsed.bundletype
anglex = parsed.anglex
reverse = parsed.reverse
surfaces_do_not_draw = parsed.do_not_draw_surfaces.split(",")
raybundles_do_not_draw = [int(s) for s in parsed.do_not_draw_raybundles.split(",") if s != '']

p = ZMXParser(file_to_read, name='ZMXParser')
lctmp = LocalCoordinates.p("tmp")

matdict = {"BK7": ConstantIndexGlass.p(lctmp, 1.5168),
           "LAFN21": ConstantIndexGlass.p(lctmp, 1.788),
           "SF53": ConstantIndexGlass.p(lctmp, 1.72)}

(s, seq) = p.create_optical_system(matdict)

if s is None:
    sys.exit()

initialbundles_dict = p.create_initial_bundle()

osa = OpticalSystemAnalysis(s, seq, name="Analysis")

ray_paths = []

if initialbundles_dict == []:
    initialbundles_dict = [{"radius": enpd*0.5}]

for d in initialbundles_dict:
    if show_spot:
        d["raster"] = raster.RectGrid()
        osa.aim(num_rays*num_rays, d, wave=standard_wavelength)
        osa.draw_spotdiagram()
    else:
        d["raster"] = raster.MeridionalFan()
        d["anglex"] = anglex
        osa.aim(num_rays, d, wave=standard_wavelength)
        ray_paths.append(osa.trace()[0])

import draw3D

if not show_spot:
    draw(s, ray_paths, do_not_draw_surfaces=surfaces_do_not_draw,
         do_not_draw_raybundles=raybundles_do_not_draw)
    #draw3D.draw3D_pyvista(s, ray_paths)
else:
    plt.show()
    #draw3D.draw3D_pyvista(s, ray_paths)

osa.prettyprint()
