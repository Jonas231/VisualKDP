import numpy as np
import scipy as sp

import pyrateoptics
from pyrateoptics import build_rotationally_symmetric_optical_system

print("ready")
#from Optalix_parser import OptalixParser  # C:\Users\Jonas\PycharmProjects\VisualOpticsPy\venv\Optalix_parser.py
a = 5

print("ready2")
#filename = r'C:\Program Files\OpTaliX-PRO\examples\Eye\EYE_ASPH.otx'
#Opt1 = OptalixParser(filename, name='OTXParser')

import warnings
warnings.filterwarnings("ignore")


import logging

import numpy as np


from pyrateoptics.raytracer.ray import RayBundle
from pyrateoptics.raytracer.analysis.optical_system_analysis import\
    OpticalSystemAnalysis

from pyrateoptics.optimize.optimize import Optimizer
from pyrateoptics.optimize.optimize_backends import ScipyBackend

from pyrateoptics import build_simple_optical_system, draw

from VisualOpticsPy import draw3D

if 1:
    logging.basicConfig(level=logging.DEBUG)

    wavelength = 0.5876e-3

    # definition of optical system
    (s, sysseq) = build_simple_optical_system(
                    [
                     ({"shape": "Conic"}, {"decz": 0.0}, None, "stop",
                      {"is_stop": True}),
                     ({"shape": "Conic"}, {"decz": 5.0}, 1.5168, "front", {}),
                     ({"shape": "Asphere", "curv": -1./50.,
                       "cc": -1., "coefficients": [0.0, 0.0, 0.0]},
                      {"decz": 20.0}, None, "back", {}),
                     ({"shape": "Conic"}, {"decz": 100.0}, None, "image", {})
                    ],
                    )

    osa = OpticalSystemAnalysis(s, sysseq, name="Analysis")

    (o, k, E0) = osa.collimated_bundle(121, {"startz": -5., "radius": 11.43},
                                       wave=wavelength)
    initialbundle = RayBundle(x0=o, k0=k, Efield0=E0, wave=wavelength)

    # initialbundle = generatebundle(openangle=10.*math.pi/180, numrays=121)


    def meritfunctionrms(my_s):
        """
        Standard rms spot radius merit function without removing centroid.
        """
        initialbundle_local = RayBundle(x0=o, k0=k, Efield0=E0, wave=wavelength)
        rpaths = my_s.seqtrace(initialbundle_local, sysseq)
        # other constructions lead to fill up of initial bundle with intersection
        # values

        # for glassy asphere only one path necessary
        x = rpaths[0].raybundles[-1].x[-1, 0, :]
        y = rpaths[0].raybundles[-1].x[-1, 1, :]

        res = np.sum(x**2 + y**2)

        return res


    backsurf = s.elements["stdelem"].surfaces["back"]
    backsurf.shape.params["curv"].to_variable()
    backsurf.shape.params["cc"].to_variable()
    # A2 not variable
    backsurf.shape.params["A4"].to_variable()
    backsurf.shape.params["A6"].to_variable()

    opt_backend = ScipyBackend(method='Nelder-Mead', tol=1e-9)
    optimi = Optimizer(s, meritfunctionrms, opt_backend,
                       name="Nelder-Mead Optimizer")
    #s = optimi.run()

    r2 = s.seqtrace(initialbundle, sysseq)

    #draw(s, r2)
    draw3D.draw3D_pyvista(s, (r2, "yellow"), vertices = 50)
