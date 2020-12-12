import warnings
warnings.filterwarnings("ignore")


import logging

import numpy as np


from VisualOpticsPy import draw3D


from pyrateoptics.sampling2d import raster
from pyrateoptics.raytracer.material.material_isotropic import\
    ConstantIndexGlass
from pyrateoptics.raytracer.surface_shape import Conic, Biconic
from pyrateoptics.raytracer.optical_element import OpticalElement
from pyrateoptics.raytracer.optical_system import OpticalSystem
from pyrateoptics.raytracer.surface import Surface
from pyrateoptics.raytracer.ray import RayBundle

from pyrateoptics.raytracer.aperture import CircularAperture
from pyrateoptics.raytracer.localcoordinates import LocalCoordinates

from pyrateoptics.raytracer.globalconstants import degree, standard_wavelength

from pyrateoptics import draw

from pyrateoptics.raytracer.analysis.optical_system_analysis import\
    OpticalSystemAnalysis
from pyrateoptics.raytracer.helpers import build_pilotbundle_complex

# definition of optical system

# Design: US patent no. 5701202 A, inventor: Koichi Takahashi
# and also: Bo Chen, Alois M. Herkommer, Opt. Express 24, 26999 (2016)

logging.basicConfig(level=logging.INFO)

s = OpticalSystem.p()

lc0 = s.addLocalCoordinateSystem(LocalCoordinates.p(name="object", decz=0.0),
                                 refname=s.rootcoordinatesystem.name)

air = ConstantIndexGlass.p(lc0, 1.0)
glass = ConstantIndexGlass.p(lc0, 1.492)
s.material_background = air

si = -1.

lcD1 = s.addLocalCoordinateSystem(LocalCoordinates.p(name="D1", decz=30.002), refname=lc0.name)
lcS1 = s.addLocalCoordinateSystem(LocalCoordinates.p(name="S1",               decy=-24.028,   decz=26.360, tiltx=-si*14.7*degree, tiltThenDecenter=False), refname=lc0.name)
lcD1prime = s.addLocalCoordinateSystem(LocalCoordinates.p(name="D1prime",     decy=0.,        decz=30.002, tiltx=-si*1.066*degree, tiltThenDecenter=False), refname=lc0.name)
lcD2 = s.addLocalCoordinateSystem(LocalCoordinates.p(name="D2",               decy=-0.251,    decz=43.485, tiltx=-si*1.066*degree, tiltThenDecenter=False), refname=lc0.name)
lcS2 = s.addLocalCoordinateSystem(LocalCoordinates.p(name="S2",               decy=19.109,    decz=33.339, tiltx=si*36.660*degree, tiltThenDecenter=False), refname=lc0.name)
lcD2prime = s.addLocalCoordinateSystem(LocalCoordinates.p(name="D2prime",     decy=-0.251,    decz=43.485, tiltx=si*38.376*degree, tiltThenDecenter=False), refname=lc0.name)
lcD3 = s.addLocalCoordinateSystem(LocalCoordinates.p(name="D3",               decy=-11.858,   decz=28.827, tiltx=si*38.376*degree, tiltThenDecenter=False), refname=lc0.name)
lcS3 = s.addLocalCoordinateSystem(LocalCoordinates.p(name="S3",               decy=-24.028,   decz=26.360, tiltx=-si*14.7*degree, tiltThenDecenter=False), refname=lc0.name)
lcD3prime = s.addLocalCoordinateSystem(LocalCoordinates.p(name="D3prime",     decy=-11.858,   decz=28.827, tiltx=-si*55.019*degree, tiltThenDecenter=False), refname=lc0.name)
lcD4 = s.addLocalCoordinateSystem(LocalCoordinates.p(name="D4",               decy=-23.067,   decz=36.667, tiltx=-si*55.019*degree, tiltThenDecenter=False), refname=lc0.name)
lcS4 = s.addLocalCoordinateSystem(LocalCoordinates.p(name="S4",               decy=-35.215,   decz=18.817, tiltx=-si*47.770*degree, tiltThenDecenter=False), refname=lc0.name)
lcD4prime = s.addLocalCoordinateSystem(LocalCoordinates.p(name="D4prime",     decy=-23.067,   decz=36.667, tiltx=-si*50.668*degree, tiltThenDecenter=False), refname=lc0.name)
lcimage = s.addLocalCoordinateSystem(LocalCoordinates.p(name="image",         decy=-30.892,   decz=43.083, tiltx=-si*50.668*degree, tiltThenDecenter=False), refname=lc0.name)

objsurf = Surface.p(lc0)
D1surf = Surface.p(lcD1)
S1surf = Surface.p(lcS1, shape=Biconic.p(lcS1, curvy=si*1./108.187,
                                         curvx=si*1./73.105,
                                         coefficients=[(0., 0.),
                                                       (-si*5.542e-7, -0.08),
                                                       (-si*8.176e-11, -1.379)]),
                   aperture=CircularAperture.p(lcS1, maxradius=40.0))
D1Psurf = Surface.p(lcD1prime)
D2surf = Surface.p(lcD2)
S2surf = Surface.p(lcS2, shape=Biconic.p(lcS2, curvy=si*1./69.871,
                                         curvx=si*1./60.374, ccy=-0.1368,
                                         ccx=-0.123,
                                         coefficients=[(0., 0.),
                                                       (si*7.233e-11, 29.075),
                                                       (si*4.529e-12, -2.085)]),
                   aperture=CircularAperture.p(lcS2, maxradius=40.0))
D2Psurf = Surface.p(lcD2prime)
D3surf = Surface.p(lcD3)
S3surf = Surface.p(lcS3, shape=Biconic.p(lcS3, curvy=si*1./108.187,
                                         curvx=si*1./73.105,
                                         coefficients=[(0., 0.),
                                                       (-si*5.542e-7, -0.08),
                                                       (-si*8.176e-11, -1.379)]),
                   aperture=CircularAperture.p(lcS3, maxradius=40.0))
D3Psurf = Surface.p(lcD3prime)
D4surf = Surface.p(lcD4)
S4surf = Surface.p(lcS4, shape=Conic.p(lcS4, curv=1./77.772),
                   aperture=CircularAperture.p(lcS4, maxradius=40.0))
D4Psurf = Surface.p(lcD4prime)
imgsurf = Surface.p(lcimage)

elem = OpticalElement.p(lc0, name="HUD")

elem.addMaterial("air", air)
elem.addMaterial("glass", glass)

elem.addSurface("object", objsurf, (None, None))
elem.addSurface("d1", D1surf, (None, None))
elem.addSurface("s1", S1surf, (None, "glass"))
elem.addSurface("d1p", D1Psurf, ("glass", "glass"))
elem.addSurface("d2", D2surf, ("glass", "glass"))
elem.addSurface("s2", S2surf, ("glass", "glass"))
elem.addSurface("d2p", D2Psurf, ("glass", "glass"))
elem.addSurface("d3", D3surf, ("glass", "glass"))
elem.addSurface("s3", S1surf, ("glass", "glass"))
elem.addSurface("d3p", D3Psurf, ("glass", "glass"))
elem.addSurface("d4", D4surf, ("glass", "glass"))
elem.addSurface("s4", S4surf, ("glass", None))
elem.addSurface("d4p", D4Psurf, (None, None))
elem.addSurface("image", imgsurf, (None, None))

s.addElement("HUD", elem)

print(s.rootcoordinatesystem.pprint())

sysseq = [("HUD",
           [
                ("object", {"is_stop": True}),
                ("d1", {}),
                ("s1", {}),
                ("d1p", {}),
                ("d2", {}),
                ("s2", {"is_mirror": True}),
                ("d2p", {}),
                ("d3", {}),
                ("s3", {"is_mirror": True}),
                ("d3p", {}),
                ("d4", {}),
                ("s4", {}),
                ("d4p", {}),
                ("image", {})
            ]
           )
          ]


osa = OpticalSystemAnalysis(s, sysseq, name="Analysis")

print("collimated bundles")
(o, k1, E1) = osa.collimated_bundle(3, {"radius": 2,
                                        "raster": raster.MeridionalFan()},
                                    wave=standard_wavelength)
(o, k2, E2) = osa.collimated_bundle(3, {"radius": 2,
                                        "raster": raster.MeridionalFan(),
                                        "anglex": 15*degree},
                                    wave=standard_wavelength)
(o, k3, E3) = osa.collimated_bundle(3, {"radius": 2,
                                        "raster": raster.MeridionalFan(),
                                        "anglex": -15*degree},
                                    wave=standard_wavelength)


initialbundle1 = RayBundle(x0=o, k0=k1, Efield0=E1, wave=standard_wavelength)
initialbundle2 = RayBundle(x0=o, k0=k2, Efield0=E2, wave=standard_wavelength)
initialbundle3 = RayBundle(x0=o, k0=k3, Efield0=E3, wave=standard_wavelength)
print("performing sequential raytrace")
r1 = s.seqtrace(initialbundle1, sysseq)
r2 = s.seqtrace(initialbundle2, sysseq)
r3 = s.seqtrace(initialbundle3, sysseq)

obj_dx = 0.1
obj_dphi = 5*degree

print("calculating pilotbundles")
pilotbundles = build_pilotbundle_complex(objsurf,
                                         air,
                                         (obj_dx, obj_dx),
                                         (obj_dphi, obj_dphi),
                                         num_sampling_points=3)

# print("seqtrace of rays_pilot")
# rays_pilot = [s.seqtrace(p, sysseq) for p in pilotbundles[2:]]
# only last two bundles hit the next surface

# print("para seqtrace")
# (pilotray, r_pilot) = s.para_seqtrace(pilotbundles[-1], initialbundle1,
#                                       sysseq)

print("calculating XYUV")
(m_obj_stop, m_stop_img) = s.extractXYUV(pilotbundles[-1], sysseq)

print(np.array_str(m_obj_stop, precision=5, suppress_small=True))
print(np.array_str(m_stop_img, precision=5, suppress_small=True))

# print(s.sequence_to_hitlist(sysseq))


#draw(s, [r1, r2, r3])
draw3D.draw3D_pyvista(s, [(r1, "blue"), (r2, "blue"), (r3, "blue")])

# pz = np.array([26.36, 33.339, 18.817])
# py = np.array([-24.028, 19.109, -35.215])