"""
Implements Gaussian Beam Decomposition:
see e.g: https://doi.org/10.1117/1.OE.54.3.035105

"""

import numpy as np

#
def create_Gaussians_mac(file, wry = 1, wrx = 1, bdy = 1, bdx = 1, TEM00 = 0):
    """
    Define a Gaussian beam
    """
    # Gaussian beam specification commands, p.69
    # WRY ... YZ-plane gaussian beam 1/e^2 diameter of surface 1
    # WRX ... XZ-plane [lens units]
    # BDY ... YZ-plane beam divergence of surface 1
    # BDX ... XZ-plane [mrad]
    # BDY TEM00 ... TEM00 mode laser beam
    # BDX TEM00

    file = open(file, "w")
    def write_line(file, L):
        file.writelines(L)
        file.write("\n")

    write_line(file, L=['u l'])
    write_line(file, L=['wry, '+str(wry)])
    write_line(file, L=['wrx, ' + str(wrx)])
    write_line(file, L=['bdy, '+str(bdy)])
    write_line(file, L=['bdx, ' + str(bdx)])
    if TEM00:
        write_line(file, L=['bdy TEM00'])
        write_line(file, L=['bdx TEM00'])
    write_line(file, L=['end'])
    file.close()

def Gaussian_prop(file, spec = 'BEAM', Nsur = 0, y_fob = 0, x_fob = 0, z_obj_shift = 0, wavelength =''):
    """
    GAUSSIAN BEAM PROPAGATION, p. 166 (KDP-Manual)
    :param file:
    :param spec:
    :return:
    """
    # GAUSSIAN BEAM PROPAGATION, p. 166
    file = open(file, "a+")
    def write_line(file, L):
        file.writelines(L)
        file.write("\n")

    write_line(file, L=['output cp'])
    #write_line(file, L=['output file GBD_results.dat'])
    write_line(file, L=['get isn'])
    write_line(file, L=['write'])
    for i in range(0, Nsur+1):
        write_line(file, L=['beam, '+str(i)+', '+str(y_fob)+', '+str(x_fob)+', '+str(z_obj_shift)+', '+str(wavelength)])
    write_line(file, L=['output tp'])
    file.close()

def multi_Gaussian_prop(file, spec = 'BEAM', Nsur = 0, y_fob = 0, x_fob = 0, z_obj_shift = 0, wavelength =''):
    # GAUSSIAN BEAM PROPAGATION, p. 166
    file = open(file, "a+")

    def write_line(file, L):
        file.writelines(L)
        file.write("\n")

    write_line(file, L=['output cp'])
    #write_line(file, L=['output file GBD_results.dat'])
    write_line(file, L=['get isn'])
    write_line(file, L=['write'])

    for i in range(0, Nsur + 1):
        #write_line(file, L=[str(Nsur)])
        for y in range(0,len(y_fob)):
            for x in range(0, len(x_fob)):

                write_line(file, L=['beam, ' + str(i) + ', ' + str(y_fob[y,x]) + ', ' + str(x_fob[y,x]) + ', ' + str(z_obj_shift) + ', ' + str(wavelength)])

    write_line(file, L=['output tp'])
    file.close()


path_KDP_mac = r'C:\D\optical design software demos (other than Zemax)\free software\KDP-2\GBD.DAT'

create_Gaussians_mac(path_KDP_mac, wry = 0.2, wrx = 0.2, bdy = 1, bdx = 1, TEM00 = 1)
Gaussian_prop(path_KDP_mac, spec = 'BEAM', Nsur = 5, y_fob = 0, x_fob = 0, z_obj_shift = 0, wavelength ='')

def fob_caob_arrays(lfob = 0, nfob = 1, lcaob = 0, ncaob = 1):
    # origin of wavefront on surface 0
    # (depends on fob definition - in [mm] or [angle])
    if lfob != 0 and nfob != 1:
        x_fob_a = np.linspace(-lfob, lfob, num = nfob)
        y_fob_a = np.linspace(-lfob, lfob, num = nfob)
    else:
        x_fob_a = 0
        y_fob_a = 0
    # center of clear aperture assigned to surface 1
    # (for each of the above defined origins)
    if lfob != 0 and nfob != 1:
        x_caob_a = np.linspace(-lcaob, lcaob, num = ncaob)
        y_caob_a = np.linspace(-lcaob, lcaob, num = ncaob)
    else:
        x_caob_a = 0
        y_caob_a = 0

    if lfob != 0 and nfob != 1:
        xx_fob, yy_fob = np.meshgrid(x_fob_a, y_fob_a)
    else:
        xx_fob = 0
        yy_fob = 0
    if lcaob != 0 and ncaob != 1:
        xx_caob, yy_caob = np.meshgrid(x_caob_a, y_caob_a)
    else:
        xx_caob = 0
        yy_caob = 0
    #return x_fob_a, y_fob_a, x_caob_a, y_caob_a
    return xx_fob, yy_fob, xx_caob, yy_caob


xx_fob, yy_fob, xx_caob, yy_caob = fob_caob_arrays(lfob = 0, nfob = 1, lcaob = 0, ncaob = 1)
xx_fob, yy_fob, xx_caob, yy_caob = fob_caob_arrays(lfob = 10, nfob = 11, lcaob = 0, ncaob = 1)

#multi_Gaussian_prop(path_KDP_mac, spec = 'BEAM', Nsur = 5, y_fob = yy_fob, x_fob = xx_fob, z_obj_shift = 0, wavelength ='')
