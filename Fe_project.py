# -*- coding: utf-8 -*-
import sys
sys.path.append("../../..")
from pyfitit import *
import os
def getProjectFolder(): return os.path.dirname(os.path.realpath(__file__))


def moleculeConstructor(project, params):
    for p in params:
        a,b = project.geometryParamRanges[p]
        v = params[p]
        if v<a or v>b: print(f'Param {p} = {v} not in range [{a}, {b}]. Molecule can be incorrect!')
    projectFolder = getProjectFolder()
    m = Molecule(join(projectFolder,'xyz/Fe.xyz'))
    t = Molecule(join(projectFolder,'xyz/tetra.xyz'))
    if project.CN == 2:
        t.rotate([0,1,0],[0,0,0],params['psi']/180*pi)
        t1 = t.copy()
        t2 = t.copy()
        t1.shift([params['r1'], 0, 0])
        m.unionWith(t1)
        t2.shift([params['r1']+params['d12'], 0, 0])
        t2.rotate([0,0,1], [0,0,0], params['phi1']/180*pi)
        m.unionWith(t2)
    elif project.CN == 3:
        t1 = t.copy()
        t2 = t.copy()
        t3 = t.copy()
        t1.shift([params['r1'], 0, 0])
        t2.shift([params['r1'], 0, 0])
        t3.shift([params['r2'], 0, 0])

        t1.rotate([0, 0, 1], [0, 0, 0], params['phi']*pi/180)
        t2.rotate([0, 0, 1], [0, 0, 0], params['phi']*pi/180)
        t3.rotate([0, 0, 1], [0, 0, 0], params['phi']*pi/180)

        t2.rotate([0, 1, 0], [0, 0, 0], 120*pi/180)
        t3.rotate([0, 1, 0], [0, 0, 0], -120*pi/180)
        m.unionWith(t1)
        m.unionWith(t2)
        m.unionWith(t3)
    elif project.CN == 4:
        t1 = t.copy()
        t2 = t.copy()
        t3 = t.copy()
        t4 = t.copy()

        t1.shift([params['r1'], 0, 0])
        t2.shift([params['r1'], 0, 0])
        t3.rotate([0,1,0], [0,0,0], pi)  # rotating so closest oxygen faces Cr
        t3.shift([-params['r2'], 0, 0])
        t4.rotate([0,1,0], [0,0,0], pi)
        t4.shift([-params['r2'], 0, 0])

        t1.rotate([0,1,0], [0,0,0], params['phi1']/360*pi)
        m.unionWith(t1)
        t2.rotate([0,-1,0], [0,0,0], params['phi1']/360*pi)
        m.unionWith(t2)
        t3.rotate([0,0,1], [0,0,0], params['phi2']/360*pi)
        m.unionWith(t3)
        t4.rotate([0,0,-1], [0,0,0], params['phi2']/360*pi)
        m.unionWith(t4)
    elif project.CN == 5:
        t1 = t.copy()
        t2 = t.copy()
        t3 = t.copy()
        t4 = t.copy()
        t5 = t.copy()
        t1.shift([params['r1'], 0, 0])
        t2.shift([params['r2'], 0, 0])
        t3.shift([params['r2']+params['d2'], 0, 0])
        t4.shift([params['r2'], 0, 0])
        t5.shift([params['r2']+params['d2'], 0, 0])
        t1.rotate([0,1,0], [0,0,0], pi/2)

        # угол с вертикальной осью
        for tt in [t2,t3,t4,t5]: tt.rotate([0, 1, 0], [0, 0, 0], -params['phi1']*pi/180)

        t2.rotate([0,0,1], [0,0,0], pi/2)
        t3.rotate([0,0,1], [0,0,0], pi)
        t4.rotate([0,0,1], [0,0,0], -pi/2)

        # ножницы в плоскости 4-х тетра
        t2.rotate([0,0,1], [0,0,0], (90-params['phi2'])*pi/180/2)
        t3.rotate([0,0,1], [0,0,0], -(90-params['phi2'])*pi/180/2)
        t4.rotate([0,0,1], [0,0,0], (90-params['phi2'])*pi/180/2)
        t5.rotate([0,0,1], [0,0,0], -(90-params['phi2'])*pi/180/2)

        m.unionWith(t1)
        m.unionWith(t2)
        m.unionWith(t3)
        m.unionWith(t4)
        m.unionWith(t5)
    elif project.CN == 6:
        t1 = t.copy()
        t2 = t.copy()
        t3 = t.copy()
        t4 = t.copy()
        t5 = t.copy()
        t6 = t.copy()
        t1.shift([params['r1'], 0, 0])
        t6.shift([params['r1'], 0, 0])
        t2.shift([params['r2'], 0, 0])
        t3.shift([params['r2']+params['d2'], 0, 0])
        t4.shift([params['r2'], 0, 0])
        t5.shift([params['r2']+params['d2'], 0, 0])
        t1.rotate([0,1,0], [0,0,0], pi/2)
        m.unionWith(t1)
        t6.rotate([0,1,0], [0,0,0], -pi/2)
        m.unionWith(t6)

        # угол с вертикальной осью
        for tt in [t2,t3,t4,t5]: tt.rotate([0, 1, 0], [0, 0, 0], -params['phi1']*pi/180)

        t2.rotate([0,0,1], [0,0,0], pi/2)
        t3.rotate([0,0,1], [0,0,0], pi)
        t4.rotate([0,0,1], [0,0,0], -pi/2)

        # ножницы в плоскости 4-х тетра
        t2.rotate([0,0,1], [0,0,0], (90-params['phi2'])*pi/180/2)
        t3.rotate([0,0,1], [0,0,0], -(90-params['phi2'])*pi/180/2)
        t4.rotate([0,0,1], [0,0,0], (90-params['phi2'])*pi/180/2)
        t5.rotate([0,0,1], [0,0,0], -(90-params['phi2'])*pi/180/2)

        m.unionWith(t2)
        m.unionWith(t3)
        m.unionWith(t4)
        m.unionWith(t5)

    # return molecule after check
    if not m.checkInteratomicDistance(minDist=0.8):
        print('Warning: there are atoms with distance < minDist')
    return m


def projectConstructor(CN, valence=2, exp='exp/a-Fe2O3.txt'):
    assert CN in [2,3,4,5,6]
    assert valence in [2,3]
    project = Project()

    project.name = f'Fe{CN}_v{valence}'
    project.CN = CN
    project.valence = valence
    
    filePath = join(getProjectFolder(), exp)
    s = readSpectrum(filePath)
    assert len(s.energy) > 0
    if s.energy[0] < 10: s.energy *= 1000
    # Number of spectrum points to use for machine learning prediction (bigger is slower)
    # project.maxSpectrumPoints = 200
    project.useFullSpectrum = True
    project.spectrum = s

    # specify part of experimental spectrum for fitting
    a = 7110; b = 7180
    project.intervals = {
      'fit_norm': [a, b],
      'fit_smooth': [a, b],
      'fit_geometry': [a, b],
      'plot': [a, b]
    }
    # specify ranges of deformations
    if CN == 2:
        project.geometryParamRanges = {
            'psi': [0,70],
            'r1': [1.7, 2.3],
            'd12': [0, 0.2],
            'phi1': [120,180]
        }
    elif CN == 3:
        project.geometryParamRanges = {
            'r1': [1.7, 2.3],
            'r2': [1.7, 2.3],
            'phi': [0, 40],
        }    
    elif CN == 4:
        project.geometryParamRanges = {
            'r1': [1.7, 2.3],
            'r2': [1.7, 2.3],
            'phi1': [65, 180],
            'phi2': [65, 180],
        }
    elif (CN == 5) or (CN == 6):
        project.geometryParamRanges = {
            'r1': [1.7, 2.3],
            'r2': [1.7, 2.3],
            'd2': [0, 0.3],
            'phi1': [0, 30],
            'phi2': [60, 90]
        } 
    # specify parameters of calculation by FDMNES
    project.FDMNES_calc = {
        'Energy range': '-5 0.1 18 0.5 30 2 54 3 130',
        'Green': False,
        'radius': 5,
    }
    # specify convolution parameters and spectrum shift
    project.FDMNES_smooth = {
        'Gamma_hole': 3.9,
        'Ecent': 48,
        'Elarg': 41,
        'Gamma_max': 25,
        'Efermi': 7115.45,
        'norm': 0.0301
    }
    shifts = {2:142, 3:142, 4:142.4, 5:143, 6:143}
    d_shift = 0 if valence == 3 else -3
    project.FDMNES_smooth['shift'] = shifts[CN]+d_shift
    project.moleculeConstructor = MethodType(moleculeConstructor, project)
    return project

