# -*- coding: utf-8 -*-
"""

Useful functions to calculate geometric properties.
Created on Tue Mar 11 17:01:43 2025

@author: hugo
"""

if 'pi' not in globals():
    from numpy import pi


def circle_area(D):
    """ Area of a circle of diameter D. 
    A = pi D²/4."""
    return pi*(D**2)/4

def circle_perimeter(D):
    """ Perimeter of a circle of diameter D. 
    L = pi D."""
    return pi*D


def hcylinder_area(D, d, L):
    """ Area of the surfaces of a hollow cylinder.
    A = pi*D*L + pi*d*L + pi(D²-d²)/2. D is the external diameter, d is the hollow diameter,
    L is the length (height). """
    return pi*(L*(D+d) +(D**2-d**2)/2)

def BATES_area(D, d, L):
    """ Area of the surfaces of a hollow cylinder  without external side surface (BATES).
    A = pi*d*L + pi(D²-d²)/2. D is the external diameter, d is the hollow diameter,
    L is the length (height). """
    return pi*(L*d +(D**2-d**2)/2)

def cyl_vol(D, d, L):
    """ Volume of a hollow cylinder.
    V = pi L (D²-d²)/4. 
    For a non-hollow cylinder, just make d = 0."""
    return pi*(D**2-d**2)*L/4

