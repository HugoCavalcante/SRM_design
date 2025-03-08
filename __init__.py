# -*- coding: utf-8 -*-
""" Init file for module my_rocket_motor.py """
from constants import *
from rocket_motor_classes import Propellant, Grain, StructuralMaterial, CombustionChamber, Bulkhead, Nozzle, RocketMotor
from propellants import *
from motors import *

import pkgutil

#__all__ = [name for _, name, _ in pkgutil.iter_modules(__path__)]
__all__ = [name for name in dir() if not name.startswith('_')]


