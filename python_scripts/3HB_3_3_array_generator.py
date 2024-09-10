# -*- coding: utf-8 -*-
"""
Created on Mon Sep  9 11:10:07 2024

@author: danie
"""

## Import python modules ##
from modules.sw_directories import *
from modules.sw_build_systems import *
import os as os

## Set up manager and builder classes ## 
main_dir = os.getcwd()
manager = PolymerSimulatorDirs(main_dir)
builder = BuildAmberSystems(manager)

## Get rescode for the base trimer ##
builder.SmilesToPDB_GenResCode('CC(CC(=O)OC(C)CC(=O)OC(C)CC(=O)O)O ', '3HB_trimer')

## Parameterize the base trimer ##
builder.parameterize_mol("3HB_trimer")

## Generate polymeric residue code units ## 
builder.GenRescode_4_PolyUnits("3HB_trimer")

## Build the decamer ##
output = builder.gen_polymer_pdb("3HB_trimer", 10)

## Load the base trimer and the decamer ##
base_trimer = manager.load_pdb_filepath("3HB_trimer")
polymer = manager.load_pdb_filepath("3HB_10_polymer")

## Build the final system ##
output = builder.build_3_3_polymer_array("3HB_trimer", "3HB_10_polymer")