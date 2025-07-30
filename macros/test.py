#!/usr/bin/python3
import sys
import importlib
sys.path.append('macros/config') 
module_name = sys.argv[1]
inputs = importlib.import_module(module_name).inputs

params = inputs()
print(params)
