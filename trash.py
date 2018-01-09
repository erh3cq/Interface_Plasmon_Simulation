# -*- coding: utf-8 -*-
"""
Created on Sat Dec  9 09:20:00 2017

@author: erhog
"""

import numpy as np

test = np.array([0,1,2])
print(type(test))
if isinstance(test,np.ndarray):
    print('True')