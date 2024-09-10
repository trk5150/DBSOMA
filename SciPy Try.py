# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 13:39:52 2024

@author: tik105
"""
from io import BytesIO
import numpy as np
import scipy
from scipy.sparse import coo_matrix
from scipy.io import mmwrite

matrix = scipy.io.mmread("C://Users//tik105//Desktop//mRNA//single_cell_read_counts//GSM5114468_S7_exp23_matrix//matrix.mtx")

print('read')
target = BytesIO()
mmwrite(target, coo_matrix(matrix), precision=1)
print(target.getvalue().decode('latin1'))