
################################################################################
# NAME : catsrc.py
# DATE STARTED : June 11, 2019
# AUTHORS : Victoria Butler & Dale Mercado
# PURPOSE : Fetches the cluster lookup table
#           This should just be some of the header material
# EXPLANATION :
# CALLING SEQUENCE :
# INPUTS :
#
#
# OUTPUTS :
# REVISION HISTORY :
################################################################################
import scipy.io
import numpy as np
from config import * #(this line will give me access to all directory variables)
import matplotlib.pyplot as plt
import math

def get_clus_params():

    verbose = 1 if not verbose else verbose
