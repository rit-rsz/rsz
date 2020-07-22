################################################################################
# Name : config.py
# Purpose : The purpose of this is to create a file that makes it easier to
# set up this code on a new computer
#Author : Benjamin Vaughan
#Start Date : Oct 11, 2019
#Additional Info
#
################################################################################
import numpy as np

IrisLookupFile = '/home/vaughan/IRISpy/info_issa_map4.txt' #please point to this textfile
# this means setting the path if it is in a different direct than the python
# scripts

IrisDir =  '/home/vaughan/New_Horizons/nh-ebl-pipeline/py/iris/IRISNOHOLES_B4/' #please point this
# to the directory where you are keeping all of your IRIS fits files.
DataDir = '/home/vaughan/IRISpy/' #please point this to where you want
# your saved data files to be stored
FieldsFile = '/home/vaughan/IRISpy/fields.txt' #please point this to the filepath that contains your
# txt file with the fields you want to look at.
# format for txt file should be:
# RA DEC
# where space is your delimeter
fields = np.loadtxt(FieldsFile, unpack=True)
# these are the potential fields to look at, at some point it may be worthwile
# to have these written in a txtfile or something of that nature and then just
# read them into the config file with a call to a function such as np.loadtxt()
# for now, they are still hardcoded into the script
# index = 0 corresponds to the first field, index = 1 to the second, etc.
