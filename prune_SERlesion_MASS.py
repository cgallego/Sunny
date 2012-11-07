#!/usr/bin/env python

import os
import os.path
import sys
import string
import time
from sys import argv, stderr, exit
import shlex, subprocess
from numpy import *
from scipy.io import loadmat, savemat
import dicom
import vtk
import Tix
import _mssql
import pymssql
from datetime import date
from operator import itemgetter, attrgetter
from vtk.util.numpy_support import numpy_to_vtk, vtk_to_numpy

import excel_writer
import Tkinter

import gdcm
import vtkgdcm
import vtkgdcmPython

# ******************************************************************************
print "usage: prune_SERLesion.py patient_list.txt"
print '''
Run on MassNonmass root folder processDicoms.py based on list of 
patient_list.txt(list of StudyIDs in folder mass or nonmass)

"usage:" 
python2.7 codeProject/prune_SERLesion.py img_folder/ batchs/bUniqueTotal_formSharmila_june21_MASS.txt
-- This script will

% Copyright (C) Cristina Gallego, University of Toronto, 2012
% Nov 07/12 -  Process previously selected VOI with seeds and handle oversegmentation or undersegmentations
		Outputs: Segmentations saved inside dynamic series/in ExamID folder as VOIlesions_idXXX (where idXXX corresponds to table field Radiology.MassLesion.LesionFindingID)
		Reports: Compilation of lesion extracted Parameters will be made available on \\labshare\amartel_data\Cristina\MassNonmass\reports
----------------------------------------------------------------------
'''

