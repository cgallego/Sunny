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
# *****************************************************************************
numberS=5

# Get Root folder ( the directory of the script being run)
path_rootFolder = os.path.dirname(os.path.abspath(__file__))
print path_rootFolder

# init vars
studyTree = []

# Get study image folder
img_folder = sys.argv[1]

# Open filename list
file_ids = open(sys.argv[2],"r")

try:
	for line in file_ids:
		# Get the line: Study#, DicomExam# 
		line = line.split()	
		StudyID = line[0]
		ExamID = line[1]
		MRN = line[2]
		LesionID = line[3]
		print "\n\nGetting studyID# DicomExam# LesionFindingID "
		print StudyID, ExamID, MRN, LesionID
		row=[]	
		
		# Execute mysqlconnect (SQL Query design in Project/SMIAL_Database/SQLQuery_Filtered_MRImagingVisits_kinetics_Foci.sql resutls 305 - 18 July 2012
		print "Executing SQL connection..."		
		conn = pymssql.connect(host='SG12-BIOMATRIX', user='cristina.gallego', password='miag_1234', database='SMIAL_BCAD')
		cur = conn.cursor()
		cur.execute('''
		SELECT     Patient.Patients.MRN, Patient.Patients.StudyNumber, Radiology.RadiologyReport.DicomExamNumber, Radiology.RadiologyReport.DateOfExam, 
                      Lesion.GeoDescription.Description, Lesion.FindingLocation.AdditionDescription, Lesion.GeoDescription.BreastSide, Lesion.FindingLocation.QuadrantID, 
                      Lesion.FindingLocation.ClockStart, Lesion.FindingLocation.ClockEnd, Lesion.FindingLocation.MRISeries, Lesion.FindingLocation.LesionCenterPixelX, 
                      Lesion.FindingLocation.LesionCenterPixelY, Lesion.LesionFindingType.Description AS Expr1, Pathology.PathoReport.HistopathWorstType, 
                      Pathology.PathoResultType.PathoResultTypeID, Pathology.PathoResultType.Description AS Expr7, Lesion.DCEDelayEnhPhaseType.Description AS Expr2, 
                      Radiology.MassLesion.LesionFindingID, Lesion.DCEInitialEnhPhaseType.Description AS Expr3, Radiology.MassLesion.LesionFindingID, Lesion.MassLesionIntEnhanceType.MassIntEnhanceTypeID, 
                      Lesion.MassLesionIntEnhanceType.Description AS Expr4, Lesion.MassLesionShapeType.MassShapeTypeID, Lesion.MassLesionShapeType.Description AS Expr5, 
                      Lesion.MassLesionMarginType.MassMarginTypeID, Lesion.MassLesionMarginType.Description AS Expr6
                FROM         Radiology.RadiologyReport INNER JOIN
                      Radiology.Lesions ON Radiology.RadiologyReport.RadReportID = Radiology.Lesions.RadReportID INNER JOIN
                      Lesion.LesionFindingType ON Radiology.Lesions.LesionFindingTypeID = Lesion.LesionFindingType.LesionFindingTypeID INNER JOIN
                      Lesion.GeoDescription ON Radiology.Lesions.GeoDescriptionID = Lesion.GeoDescription.GeoDescriptionID INNER JOIN
                      Lesion.FindingLocation ON Lesion.GeoDescription.FindingLocationID = Lesion.FindingLocation.FindingLocationID AND 
                      Lesion.GeoDescription.FindingLocationID = Lesion.FindingLocation.FindingLocationID INNER JOIN
                      Visit.ImagingVisits ON Radiology.RadiologyReport.ImagingVisitID = Visit.ImagingVisits.ImagingVisitID INNER JOIN
                      Patient.Patients ON Visit.ImagingVisits.StudyPatID = Patient.Patients.StudyPatID AND Visit.ImagingVisits.StudyPatID = Patient.Patients.StudyPatID AND 
                      Visit.ImagingVisits.StudyPatID = Patient.Patients.StudyPatID AND Visit.ImagingVisits.StudyPatID = Patient.Patients.StudyPatID AND 
                      Visit.ImagingVisits.StudyPatID = Patient.Patients.StudyPatID INNER JOIN
                      Pathology.PathoReport ON Visit.ImagingVisits.ImagingVisitID = Pathology.PathoReport.ImagingVisitID AND 
                      Visit.ImagingVisits.ImagingVisitID = Pathology.PathoReport.ImagingVisitID AND Visit.ImagingVisits.ImagingVisitID = Pathology.PathoReport.ImagingVisitID AND 
                      Visit.ImagingVisits.ImagingVisitID = Pathology.PathoReport.ImagingVisitID AND Visit.ImagingVisits.ImagingVisitID = Pathology.PathoReport.ImagingVisitID INNER JOIN
                      Pathology.PathoResultType ON Pathology.PathoReport.PathoResultTypeID = Pathology.PathoResultType.PathoResultTypeID AND 
                      Pathology.PathoReport.PathoResultTypeID = Pathology.PathoResultType.PathoResultTypeID AND 
                      Pathology.PathoReport.PathoResultTypeID = Pathology.PathoResultType.PathoResultTypeID AND 
                      Pathology.PathoReport.PathoResultTypeID = Pathology.PathoResultType.PathoResultTypeID AND 
                      Pathology.PathoReport.PathoResultTypeID = Pathology.PathoResultType.PathoResultTypeID AND 
                      Pathology.PathoReport.PathoResultTypeID = Pathology.PathoResultType.PathoResultTypeID INNER JOIN
                      Visit.ImagingVisitType ON Visit.ImagingVisits.ImagingVisitTypeID = Visit.ImagingVisitType.ImagingVisitTypeID INNER JOIN
                      Visit.VisitReasonType ON Visit.ImagingVisits.VisitReasonTypeID = Visit.VisitReasonType.VisitReasonTypeID INNER JOIN
                      Radiology.MassLesion ON Radiology.Lesions.LesionFindingID = Radiology.MassLesion.LesionFindingID INNER JOIN
                      Lesion.DCEDelayEnhPhaseType ON Radiology.MassLesion.DCEDelayedPhaseEnhID = Lesion.DCEDelayEnhPhaseType.DCEDelayedPhaseEnhID AND 
                      Radiology.MassLesion.DCEDelayedPhaseEnhID = Lesion.DCEDelayEnhPhaseType.DCEDelayedPhaseEnhID AND 
                      Radiology.MassLesion.DCEDelayedPhaseEnhID = Lesion.DCEDelayEnhPhaseType.DCEDelayedPhaseEnhID AND 
                      Radiology.MassLesion.DCEDelayedPhaseEnhID = Lesion.DCEDelayEnhPhaseType.DCEDelayedPhaseEnhID AND 
                      Radiology.MassLesion.DCEDelayedPhaseEnhID = Lesion.DCEDelayEnhPhaseType.DCEDelayedPhaseEnhID AND 
                      Radiology.MassLesion.DCEDelayedPhaseEnhID = Lesion.DCEDelayEnhPhaseType.DCEDelayedPhaseEnhID AND 
                      Radiology.MassLesion.DCEDelayedPhaseEnhID = Lesion.DCEDelayEnhPhaseType.DCEDelayedPhaseEnhID AND 
                      Radiology.MassLesion.DCEDelayedPhaseEnhID = Lesion.DCEDelayEnhPhaseType.DCEDelayedPhaseEnhID AND 
                      Radiology.MassLesion.DCEDelayedPhaseEnhID = Lesion.DCEDelayEnhPhaseType.DCEDelayedPhaseEnhID AND 
                      Radiology.MassLesion.DCEDelayedPhaseEnhID = Lesion.DCEDelayEnhPhaseType.DCEDelayedPhaseEnhID INNER JOIN
                      Lesion.DCEInitialEnhPhaseType ON Radiology.MassLesion.DCEInitialPhaseEnhTypeID = Lesion.DCEInitialEnhPhaseType.DCEInitialPhaseEnhTypeID AND 
                      Radiology.MassLesion.DCEInitialPhaseEnhTypeID = Lesion.DCEInitialEnhPhaseType.DCEInitialPhaseEnhTypeID AND 
                      Radiology.MassLesion.DCEInitialPhaseEnhTypeID = Lesion.DCEInitialEnhPhaseType.DCEInitialPhaseEnhTypeID AND 
                      Radiology.MassLesion.DCEInitialPhaseEnhTypeID = Lesion.DCEInitialEnhPhaseType.DCEInitialPhaseEnhTypeID AND 
                      Radiology.MassLesion.DCEInitialPhaseEnhTypeID = Lesion.DCEInitialEnhPhaseType.DCEInitialPhaseEnhTypeID AND 
                      Radiology.MassLesion.DCEInitialPhaseEnhTypeID = Lesion.DCEInitialEnhPhaseType.DCEInitialPhaseEnhTypeID AND 
                      Radiology.MassLesion.DCEInitialPhaseEnhTypeID = Lesion.DCEInitialEnhPhaseType.DCEInitialPhaseEnhTypeID AND 
                      Radiology.MassLesion.DCEInitialPhaseEnhTypeID = Lesion.DCEInitialEnhPhaseType.DCEInitialPhaseEnhTypeID AND 
                      Radiology.MassLesion.DCEInitialPhaseEnhTypeID = Lesion.DCEInitialEnhPhaseType.DCEInitialPhaseEnhTypeID AND 
                      Radiology.MassLesion.DCEInitialPhaseEnhTypeID = Lesion.DCEInitialEnhPhaseType.DCEInitialPhaseEnhTypeID INNER JOIN
                      Lesion.MassLesionIntEnhanceType ON Radiology.MassLesion.MassIntEnhanceTypeID = Lesion.MassLesionIntEnhanceType.MassIntEnhanceTypeID AND 
                      Radiology.MassLesion.MassIntEnhanceTypeID = Lesion.MassLesionIntEnhanceType.MassIntEnhanceTypeID AND 
                      Radiology.MassLesion.MassIntEnhanceTypeID = Lesion.MassLesionIntEnhanceType.MassIntEnhanceTypeID AND 
                      Radiology.MassLesion.MassIntEnhanceTypeID = Lesion.MassLesionIntEnhanceType.MassIntEnhanceTypeID AND 
                      Radiology.MassLesion.MassIntEnhanceTypeID = Lesion.MassLesionIntEnhanceType.MassIntEnhanceTypeID AND 
                      Radiology.MassLesion.MassIntEnhanceTypeID = Lesion.MassLesionIntEnhanceType.MassIntEnhanceTypeID AND 
                      Radiology.MassLesion.MassIntEnhanceTypeID = Lesion.MassLesionIntEnhanceType.MassIntEnhanceTypeID AND 
                      Radiology.MassLesion.MassIntEnhanceTypeID = Lesion.MassLesionIntEnhanceType.MassIntEnhanceTypeID AND 
                      Radiology.MassLesion.MassIntEnhanceTypeID = Lesion.MassLesionIntEnhanceType.MassIntEnhanceTypeID AND 
                      Radiology.MassLesion.MassIntEnhanceTypeID = Lesion.MassLesionIntEnhanceType.MassIntEnhanceTypeID INNER JOIN
                      Lesion.MassLesionMarginType ON Radiology.MassLesion.MassMarginTypeID = Lesion.MassLesionMarginType.MassMarginTypeID AND 
                      Radiology.MassLesion.MassMarginTypeID = Lesion.MassLesionMarginType.MassMarginTypeID AND 
                      Radiology.MassLesion.MassMarginTypeID = Lesion.MassLesionMarginType.MassMarginTypeID AND 
                      Radiology.MassLesion.MassMarginTypeID = Lesion.MassLesionMarginType.MassMarginTypeID AND 
                      Radiology.MassLesion.MassMarginTypeID = Lesion.MassLesionMarginType.MassMarginTypeID AND 
                      Radiology.MassLesion.MassMarginTypeID = Lesion.MassLesionMarginType.MassMarginTypeID AND 
                      Radiology.MassLesion.MassMarginTypeID = Lesion.MassLesionMarginType.MassMarginTypeID AND 
                      Radiology.MassLesion.MassMarginTypeID = Lesion.MassLesionMarginType.MassMarginTypeID AND 
                      Radiology.MassLesion.MassMarginTypeID = Lesion.MassLesionMarginType.MassMarginTypeID AND 
                      Radiology.MassLesion.MassMarginTypeID = Lesion.MassLesionMarginType.MassMarginTypeID INNER JOIN
                      Lesion.MassLesionShapeType ON Radiology.MassLesion.MassShapeTypeID = Lesion.MassLesionShapeType.MassShapeTypeID AND 
                      Radiology.MassLesion.MassShapeTypeID = Lesion.MassLesionShapeType.MassShapeTypeID AND 
                      Radiology.MassLesion.MassShapeTypeID = Lesion.MassLesionShapeType.MassShapeTypeID AND 
                      Radiology.MassLesion.MassShapeTypeID = Lesion.MassLesionShapeType.MassShapeTypeID AND 
                      Radiology.MassLesion.MassShapeTypeID = Lesion.MassLesionShapeType.MassShapeTypeID AND 
                      Radiology.MassLesion.MassShapeTypeID = Lesion.MassLesionShapeType.MassShapeTypeID AND 
                      Radiology.MassLesion.MassShapeTypeID = Lesion.MassLesionShapeType.MassShapeTypeID AND 
                      Radiology.MassLesion.MassShapeTypeID = Lesion.MassLesionShapeType.MassShapeTypeID AND 
                      Radiology.MassLesion.MassShapeTypeID = Lesion.MassLesionShapeType.MassShapeTypeID AND 
                      Radiology.MassLesion.MassShapeTypeID = Lesion.MassLesionShapeType.MassShapeTypeID
                WHERE Radiology.RadiologyReport.DicomExamNumber=%d    
			''', int(ExamID))
		lesions_id = []
		for row in cur:
			print  row # "\tStudyNumber %d\t DicomExamNumber %d\t Description: %s\n AdditonalDescription: %s\n BreastSide: %s\t ClockStart: %s\t ClockEnd: %s\t MRISeries %s\t\n LesionCenterPixelX: %s\t LesionCenterPixelY: %s\t LesionFindingType.Description: %s\t\n\n" % (row[0], row[1], row[3], row[4], row[5], row[7], row[8], row[9], row[10], row[11], row[12])
			lesions_id.append( row[18] )	
			print "Lesions found:"
			print lesions_id
		conn.close()
			
		# *****************************************************************************
		# get_series and print which series to load
		print " "
		chdirname='Z:/Cristina/MassNonmass/'
		os.chdir(chdirname)  
		print os.getcwd()
		
		[abspath_ExamID, eID, SeriesID, studyFolder, dicomInfo] = get_series(StudyID, img_folder)

		## exit ChoseTodo								
								
finally:
	file_ids.close()		
		
			
