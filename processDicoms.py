#!/usr/bin/env python

import os
import os.path
import sys
from sys import argv, stderr, exit
import shlex, subprocess
import numpy
from scipy.io import loadmat, savemat
import dicom
import vtk
import Tix
import _mssql
import pymssql
from datetime import date
from operator import itemgetter, attrgetter
from dicom.contrib.pydicom_PIL import show_PIL

import gdcm
import vtkgdcm
import vtkgdcmPython
from matplotlib.pylab import *

# ******************************************************************************
print "usage: processDicoms.py patient_list.txt"
print '''
Run on mass/nonmass/foci root folder processDicoms.py based on list of 
patient_list.txt(list of StudyIDs in folder mass or nonmass)

"usage:" 
python ../codeProject/processDicoms.py ../bUniqueTotal_formSharmila_june21_MASS.txt
-- This script will 

% Copyright (C) Cristina Gallego, University of Toronto, 2012
% June 12/12 - 	Added dicom_dict python comparison
% June 27/12 - 	Now works with more added columns afer ExamsNo on TotalFromSharmila.txt such as MRN
		SQL QUERY with PYTHON to identify more information about the lesion present	
% July 12 -     Now works with 4D image sequence viewer (extracting series info/phase/time from dicom) 
		Sorts and generates folders/links of all locations of the T13D Gad Dynamic sequences by:
		1) 5time-points per Location folders
		2) L/R breast Volumes and

-----------------------------------------
'''
def get_gdcm_to_numpy_typemap():
    """Returns the GDCM Pixel Format to numpy array type mapping."""
    _gdcm_np = {gdcm.PixelFormat.UINT8  :numpy.int8,
                gdcm.PixelFormat.INT8   :numpy.uint8,
                gdcm.PixelFormat.UINT16 :numpy.uint16,
                gdcm.PixelFormat.INT16  :numpy.int16,
                gdcm.PixelFormat.UINT32 :numpy.uint32,
                gdcm.PixelFormat.INT32  :numpy.int32,
                gdcm.PixelFormat.FLOAT32:numpy.float32,
                gdcm.PixelFormat.FLOAT64:numpy.float64 }
    return _gdcm_np

def get_numpy_array_type(gdcm_pixel_format):
    """Returns a numpy array typecode given a GDCM Pixel Format."""
    return get_gdcm_to_numpy_typemap()[gdcm_pixel_format]

def gdcm_to_numpy(image):
    """Converts a GDCM image to a numpy array.
    """
    pf = image.GetPixelFormat().GetScalarType()
    print 'pf', pf
    print image.GetPixelFormat().GetScalarTypeAsString()
    assert pf in get_gdcm_to_numpy_typemap().keys(), \
           "Unsupported array type %s"%pf
    d = image.GetDimension(0), image.GetDimension(1)
    print 'Image Size: %d x %d' % (d[0], d[1])
    dtype = get_numpy_array_type(pf)
    gdcm_array = image.GetBuffer()
    result = numpy.frombuffer(gdcm_array, dtype=dtype)
    maxV = float(result[result.argmax()])
    ## linear gamma adjust
    #result = result + .5*(maxV-result)
    ## log gamma
    result = numpy.log(result+50) ## 50 is apprx background level
    maxV = float(result[result.argmax()])
    result = result*(2.**8/maxV) ## histogram stretch
    result.shape = d
    return result
    
def readlink(path): 
    return path 

def get_immediate_subdirectories(dir):
    return [name for name in os.listdir(dir) 
            if os.path.isdir(os.path.join(dir, name))]

def get_only_linksindirectory(dir):
    return [name for name in os.listdir(dir) 
            if os.path.islink(os.path.join(dir, name))]
 
def get_only_filesindirectory(dir):
    return [name for name in os.listdir(dir) 
            if os.path.isfile(os.path.join(dir, name))]
            
def find(strng, ch):
    index = 0
    while index < len(strng):
        if strng[index] == ch:
            return index
        index += 1                  
    return -1

def get_display(abspath_SeriesID):
	arranged_folders = get_immediate_subdirectories(abspath_SeriesID);
	
	# Initialize series count
	s=0;
	print " "
	print "%s \t %s \t\t %s " % ('n', 'Series#', '#Images')
	
	# Find all dicoms Series in sudyfolder (will go into including subfolders)
	# Iterate
	for arrangedF in arranged_folders:
		#print "arrangedFolder: %s" % arrangedF
		path_arrangedF_ID = abspath_SeriesID+'/'+arrangedF
		#print path_arrangedF_ID
		
		# Enter Studyfolder to process all Dicom Series of ExamID
		path_arrangedF_images = get_only_filesindirectory(path_arrangedF_ID);
		#print path_arrangedF_images

		print "%s \t %s \t\t %s " % (s, arrangedF, len(path_arrangedF_images))
		
		# increment series number
		s=s+1;	
		
	# Go back to rootfolder
	os.chdir(path_rootFolder)
    
    	return arranged_folders		
	
def get_series(StudyID):
	# obtain subdires in the StudyID directory, correspond to ExamsID 
	# check for one in subdirectory tree, (e.g Mass or NonMass) 
	path_studyID = '../mass/'+StudyID
	ExamsID = get_immediate_subdirectories(path_studyID);
	#print ExamsID
	
	studyFolder = os.path.abspath(path_studyID)
	#print studyFolder
	
	# Find all dicoms Series in sudyfolder (will go into including subfolders)
	# Iterate
	for eID in ExamsID:
		print "ExamID: %s" % eID
		path_ExamID = '../mass/'+StudyID+'/'+eID
		abspath_ExamID = os.path.abspath(path_ExamID)
		#print abspath_ExamID
		
		# Enter Studyfolder to process all Dicom Series of ExamID
		SeriesID = get_immediate_subdirectories(path_ExamID);
		#print SeriesID
		
		# Initialize series count
		s=0;
		print "Total number of series: %d" % len(SeriesID)
		print " "
		print "%s     %s		%s 	%s " % ('n', 'Series#', '#Images', 'SeriesDescription')
				
		# Iterate for each series in ExamID
		for sID in SeriesID:

			path_SeriesID = '../mass/'+StudyID+'/'+eID+'/'+sID
			#print path_SeriesID
							
			abspath_SeriesID = os.path.abspath(path_SeriesID)

			# Get total number of files
			listSeries_files = get_only_filesindirectory(abspath_SeriesID)
						
			# Use only the first slice one file and get DICOM DICTIONARY
			path_filenameID = '../mass/'+StudyID+'/'+eID+'/'+sID+'/'+listSeries_files[0]
						
			dicomInfo = dicom.read_file(os.path.abspath(path_filenameID)) 
					
			# Get the main dataset (they are in fact two separate datasets in the DICOM standard). 
			# That dicom dataset is now stored in the file_meta attribute of the dataset
			PatientID = dicomInfo.PatientID
			SeriesNumber = dicomInfo.SeriesNumber
									
			# Get structure of study (all files in directory consistent with studyID and patientID
			studyTree = []
			SeriesNumber = dicomInfo.SeriesNumber
			FileNames=listSeries_files;
			SeriesDescription = dicomInfo.SeriesDescription;
			MinSliceLocation = dicomInfo.SliceLocation; # Infos to identify number of slices/volumes
			ImageOrientationPatient = dicomInfo.ImageOrientationPatient;
			NumberOfVolumes = 1 # default
			
			# iterate number of slices and get full volume (#slices, each slice loc)
			slices = []
			num_images=0;
			for filename in listSeries_files:
				num_images = num_images+1
					
			# Print series info						
			print "%d	%s		%d		%s" % (s, dicomInfo.SeriesNumber, num_images, dicomInfo.SeriesDescription) 
			# increment series number
			s=s+1;	
			
			# Go back to rootfolder
			os.chdir(path_rootFolder)	
			
	return abspath_ExamID, eID, SeriesID
		
			
def get_slices_at_all_locs(abspath_SeriesID, len_listSeries_files, StudyID, eID, SeriesID, choseSerie, listSeries_files, num_series, abspath_ExamID ):
	# enter folder of Series selection /examID/S_XXX 
	print "entering folder abspath_SeriesID"
	os.chdir(abspath_SeriesID)
	print abspath_SeriesID
				
	# Obtain all datasets in current series
	# iterate number of slices and get full volume (#slices, each slice loc)
	slices = []
	FileNms_slices =  []
				
	for n in range(len_listSeries_files):
		# Use all DICOM slices on series
		''' EXTRACT DICOM SLICE LOCATION '''
		absp_fsID = 'Z:/Cristina/MassNonMass/mass/'+StudyID+'/'+eID+'/'+SeriesID[int(choseSerie)]+'/'+listSeries_files[n]
							
		dInfo = dicom.read_file(absp_fsID)
		slices.append(dInfo.SliceLocation)
		FileNms_slices.append(listSeries_files[n])
		FileNms_slices.append(dInfo.SliceLocation)
		
	print '''\nPROCESS STACKS BY SLICE LOCATIONS '''
	FileNms_slices_stack = reshape(FileNms_slices, [len_listSeries_files,2])
	FileNms_slices_sorted = sorted(FileNms_slices_stack, key=itemgetter(1))
	FileNms_slices_sorted_stack = reshape(FileNms_slices_sorted, [len_listSeries_files,2])
	
	current_slice = FileNms_slices_sorted_stack[0,1]
	print current_slice
	stack_byLocation = []
	name_byLocation = []
	scount = 0
	
	for sliceS in FileNms_slices_sorted_stack:
		# Get the num_series = 5 name_byLocations for a given Location
		if( scount < num_series):
			print "Name: %s" % sliceS[0]
			#print "Slice_loc: %s" % sliceS[1]
			stack_byLocation.append(sliceS[1])
			name_byLocation.append(sliceS[0])
			scount = scount+1
			
		# Finish getting all series for a given Location
		else:
			'''-----\t NOW HAVE ALL SLICES IN A stack_byLocation '''
			# Get the new Location folder
			# To verify "Above lengths should be = num_series + 1" len(stack_byLocation) len(name_byLocation)'''
								
			current_loc = stack_byLocation[0]
									
			# Makedir of current loc
			# if current loc folder doesn't exist create it		
			if not os.path.exists(str(current_loc)):
				os.makedirs(str(current_loc))
			
			# Get inside location directory
			os.chdir(str(current_loc))
			#print os.getcwd()
			
			# Now link slices at location to folder
			filename = str(name_byLocation[0])
			filename = filename[0:-10]
			file_ending = '.MR.dcm'
			
			# Save the file list to read as series later
			filename_series = 'DIRCONTENTS.txt'
			file_series = open(filename_series, 'w')
			
			for j in range(num_series):
				link_to = '../'+name_byLocation[j]
				name4link_to = filename+'00'+str(j)+file_ending
				print "linking file: %s to: %s" % (link_to, name4link_to)
				ln_subp = subprocess.Popen(['ln', link_to, name4link_to], stdout=subprocess.PIPE)
				ln_subp.wait()
				file_series.write(str(name4link_to)+'\n')
				
			file_series.close()
				
			# Get back inside the  Series directory
			os.chdir(abspath_SeriesID)
			
			print '''\n-----\tchdir out GET NEXT LOCATIONS SLICES'''
			current_slice = sliceS[1]
			scount = 0
			stack_byLocation = []
			name_byLocation = []
			print current_slice
			
			print "Name: %s" % sliceS[0]
			#print "Slice_loc: %s" % sliceS[1]
			stack_byLocation.append(sliceS[1])
			name_byLocation.append(sliceS[0])
			scount = scount+1
								
	print "\n------------------------------------------------------------------"						
	print '''FINISH get_slices_at_all_locs \n'''	
	os.chdir(abspath_ExamID)
	return

def get_slices_for_volumes(abspath_SeriesID, len_listSeries_files, StudyID, eID, SeriesID, choseSerie, listSeries_files, num_series, abspath_ExamID ):
	print "entering folder abspath_SeriesID"
	os.chdir(abspath_SeriesID)
	print abspath_SeriesID
				
	# Obtain all datasets in current series
	# iterate number of slices and get full volume (#slices, each slice loc)
	slices = []
	FileNms_slices =  []
	
	for n in range(len_listSeries_files):
		# Use all DICOM slices on series
		''' EXTRACT DICOM SLICE LOCATION '''
		absp_fsID = 'Z:/Cristina/MassNonMass/mass/'+StudyID+'/'+eID+'/'+SeriesID[int(choseSerie)]+'/'+listSeries_files[n]
							
		dInfo = dicom.read_file(absp_fsID)
		slices.append(dInfo.SliceLocation)
		FileNms_slices.append(listSeries_files[n])
		FileNms_slices.append(dInfo.SliceLocation)
		
	print '''\nPROCESS STACKS BY SLICE LOCATIONS '''
	FileNms_slices_stack = reshape(FileNms_slices, [len_listSeries_files,2])
	FileNms_slices_sorted = sorted(FileNms_slices_stack, key=itemgetter(0))
	FileNms_slices_sorted_stack = reshape(FileNms_slices_sorted, [len_listSeries_files,2])
	
	stack_byLocation = []
	name_byLocation = []
	scount = 0
	
	# Extracting the number of slices per bilateral 3D volumes based on num_series for 280/5 = 56
	numSlices_perBilat_vol = len_listSeries_files/num_series
	print "------\tNumber of Slices per Bilateral Volume: %d " % numSlices_perBilat_vol
				
	# Get the folder names based on num_series
	stack_byBilatVol = [] 
	stack_byBilatVol.append('pre-Contrast')
	
	for k in range(1,num_series):
		# Now link slices at location to folder
		stack_byBilatVol.append('post_Contrast-'+str(k))
		
	# Initialized BilatVol
	print stack_byBilatVol
	BilatVol_folder = 0
		
	for sliceS in FileNms_slices_sorted_stack:
		# Get the num_series = 5 name_byLocations for a given Location
		if( scount < numSlices_perBilat_vol):
			print "--\tName: %s" % sliceS[0]
			print "--\tSLoc: %s" % sliceS[1]
			stack_byLocation.append(sliceS[1])
			name_byLocation.append(sliceS[0])
			scount = scount+1
			
		# Finish getting all series for a given Location
		else:
			'''-----\t NOW HAVE ALL SLICES IN A bilateral vol '''
			# Get the new Location folder
			print "bit vol %d" % BilatVol_folder
			current_vol = stack_byBilatVol[BilatVol_folder]
									
			# Makedir of current loc
			# if current loc folder doesn't exist create it		
			if not os.path.exists(str(current_vol)):
				os.makedirs(str(current_vol))
			
			# Get inside location directory
			os.chdir(str(current_vol))
			print os.getcwd()
			
			# Now link slices at location to folder
			filename = str(name_byLocation[0])
			filename = filename[0:-10]
			file_ending = '.MR.dcm'
			
			# Save the file list to read as series later
			filename_series = 'DIRCONTENTS.txt'
			file_series = open(filename_series, 'w')
			
			for j in range(numSlices_perBilat_vol):
				link_to = '../'+name_byLocation[j]
				name4link_to = filename+'00'+str(j+1)+file_ending
				print "linking file: %s to: %s" % (link_to, name4link_to)
				ln_subp = subprocess.Popen(['ln', link_to, name4link_to], stdout=subprocess.PIPE)
				ln_subp.wait()
				file_series.write(str(name4link_to)+'\n')
				
			file_series.close()
				
			# Get back inside the  Series directory
			os.chdir(abspath_SeriesID)
			
			print '''\n-----\tchdir out GET NEXT LOCATIONS SLICES'''
			scount = 0
			stack_byLocation = []
			name_byLocation = []
			
			print "--\tName: %s" % sliceS[0]
			print "--\tSLoc: %s" % sliceS[1]
			stack_byLocation.append(sliceS[1])
			name_byLocation.append(sliceS[0])
			scount = scount+1
			BilatVol_folder = BilatVol_folder+1
			
	'''-----\tFINISH THE LAST bilateral vol '''
	# Get the new Location folder
	print "bit vol %d" % BilatVol_folder
	current_vol = stack_byBilatVol[BilatVol_folder]
							
	# Makedir of current loc
	# if current loc folder doesn't exist create it		
	if not os.path.exists(str(current_vol)):
		os.makedirs(str(current_vol))
	
	# Get inside location directory
	os.chdir(str(current_vol))
	print os.getcwd()
	
	# Now link slices at location to folder
	filename = str(name_byLocation[0])
	filename = filename[0:-10]
	file_ending = '.MR.dcm'
		
	# Save the file list to read as series later
	filename_series = 'DIRCONTENTS.txt'
	file_series = open(filename_series, 'w')

	for j in range(numSlices_perBilat_vol):
		link_to = '../'+name_byLocation[j]
		name4link_to = filename+'00'+str(j+1)+file_ending
		print "linking file: %s to: %s" % (link_to, name4link_to)
		ln_subp = subprocess.Popen(['ln', link_to, name4link_to], stdout=subprocess.PIPE)
		ln_subp.wait()
		file_series.write(str(name4link_to)+'\n')
				
	file_series.close()
		
	# Get back inside the  Series directory
	os.chdir(abspath_SeriesID)
	
	print "\n------------------------------------------------------------------"						
	print '''FINISH get_slices_for_volumes \n'''	
	os.chdir(abspath_ExamID)
	
	return
	
				

# ******************************************************************************
epsilon = .0001         

# use cell picker for interacting with the image orthogonal views.
picker = vtk.vtkCellPicker()
picker.SetTolerance(0.005)

# Create 3 orthogonal view using the ImagePlaneWidget
xImagePlaneWidget = vtk.vtkImagePlaneWidget()
yImagePlaneWidget = vtk.vtkImagePlaneWidget()
zImagePlaneWidget = vtk.vtkImagePlaneWidget()

#  The 3 image plane widgets
xImagePlaneWidget.DisplayTextOn();
xImagePlaneWidget.SetPicker(picker);
xImagePlaneWidget.RestrictPlaneToVolumeOn();
xImagePlaneWidget.SetKeyPressActivationValue('x');
xImagePlaneWidget.GetPlaneProperty().SetColor(1, 0, 0);
xImagePlaneWidget.SetResliceInterpolateToNearestNeighbour();

yImagePlaneWidget.DisplayTextOn();
yImagePlaneWidget.SetPicker(picker);
yImagePlaneWidget.RestrictPlaneToVolumeOn();
yImagePlaneWidget.SetKeyPressActivationValue('y');
yImagePlaneWidget.GetPlaneProperty().SetColor(0, 1, 0);
yImagePlaneWidget.SetLookupTable(xImagePlaneWidget.GetLookupTable());

zImagePlaneWidget.DisplayTextOn();
zImagePlaneWidget.SetPicker(picker);
zImagePlaneWidget.SetKeyPressActivationValue('z');
zImagePlaneWidget.GetPlaneProperty().SetColor(0, 0, 1);
zImagePlaneWidget.SetLookupTable(xImagePlaneWidget.GetLookupTable());
zImagePlaneWidget.SetRightButtonAutoModifier(1);

# Create a renderer, render window, and render window interactor to
# display the results.
renderer1 = vtk.vtkRenderer()
renWin1 = vtk.vtkRenderWindow()
iren1 = vtk.vtkRenderWindowInteractor()

renWin1.SetSize(1000, 800);
renWin1.AddRenderer(renderer1);
iren1.SetRenderWindow(renWin1);

xImagePlaneWidget.SetInteractor( iren1 )
xImagePlaneWidget.SetSliceIndex(20)
xImagePlaneWidget.On()

yImagePlaneWidget.SetInteractor( iren1 )
yImagePlaneWidget.SetSliceIndex(20)
yImagePlaneWidget.On()

zImagePlaneWidget.SetInteractor( iren1 )
zImagePlaneWidget.SetSliceIndex(10)
zImagePlaneWidget.On()

# *****************************************************************************
# Get Root folder ( the directory of the script being run)
path_rootFolder = os.path.dirname(os.path.abspath(__file__))
#print path_rootFolder

# init vars
studyTree = []

# Open filename list
file_ids = open(sys.argv[1],"r")

try:
	for line in file_ids:
		# Get the line: Study#, DicomExam# 
		line = line.split()	
		StudyID = line[0]
		ExamID = line[1]
		print "Getting studyID# DicomExam# "
		print StudyID, ExamID
		
		# Execute mysqlconnect
		print "Executing SQL connection..."
		conn = pymssql.connect(host='SG12-BIOMATRIX', user='cristina.gallego', password='miag_1234', database='SMIAL_BCAD')
		cur = conn.cursor()
		cur.execute('''
			SELECT  Patient.Patients.StudyNumber, Radiology.RadiologyReport.DicomExamNumber, Radiology.RadiologyReport.DateOfExam, Lesion.GeoDescription.Description AS Expr1, Lesion.FindingLocation.AdditionDescription, Lesion.GeoDescription.BreastSide, Lesion.FindingLocation.ClockStart, 
                      Lesion.FindingLocation.ClockEnd, Lesion.FindingLocation.MRISeries, Lesion.FindingLocation.LesionCenterPixelX, Lesion.FindingLocation.LesionCenterPixelY,Lesion.LesionFindingType.Description, 
                      Lesion.LesionFindingType.LesionFindingTypeID,  Lesion.GeoDescription.CurveTypeID, Lesion.GeoDescription.T2SignalTypeID 
			FROM Radiology.RadiologyReport INNER JOIN
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
                      Visit.VisitReasonType ON Visit.ImagingVisits.VisitReasonTypeID = Visit.VisitReasonType.VisitReasonTypeID 
			WHERE DicomExamNumber=%d  
			''', int(ExamID))
		for row in cur:
			print  "\tStudyNumber %d\t DicomExamNumber %d\t Description: %s\n AdditonalDescription: %s\n BreastSide: %s\t ClockStart: %s\t ClockEnd: %s\t MRISeries %s\t\n LesionCenterPixelX: %s\t LesionCenterPixelY: %s\t LesionFindingType.Description: %s\t\n\n" % (row[0], row[1], row[3], row[4], row[5], row[7], row[8], row[9], row[10], row[11], row[12])
		
		conn.close()
		
		# *****************************************************************************
		# get_series and print which series to load
		print "first time"
		print os.getcwd()
		
		[abspath_ExamID, eID, SeriesID] = get_series(StudyID)
		
		# *****************************************************************************
		# ask for which series to load
		print "\n----------------------------------------------------------"
		choseSerie = raw_input('Enter n Series to load (0-n), or x to exit: ')
		
		while(choseSerie!='x'):
			# Enter folder
			os.chdir(abspath_ExamID)
			
			path_SeriesID = 'mass/'+StudyID+'/'+eID+'/'+SeriesID[int(choseSerie)]
			print path_SeriesID
							
			# abspath_SeriesID = os.path.abspath(path_SeriesID) -- doesn't work bug
			abspath_SeriesID = 'Z:/Cristina/MassNonMass/'+path_SeriesID
			
			# Get total number of files
			listSeries_files = sorted(list(get_only_filesindirectory(abspath_SeriesID)))
			
			print "Total images in series: %d " % len(listSeries_files)
			len_listSeries_files = len(listSeries_files)
			
			# Use only the first slice one file and get DICOM DICTIONARY
			path_filenameID = 'mass/'+StudyID+'/'+eID+'/'+SeriesID[int(choseSerie)]+'/'+listSeries_files[0]
						
			# os.path.abspath(path_filenameID) not working
			abspath_path_filenameID = 'Z:/Cristina/MassNonMass/'+path_filenameID
											
			''' EXTRACT DICOM HEADER '''
			dicomInfo = dicom.read_file(abspath_path_filenameID) 
			
			print("Filename.........:", listSeries_files[0])
			print dicomInfo[0x0008,0x0020]	#print("StudyDate...:", dicomInfo.StudyDate)
			print dicomInfo[0x0010,0x0020]	#print("Patient id (MRN).:", dicomInfo.PatientID)
			print dicomInfo[0x0012,0x0040]	#print("ClinicalTrialID.......:", dicomInfo.ClinicalTrialID)
			print dicomInfo[0x0020,0x0010]	# StudyID				
			print dicomInfo[0x0010,0x1010]	#print("PatientAge.:", dicomInfo.PatientAge)
			
			print dicomInfo[0x0019,0x109e] 	# Internal Pulse sequence Name
			print dicomInfo[0x0018,0x5100]	# print("Patient Position.......:", dicomInfo.PatientPosition)
			print dicomInfo[0x0020,0x0032]	# Image Position (patient)
			print dicomInfo[0x0020,0x0037]	# Image Orientation (patient)
			print dicomInfo[0x0019,0x1019]	# First Scan location
			print dicomInfo[0x0019,0x109c] 	# Pulse sequence Nmae
			print dicomInfo[0x0019,0x10a4] 	# SAT Fat/water/bone
			
			
			print dicomInfo[0x0008,0x0060]	#print("Modality.........:", dicomInfo.Modality)
			print dicomInfo[0x0008,0x0030]	#print("Study Time.......:", dicomInfo.StudyTime)
			print dicomInfo[0x0008,0x0031]	#print("Series Time.......:", dicomInfo.SeriesTime)
			print dicomInfo[0x0008,0x0033]	# Imagetime
			print dicomInfo[0x0020,0x0011] 	# SeriesNumber
			print dicomInfo[0x0008,0x103e]	#print("SeriesDescription.......:", dicomInfo.SeriesDescription)
			print dicomInfo[0x0020,0x000E]	#Series Instance UID
							
			#print dicomInfo[0x0008,0x1140]	# Reference Image sequence
			#print dicomInfo[0x0020,0x1002]	#Images in Acquisition
			#print dicomInfo[0x0020,0x9056]	# In-Stack Position Number. Filled in for some applications only. Slice number within the stack that this image belongs to.
			#print dicomInfo[0x0020,0x9057]	# Stack ID. Filled in for some applications only. Number (starting at 1) of the graphic prescription slice group the image belongs to.
			
			print dicomInfo[0x0018,0x1050]	#print("SpatialResolution.......:", dicomInfo.SpatialResolution)
			print dicomInfo[0x0020,0x0110]	#print("TemporalResolution.......:", dicomInfo.TemporalResolution)
			#print dicomInfo[0x0020,0x0100]	# Temporal Position Identifier
			#print dicomInfo[0x0020,0x0105]	# Number of Temporal Positions
			print dicomInfo[0x0018,0x0080]	# TR #print("Repetition Time.......:", dicomInfo.AcquisitionTime)
			print dicomInfo[0x0018,0x0081]	# TE #print("Echo Time.......:", dicomInfo.AcquisitionTime)
			print dicomInfo[0x0018,0x1314]	# Flip Angle
			#print dicomInfo[0x0018,0x0082]	# Inversion Time
			print dicomInfo[0x0018,0x0091]	# Echo Train Length
			
			print dicomInfo[0x0018,0x0050] 	#print("SliceThickness.......:", dicomInfo.SliceThickness)
			print dicomInfo[0x0018,0x0088] 	#print("SpacingBetweenSlices.......:", dicomInfo.SpacingBetweenSlices)
			print dicomInfo[0x7fe0,0x0010]	#array
			print dicomInfo[0x0019,0x1017]
			
			#///////////////////////////////////////////////////////////////////////////////////////////#
			## Get number of series, + 1 for precontrast
			num_series = dicomInfo[0x0019,0x1017].value  + 1
			print "\n------------------------------------------------------------------\n"
			print "Number of series: %d" % num_series
			
			# See what kind of Series are dealing with
			# TO DO if(len_listSeries_files>50):
			
			print "choseto firs time"
			print os.getcwd()
			
			choseTodo = raw_input('Enter option to load: \n\t1) get_slices_at_all_locs or \n\t2) get_slices_for_volumes or \n\tx) to exit \n\t input: ')
			while(choseTodo!='x'):
				if(choseTodo == '1'):
					print "\n------------------------------------------------------------------"						
					print '''GETTING SLICE AT ALL LOCS \n'''	
					get_slices_at_all_locs(abspath_SeriesID, len_listSeries_files, StudyID, eID, SeriesID, choseSerie, listSeries_files, num_series, abspath_ExamID )
				if(choseTodo == '2'):	
					print "\n------------------------------------------------------------------"						
					print '''NOW GET SLICES FOR VOLUMES \n'''	
					get_slices_for_volumes(abspath_SeriesID, len_listSeries_files, StudyID, eID, SeriesID, choseSerie, listSeries_files, num_series, abspath_ExamID )
								
				choseTodo = raw_input('Enter option to load: \n\t1) get_slices_at_all_locs or \n\t2) get_slices_for_volumes or \n\tx) to exit \n\t input: ')
						
			# *****************************************************************************
			# ask for which series to display
			print "\n----------------------------------------------------------"
			os.chdir(abspath_SeriesID)
			arranged_folders = get_display(abspath_SeriesID)
					
			choseDisplay = raw_input('Enter n Folder to load (0-n), or x to exit: ')
						
			while(choseDisplay!='x'):
				print "\n------------------------------------------------------------------"
				print ''' LOAD IMAGES '''
				chosen_folder = arranged_folders[int(choseDisplay)]
				
				abspath_FolderID = abspath_SeriesID+'/'+chosen_folder
				os.chdir(abspath_FolderID)
				
				#proc = subprocess.Popen('ls', stdout=subprocess.PIPE)
				#listchosen_images = proc.stdout.read().split()
				#listchosen_images = list(listchosen_images)
				#print listchosen_images
				
				'''-----------------------------'''
				vtkStringArray = vtk.vtkStringArray()
				filename_series = 'DIRCONTENTS.txt'
				file_series = open(filename_series,'r')
				files_toRead = file_series.read()
				i = 0
				try:
					for slicename in files_toRead.split():
						#print "adding : %s" %  slicename
						vtkStringArray.InsertNextValue( slicename )
				finally:
					file_series.close()
					
				# go back
				#vtkStringArray = vtk.vtkStringArray()
				#for slicename in listchosen_images:
				#	print "adding : %s" %  slicename
				#	vtkStringArray.InsertNextValue( slicename )
									
				#ds = dicom.read_file(listchosen_images[0:19])
				#show_PIL(ds)
				
				#r  = gdcm.ImageReader()
				#r.SetFileName( listchosen_images[0:19])
				''' '''
				
				''' '''
				dicomReader  = vtkgdcmPython.vtkGDCMImageReader()
				dicomReader.SetFileNames( vtkStringArray )
				dicomReader.Update()
				
				# Get image from reader
				im = dicomReader.GetOutput()
				im_scalars = im.GetPointData().GetScalars()
				dims = im.GetDimensions()
					
				# check wether DicomTag corresponds to Image loaded in vtk
				print dims
				spacing = im.GetSpacing()
				print spacing
					
				# Visualize results
				xImagePlaneWidget.SetInput( dicomReader.GetOutput() )
				xImagePlaneWidget.SetPlaneOrientationToXAxes();
				yImagePlaneWidget.SetInput( dicomReader.GetOutput() )
				yImagePlaneWidget.SetPlaneOrientationToYAxes();
				zImagePlaneWidget.SetInput( dicomReader.GetOutput() )
				zImagePlaneWidget.SetPlaneOrientationToZAxes();
				
				xImagePlaneWidget.SetInteractor( iren1 )
				xImagePlaneWidget.SetSliceIndex(int(dims[0]))
				xImagePlaneWidget.On()
				
				yImagePlaneWidget.SetInteractor( iren1 )
				yImagePlaneWidget.SetSliceIndex(int(dims[1]/2.0))
				yImagePlaneWidget.On()
		
				zImagePlaneWidget.SetInteractor( iren1 )
				zImagePlaneWidget.SetSliceIndex(int(dims[2]/2.0))
				zImagePlaneWidget.On()
				
				# Set Up Camera view
				renderer1.SetBackground(0.0, 0.0, 0.0)
				camera = renderer1.GetActiveCamera()
		
				#bounds and initialize camera
				b = im.GetBounds()
				renderer1.ResetCamera(b)	
				renderer1.ResetCameraClippingRange()
				camera.SetViewUp(0.0,1.0,0.0)
				camera.Azimuth(315)
					
				# Initizalize
				renWin1.Render()
				renderer1.Render()
				iren1.Initialize()
				iren1.Start()
				#
				#
				choseDisplay = raw_input('Enter n Series to load (0-n), or x to exit: ')
	
			print "\nFINISH DISPLAYING----------------------------------------------------------------"
			print '''Ask for another series or continue \n'''
			os.chdir("Z:\Cristina\MassNonmass\mass")
			
			print "Current folder time"
			print os.getcwd()

			# Save the file list to read as series later
			summary_cases = 'SUMMARY_CASES.txt'
			file_summary_cases = open(summary_cases, 'a')
			
			file_summary_cases.write(str(name4link_to)+'\n')
			file_summary_cases.close()
			
			[abspath_ExamID, eID, SeriesID] = get_series(StudyID)
			choseSerie = raw_input('Enter n Series to load (0-n), or x to exit: ')
	
						
						
finally:
	file_ids.close()		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
