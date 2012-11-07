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
print "usage: processDicoms.py patient_list.txt"
print '''
Run on mass/nonmass/foci root folder processDicoms.py based on list of 
patient_list.txt(list of StudyIDs in folder mass or nonmass)

"usage:" 
python2.7  codeProject/processDicoms.py img_folder batchs/bUniqueTotal_formSharmila_june21_MASS.txt
-- This script will 

% Copyright (C) Cristina Gallego, University of Toronto, 2012
% June 12/12 - 	Added dicom_dict python comparison
% June 27/12 - 	Now works with more added columns afer ExamsNo on TotalFromSharmila.txt such as MRN
		SQL QUERY with PYTHON to identify more information about the lesion present	
% July 12 -     Now works with 4D image sequence viewer (extracting series info/phase/time from dicom) 
		Sorts and generates folders/links of all locations of the T13D Gad Dynamic sequences by:
		1) nth-time-points per Location folders
		2) bilateral post-contrast Volumes when phases are not separated. (i.e when Ph1/Ph2/Ph3.. is unavailable)
% Aug 22 -      Added function get_locs_for_all_phases() 
		3) Loopps through each one of the Ph1/Ph2/Ph3/Ph4 series and finds equal locations at 4 post-contrast timepoints in the T13D Gad Dynamic sequences by
% Aug 24 -      Initial lesion detection and lesion center identification:  
		- Translated DICOM patientcoordinate system to VTK world coordinate system to match lesion location (from SMIAL database) More info: \\labshare\amartel_data\Cristina\MassNonmass\codeProject\worldCoordinates.txt
		- Initial lesion detection and lesion center identification:  Based on vtkImageThresholdConnectivity to select VOI
		- Create a Signal enhancement Map on Volumetric region of interest (VOI), saved as a vtkPolyData and xml (compatible with BMRIVIewer)
		- Extract lesion parameters of interest from VOI. More info refer to white_paper.doc
		 
		Outputs: Saved inside dynamic series/in ExamID folder as VOIlesions_idXXX (where idXXX corresponds to table field Radiology.MassLesion.LesionFindingID)
		Reports: Compilation of lesion extracted Parameters will be made available on \\labshare\amartel_data\Cristina\MassNonmass\reports\ALLCASES_lesion_params.xml
-----------------------------------------------------------------------------------------------------
'''
def FileCheck(filename):       
       try:
           fn=open(filename,"r") 
           fn.close()
           return True
       except IOError: 
           print "Error: DICOMDIR.txt doesn't exit, loading Series using vtkDICOMImageReader by setDirectoryName."
       return False
       
def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False
   
def readlink(path): 
    return path 

def get_immediate_subdirectories(mydir):
    return [name for name in os.listdir(mydir) 
            if os.path.isdir(os.path.join(mydir, name))]

def get_only_linksindirectory(mydir):
    return [name for name in os.listdir(mydir) 
            if os.path.islink(os.path.join(mydir, name))]
 
def get_only_filesindirectory(mydir):
     return [name for name in os.listdir(mydir) 
            if os.path.isfile(os.path.join(mydir, name))]
            
def find(strng, ch):
    index = 0
    while index < len(strng):
        if strng[index] == ch:
            return index
        index += 1                  
    return -1

def get_display_series(abspath_SeriesID):
	arranged_folders = get_immediate_subdirectories(abspath_SeriesID);

	# Initialize series count
	s=0;
	print "Total number of series: %d" % len(arranged_folders)
	print " "
	print "%s     %s		%s 	%s " % ('n', 'Series#', '#Images', 'SeriesDescription')
			
	# Iterate for each series in ExamID
	for arrangedF in arranged_folders:
		path_arrangedF_ID = abspath_SeriesID+os.sep+arrangedF
		#print path_SeriesID
		
		# Get total number of files
		listSeries_files = get_only_filesindirectory(path_arrangedF_ID)
					
		# Use only the first slice one file and get DICOM DICTIONARY
		if(listSeries_files != []):
			path_filenameID = abspath_SeriesID+os.sep+arrangedF+os.sep+listSeries_files[0]
						
			dicomInfo = dicom.read_file(os.path.abspath(path_filenameID)) 
					
			# Get the main dataset (they are in fact two separate datasets in the DICOM standard). 
			# That dicom dataset is now stored in the file_meta attribute of the dataset
					
			# Get structure of study (all files in directory consistent with studyID and patientID
			studyTree = []
			FileNames=listSeries_files;
			
			if("PatientID" in dicomInfo):
				PatientID = dicomInfo.PatientID
			else:	PatientID=''
			if("SeriesNumber" in dicomInfo):
				SeriesNumber = dicomInfo.SeriesNumber
			else:	SeriesNumber=''
			if("SeriesDescription" in dicomInfo):
				SeriesDescription = dicomInfo.SeriesDescription;
			else:	SeriesDescription=''
			if("SliceLocation" in dicomInfo):
				MinSliceLocation = dicomInfo.SliceLocation; # Infos to identify number of slices/volumes
			else:	MinSliceLocation=''
			if('ImageOrientationPatient' in dicomInfo):
				ImageOrientationPatient = dicomInfo.ImageOrientationPatient; # Infos to identify number of slices/volumes
			else:	ImageOrientationPatient=''
			
			NumberOfVolumes = 1 # default
			
			# iterate number of slices and get full volume (#slices, each slice loc)
			slices = []
			num_images=0;
			for filename in listSeries_files:
				num_images = num_images+1
					
			# Print series info						
			print "%d	%s		%d		%s" % (s, SeriesNumber, num_images, SeriesDescription) 
			# increment series number
			s=s+1;	
		else:
			print "%d	%s		%d		%s" % (s, "NONE", 0, "NULL") 
			# increment series number
			s=s+1;	
		
	# Go back to rootfolder
	os.chdir(path_rootFolder)
    
    	return arranged_folders
    	
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
		path_arrangedF_ID = abspath_SeriesID+os.sep+arrangedF
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
	
def get_series(StudyID,img_folder):
	# obtain subdires in the StudyID directory, correspond to ExamsID 
	# check for one in subdirectory tree, (e.g Mass or NonMass) 
	global abspath_SeriesID
	
	path_studyID = img_folder+StudyID
	studyFolder = os.path.abspath(path_studyID)
	print studyFolder
	
	ExamsID = get_immediate_subdirectories(path_studyID);
	#print ExamsID
	c = 0
	if(len(ExamsID)>1):
		
		print "%s     %s	" % ('n', 'Series#')
		for iexam in ExamsID:
			print "%d	%s " % (c, str(iexam)) 
			c=c+1
					
		choseSerie = raw_input('Enter n Series to load (0-n), or x to exit: ')
		if(choseSerie != 'x'):
			c = 0
			for iexam in ExamsID:
				if(int(choseSerie) == c):
					eID = iexam
				c=c+1	
		else:
			return '', 0, '', studyFolder
			
		
		print "ExamID: %s" % eID
		path_ExamID = img_folder+StudyID+os.sep+eID
		abspath_ExamID = os.path.abspath(path_ExamID)
		print abspath_ExamID
		
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

			path_SeriesID = img_folder+StudyID+os.sep+eID+os.sep+sID
							
			abspath_SeriesID = os.path.abspath(path_SeriesID)
			
			# Get total number of files
			listSeries_files = get_only_filesindirectory(abspath_SeriesID)
						
			# Use only the first slice one file and get DICOM DICTIONARY
			path_filenameID = img_folder+StudyID+os.sep+eID+os.sep+sID+os.sep+listSeries_files[0]
			
			try:
				dicomInfo = dicom.read_file(os.path.abspath(path_filenameID))
			except ValueError:
				dicomInfo = []
				
			# Get the main dataset (they are in fact two separate datasets in the DICOM standard). 
			# That dicom dataset is now stored in the file_meta attribute of the dataset
											
			# Get structure of study (all files in directory consistent with studyID and patientID
			studyTree = []
			FileNames=listSeries_files;
			if("PatientID" in dicomInfo):
				PatientID = dicomInfo.PatientID#
			else:	PatientID=''
			if("SeriesNumber" in dicomInfo):
				SeriesNumber = dicomInfo.SeriesNumber#
			else:	SeriesNumber=''
			if("SeriesDescription" in dicomInfo):
				SeriesDescription = dicomInfo.SeriesDescription; #
			else:	SeriesDescription=''
			if("SliceLocation" in dicomInfo):
				MinSliceLocation = dicomInfo.SliceLocation; # Infos to identify number of slices/volumes
			else:	MinSliceLocation=''
			if('ImageOrientationPatient' in dicomInfo):
				ImageOrientationPatient = dicomInfo.ImageOrientationPatient; # Infos to identify number of slices/volumes
			else:	ImageOrientationPatient=''
			
			NumberOfVolumes = 1 # default
			
			# iterate number of slices and get full volume (#slices, each slice loc)
			slices = []
			num_images=0;
			for filename in listSeries_files:
				num_images = num_images+1
					
			# Print series info						
			print "%d	%s		%d		%s" % (s, SeriesNumber, num_images, SeriesDescription) 
			# increment series number
			s=s+1;	
			
			# Go back to rootfolder
			chdirname='Z:/Cristina/MassNonmass/'
			os.chdir(chdirname)  
			#print os.getcwd()
	else:	
		studyFolder = os.path.abspath(path_studyID)
		#print studyFolder
		
		# Find all dicoms Series in sudyfolder (will go into including subfolders)
		# Iterate
		for eID in ExamsID:
			print "ExamID: %s" % eID
			path_ExamID = img_folder+StudyID+os.sep+eID
			abspath_ExamID = os.path.abspath(path_ExamID)
			print abspath_ExamID
			
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
	
				path_SeriesID = img_folder+StudyID+os.sep+eID+os.sep+sID
				
				abspath_SeriesID = os.path.abspath(path_SeriesID)
				#print abspath_SeriesID
				
				# Get total number of files
				listSeries_files = get_only_filesindirectory(abspath_SeriesID)
			
				# Use only the first slice one file and get DICOM DICTIONARY
				if(listSeries_files != []):
					path_filenameID = img_folder+StudyID+os.sep+eID+os.sep+sID+os.sep+listSeries_files[0]
					try:
						dicomInfo = dicom.read_file(os.path.abspath(path_filenameID)) 
					except ValueError:
						dicomInfo = []
										
					# Get the main dataset (they are in fact two separate datasets in the DICOM standard). 
					# That dicom dataset is now stored in the file_meta attribute of the dataset
					PatientID = dicomInfo.PatientID
					SeriesNumber = dicomInfo.SeriesNumber
											
					# Get structure of study (all files in directory consistent with studyID and patientID
					studyTree = []
					FileNames=listSeries_files;
						
					if("SeriesNumber" in dicomInfo):
						SeriesNumber = dicomInfo.SeriesNumber
					else:	SeriesNumber=''
					if("SeriesDescription" in dicomInfo):
						SeriesDescription = dicomInfo.SeriesDescription;
					else:	SeriesDescription=''
					if("SliceLocation" in dicomInfo):
						MinSliceLocation = dicomInfo.SliceLocation; # Infos to identify number of slices/volumes
					else:	MinSliceLocation=''
					if('ImageOrientationPatient' in dicomInfo):
						ImageOrientationPatient = dicomInfo.ImageOrientationPatient; # Infos to identify number of slices/volumes
					else:	ImageOrientationPatient=''
					
					NumberOfVolumes = 1 # default
					
					# iterate number of slices and get full volume (#slices, each slice loc)
					slices = []
					num_images=0;
					for filename in listSeries_files:
						num_images = num_images+1
							
					# Print series info						
					print "%d	%s		%d		%s" % (s, SeriesNumber, num_images, SeriesDescription) 
					# increment series number
					s=s+1;	
				else:
					print "%d	%s		%d		%s" % (s, "NONE", 0, "NULL") 
					# increment series number
					s=s+1;	
				
		# Go back to rootfolder
		os.chdir(path_rootFolder)	
				
	return abspath_ExamID, eID, SeriesID, studyFolder, dicomInfo


def display(image, image_pos_pat, image_ori_pat, xImagePlaneWidget, yImagePlaneWidget, zImagePlaneWidget, renWin1, renderer1, iren1):
	global xMin, xMax, yMin, yMax, zMin, zMax, xSpacing, ySpacing, zSpacing
	# The box widget observes the events invoked by the render window
	# interactor.  These events come from user interaction in the render
	# window.
	# boxWidget = vtk.vtkBoxWidget()
	# boxWidget.SetInteractor(iren1)
	# boxWidget.SetPlaceFactor(1)
	
	# Initialize Image orienation
	IO = matrix(	[[0, 0,-1, 0],
			[1, 0, 0, 0],
			[0,-1, 0, 0],
			[0, 0, 0, 1]])
	# Assign the 6-Image orientation patient coordinates (from Dicomtags)
	IO[0,0] = image_ori_pat[0]; IO[1,0] = image_ori_pat[1]; IO[2,0] = image_ori_pat[2]; 
	IO[0,1] = image_ori_pat[3]; IO[1,1] = image_ori_pat[4]; IO[2,1] = image_ori_pat[5]; 
	
	# obtain thrid column as the cross product of column 1 y 2
	IO_col1 = [image_ori_pat[0], image_ori_pat[1], image_ori_pat[2]]
	IO_col2 = [image_ori_pat[3], image_ori_pat[4], image_ori_pat[5]]
	IO_col3 = cross(IO_col1, IO_col2)
	
	# assign column 3	
	IO[0,2] = IO_col3[0]; IO[1,2] = IO_col3[1]; IO[2,2] = IO_col3[2]; 
	
	IP =  [0, 0, 0, 1] # Initialization Image Position
	IP[0] = image_pos_pat[0]; IP[1] = image_pos_pat[1]; IP[2] = image_pos_pat[2];  
	
	# Calculate the center of the volume
	image.UpdateInformation()
	dims = image.GetDimensions()
	(xMin, xMax, yMin, yMax, zMin, zMax) = image.GetWholeExtent()
	(xSpacing, ySpacing, zSpacing) =image.GetSpacing()
	(x0, y0, z0) = image.GetOrigin()
	zcen = (zMax - zMin)/2
	ycen = (yMax - yMin)/2
	xcen = (xMax - xMin)/2
	
	print "Image center"
	print xcen, ycen, zcen
	print "Image Extension"
	print xMin, xMax, yMin, yMax, zMin, zMax

	# Visualize results
	xImagePlaneWidget.SetInput( image )
	xImagePlaneWidget.SetPlaneOrientationToXAxes();
	yImagePlaneWidget.SetInput( image )
	yImagePlaneWidget.SetPlaneOrientationToYAxes();
	zImagePlaneWidget.SetInput( image )
	zImagePlaneWidget.SetPlaneOrientationToZAxes();
	
	xImagePlaneWidget.SetInteractor( iren1 )
	xImagePlaneWidget.SetSliceIndex(int(dims[0]))
	xImagePlaneWidget.Off()
	
	yImagePlaneWidget.SetInteractor( iren1 )
	yImagePlaneWidget.SetSliceIndex(int(dims[1]/2.0))
	yImagePlaneWidget.Off()

	zImagePlaneWidget.SetInteractor( iren1 )
	zImagePlaneWidget.SetSliceIndex(int(dims[2]/2.0))
	zImagePlaneWidget.On()
	
	cube = vtk.vtkAnnotatedCubeActor()
	cube.SetXPlusFaceText( "R" );
	cube.SetXMinusFaceText( "L" );
	cube.SetYPlusFaceText( "A" );
	cube.SetYMinusFaceText( "P" );
	cube.SetZPlusFaceText( "H" );
	cube.SetZMinusFaceText( "F" );
	cube.SetFaceTextScale( 0.666667 );
		
	invert = vtk.vtkMatrix4x4();
	invert.SetElement(0, 0, IO[0,0])
	invert.SetElement(0, 1, IO[0,1])
	invert.SetElement(0, 2, IO[0,2])
	invert.SetElement(0, 3, IO[0,3])
	
	invert.SetElement(1, 0, IO[1,0])
	invert.SetElement(1, 1, IO[1,1])
	invert.SetElement(1, 2, IO[1,2])
	invert.SetElement(1, 3, IO[1,3])
	
	invert.SetElement(2, 0, IO[2,0])
	invert.SetElement(2, 1, IO[2,1])
	invert.SetElement(2, 2, IO[2,2])
	invert.SetElement(2, 3, IO[2,3])
	
	invert.SetElement(3, 0, IO[3,0])
	invert.SetElement(3, 1, IO[3,1])
	invert.SetElement(3, 2, IO[3,2])
	invert.SetElement(3, 3, IO[3,3])
	invert.Invert()
		
	transform = vtk.vtkTransform()
	transform.Identity();
	transform.Concatenate(invert);
	cube.GetAssembly().SetUserTransform( transform );
	
	axes2 = vtk.vtkAxesActor()
	axes2.SetShaftTypeToCylinder();
	axes2.SetUserTransform( transform );		 
	axes2.SetTotalLength( 1.5, 1.5, 1.5 );
	axes2.SetCylinderRadius( 0.500 * axes2.GetCylinderRadius() );
	axes2.SetConeRadius( 1.025 * axes2.GetConeRadius() );
	axes2.SetSphereRadius( 1.500 * axes2.GetSphereRadius() );

	tprop2 = axes2.GetXAxisCaptionActor2D()
	tprop2.GetCaptionTextProperty();

	assembly = vtk.vtkPropAssembly();
	assembly.AddPart( axes2 );
	assembly.AddPart( cube );

	widget = vtk.vtkOrientationMarkerWidget();
	widget.SetOutlineColor( 0.9300, 0.5700, 0.1300 );
	widget.SetOrientationMarker( assembly );
	widget.SetInteractor( iren1 );
	widget.SetViewport( 0.0, 0.0, 0.4, 0.4 );
	widget.SetEnabled( 1 );
	widget.InteractiveOff();
			
	# Set Up Camera view
	renderer1.SetBackground(0.0, 0.0, 0.0)
	camera = renderer1.GetActiveCamera()

	#bounds and initialize camera
	b = image.GetBounds()
	renderer1.ResetCamera(b)	
	renderer1.ResetCameraClippingRange()
	camera.SetViewUp(0.0,-1.0,0.0)
	camera.Azimuth(315)
	
	# Create a text property for both cube axes
	tprop = vtk.vtkTextProperty()
	tprop.SetColor(1, 1, 1)
	tprop.ShadowOff()
	
	# Create a vtkCubeAxesActor2D.  Use the outer edges of the bounding box to
	# draw the axes.  Add the actor to the renderer.
	axes = vtk.vtkCubeAxesActor2D()
	axes.SetInput(image)
	axes.SetCamera(renderer1.GetActiveCamera())
	axes.SetLabelFormat("%6.4g")
	axes.SetFlyModeToOuterEdges()
	axes.SetFontFactor(1.2)
	axes.SetAxisTitleTextProperty(tprop)
	axes.SetAxisLabelTextProperty(tprop)      
	renderer1.AddViewProp(axes)
	
	############
	# Place the interactor initially. The input to a 3D widget is used to
	# initially position and scale the widget. The "EndInteractionEvent" is
  	# observed which invokes the SelectPolygons callback.
    	# boxWidget.SetInput(image)
	# boxWidget.PlaceWidget()
	# boxWidget.AddObserver("InteractionEvent", SelectPolygons)
	# boxWidget.On()

	# Initizalize
	renWin1.Render()
	renderer1.Render()
	iren1.Initialize()
	iren1.Start()
	renderer1.RemoveViewProp(axes)
	
				
	return zImagePlaneWidget.GetSliceIndex()

def display_pick(image, spatial_res, slice_thickn, image_pos_pat, image_ori_pat, dicomReader, xImagePlaneWidget, yImagePlaneWidget, zImagePlaneWidget, renWin1, renderer1, iren1, picker, zplane):
	global origExt, xMin, xMax, yMin, yMax, zMin, zMax, xSpacing, ySpacing, zSpacing, origin, seeds
	
	# Initialize Image orienation
	IO = matrix(	[[0, 0,-1, 0],
			[1, 0, 0, 0],
			[0,-1, 0, 0],
			[0, 0, 0, 1]])
	# Assign the 6-Image orientation patient coordinates (from Dicomtags)
	IO[0,0] = image_ori_pat[0]; IO[1,0] = image_ori_pat[1]; IO[2,0] = image_ori_pat[2]; 
	IO[0,1] = image_ori_pat[3]; IO[1,1] = image_ori_pat[4]; IO[2,1] = image_ori_pat[5]; 
	
	# obtain thrid column as the cross product of column 1 y 2
	IO_col1 = [image_ori_pat[0], image_ori_pat[1], image_ori_pat[2]]
	IO_col2 = [image_ori_pat[3], image_ori_pat[4], image_ori_pat[5]]
	IO_col3 = cross(IO_col1, IO_col2)
	
	# assign column 3	
	IO[0,2] = IO_col3[0]; IO[1,2] = IO_col3[1]; IO[2,2] = IO_col3[2]; 
	
	IP =  [0, 0, 0, 1] # Initialization Image Position
	IP[0] = image_pos_pat[0]; IP[1] = image_pos_pat[1]; IP[2] = image_pos_pat[2];  
	
	print "image_pos_pat :"
	print image_pos_pat
	print "image_ori_pat:"
	print image_ori_pat
	
	origin = IP*IO.I
	print "Volume Origin:"
	print origin[0,0], origin[0,1], origin[0,2]
	
	# Calculate the center of the volume
	image.UpdateInformation()
	origExt = image.GetWholeExtent()
	(xMin, xMax, yMin, yMax, zMin, zMax) = image.GetWholeExtent()
	(xSpacing, ySpacing, zSpacing) = image.GetSpacing()
	(x0, y0, z0) = image.GetOrigin()
	zcen = (zMax - zMin)/2
	ycen = (yMax - yMin)/2
	xcen = (xMax - xMin)/2
	
	print "Image center"
	print xcen, ycen, zcen
	print "Image Extension"
	print xMin, xMax, yMin, yMax, zMin, zMax

	# Visualize results
	xImagePlaneWidget.SetInput( image )
	xImagePlaneWidget.SetPlaneOrientationToXAxes();
	yImagePlaneWidget.SetInput( image )
	yImagePlaneWidget.SetPlaneOrientationToYAxes();
	zImagePlaneWidget.SetInput( image )
	zImagePlaneWidget.SetPlaneOrientationToZAxes();
	zImagePlaneWidget.SetSliceIndex(zplane)
	
	#xImagePlaneWidget.SetInteractor( iren1 )
	xImagePlaneWidget.On()
	
	#yImagePlaneWidget.SetInteractor( iren1 )
	yImagePlaneWidget.On()

	#zImagePlaneWidget.SetInteractor( iren1 )
	zImagePlaneWidget.On()
			
	# Initizalize
	iren1.SetPicker(picker)
	renderer1.AddActor2D(textActor)
	picker.AddObserver("EndPickEvent", annotatePick)
	renWin1.Render()
	renderer1.Render()
	iren1.Start()
				
	return seeds
	
# Create a Python function to create the text for the text mapper used
# to display the results of picking.
def annotatePick(object, event):
	global picker, textActor, textMapper, seeds
	print "pick"
	seeds = vtk.vtkPoints()
	
	if(picker.GetCellId() < 0):
		textActor.VisibilityOff()     
	else:
		selPt = picker.GetSelectionPoint()
		pickPos = picker.GetPickPosition()
		seeds.InsertNextPoint(float(pickPos[0]), float(pickPos[1]), float(pickPos[2]) )
		print pickPos
		
		textMapper.SetInput("(%.6f, %.6f, %.6f)"%pickPos)
		textActor.SetPosition(selPt[:2])
		textActor.VisibilityOn()
	return
	
def SelectPolygons(object, event):
	# This callback funciton does the actual work: updates the vtkPlanes
	# implicit function.  This in turn causes the pipeline to update.
	# object will be the boxWidget
	global selectActor, planes
	planes = vtk.vtkPlanes()
	object.GetPlanes(planes)
		
	
def create_SER(img_folder, SeriesID, eID, StudyID, abspath_ExamID, row, dicomInfo, tablefile_name, tablewriter, irow_DynSeries):
	global chosen_lesions_id, dicomInfo_series, abspath_SeriesID, lThre, uThre
	
	# Makedir of current loc
	os.chdir(abspath_ExamID)
	DynVolstack =  vtk.vtkImageAppendComponents() 
	
	if os.path.exists('DynPhases'):
		print '''DynPhases'''
		
		abspath_SeriesID = 'Z:/Cristina/MassNonMass/'+img_folder+StudyID+os.sep+eID
						
		chosen_folderID = raw_input('Enter Series# corresponding to the Ph1-SCAN: ')
		phases_series=[]
		testSID = SeriesID[int(chosen_folderID)]
		if 'S' in str(testSID):
			print testSID[1:]
			chosen_phase = int(testSID[1:])
		else:
			chosen_phase = int(testSID)
			
		phases_series.append(chosen_phase)
		
		print "Printing available SeriesID"
		print "%s     %s " % ('i', 'Series#')
		
		i=1
				
		for chSer in [chosen_phase+1, chosen_phase+2, chosen_phase+3, chosen_phase+4]:
			for SerID in SeriesID:
				if 'S' in str(SerID):
					SerID = SerID[1:]
				try:	
					DynphasesFolder = int(SerID)
				except ValueError:
					DynphasesFolder = 0
				if ( DynphasesFolder == chSer):
					print "%d     %s " % (i, SerID)
					phases_series.append( int(SerID) )
					i=i+1
		# *****************************************************************************
		print '''\nLOADING IMAGES '''
		print '''LESION INFO FROM DATABASE '''
		print "%s" % str("Description: ")+str( row[4])+'\n'
		print "%s" % str("AdditonalDescription: ")+str( row[5])+'\n'
		print "%s" % str("BreastSide: ")+str(row[6])+'\n'+str("QuadrantIQ: ")+str(row[7])+'\n'+str("ClockStart: ")+str(row[8])+'\n'+str("ClockEnd: ")+str(row[9])+'\n'
		print "%s" % str("MRISeries: ")+str(row[10])+'\n'+str("LesionCenterPixelX: ")+str( row[11])+'\n'+str("LesionCenterPixelY: ")+str( row[12])+'\n'+str("LesionFindingType.Description ")+str( row[13])+'\n'
		print "%s" % str("WorstHistoType: ")+str( row[14])+'\n'+str("PathoResultTypeID: ")+str( row[15])+'\n'+str("PathoResultTypeID.Des: ")+str( row[16])+'\n'
		print "%s" % str("DCEDelayEnhPhaseType: ")+str( row[17])+'\n'+ str("LesionFindingID: ")+str( row[18])+str("DCEInitialEnhPhaseType.Desc ")+str( row[19])+'\n'+str("MassIntEnhanceTypeID: ")+str( row[21])
		print "%s" % str("MassLesionIntEnhanceType.Desc: ")+str( row[22])+'\n'+str("MassShapeTypeID.Desc: ")+str( row[23])+'\n'+str("MassLesionShapeType.Desc: ")+str( row[24])+'\n'+str("MassMarginTypeID: ")+str( row[25])+'\n'+str("MassLesionMarginType.Desc: ")+str( row[26])+'\n'		
		'''-----------------------------'''
				
		############ Get number of lesions and their IDs
		print lesions_id
		chosen_id = raw_input('Enter which lesion to segment: ')
		for lesion_j in lesions_id:
			if(int(chosen_id) == int(lesion_j)):
				chosen_lesions_id = lesion_j
				print "Chosen id: %d" % chosen_lesions_id
					
		###### OBTAIN INFO ABOUT LOCATION
		# Make sure you're in root ExamId dir
		os.chdir(abspath_ExamID)  
		
		# Obtain all locations in all phases_series
		for phaseS in phases_series:
			if 'S' in str(testSID):
				phaseS = 'S'+str(phaseS)
				
			abspath_PhaseID = abspath_ExamID+os.sep+str(phaseS) # this will return last element on list (so phase1)
			print abspath_PhaseID
			
			# Get total number of files
			list_files = get_only_filesindirectory(abspath_PhaseID)
			path_filenameID = abspath_PhaseID+os.sep+list_files[0]										
			dicomInfo_series = dicom.read_file(os.path.abspath(path_filenameID)) 
			
			dicomReader  = vtk.vtkDICOMImageReader()
			dicomReader.SetDirectoryName( abspath_PhaseID )
			dicomReader.Update()
			
			# Get image from reader
			im = dicomReader.GetOutput()
			im_scalars = im.GetPointData().GetScalars()
			dims = im.GetDimensions()
				
			# check wether DicomTag corresponds to Image loaded in vtk
			print "\n... Reading %s" % abspath_PhaseID
			print "VTK Dimensions im.GetDimensions(): %d %d %d" % dims
			spacing = im.GetSpacing()
			print "VTK Spacing im.GetSpacing(): %f %f %f\n" % spacing
			
			# vtkImageAppendComponents of the 5 dyn volumnes
			DynVolstack.AddInput(dicomReader.GetOutput())
			DynVolstack.Update()
			
			print ''' DICOM INFO HEADER '''
			print "%s" % str("Patient Position: ")+str(dicomInfo_series[0x0018,0x5100].value)+'\n'+str("Image Position (patient): ")+str(dicomInfo_series[0x0020,0x0032].value)
			print "%s" % str("Image Orientation (patient): ")+str(dicomInfo_series[0x0020,0x0037].value)+'\n'
			print "%s" % str("SpatialResolution: ")+str(spacing[0])+'\n'
			print "%s" % str("SliceThickness: ")+str(dicomInfo_series[0x0018,0x0050].value)+'\n'+str("SpacingBetweenSlices: ")+str(dicomInfo_series[0x0018,0x0088].value)
			print "%s" % str("SliceLocation: ")+str(dicomInfo_series.SliceLocation)+'\n'
			spatial_res = float(spacing[0])
			slice_thickn = float(dicomInfo_series[0x0018,0x0050].value)  
			image_pos_pat = list(dicomInfo_series[0x0020,0x0032].value)
			image_ori_pat = list(dicomInfo_series[0x0020,0x0037].value)		
	else:
		print '''SeriesPhases'''
				
		choseSerie = raw_input('Enter Series# corresponding to the ALL-PHASES-SCAN: ')
		path_SeriesID = img_folder+StudyID+os.sep+eID+os.sep+SeriesID[int(choseSerie)]
		print path_SeriesID
		
		# abspath_SeriesID = os.path.abspath(path_SeriesID) -- doesn't work bug
		abspath_SeriesID = 'Z:'+os.sep+'Cristina'+os.sep+'MassNonMass'+os.sep+path_SeriesID
				
		listLocations = sorted(list(get_immediate_subdirectories(abspath_SeriesID)))
		print "Total:  %d " % len(listLocations)
		
		#print ExamsID
		c = 0
		if(len(listLocations)>1):
			print "%s \t %s \t\t\t %s " % ('n', 'Series#', '#Images')
			for iloc in listLocations:
				path_loc_ID = abspath_SeriesID+os.sep+iloc
				path_loc_images = get_only_filesindirectory(path_loc_ID);
				
				print "%s \t %s \t\t\t %s " % (c, str(iloc), len(path_loc_images))
				c=c+1	 
		
		# *****************************************************************************
		print '''\nLOADING IMAGES '''
		print '''LESION INFO FROM DATABASE '''
		print "%s" % str("Description: ")+str( row[4])+'\n'
		print "%s" % str("AdditonalDescription: ")+str( row[5])+'\n'
		print "%s" % str("BreastSide: ")+str(row[6])+'\n'+str("QuadrantIQ: ")+str(row[7])+'\n'+str("ClockStart: ")+str(row[8])+'\n'+str("ClockEnd: ")+str(row[9])+'\n'
		print "%s" % str("MRISeries: ")+str(row[10])+'\n'+str("LesionCenterPixelX: ")+str( row[11])+'\n'+str("LesionCenterPixelY: ")+str( row[12])+'\n'+str("LesionFindingType.Description ")+str( row[13])+'\n'
		print "%s" % str("WorstHistoType: ")+str( row[14])+'\n'+str("PathoResultTypeID: ")+str( row[15])+'\n'+str("PathoResultTypeID.Des: ")+str( row[16])+'\n'
		print "%s" % str("DCEDelayEnhPhaseType: ")+str( row[17])+'\n'+ str("LesionFindingID: ")+str( row[18])+str("DCEInitialEnhPhaseType.Desc ")+str( row[19])+'\n'+str("MassIntEnhanceTypeID: ")+str( row[21])
		print "%s" % str("MassLesionIntEnhanceType.Desc: ")+str( row[22])+'\n'+str("MassShapeTypeID.Desc: ")+str( row[23])+'\n'+str("MassLesionShapeType.Desc: ")+str( row[24])+'\n'+str("MassMarginTypeID: ")+str( row[25])+'\n'+str("MassLesionMarginType.Desc: ")+str( row[26])+'\n'		
		'''-----------------------------'''
				
		############ Get number of lesions and their IDs
		print lesions_id
		chosen_id = raw_input('Enter which lesion to segment: ')
		for lesion_j in lesions_id:
			if(int(chosen_id) == int(lesion_j)):
				chosen_lesions_id = lesion_j
				print "Chosen id: %d" % chosen_lesions_id
					
		# process all Volumes when in stacks of Dyn Volumes
		# listLocations[int(choseDisplay)] == 'pre-Contrast' ):
		filename_series = 'DIRCONTENTS.txt'
		for i in range(5):
			print i
			if (i==0):
				chosen_folderID = abspath_SeriesID+os.sep+'pre-Contrast'
				os.chdir(chosen_folderID)
				 
			else:
				chosen_folderID = abspath_SeriesID+os.sep+'post_Contrast-'+str(i)
				os.chdir(chosen_folderID)
					
			# Get total number of files
			list_files = get_only_filesindirectory(chosen_folderID)
			path_filenameID = chosen_folderID+os.sep+list_files[0]										
			dicomInfo_series = dicom.read_file(os.path.abspath(path_filenameID)) 
					
			#check firs if it exists
			if(FileCheck(filename_series) ):
				vtkStringArray = vtk.vtkStringArray()
				file_series = open(filename_series,'r')
				files_toRead = file_series.read()
				try:
					for slicename in files_toRead.split():
						#print "adding : %s" %  slicename
						vtkStringArray.InsertNextValue( slicename )
				finally:
					file_series.close()
				
			# Read dicom Vol from DIRCONTENTS.txt
			dicomReader  = vtkgdcmPython.vtkGDCMImageReader()
			dicomReader.SetFileNames( vtkStringArray )
			dicomReader.FileLowerLeftOn()
			dicomReader.Update()
			
			# Get image from reader
			im_scalars = dicomReader.GetOutput().GetPointData().GetScalars()
			dims = dicomReader.GetOutput().GetDimensions()
				
			# check wether DicomTag corresponds to Image loaded in vtk
			print "\n... Reading %s" % chosen_folderID
			print "VTK Dimensions im.GetDimensions(): %d %d %d" % dims
			spacing =  dicomReader.GetOutput().GetSpacing()
			print "VTK Spacing im.GetSpacing(): %f %f %f\n" % spacing
			
			# vtkImageAppendComponents of the 5 dyn volumnes
			DynVolstack.AddInput(dicomReader.GetOutput())
			DynVolstack.Update()
			
			os.chdir(abspath_SeriesID)
										
			print ''' DICOM INFO HEADER '''
			print "%s" % str("Patient Position: ")+str(dicomInfo_series[0x0018,0x5100].value)+'\n'+str("Image Position (patient): ")+str(dicomInfo_series[0x0020,0x0032].value)
			print "%s" % str("Image Orientation (patient): ")+str(dicomInfo_series[0x0020,0x0037].value)+'\n'
			print "%s" % str("SpatialResolution: ")+str(dicomInfo_series.SpatialResolution)+'\n'+str("TemporalResolution: ")+str(dicomInfo_series.TemporalResolution)
			print "%s" % str("SliceThickness: ")+str(dicomInfo_series[0x0018,0x0050].value)+'\n'+str("SpacingBetweenSlices: ")+str(dicomInfo_series[0x0018,0x0088].value)
			print "%s" % str("SliceLocation: ")+str(dicomInfo_series.SliceLocation)+'\n'
			spatial_res = float(dicomInfo_series.SpatialResolution)
			slice_thickn = float(dicomInfo_series[0x0018,0x0050].value)  
			image_pos_pat = list(dicomInfo_series[0x0020,0x0032].value)
			image_ori_pat = list(dicomInfo_series[0x0020,0x0037].value)
				
	##############################################################################
	DynVol_pre = vtk.vtkImageCast()
	DynVol_early_pos = vtk.vtkImageCast()
	DynVol_late_pos = vtk.vtkImageCast()
	DynVol_time_pos = vtk.vtkImageCast()
	
	# Obtain pre contrast volume
	DynVol_pre.SetInput(DynVolstack.GetInput(0))
	DynVol_pre.SetOutputScalarTypeToFloat()
	DynVol_pre.Update()
	
	# Obtain early post-contrast volume
	DynVol_early_pos.SetInput(DynVolstack.GetInput(1))
	DynVol_early_pos.SetOutputScalarTypeToFloat()
	DynVol_early_pos.Update()
	
	# substratc
	sub_pre_early = vtk.vtkImageMathematics()
	sub_pre_early.SetOperationToSubtract()
	sub_pre_early.SetInput1(DynVol_early_pos.GetOutput())
	sub_pre_early.SetInput2(DynVol_pre.GetOutput())
	sub_pre_early.Update()
	
	#send to display
	#image = vtk.vtkImageData()
	#image.DeepCopy(sub_pre_early.GetOutput())
		
	zplane = display(sub_pre_early.GetOutput(), image_pos_pat, image_ori_pat, xImagePlaneWidget, yImagePlaneWidget, zImagePlaneWidget, renWin1, renderer1, iren1)
		
	#initialize vcars
	num_voxels = []
	volume_lesion = []
	
		
	# ############ create slider widget
	# root = Tkinter.Tk()
	# root.title("Lesion Segment")
	# root.iconname("LesS")
	# root.withdraw()
	# top = Tkinter.Toplevel(root)
	# 
	# # create the main frame and the main panel for the render widget
	# mainPanel = Tkinter.Frame(top, relief='raised', width=150, borderwidth=2)
	# mainPanel.pack(side='left', fill='y', expand='true')
	# 
	# f4 = Tkinter.Frame(mainPanel, relief='groove', borderwidth=2, padx=20, pady=20)
	# l4 = Tkinter.Label(mainPanel, text='Lower Level', bg='black', fg='white')
	# 
	# slice_lower = Tkinter.IntVar()
	# slice_lower.set(50)
	# lThre=50
	# uThre=1000
	# lower = Tkinter.Scale(f4, from_=lThre, to=uThre, orient="horizontal", command=SetSlice, variable=slice_lower)
	# lower.pack(side='top', fill='x', expand='true')
	# 
	# f5 = Tkinter.Frame(mainPanel, relief='groove', borderwidth=2, padx=20, pady=20)
	# l5 = Tkinter.Label(mainPanel, text='Upper Level', bg='black', fg='white')
	# slice_upper = Tkinter.IntVar()
	# slice_upper.set(1000)
	# upper = Tkinter.Scale(f4, from_=lThre, to=uThre, orient="horizontal", command=SetSlice, variable=slice_upper)
	# upper.pack(side='top', fill='x', expand='true')
	# 
	# for i in (l4, f4, l5, f5):
	    # i.pack(side='top', expand='false', fill='x')
	
	for i in range(1,5):
		print "Procesing Substraction image %d " % i
		#DynVolstack.Update()		
		DynVol_time_pos.SetInput(DynVolstack.GetInput(i))
		DynVol_time_pos.SetOutputScalarTypeToFloat()
		DynVol_time_pos.Update()
		
		# substratc
		sub_pre = vtk.vtkImageMathematics()
		sub_pre.SetOperationToSubtract()
		sub_pre.SetInput1(DynVol_time_pos.GetOutput())
		sub_pre.SetInput2(DynVol_pre.GetOutput())
		sub_pre.Update()
		
		#send to display
		#image = vtk.vtkImageData()
		#image.DeepCopy(sub_pre.GetOutput())
		
		print "Displaying binary file and lesion segmentation: %d" % i
		seeds = vtk.vtkPoints()
		seeds = display_pick(sub_pre.GetOutput(), spatial_res, slice_thickn, image_pos_pat, image_ori_pat, dicomReader, xImagePlaneWidget, yImagePlaneWidget, zImagePlaneWidget, renWin1, renderer1, iren1, picker, zplane)
		
		# vtkImageThresholdConnectivity will perform a flood fill on an image, given upper and lower pixel intensity 
		# thresholds. It works similarly to vtkImageThreshold, but also allows the user to set seed points to limit
		# the threshold operation to contiguous regions of the image. The filled region, or the "inside", will be passed 
		# through to the output by default, while the "outside" will be replaced with zeros. The scalar type of the output is the same as the input.
		
		image_scalar_range = sub_pre.GetOutput().GetScalarRange() 
		print "Image Scalar Range:"
		print image_scalar_range[0], image_scalar_range[1]
		lThre = image_scalar_range[0]
		uThre = image_scalar_range[1]
				
			
		# Extract each time point substraction to segment and finally extract SER from first time point subs-ROI
		thresh_sub = vtk.vtkImageThresholdConnectivity() 
		thresh_sub.SetSeedPoints(seeds)
		# #thresh_sub.SetSliceRangeX( xExt )
		# #thresh_sub.SetSliceRangeY( yExt )
		#thresh_sub.SetSliceRangeZ( zExt )
		thresh_sub.SetInValue(1) 
		thresh_sub.SetOutValue(0)
		thresh_sub.ReplaceInOff()
		thresh_sub.ReplaceOutOn()
		# thresh_sub.SetNeighborhoodRadius(3, 3, 3) #The radius of the neighborhood that must be within the threshold values in order for the voxel to be included in the mask. The default radius is zero (one single voxel). The radius is measured in voxels
		#thresh_sub.SetNeighborhoodFraction(0.10) #The fraction of the neighborhood that must be within the thresholds. The default value is 0.5.
		thresh_sub.ThresholdBetween(0.35*uThre, uThre); 
		thresh_sub.SetInput(sub_pre.GetOutput())
		thresh_sub.Update()
		
			
		#find out how big lesion is GetNumberOfInVoxels * VoxelVolume = lesion Volume
		num_voxels.append(thresh_sub.GetNumberOfInVoxels())
		print "Number of Voxels #: %d" % num_voxels[i-1] # After the filter has executed, use GetNumberOfVoxels() to find out how many voxels were filled.
		
		# Nothing found
		if( num_voxels[i-1] == 0):
			xExt = []
			yExt = []
			zExt = []
			# chosenExtX = int(raw_input('Enter min int for updated X extent: '))
			# xExt.append( int(chosenExtX*xSpacing) )
			# chosenExtX = int(raw_input('Enter max int for updated X extent: '))
			# xExt.append( int(chosenExtX*xSpacing) )
			# 
			# chosenExtY = int(raw_input('Enter min int for updated Y extent: '))
			# yExt.append( int(chosenExtY*ySpacing) )
			# chosenExtY = int(raw_input('Enter max int for updated Y extent: '))
			# yExt.append( int(chosenExtY*ySpacing) )
			# 
			# chosenExtZ = int(raw_input('Enter min int for updated Z extent: '))
			# zExt.append( int(chosenExtZ*zSpacing) )
			# chosenExtZ = int(raw_input('Enter max int for updated Z extent: '))
			# zExt.append( int(chosenExtZ*zSpacing) )
			
			thresh_sub = vtk.vtkImageThreshold()
			thresh_sub.SetInput(sub_pre.GetOutput())
			thresh_sub.SetOutValue(0)
			thresh_sub.SetInValue(1)
			thresh_sub.ThresholdBetween(0.25*uThre, uThre)
			
		# After the filter has executed, use GetNumberOfVoxels() to find out how many voxels were filled.
		volume_lesion.append(num_voxels[i-1]*(spatial_res)*(spatial_res)*(dicomInfo_series[0x0018,0x0050].value)) 
		 
		print "Volume of Lesion mm3: %d" % volume_lesion[i-1]
		
		# Makedir of current loc
		# if current loc folder doesn't exist create it		
		os.chdir(abspath_SeriesID)
			
		VOIlesion_ID = 'VOIlesions_id'+str(chosen_lesions_id)
		if not os.path.exists(VOIlesion_ID):
			os.makedirs(VOIlesion_ID)
		
		os.chdir(VOIlesion_ID)
		image_VOIlesion = vtk.vtkImageThreshold()
		image_VOIlesion.ThresholdByUpper(1)
		image_VOIlesion.SetInValue(255)
		image_VOIlesion.SetOutValue(0)
		image_VOIlesion.SetInputConnection(thresh_sub.GetOutputPort())
		image_VOIlesion.Update()
		
		VOIlesion_name = 'VOIlesion_subs_'+str(i)+'.mhd'
		print "Writing VOIlesion to binary file: %s" % VOIlesion_name 
		VOIlesion_writer = vtk.vtkMetaImageWriter()
		VOIlesion_writer.SetFileName(VOIlesion_name)
		VOIlesion_writer.SetInput(thresh_sub.GetOutput())
		VOIlesion_writer.Write() 
		
		# Convert VOIlesion into polygonal struct
		VOIlesion_poly = vtk.vtkMarchingCubes() 
		VOIlesion_poly.SetValue(0,255)
		VOIlesion_poly.SetInput(image_VOIlesion.GetOutput())
		VOIlesion_poly.ComputeNormalsOff()
		VOIlesion_poly.Update()
		 
		# Visualize iso surface
		# Store Marching cubes output as PolyData 
		VOIlesion_mapper = vtk.vtkPolyDataMapper()
		VOIlesion_actor  = vtk.vtkActor()
		
		VOIlesion_mapper.SetInput( VOIlesion_poly.GetOutput() )
		VOIlesion_mapper.ScalarVisibilityOff()					
		VOIlesion_actor.SetMapper( VOIlesion_mapper )
		VOIlesion_actor.GetProperty().SetColor( 0,1,0 )
		VOIlesion_actor.GetProperty().SetOpacity( 0.6 )
		VOIlesion_actor.GetProperty().SetLineWidth(2.0)
		VOIlesion_actor.GetProperty().SetRepresentationToWireframe()
		renderer1.AddActor(VOIlesion_actor)
		
		# Write OIlesion_name_poly to file
		VOIlesion_name_poly = VOIlesion_name[0:-3]+'vtk'
		VOIlesion_poly_writer = vtk.vtkPolyDataWriter()
		VOIlesion_poly_writer.SetInput(VOIlesion_poly.GetOutput())
		VOIlesion_poly_writer.SetFileName(  VOIlesion_name_poly )
		VOIlesion_poly_writer.Write()
		VOIlesion_poly_writer.Update()
		
		# Write VTK XML PolyData files.
		# vtkXMLPolyDataWriter writes the VTK XML PolyData file format. One polygonal data input 
		# can be written into one file in any number of streamed pieces (if supported by the rest of the
		# pipeline). The standard extension for this writer's file format is "vtp". 
		# This writer is also used to write a single piece of the parallel file format.
		VOIlesion_name_polyXML = VOIlesion_name[0:-3]+'xml'
		VOIlesion_polyXML_writer = vtk.vtkXMLPolyDataWriter()
		VOIlesion_polyXML_writer.SetInput(VOIlesion_poly.GetOutput())
		VOIlesion_polyXML_writer.SetFileName(  VOIlesion_name_polyXML )
		VOIlesion_polyXML_writer.SetDataModeToAscii()
		VOIlesion_polyXML_writer.Write()
		
		# Initizalize
		renWin1.Render()
		renderer1.Render()
		iren1.Start()
		renderer1.RemoveViewProp(VOIlesion_actor)
				
		
	# Finish stacking 5 dynamic volumes	
	############################# Apply SER
	# Extract S0,S1 and S2 are the precontrast (baseline), early postcontrast 
	# (2.8 minutes after injection) and late postcontrast 
	# (11.2 minutes after injection) signal intensities.
	DynVol_late_pos.SetInput(DynVolstack.GetInput(4))
	DynVol_late_pos.SetOutputScalarTypeToFloat()
	DynVol_late_pos.Update()

	sub_pre_late = vtk.vtkImageMathematics()
	sub_pre_late.SetOperationToSubtract()
	sub_pre_late.SetInput1(DynVol_late_pos.GetOutput())
	sub_pre_late.SetInput2(DynVol_pre.GetOutput())
	sub_pre_late.Update()
	
	# divide to obtain SER
	divide_SER = vtk.vtkImageMathematics()
	divide_SER.SetOperationToDivide()
	divide_SER.SetInput1(sub_pre_early.GetOutput())
	divide_SER.SetInput2(sub_pre_late.GetOutput())
	divide_SER.Update()
	

	# Create an Image of model_Updated
	white_image = vtk.vtkImageData()
	white_image.SetDimensions(dims[0], dims[1], dims[2])
	white_image.SetWholeExtent(xMin, xMax, yMin, yMax, zMin, zMax)
	white_image.SetSpacing(xSpacing, ySpacing, zSpacing)
	white_image.SetOrigin(origin[0,0], origin[0,1], origin[0,2])
	white_image.SetScalarTypeToUnsignedChar()
	white_image.AllocateScalars()
	white_image.DeepCopy(divide_SER.GetOutput())
		
	############################# finished SER
	# polygonal data --> image stencil:
	SER_pol2stenc = vtk.vtkPolyDataToImageStencil()
	SER_pol2stenc.SetInput(VOIlesion_poly.GetOutput())
	SER_pol2stenc.SetOutputOrigin(origin[0,0], origin[0,1], origin[0,2])
	SER_pol2stenc.SetOutputSpacing(xSpacing, ySpacing, zSpacing)
	SER_pol2stenc.SetOutputWholeExtent(xMin, xMax, yMin, yMax, zMin, zMax)
	SER_pol2stenc.Update()

	# cut the corresponding SER image and set the background:
	SER_imgstenc = vtk.vtkImageStencil()
	SER_imgstenc.SetInput(white_image)
	SER_imgstenc.SetStencil(SER_pol2stenc.GetOutput())
	SER_imgstenc.ReverseStencilOff()
	SER_imgstenc.SetBackgroundValue(0.0)
	SER_imgstenc.Update()
	
	SERlesion_name = 'SERlesion.mhd'
	print "Writing VOIlesion to binary file: %s" % SERlesion_name 
	SERlesion_writer = vtk.vtkMetaImageWriter()
	SERlesion_writer.SetFileName(SERlesion_name)
	SERlesion_writer.SetInput(SER_imgstenc.GetOutput())
	SERlesion_writer.Write() 
	
	######################### Add to table
	# within a loop of rows returned from a MySQL database (r refers to row number, c to column number, and content is what you want to put in the cell):
	tablewriter.activateSheet("LesionSegment") 
	datarow_DynSeries = [row[1], row[2], row[0], str(dicomInfo[0x0012,0x0040].value), str(dicomInfo[0x0008,0x0020].value), str(dicomInfo[0x0008,0x103e].value), str(dicomInfo.SpatialResolution), str(dicomInfo[0x0018,0x0080].value), str(dicomInfo[0x0018,0x0081].value), str(dicomInfo[0x0018,0x1314].value), str(dicomInfo[0x0018,0x0050].value), str(dicomInfo[0x0018,0x0088].value), row[6], row[7], row[8], row[9], row[11], row[12], row[13], chosen_lesions_id, num_voxels[0], volume_lesion[0], num_voxels[1], volume_lesion[1], num_voxels[2], volume_lesion[2], num_voxels[3], volume_lesion[3] ]      
	tablewriter.addRow(irow_DynSeries, datarow_DynSeries)
	irow_DynSeries = irow_DynSeries+1
	
	#send to display
	#image = vtk.vtkImageData()
	#image.DeepCopy(SER_imgstenc.GetOutput())
	#print image
	
	#display(image, image_pos_pat, image_ori_pat, xImagePlaneWidget, yImagePlaneWidget, zImagePlaneWidget, renWin1, renderer1, iren1)
	DynVolstack =  [] 
	dicomReader = []			
	return

def SetSlice(current_widget):
	xyz = current_widget.GetCurrentCursorPosition()
	x0, y0, z0 = xyz
	pixVal = image.GetScalarComponentAsFloat( x0, y0, z0, 0 )
	
   	return pixVal
    
# ******************************************************************************
# initialize some variables
irow_DynSeries = 1
irow_lesionSegment = 1

########## Write to table
#initialize Table writer
tablefile_name = "table_SERlesion_MASS.xls" 
tablewriter = excel_writer.ExcelWriter(tablefile_name, "LesionSegment", make_visible=True) 
title_range = tablewriter.getRangeByCells((1, 1), (1, 19)) #where j = number of columns 
tablewriter.formatRange(title_range, excel_writer.STYLE_GREY_CELL)

# start by writing headings 
headings_DynSeries = ['StudyNumber', 'DicomExamNumber', 'MRN', 'ClinicalTrialID', 'StudyDate', 'SeriesDescription', 'SpatialResolution', 'TR', 'TE', 'FlipAngle', 'SliceThickness', 'SpacingBetweenSlices', 'BreastSide', 'QuadrantIQ', 'ClockStart', 'ClockEnd', 'LesionCenterPixelX', 'LesionCenterPixelY', 'LesionFindingType.Description', 'chosen_lesions_id', 'num_vox_sub1', 'volume_lesion_sub1', 'num_vox_sub2', 'volume_lesion_sub2', 'num_vox_sub3', 'volume_lesion_sub3','num_vo_sub4', 'volume_lesion_sub4']
tablewriter.addRow(irow_DynSeries, headings_DynSeries)
irow_DynSeries = irow_DynSeries+1
tablewriter.save() 


############
epsilon = .0001         

# use cell picker for interacting with the image orthogonal views.
picker = vtk.vtkCellPicker()
picker.SetTolerance(0.005)

textMapper = vtk.vtkTextMapper()
tprop = textMapper.GetTextProperty()
tprop.SetFontFamilyToArial()
tprop.SetFontSize(10)
tprop.BoldOn()
tprop.ShadowOn()
tprop.SetColor(1, 0, 0)

textActor = vtk.vtkActor2D()
textActor.VisibilityOff() 
textActor.SetMapper(textMapper)

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
		print "Getting studyID# DicomExam# "
		print StudyID, ExamID, MRN
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
						
		# *****************************************************************************
		# ask for which series to load
		print "\n----------------------------------------------------------"
		choseTodo = raw_input('Enter option to run: \n\t1) create_SERmap	\n\tx) to exit \n\t input: ')
		while(choseTodo!='x'):
			if(choseTodo == '1'):
				print '''\nCREATE SER MAP \n'''
				create_SER(img_folder, SeriesID, eID, StudyID, abspath_ExamID, row, dicomInfo, tablefile_name, tablewriter, irow_DynSeries)
										
						
			# *****************************************************************************
			print "Writing output log to file"
			chdirname='Z:/Cristina/MassNonmass/'+img_folder
			os.chdir(chdirname)  
			print os.getcwd()

			# Already run header line - now just append infor
			#file_summary_cases.write(str("MRN\t")+' '+str("StudyNumber\t")+' '+str("DicomExamNumber \t")+' '+str("Description \t")+' '+str("AdditonalDescription \t")+' '+str("AllFindings \t")+' '+str("BreastSide \t")+' '+str("QuadrantIQ \t")+' '+str("ClockStart \t")+' '+str("ClockEnd \t")+' '+str("MRISeries \t")+' '+str("LesionCenterPixelX \t")+' '+str("LesionCenterPixelY \t")+' '+str("LesionFindingType.Description \t")+' '+' '+str("PathoResultTypeID.Des \t")+' '+str("ACRScoreID \t")+' '+str("DCEInitialPhaseEnhTypeID \t")+' '+str("DCEInitialPhaseEnhTypeID.Desc \t")+' '+str("DCEDelayedPhaseEnhID \t")+' '+str("DCEDelayedPhaseEnhID.Desc \t")+' '+str("MassLesionIntEnhanceType \t")+' '+str("MassLesionIntEnhanceType.Desc \t")+' '+str("MassLesionShapeType \t")+' '+str("MassLesionShapeType.Desc \t")+' '+str("MassMarginTypeID \t")+' '+str("MassMarginTypeID.Desc \t")+' '+str("StudyDate\t")+' '+str("Patient id (MRN)\t")+' '+str("ClinicalTrialID\t")+' '+str("StudyID\t")+' '+str("PatientAge\t")+' '+str("Pulse sequence Name\t")+' '+str("Patient Position\t")+' '+str("Image Position (patient)\t")+' '+str("Image Orientation (patient)\t")+' '+str("First Scan location\t")+' '+str("Pulse sequence\t")+' '+str("SAT Fat/water/bone\t")+' '+str("Modality\t")+' '+str("SeriesNumber\t")+' '+str("SeriesDescription\t")+' '+str("SpatialResolution\t")+' '+str("TemporalResolution\t")+' '+str("TR\t")+' '+str("TE\t")+' '+str("Flip Angle\t")+' '+str("ETL\t")+' '+str("SliceThickness\t")+' '+str("SpacingBetweenSlices\t")+' '+str("Series\t")+'\n')
			if len(row) != 0:
				print "%s" % str("MRN: ")+str(row[0])+'\n'+str("StudyNumber: ")+str( row[1])+'\n'+str("DicomExamNumber:  ")+str(row[2])+'\n'
				print "%s" % str("Description: ")+str(row[4])+'\n'
				print "%s" % str("AdditonalDescription: ")+str( row[5])+'\n'
				print "%s" % str("BreastSide: ")+str(row[6])+'\n'+str("QuadrantIQ: ")+str(row[7])+'\n'+str("ClockStart: ")+str(row[8])+'\n'+str("ClockEnd: ")+str(row[9])+'\n'
				print "%s" % str("MRISeries: ")+str(row[10])+'\n'+str("LesionCenterPixelX: ")+str( row[11])+'\n'+str("LesionCenterPixelY: ")+str( row[12])+'\n'+str("LesionFindingType.Description ")+str( row[13])+'\n'
				print "%s" % str("WorstHistoType: ")+str( row[14])+'\n'+str("PathoResultTypeID: ")+str( row[15])+'\n'+str("PathoResultTypeID.Des: ")+str( row[16])+'\n'
				print "%s" % str("DCEDelayEnhPhaseType: ")+str( row[17])+ str("LesionFindingID: ")+str( row[18])+'\n'+str("DCEInitialEnhPhaseType.Desc ")+str( row[19])+'\n'+str("MassIntEnhanceTypeID: ")+str( row[21])
				print "%s" % str("MassLesionIntEnhanceType.Desc: ")+str( row[22])+'\n'+str("MassShapeTypeID.Desc: ")+str( row[23])+'\n'+str("MassLesionShapeType.Desc: ")+str( row[24])+'\n'+str("MassMarginTypeID: ")+str( row[25])+'\n'+str("MassLesionMarginType.Desc: ")+str( row[26])+'\n'
			print "%s" % str("StudyDate: ")+str(dicomInfo[0x0008,0x0020].value)+'\n'+str("Patient id (MRN): ")+str(dicomInfo[0x0010,0x0020].value)+'\n'+str("ClinicalTrialID: ")+str(dicomInfo[0x0012,0x0040].value)+'\n'
			print "%s" % str("StudyID: ")+str(dicomInfo[0x0020,0x0010].value)+'\n'+str("PatientAge: ")+str(dicomInfo[0x0010,0x1010].value)+'\n'+str("Pulse sequence Name: ")+str(dicomInfo[0x0019,0x109e].value)+'\n'
			print "%s" % str("Patient Position: ")+str(dicomInfo[0x0018,0x5100].value)+'\n'+str("Image Position (patient): ")+str(dicomInfo[0x0020,0x0032].value)+'\n'
			print "%s" % str("Image Orientation (patient): ")+str(dicomInfo[0x0020,0x0037].value)+'\n'+str("First Scan location: ")+str(dicomInfo[0x0019,0x1019].value)+'\n'+str("Pulse sequence: ")+str(dicomInfo[0x0019,0x109c].value)+'\n'
			print "%s" % str("Modality: ")+str(dicomInfo[0x0008,0x0060].value)+'\n'+str("SeriesNumber: ")+str(dicomInfo[0x0020,0x0011].value)+'\n'+str("SeriesDescription: ")+str(dicomInfo[0x0008,0x103e].value)+'\n'
			print "%s" % str("SpatialResolution: ")+str(dicomInfo_series.SpatialResolution)+'\n'+str("TemporalResolution: ")+str(dicomInfo_series.TemporalResolution)+'\n'+str("TR: ")+str(dicomInfo[0x0018,0x0080].value)+'\n'+str("TE: ")+str(dicomInfo[0x0018,0x0081].value)+'\n'+str("Flip Angle: ")+str(dicomInfo[0x0018,0x1314].value)+'\n'+str("ETL: ")+str(dicomInfo[0x0018,0x0091].value)+'\n'
			print "%s" % str("SliceThickness: ")+str(dicomInfo[0x0018,0x0050].value)+'\n'+str("SpacingBetweenSlices: ")+str(dicomInfo[0x0018,0x0088].value)+'\n'+str("Series: ")+str(dicomInfo[0x0019,0x1017].value)+'\n'
				
			choseTodo = raw_input('Enter option to load: \n\t1) create_SERmap	\n\tx) to exit \n\t input: ')
								
						
finally:
	file_ids.close()		
		
			
# Save the file list to read as series later
# Finally, save the file and clean up:
tablewriter.save()		
tablewriter.close()		
		
		
