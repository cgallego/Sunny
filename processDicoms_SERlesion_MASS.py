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
Run on MassNonmass root folder processDicoms.py based on list of 
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
----------------------------------------------------------------------
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
			
			NumberOfVolumes = 1 # default
			
			# iterate number of slices and get full volume (#slices, each slice loc)
			slices = []
			num_images=0;
			for filename in listSeries_files:
				num_images = num_images+1
					
			# Print series info						
			print "%d	%d		%s" % (s, num_images, arrangedF) 
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
			if(listSeries_files != []):
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
						
					# Get structure of study (all files in directory consistent with studyID and patientID
					studyTree = []
					FileNames=listSeries_files;
					
					if("PatientID" in dicomInfo):
						PatientID = dicomInfo.PatientID#
					else:	PatientID=''
					if("SeriesNumber" in dicomInfo):
						SeriesNumber = dicomInfo.SeriesNumber#
					else:	SeriesNumber=''
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
		#os.chdir(path_rootFolder)	
				
	return abspath_ExamID, eID, SeriesID, studyFolder, dicomInfo


def display(image, image_pos_pat, image_ori_pat):
	global xMin, xMax, yMin, yMax, zMin, zMax, xSpacing, ySpacing, zSpacing, origin
	# The box widget observes the events invoked by the render window
	# interactor.  These events come from user interaction in the render
	# window.
	# boxWidget = vtk.vtkBoxWidget()
	# boxWidget.SetInteractor(iren1)
	# boxWidget.SetPlaceFactor(1)
	# If one considers the localizer plane as a "viewport" onto the DICOM 3D coordinate space, then that viewport is described by its origin, its row unit vector, column unit vector and a normal unit vector (derived from the row and column vectors by taking the cross product). Now if one moves the origin to 0,0,0 and rotates this viewing plane such that the row vector is in the +X direction, the column vector the +Y direction, and the normal in the +Z direction, then one has a situation where the X coordinate now represents a column offset in mm from the localizer's top left hand corner, and the Y coordinate now represents a row offset in mm from the localizer's top left hand corner, and the Z coordinate can be ignored. One can then convert the X and Y mm offsets into pixel offsets using the pixel spacing of the localizer imag
	# Initialize Image orienation
	IO = matrix(	[[0, 0,-1, 0],
			[1, 0, 0, 0],
			[0,-1, 0, 0],
			[0, 0, 0, 1]])
	# Assign the 6-Image orientation patient coordinates (from Dicomtags)
	IO[0,0] = image_ori_pat[0]; IO[0,1] = image_ori_pat[1]; IO[0,2] = image_ori_pat[2]; 
	IO[1,0] = image_ori_pat[3]; IO[1,1] = image_ori_pat[4]; IO[1,2] = image_ori_pat[5]; 
	
	# obtain thrid column as the cross product of column 1 y 2
	IO_col1 = [image_ori_pat[0], image_ori_pat[1], image_ori_pat[2]]
	IO_col2 = [image_ori_pat[3], image_ori_pat[4], image_ori_pat[5]]
	IO_col3 = cross(IO_col1, IO_col2)
	
	# assign column 3	
	IO[2,0] = IO_col3[0]; IO[2,1] = IO_col3[1]; IO[2,2] = IO_col3[2]; 
	
	IP =  array([0, 0, 0, 1]) # Initialization Image Position
	IP[0] = image_pos_pat[0]; IP[1] = image_pos_pat[1]; IP[2] = image_pos_pat[2];  
	IO[0,3] = -image_pos_pat[0]; IO[1,3] = -image_pos_pat[1]; IO[2,3] = -image_pos_pat[2]
	
	print "printing Image Orientation Patient Matrix"
	print IO
	
	print "image_pos_pat :"
	print image_pos_pat
	print "image_ori_pat:"
	print image_ori_pat
	
	origin = IP*IO.I
	print "Volume Origin:"
	print origin[0,0], origin[0,1], origin[0,2]
	
	# Create matrix 4x4
	DICOM_mat = vtk.vtkMatrix4x4();
	DICOM_mat.SetElement(0, 0, IO[0,0])
	DICOM_mat.SetElement(0, 1, IO[0,1])
	DICOM_mat.SetElement(0, 2, IO[0,2])
	DICOM_mat.SetElement(0, 3, IO[0,3])
	
	DICOM_mat.SetElement(1, 0, IO[1,0])
	DICOM_mat.SetElement(1, 1, IO[1,1])
	DICOM_mat.SetElement(1, 2, IO[1,2])
	DICOM_mat.SetElement(1, 3, IO[1,3])
	
	DICOM_mat.SetElement(2, 0, IO[2,0])
	DICOM_mat.SetElement(2, 1, IO[2,1])
	DICOM_mat.SetElement(2, 2, IO[2,2])
	DICOM_mat.SetElement(2, 3, IO[2,3])
	
	DICOM_mat.SetElement(3, 0, IO[3,0])
	DICOM_mat.SetElement(3, 1, IO[3,1])
	DICOM_mat.SetElement(3, 2, IO[3,2])
	DICOM_mat.SetElement(3, 3, IO[3,3])
	#DICOM_mat.Invert()
	
	# Set up the axes	
	transform = vtk.vtkTransform()
	transform.Concatenate(DICOM_mat)
	transform.Update()
	
	# Set up the cube (set up the translation back to zero	
	DICOM_mat_cube = vtk.vtkMatrix4x4();
	DICOM_mat_cube.DeepCopy(DICOM_mat)
	DICOM_mat_cube.SetElement(0, 3, 0)
	DICOM_mat_cube.SetElement(1, 3, 0)
	DICOM_mat_cube.SetElement(2, 3, 0)
		
	transform_cube = vtk.vtkTransform()
	transform_cube.Concatenate(DICOM_mat_cube)
	transform_cube.Update()
	
	# Change info
	change_image = vtk.vtkImageChangeInformation()
	change_image.SetInput( image );
	change_image.SetOutputOrigin(origin[0,0], origin[0,1], origin[0,2])
	change_image.Update()
    
    	# Calculate the center of the volume
	change_image.GetOutput().UpdateInformation()
	dims = change_image.GetOutput().GetDimensions()
	print "Image Dimensions"
	print dims
	(xMin, xMax, yMin, yMax, zMin, zMax) = change_image.GetOutput().GetWholeExtent()
	print "Image Extension"
	print xMin, xMax, yMin, yMax, zMin, zMax
	(xSpacing, ySpacing, zSpacing) = change_image.GetOutput().GetSpacing()
	print "Image Spacing"
	print xSpacing, ySpacing, zSpacing
	(x0, y0, z0) = change_image.GetOutput().GetOrigin()
	print "Image Origin"
	print x0, y0, z0
		
	# Set up ortogonal planes
	xImagePlaneWidget.SetInput( change_image.GetOutput() )
	xImagePlaneWidget.SetPlaneOrientationToXAxes()
	yImagePlaneWidget.SetInput( change_image.GetOutput() )
	yImagePlaneWidget.SetPlaneOrientationToYAxes()
	zImagePlaneWidget.SetInput( change_image.GetOutput() )
	zImagePlaneWidget.SetPlaneOrientationToZAxes();
		
	xImagePlaneWidget.SetInteractor( iren1 )
	xImagePlaneWidget.On()
	
	yImagePlaneWidget.SetInteractor( iren1 )
	yImagePlaneWidget.On()

	zImagePlaneWidget.SetInteractor( iren1 )
	zImagePlaneWidget.On()
	
	# set up cube actor with Orientation(R-L, A-P, S-O) using transform_cube
	# Set up to ALS (+X=A, +Y=S, +Z=L) source:
	cube = vtk.vtkAnnotatedCubeActor()
	cube.SetXPlusFaceText( "R" );
	cube.SetXMinusFaceText( "L" );
	cube.SetYPlusFaceText( "P" );
	cube.SetYMinusFaceText( "A" );
	cube.SetZPlusFaceText( "I" );
	cube.SetZMinusFaceText( "S" );
	cube.SetFaceTextScale( 0.5 );
	cube.GetAssembly().SetUserTransform( transform_cube );
		
	# Set UP the axes
	axes2 = vtk.vtkAxesActor()
	axes2.SetShaftTypeToCylinder();
	#axes2.SetUserTransform( transform_cube );		 
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

	# bounds and initialize camera
	b = change_image.GetOutput().GetBounds()
	renderer1.ResetCamera(b)	
	renderer1.ResetCameraClippingRange()
	camera.SetViewUp(0.0,1.0,0.0)
	camera.Azimuth(315)
	
	# Create a text property for both cube axes
	tprop = vtk.vtkTextProperty()
	tprop.SetColor(1, 1, 1)
	tprop.ShadowOff()
	
	# Create a vtkCubeAxesActor2D.  Use the outer edges of the bounding box to
	# draw the axes.  Add the actor to the renderer.
	axes = vtk.vtkCubeAxesActor2D()
	axes.SetInput(change_image.GetOutput())
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
	iren1.Start()
	renderer1.RemoveViewProp(axes)
						
	return  origin, DICOM_mat, zImagePlaneWidget.GetSliceIndex()
		 
def display_pick(DICOM_mat, image, spatial_res, slice_thickn, image_pos_pat, image_ori_pat, dicomReader, zplane):
	global origExt, xMin, xMax, yMin, yMax, zMin, zMax, xSpacing, ySpacing, zSpacing
		
	# Set up the axes	
	transform = vtk.vtkTransform()
	transform.Concatenate(DICOM_mat)
	transform.Update()
	
	# Change info
	change_image = vtk.vtkImageChangeInformation()
	change_image.SetInput( image );
	change_image.SetOutputOrigin(origin[0,0], origin[0,1], origin[0,2])
	change_image.Update()
	
	# Calculate the center of the volume
	change_image.GetOutput().UpdateInformation()
	origExt = change_image.GetOutput().GetWholeExtent()
	(xMin, xMax, yMin, yMax, zMin, zMax) = change_image.GetOutput().GetWholeExtent()
	(xSpacing, ySpacing, zSpacing) = change_image.GetOutput().GetSpacing()
	(x0, y0, z0) = change_image.GetOutput().GetOrigin()
		
	print "Image Origin"
	print x0, y0, z0
	print "Image Extension"
	print xMin, xMax, yMin, yMax, zMin, zMax

	# Visualize results
	xImagePlaneWidget.SetInput( change_image.GetOutput() )
	xImagePlaneWidget.SetPlaneOrientationToXAxes();
	yImagePlaneWidget.SetInput( change_image.GetOutput() )
	yImagePlaneWidget.SetPlaneOrientationToYAxes();
	zImagePlaneWidget.SetInput( change_image.GetOutput() )
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
		file_Seeds.write(str(float(pickPos[0]))+'\t'+str(float(pickPos[1]))+'\t'+str(float(pickPos[2]))+'\n')
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
		
	
def create_SER(selectMethod, numberS, img_folder, LesionID, SeriesID, eID, StudyID, abspath_ExamID, row, dicomInfo, tablefile_name, tablewriter, irow_lesionSegment):
	global file_Seeds, DynVolstack, chosen_lesions_id, dicomInfo_series, abspath_SeriesID, lThre, uThre, spatial_res
	
	# Makedir of current loc
	os.chdir(abspath_ExamID)
	#DynVolstack =  vtk.vtkImageAppendComponents() 
	#DynVolstack.SetNumberOfThreads(4)
	
	############ Get number of lesions and their IDs
	chosen_lesions_id = int(LesionID)
	print chosen_lesions_id
	print "Number of dynamic time points: %i" % numberS
	
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
		
	if(selectMethod == '1'):
		abspath_SeriesID = 'Z:/Cristina/MassNonMass/'+img_folder+StudyID
		# abspath_SeriesID = os.path.abspath(path_SeriesID) -- doesn't work bug
		abspath_lesionID = abspath_SeriesID
		
		choseSerie = raw_input('Enter Series# corresponding to the Ph1-SCAN: ')
		phases_series=[]
		testSID = SeriesID[int(choseSerie)]
		if 'S' in str(testSID):
			print testSID[1:]
			chosen_phase = int(testSID[1:])
		else:
			chosen_phase = int(testSID)
			
		phases_series.append(chosen_phase)
		
		print "Printing available SeriesID"
		print "%s     %s " % ('i', 'Series#')
		
		i=1		
		for chSer in [chosen_phase+1, chosen_phase+2, chosen_phase+3]:
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
			
		###### OBTAIN INFO ABOUT LOCATION
		# Make sure you're in root ExamId dir
		os.chdir(abspath_ExamID)  
		
		# Obtain all locations in all phases_series
		k=0
		for phaseS in phases_series:
			if 'S' in str(testSID):
				phaseS = 'S'+str(phaseS)
				phases_series[k] = phaseS
				k=k+1
					
	if os.path.exists('DynPhases') and selectMethod == '0':
		print '''DynPhases'''
		abspath_SeriesID = 'Z:/Cristina/MassNonMass/'+img_folder+StudyID+os.sep+eID
		abspath_lesionID = abspath_SeriesID+os.sep+'DynPhases'

		choseSerie = raw_input('Enter Series# corresponding to what should be preCont-SCAN: ')
		phases_series=[]
		testSID = SeriesID[int(choseSerie)]
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
					
		###### OBTAIN INFO ABOUT LOCATION
		# Make sure you're in root ExamId dir
		os.chdir(abspath_ExamID)  
		
		# Obtain all locations in all phases_series
		k=0
		for phaseS in phases_series:
			if 'S' in str(testSID):
				phaseS = 'S'+str(phaseS)
				phases_series[k] = phaseS
				k=k+1
			
	if not os.path.exists('DynPhases') and selectMethod == '0':
		print '''SeriesPhases'''
				
		choseSerie = raw_input('Enter Series# corresponding to the ALL-PHASES-SCAN: ')
		path_SeriesID = img_folder+StudyID+os.sep+eID+os.sep+SeriesID[int(choseSerie)]
		print path_SeriesID
		
		# abspath_SeriesID = os.path.abspath(path_SeriesID) -- doesn't work bug
		abspath_SeriesID = 'Z:'+os.sep+'Cristina'+os.sep+'MassNonMass'+os.sep+path_SeriesID
		abspath_lesionID = abspath_SeriesID
		
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
								
		# process all Volumes when in stacks of Dyn Volumes
		# listLocations[int(choseDisplay)] == 'pre-Contrast' ):
		if os.path.exists(abspath_SeriesID+os.sep+'pre-Contrast'):
			phases_series = []
			chosen_folderID = abspath_SeriesID+os.sep+'pre-Contrast'
			os.chdir(chosen_folderID)
			phases_series.append(SeriesID[int(choseSerie)]+os.sep+'pre-Contrast')
					
			print "Arranging series scans"
			for i in range(1,numberS):
				#print i
				chosen_folderID = abspath_SeriesID+os.sep+'post_Contrast-'+str(i)
				os.chdir(chosen_folderID)
				phases_series.append(SeriesID[int(choseSerie)]+os.sep+'post_Contrast-'+str(i))
							
	##############################################################################
	if(selectMethod == '0'):
		###### OBTAIN INFO ABOUT LOCATION
		# Make sure you're in root ExamId dir
		os.chdir(abspath_ExamID)  
		
		abspath_PhaseID = abspath_ExamID+os.sep+str(phases_series[0]) # this will return last element on list (so phase1)
		print abspath_PhaseID
		os.chdir(abspath_PhaseID)
	
		dicomReader  = vtk.vtkDICOMImageReader()
		dicomReader.SetDirectoryName( abspath_PhaseID )
		dicomReader.FileLowerLeftOn()
		dicomReader.Update()
				
		# Get total number of files
		list_files = get_only_filesindirectory(abspath_PhaseID)
		kfile = 0
		path_filenameID = []
		for nameFile in list_files:
			if( nameFile=='DIRCONTENTS.txt'):
				path_filenameID = abspath_PhaseID+os.sep+list_files[kfile+1]			
				path_filenameID = path_filenameID[0:-10]+'001.MR.dcm'
				kfile+=1
		if(path_filenameID==[]):
			path_filenameID = abspath_PhaseID+os.sep+list_files[0]
			path_filenameID = path_filenameID[0:-10]+'001.MR.dcm'
			print path_filenameID
			
		dicomInfo_series = dicom.read_file(os.path.abspath(path_filenameID)) 
				
		slice_thickn = float(dicomInfo_series[0x0018,0x0050].value)  
		
		# Image Position (0020,0032): specifies the x, y, and z coordinates 
		# of the upper left hand corner of the image. This tag specifies the coordinates 
		# of the the first voxel transmitted.
		image_pos_pat = list(dicomInfo_series[0x0020,0x0032].value)
		# Image Orientation (0020,0037): specifies the direction cosines 
		# of the first row and the first column with respect to the patient. 
		# The direction of the axes are defined by the patients orientation 
		# to ensure LPS system ( x-axis increasing to the left hand side of the patient, 
		# y-axis increasing to the posterior side of the patient and z-axis increasing toward
		# the head of the patient )
		image_ori_pat = list(dicomInfo_series[0x0020,0x0037].value)	
		
		#send to display
		[origin, DICOM_mat, zplane] = display(dicomReader.GetOutput(), image_pos_pat, image_ori_pat)
		
		# Change info
		dicomReader  = vtk.vtkDICOMImageReader()
		dicomReader.SetDirectoryName( abspath_PhaseID )
		dicomReader.FileLowerLeftOn()
		dicomReader.Update()
			
		# Get image from reader
		pre_im = dicomReader.GetOutput()
		im_scalars = pre_im.GetPointData().GetScalars()
		dims = pre_im.GetDimensions()
				
		# check wether DicomTag corresponds to Image loaded in vtk
		print "\n... Reading %s" % abspath_PhaseID
		print "VTK Dimensions im.GetDimensions(): %d %d %d" % dims
		spacing = pre_im.GetSpacing()
		print "VTK Spacing im.GetSpacing(): %f %f %f\n" % spacing
		
		print ''' DICOM INFO HEADER '''
		print "%s" % str("Patient Position: ")+str(dicomInfo_series[0x0018,0x5100].value)+'\n'+str("Image Position (patient): ")+str(dicomInfo_series[0x0020,0x0032].value)
		print "%s" % str("Image Orientation (patient): ")+str(dicomInfo_series[0x0020,0x0037].value)+'\n'
		print "%s" % str("SpatialResolution: ")+str(spacing[0])+'\n'
		print "%s" % str("SliceThickness: ")+str(dicomInfo_series[0x0018,0x0050].value)+'\n'+str("SpacingBetweenSlices: ")+str(dicomInfo_series[0x0018,0x0088].value)
		print "%s" % str("SliceLocation: ")+str(dicomInfo_series.SliceLocation)+'\n'
		spatial_res = float(spacing[0])
		
	if(selectMethod == '1'):
		
		abspath_PhaseID = abspath_ExamID+os.sep+str(phases_series[0]) # this will return last element on list (so phase1)
		print abspath_PhaseID
		
		# Get total number of files
		list_files = get_only_filesindirectory(abspath_PhaseID)
		kfile = 0
		path_filenameID = []
		for nameFile in list_files:
			if( nameFile=='DIRCONTENTS.txt'):
				path_filenameID = abspath_PhaseID+os.sep+list_files[kfile+1]			
				path_filenameID = path_filenameID[0:-10]+'001.MR.dcm'
				kfile+=1
		if(path_filenameID==[]):
			path_filenameID = abspath_PhaseID+os.sep+list_files[1]
			path_filenameID = path_filenameID[0:-10]+'001.MR.dcm'
			print path_filenameID
		
		# Get DICOM header information
		dicomInfo_series = dicom.read_file(os.path.abspath(path_filenameID)) 
		
		dicomReader  = vtk.vtkDICOMImageReader()
		dicomReader.SetDirectoryName( abspath_PhaseID )
		dicomReader.FileLowerLeftOn()
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
		#DynVolstack.AddInputConnection(dicomReader.GetOutputPort())
		#DynVolstack.Update()
		
		print ''' DICOM INFO HEADER '''
		print "%s" % str("Patient Position: ")+str(dicomInfo_series[0x0018,0x5100].value)+'\n'+str("Image Position (patient): ")+str(dicomInfo_series[0x0020,0x0032].value)
		print "%s" % str("Image Orientation (patient): ")+str(dicomInfo_series[0x0020,0x0037].value)+'\n'
		print "%s" % str("SpatialResolution: ")+str(spacing[0])+'\n'
		print "%s" % str("SliceThickness: ")+str(dicomInfo_series[0x0018,0x0050].value)+'\n'+str("SpacingBetweenSlices: ")+str(dicomInfo_series[0x0018,0x0088].value)
		print "%s" % str("SliceLocation: ")+str(dicomInfo_series.SliceLocation)+'\n'
		spatial_res = float(spacing[0])
		
		slice_thickn = float(dicomInfo_series[0x0018,0x0050].value)  
		# Image Position (0020,0032): specifies the x, y, and z coordinates 
		# of the upper left hand corner of the image. This tag specifies the coordinates 
		# of the the first voxel transmitted.
		image_pos_pat = list(dicomInfo_series[0x0020,0x0032].value)
		# Image Orientation (0020,0037): specifies the direction cosines 
		# of the first row and the first column with respect to the patient. 
		# The direction of the axes are defined by the patients orientation 
		# to ensure LPS system ( x-axis increasing to the left hand side of the patient, 
		# y-axis increasing to the posterior side of the patient and z-axis increasing toward
		# the head of the patient )
		image_ori_pat = list(dicomInfo_series[0x0020,0x0037].value)
		
		# select z plane
		[origin, DICOM_mat, zplane] = display(im, image_pos_pat, image_ori_pat)
		
	#initialize vcars
	num_voxels = []
	volume_lesion = []
	## create input to table
	datarow_DynSeries = [row[1], row[2], row[0], str(dicomInfo_series[0x0012,0x0040].value), str(dicomInfo_series[0x0008,0x0020].value), str(dicomInfo_series[0x0008,0x103e].value), str(SeriesID[int(choseSerie)]), str(spatial_res), str(dicomInfo_series[0x0018,0x0080].value), str(dicomInfo_series[0x0018,0x0081].value), str(dicomInfo_series[0x0018,0x1314].value), str(dicomInfo_series[0x0018,0x0050].value), str(image_pos_pat), str(image_ori_pat), row[6], row[7], row[8], row[9], row[11], row[12], row[13], chosen_lesions_id]
	
	# Remove pre-contrast from phase sieries
	if(selectMethod == '0'):
		phases_series.pop(0)
	
	for i in range(1,numberS):
		print "Procesing Substraction image %d " % i
		abspath_PhaseID = abspath_ExamID+os.sep+str(phases_series[i-1]) # this will return last element on list (so phase1)
		print abspath_PhaseID

		# Get DICOM header information
		list_files = get_only_filesindirectory(abspath_PhaseID)
		kfile = 0
		path_filenameID = []
		for nameFile in list_files:
			if( nameFile=='DIRCONTENTS.txt'):
				path_filenameID = abspath_PhaseID+os.sep+list_files[kfile+1]			
				path_filenameID = path_filenameID[0:-10]+'001.MR.dcm'
				kfile+=1
		if(path_filenameID==[]):
			path_filenameID = abspath_PhaseID+os.sep+list_files[1]
			path_filenameID = path_filenameID[0:-10]+'001.MR.dcm'
			print path_filenameID
						
		# Get dicom header			
		dicomInfo_series = dicom.read_file(os.path.abspath(path_filenameID)) 
		
		# Change info
		dicomReader  = vtk.vtkDICOMImageReader()
		dicomReader.SetDirectoryName( abspath_PhaseID )
		dicomReader.FileLowerLeftOn()
		dicomReader.Update()
		
		# Change info
		change_image = vtk.vtkImageChangeInformation()
		change_image.SetInput( dicomReader.GetOutput() );
		change_image.SetOutputOrigin(origin[0,0], origin[0,1], origin[0,2])
		change_image.Update()
				
		if(selectMethod == '1'):
			# Get image from reader
			sub_pre = vtk.vtkImageData()
			sub_pre = change_image.GetOutput()
					
		if(selectMethod == '0'):
			sub_preMat = vtk.vtkImageMathematics()
			sub_preMat.SetOperationToSubtract()
			sub_preMat.SetInput1(dicomReader.GetOutput())
			sub_preMat.SetInput2(pre_im)
			sub_preMat.Update()
			
			# Change info
			sub_preChange = vtk.vtkImageChangeInformation()
			sub_preChange.SetInput( sub_preMat.GetOutput() );
			sub_preChange.SetOutputOrigin(origin[0,0], origin[0,1], origin[0,2])
			sub_preChange.Update()
			
			sub_pre = vtk.vtkImageData()
			sub_pre = sub_preChange.GetOutput()
		
		im_scalars = sub_pre.GetPointData().GetScalars()
		dims = sub_pre.GetDimensions()
		
		# check wether DicomTag corresponds to Image loaded in vtk
		print "\n... Reading %s" % abspath_PhaseID
		print "VTK Dimensions im.GetDimensions(): %d %d %d" % dims
		spacing = sub_pre.GetSpacing()
		print "VTK Spacing im.GetSpacing(): %f %f %f\n" % spacing
		
		if(i == 1 and selectMethod == '1'):
			sub_pre_early = vtk.vtkImageData()
			sub_pre_early = change_image.GetOutput()
		
		if(i == 1 and selectMethod == '0'):
			sub_pre_early = vtk.vtkImageData()
			sub_pre_early = sub_preChange.GetOutput()
			
		if(i == numberS-1 and selectMethod == '1'):
			sub_pre_late = vtk.vtkImageData()
			sub_pre_late = change_image.GetOutput()
			
		if(i == numberS-1 and selectMethod == '0'):
			sub_pre_late = vtk.vtkImageData()
			sub_pre_late = sub_preChange.GetOutput()
		
		# Makedir of current loc
		# if current loc folder doesn't exist create it		
		os.chdir(abspath_lesionID)
			
		VOIlesion_ID = 'VOIlesions_id'+str(chosen_lesions_id)
		if not os.path.exists(VOIlesion_ID):
			os.makedirs(VOIlesion_ID)
		
		os.chdir(VOIlesion_ID)
		
		# save seed points into a file for further reference
		file_Seeds = open('Seeds'+str(i)+'_'+str(chosen_lesions_id)+'.txt',"a")
		file_Seeds.write('origin:\t')
		file_Seeds.write(str(float(origin[0,0]))+'\t'+str(float(origin[0,1]))+'\t'+str(float(origin[0,2]))+'\n')
		
		print "Displaying picker for lesion segmentation: %d" % i
		seeds = display_pick(DICOM_mat, sub_pre, spatial_res, slice_thickn, image_pos_pat, image_ori_pat, dicomReader, zplane)
		file_Seeds.close()
		
		# vtkImageThresholdConnectivity will perform a flood fill on an image, given upper and lower pixel intensity 
		# thresholds. It works similarly to vtkImageThreshold, but also allows the user to set seed points to limit
		# the threshold operation to contiguous regions of the image. The filled region, or the "inside", will be passed 
		# through to the output by default, while the "outside" will be replaced with zeros. The scalar type of the output is the same as the input.
		
		image_scalar_range = sub_pre.GetScalarRange() 
		print "Image Scalar Range:"
		print image_scalar_range[0], image_scalar_range[1]
		lThre = image_scalar_range[0]
		uThre = image_scalar_range[1]
		
		# Extract each time point substraction to segment and finally extract SER from first time point subs-ROI
		thresh_sub = vtk.vtkImageThresholdConnectivity() 
		thresh_sub.SetSeedPoints(seeds)
		thresh_sub.SetInValue(1) 
		thresh_sub.SetOutValue(0)
		thresh_sub.ReplaceInOff()
		thresh_sub.ReplaceOutOn()
		thresh_sub.SetNeighborhoodRadius(3, 3, 2) #The radius of the neighborhood that must be within the threshold values in order for the voxel to be included in the mask. The default radius is zero (one single voxel). The radius is measured in voxels
		thresh_sub.SetNeighborhoodFraction(0.10) #The fraction of the neighborhood that must be within the thresholds. The default value is 0.5.
		thresh_sub.ThresholdBetween(0.25*uThre, uThre); 
		thresh_sub.SetInput(sub_pre)
		thresh_sub.Update()
		
		#find out how big lesion is GetNumberOfInVoxels * VoxelVolume = lesion Volume
		num_voxels.append(thresh_sub.GetNumberOfInVoxels())
		print "Number of Voxels #: %d" % num_voxels[i-1] # After the filter has executed, use GetNumberOfVoxels() to find out how many voxels were filled.
				
		# Nothing found
		if( num_voxels[i-1] == 0):
			thresh_sub = vtk.vtkImageThreshold()
			thresh_sub.SetInput(sub_pre)
			thresh_sub.SetOutValue(0)
			thresh_sub.SetInValue(1)
			thresh_sub.ThresholdBetween(0.10*uThre, uThre)
			#num_voxels.append(thresh_sub.GetOutput().GetNumberOfInVoxels
			vol_lesion = float(num_voxels[i-1])*float(spatial_res)*float(spatial_res)*float(dicomInfo_series[0x0018,0x0050].value)
		else:
			# After the filter has executed, use GetNumberOfVoxels() to find out how many voxels were filled.
			vol_lesion = float(num_voxels[i-1])*float(spatial_res)*float(spatial_res)*float(dicomInfo_series[0x0018,0x0050].value)
			
		volume_lesion.append(vol_lesion)
		print "Volume of Lesion mm3: %d" % volume_lesion[i-1]
		
		# Append Num_voxels and vol_lesion to table row
		datarow_DynSeries.append(num_voxels[i-1])
		datarow_DynSeries.append(volume_lesion[i-1])
		
		# Save		
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
			
	print "Doing SER:"
	divide_SER = vtk.vtkImageMathematics()
	divide_SER.SetOperationToDivide()
	divide_SER.SetInput1(sub_pre_early)
	divide_SER.SetInput2(sub_pre_late)
	divide_SER.Update()
			
	# Create an Image of model_Updated
	white_image = vtk.vtkImageData()
	white_image.SetDimensions(dims[0], dims[1], dims[2])
	white_image.SetWholeExtent(xMin, xMax, yMin, yMax, zMin, zMax)
	white_image.SetSpacing(xSpacing, ySpacing, zSpacing)
	white_image.SetOrigin(origin[0,0], origin[0,1], origin[0,2])
	white_image.SetScalarTypeToFloat()
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
	# within a loop of rows returned from a MySQL database (r refers to row number, c to column number, and content is what you want to put in the cell):																									'SliceThickness', 'image_pos_pat', 'image_ori_pat'			
	tablewriter.addRow(irow_lesionSegment, datarow_DynSeries)
		
	# Dealing with  DynVolstack"
	#DynVolstack.RemoveAllInputs()
	
	return irow_lesionSegment

def SetSlice(current_widget):
	xyz = current_widget.GetCurrentCursorPosition()
	x0, y0, z0 = xyz
	pixVal = image.GetScalarComponentAsFloat( x0, y0, z0, 0 )
	
   	return pixVal
  
def get_slices_at_all_locs(img_folder, abspath_SeriesID, len_listSeries_files, StudyID, eID, SeriesID, choseSerie, listSeries_files, abspath_ExamID):
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
		absp_fsID = 'Z:/Cristina/MassNonMass/'+img_folder+StudyID+'/'+eID+'/'+SeriesID[int(choseSerie)]+'/'+listSeries_files[n]
							
		dInfo = dicom.read_file(absp_fsID)
		slices.append(dInfo.SliceLocation)
		FileNms_slices.append(listSeries_files[n])
		FileNms_slices.append(dInfo.SliceLocation)
		
	print '''\nPROCESS STACKS BY SLICE LOCATIONS '''
	FileNms_slices_stack = reshape(FileNms_slices, [len_listSeries_files,2])
	
	FileNms_slices_sorted = sorted(FileNms_slices_stack, key=itemgetter(1))
	#print FileNms_slices_sorted
	#time.sleep(5) #will sleep for 5 seconds
	
	FileNms_slices_sorted_stack = reshape(FileNms_slices_sorted, [len_listSeries_files,2])
	#print FileNms_slices_sorted_stack
	#time.sleep(5) #will sleep for 5 seconds
	
	current_slice = FileNms_slices_sorted_stack[0,1]
	print current_slice
	stack_byLocation = []
	name_byLocation = []
	scount = 0
	
	
	for sliceS in FileNms_slices_sorted_stack:
		# Get the num_series = 5 name_byLocations for a given Location
		if( current_slice == sliceS[1]):
			print "Name: %s" % sliceS[0]
			#print "Slice_loc: %s" % sliceS[1]
			stack_byLocation.append(sliceS[1])
			name_byLocation.append(sliceS[0])
			scount = scount+1
	
		# Finish getting all series for a given Location
		else:
			print scount
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
			file_series = open(filename_series, 'a')
						
			for j in range(scount):
				# list to read as series later
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
	
	# finalize the last location move
	if( current_slice ==  FileNms_slices_sorted_stack[-1,1]):
		print scount
		print '''-----\t NOW HAVE ALL SLICES IN A stack_byLocation '''
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
		
		for j in range(scount):
			# list to read as series later
			link_to = '../'+name_byLocation[j]
			name4link_to = filename+'00'+str(j)+file_ending
			print "linking file: %s to: %s" % (link_to, name4link_to)
			ln_subp = subprocess.Popen(['ln', link_to, name4link_to], stdout=subprocess.PIPE)
			ln_subp.wait()
			file_series.write(str(name4link_to)+'\n')
			
		file_series.close()
			
		# Get back inside the  Series directory
		os.chdir(abspath_SeriesID)
	
	print "\n------------------------------------------------------------------"						
	print '''FINISH get_slices_at_all_locs \n'''	
	os.chdir(abspath_ExamID)
	
	return scount


def get_slices_for_volumes(scount, img_folder, abspath_SeriesID, len_listSeries_files, StudyID, eID, SeriesID, choseSerie, listSeries_files, abspath_ExamID ):
	global num_locs_per_Vol
	
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
		absp_fsID = 'Z:/Cristina/MassNonMass/'+img_folder+StudyID+'/'+eID+'/'+SeriesID[int(choseSerie)]+'/'+listSeries_files[n]
							
		dInfo = dicom.read_file(listSeries_files[n])
		slices.append(dInfo.SliceLocation)
		FileNms_slices.append(listSeries_files[n])
		FileNms_slices.append(dInfo.SliceLocation)
	
	print '''\nPROCESS STACKS BY SLICE LOCATIONS '''
	print "Total number of locations found: %d" % scount
	FileNms_slices_stack = reshape(FileNms_slices, [len_listSeries_files,2])
	#print FileNms_slices_stack
	FileNms_slices_sorted = sorted(FileNms_slices_stack, key=itemgetter(0))
	FileNms_slices_sorted_stack = reshape(FileNms_slices_sorted, [len_listSeries_files,2])
	print   FileNms_slices_sorted_stack

	stack_byLocation = []
	name_byLocation = []
	
	# Extracting the number of slices per bilateral 3D volumes based on num_series for 280/5 = 56
	num_locs_per_Vol = int(float(len_listSeries_files)/float(scount))
	scount = int(scount)
	
	print "------\tNumber of Dynamic Volumes (series time points: including pre-contrast): %d" % scount 
	print "------\tNumber of Locations per Bilateral Volume: %d " % num_locs_per_Vol
	stack_byBilatVol = []
		
	k=0
	# Get the folder names based on num_series
	for numlocs in range(scount):
		if ( numlocs == 0):	
			stack_byBilatVol.append('pre-Contrast')
		else:# Now link slices at location to folder
			stack_byBilatVol.append('post_Contrast-'+str(k))
		k=k+1
			
	# Initialized BilatVol
	print stack_byBilatVol
	slice_i = 0
	
	for Vol_i in range(scount):	
		print'''-----\t NOW HAVE ALL SLICES IN A bilateral vol '''
		# Get the new Location folder
		print "bit vol %d" % Vol_i
		current_vol = stack_byBilatVol[Vol_i]
								
		# Makedir of current loc
		# if current loc folder doesn't exist create it		
		if not os.path.exists(str(current_vol)):
			os.makedirs(str(current_vol))
		
		# Get inside location directory
		os.chdir(str(current_vol))
		print os.getcwd()
		
		#FileNms_slices_sorted_stack[slice_i,0]
		# Now link slices at location to folder
		filename = str(FileNms_slices_sorted_stack[slice_i,0])
		filename = filename[0:-10]
		file_ending = '.MR.dcm'
		
		# Save the file list to read as series later
		filename_series = 'DIRCONTENTS.txt'
		file_series = open(filename_series, 'a')
		
		for j in range(num_locs_per_Vol):
			link_to = '../'+FileNms_slices_sorted_stack[slice_i,0]
			if ( j < 9):
				name4link_to = filename+'00'+str(j+1)+file_ending
			else:
				name4link_to = filename+'0'+str(j+1)+file_ending
			print "linking file: %s to: %s" % (link_to, name4link_to)
			ln_subp = subprocess.Popen(['ln', link_to, name4link_to], stdout=subprocess.PIPE)
			ln_subp.wait()
			file_series.write(str(name4link_to)+'\n')
			slice_i = slice_i+1
			
		file_series.close()
			
		# Get back inside the  Series directory
		os.chdir(abspath_SeriesID)
						
	print "\n------------------------------------------------------------------"						
	print '''FINISH get_slices_for_volumes \n'''	
	os.chdir(abspath_ExamID)
	
	return scount
	
	
# ******************************************************************************
# initialize some variables
irow_lesionSegment = 1

########## Write to table
#initialize Table writer
tablefile_name = "table_SERlesion_MASS.xls" 
tablewriter = excel_writer.ExcelWriter(tablefile_name, "LesionSegment", make_visible=True) 
title_range = tablewriter.getRangeByCells((1, 1), (1, 30)) #where j = number of columns 
tablewriter.formatRange(title_range, excel_writer.STYLE_GREY_CELL)

# start by writing headings 
headings_DynSeries = ['StudyNumber', 'DicomExamNumber', 'MRN', 'ClinicalTrialID', 'StudyDate', 'SeriesDescription', 'SeriesID', 'SpatialResolution', 'TR', 'TE', 'FlipAngle', 'SliceThickness', 'image_pos_pat', 'image_ori_pat', 'BreastSide', 'QuadrantIQ', 'ClockStart', 'ClockEnd', 'LesionCenterPixelX', 'LesionCenterPixelY', 'LesionFindingType.Description', 'chosen_lesions_id', 'num_vox_sub1', 'volume_lesion_sub1', 'num_vox_sub2', 'volume_lesion_sub2', 'num_vox_sub3', 'volume_lesion_sub3','num_vo_sub4', 'volume_lesion_sub4']
tablewriter.addRow(irow_lesionSegment, headings_DynSeries)

############
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
						
		# *****************************************************************************
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
		
		# Set Up Camera view
		renderer1.SetBackground(0.0, 0.0, 0.0)
		camera = renderer1.GetActiveCamera()
		iren1.SetPicker(picker)

		# ask for which series to load
		print "\n----------------------------------------------------------"
		choseTodo = raw_input('Enter option to run: \n\t1) Visualize_location \n\t2) get_slices_at_all_locs \n\t3) get_slices_for_volumes  \n\t4) create_SERmap 	\n\tx) to exit \n\t input: ')
		while(choseTodo != 'x'):
			if(choseTodo == '1'):	
				[abspath_ExamID, eID, SeriesID, studyFolder, dicomInfo] = get_series(StudyID, img_folder)
 				
				visual_SerieID = raw_input('Enter Series# to load or x to exit \n\t input: ')
				chosen_folderID = abspath_ExamID+os.sep+str(SeriesID[int(visual_SerieID)])
				print chosen_folderID
				os.chdir(chosen_folderID)
				
				# Get DICOM header information
				list_files = get_only_filesindirectory(chosen_folderID)
				kfile = 0
				path_filenameID = []
				for nameFile in list_files:
					if( nameFile=='DIRCONTENTS.txt'):
						path_filenameID = chosen_folderID+os.sep+list_files[kfile+1]			
						path_filenameID = path_filenameID[0:-10]+'001.MR.dcm'
						kfile+=1
				if(path_filenameID==[]):
					path_filenameID = chosen_folderID+os.sep+list_files[1]
					path_filenameID = path_filenameID[0:-10]+'001.MR.dcm'
					print path_filenameID
				
				dicomInfo_visual_SerieID = dicom.read_file(os.path.abspath(path_filenameID)) 
				
				# Get position and orientation				
				image_pos_pat_visual_SerieID = list(dicomInfo_visual_SerieID[0x0020,0x0032].value) 
				image_ori_pat_visual_SerieID = list(dicomInfo_visual_SerieID[0x0020,0x0037].value)
				
				filename_series = 'DIRCONTENTS.txt'
				#check firs if it exists
				if(FileCheck(filename_series) ):
					vtkStringArray = vtk.vtkStringArray()
					file_series = open(filename_series,'r')
					files_toRead = file_series.read()
					try:
						for slicename in files_toRead.split():
							print "adding : %s" %  slicename
							vtkStringArray.InsertNextValue( slicename )
					finally:
						file_series.close()
					
					# Read dicom Vol from DIRCONTENTS.txt
					dicomReader  = vtkgdcmPython.vtkGDCMImageReader()
					dicomReader.SetFileNames( vtkStringArray )
					dicomReader.FileLowerLeftOn()
					dicomReader.Update()
				else:
					dicomReader  = vtk.vtkDICOMImageReader()
					dicomReader.SetDirectoryName( chosen_folderID )
					dicomReader.FileLowerLeftOn()
					dicomReader.Update()
			
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
							
				[origin, DICOM_mat, z_visualSerie] = display(dicomReader.GetOutput(), image_pos_pat_visual_SerieID, image_ori_pat_visual_SerieID)
				
				print " "
				chdirname='Z:/Cristina/MassNonmass/'
				os.chdir(chdirname)  
				print os.getcwd()
				
			if(choseTodo == '2'):	
				print '''\n get_slices_at_all_locs\n''' 
				choseSerie = raw_input('Enter n Series to load (0-n), or x to exit: ')
				# Enter folder
				os.chdir(abspath_ExamID)
				
				path_SeriesID = img_folder+StudyID+'/'+eID+'/'+SeriesID[int(choseSerie)]
				print path_SeriesID
								
				# abspath_SeriesID = os.path.abspath(path_SeriesID) -- doesn't work bug
				abspath_SeriesID = 'Z:/Cristina/MassNonMass/'+path_SeriesID
				
				# Get total number of files
				listSeries_files = sorted(list(get_only_filesindirectory(abspath_SeriesID)))
				
				print "Total images in series: %d " % len(listSeries_files)
				len_listSeries_files = len(listSeries_files)
				
				#!!!!!!!!!!!!!!!!!!!!! IMPORTANT: scount CONTAINS THE NUMBER OF SLICE LOCATIONS (e.g 56)
				scount = get_slices_at_all_locs(img_folder, abspath_SeriesID, len_listSeries_files, StudyID, eID, SeriesID, choseSerie, listSeries_files, abspath_ExamID )	
			
			if(choseTodo == '3'):	
				print '''\nget_slices_for_volumes\n'''  
				
				choseSerie = raw_input('Enter n Series to load (0-n), or x to exit: ')
				# Enter folder
				os.chdir(abspath_ExamID)
				
				path_SeriesID = img_folder+StudyID+'/'+eID+'/'+SeriesID[int(choseSerie)]
				print path_SeriesID
								
				# abspath_SeriesID = os.path.abspath(path_SeriesID) -- doesn't work bug
				abspath_SeriesID = 'Z:/Cristina/MassNonMass/'+path_SeriesID
				
				# Get total number of files
				listSeries_files = sorted(list(get_only_filesindirectory(abspath_SeriesID)))
								
				print "Total images in series: %d " % len(listSeries_files)
				len_listSeries_files = len(listSeries_files)
				
				numberS = get_slices_for_volumes(scount, img_folder, abspath_SeriesID, len_listSeries_files, StudyID, eID, SeriesID, choseSerie, listSeries_files, abspath_ExamID )
			
			if(choseTodo == '4'):
				print '''\nCREATE SER MAP \n''' 
				irow_lesionSegment = irow_lesionSegment+1
				#print irow_lesionSegment
				
				selectMethod = raw_input('Enter 1 for Manual selection or 0 for Automatic(ALL-PHASES vs DynPhases): \t')
				irow_lesionSegment = create_SER(selectMethod, numberS, img_folder, LesionID, SeriesID, eID, StudyID, abspath_ExamID, row, dicomInfo, tablefile_name, tablewriter, irow_lesionSegment)
				
				# *****************************************************************************
				chdirname='Z:/Cristina/MassNonmass/'
				os.chdir(chdirname)  
				print os.getcwd()
				numberS=5				
			
			choseTodo = raw_input('Enter option to run: \n\t1) Visualize_location \n\t2) get_slices_at_all_locs \n\t3) get_slices_for_volumes  \n\t4) create_SERmap 	\n\tx) to exit \n\t input: ')
			print " "
			chdirname='Z:/Cristina/MassNonmass/'
			os.chdir(chdirname)  
			print os.getcwd()
			
			
		## exit ChoseTodo								
								
finally:
	file_ids.close()		
		
			
# Save the file list to read as series later
# Finally, save the file and clean up:
tablewriter.save()		
tablewriter.close()		
		
		
