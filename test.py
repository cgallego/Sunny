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

import gdcm
import vtkgdcm
import vtkgdcmPython

def FileCheck(filename):       
       try:
           fn=open(filename,"r") 
           fn.close()
           return True
       except IOError: 
           print "Error: DICOMDIR.txt doesn't exit, loading Series using vtkDICOMImageReader by setDirectoryName."
       return False
       
def get_only_filesindirectory(mydir):
     return [name for name in os.listdir(mydir) 
            if os.path.isfile(os.path.join(mydir, name))]
 
def get_immediate_subdirectories(mydir):
    return [name for name in os.listdir(mydir) 
            if os.path.isdir(os.path.join(mydir, name))]

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
	global xMin, xMax, yMin, yMax, zMin, zMax, xSpacing, ySpacing, zSpacing, interactor, actions, reslice, interactorStyle
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
	
	IP =  array([0, 0, 0, 1]) # Initialization Image Position
	IP[0] = image_pos_pat[0]; IP[1] = image_pos_pat[1]; IP[2] = image_pos_pat[2];  
	IO[0,3] = image_pos_pat[0]; IO[1,3] = image_pos_pat[1]; IO[2,3] = image_pos_pat[2]
	
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
	
	##########################
	# Calculate the center of the volume
	(xMin, xMax, yMin, yMax, zMin, zMax) = image.GetWholeExtent()
	(xSpacing, ySpacing, zSpacing) = image.GetSpacing()
	(x0, y0, z0) = image.GetOrigin()
	
	center = [x0 + xSpacing * 0.5 * (xMin + xMax),
		  y0 + ySpacing * 0.5 * (yMin + yMax),
		  z0 + zSpacing * 0.5 * (zMin + zMax)]
	
	# Matrices for axial, coronal, sagittal, oblique view orientations
	axial = vtk.vtkMatrix4x4()
	axial.DeepCopy((1, 0, 0, center[0],
			0, 1, 0, center[1],
			0, 0, 1, center[2],
			0, 0, 0, 1))
	
	coronal = vtk.vtkMatrix4x4()
	coronal.DeepCopy((1, 0, 0, center[0],
			  0, 0, 1, center[1],
			  0,-1, 0, center[2],
			  0, 0, 0, 1))
	
	sagittal = vtk.vtkMatrix4x4()
	sagittal.DeepCopy((0, 0,-1, center[0],
			   1, 0, 0, center[1],
			   0,-1, 0, center[2],
			   0, 0, 0, 1))
	
	oblique = vtk.vtkMatrix4x4()
	oblique.DeepCopy((1, 0, 0, center[0],
			  0, 0.866025, -0.5, center[1],
			  0, 0.5, 0.866025, center[2],
			  0, 0, 0, 1))
	
	# Extract a slice in the desired orientation
	reslice = vtk.vtkImageReslice()
	reslice.SetInput(image)
	reslice.SetOutputDimensionality(2)
	reslice.SetResliceAxes(sagittal)
	reslice.SetInterpolationModeToLinear()
	
	# Create a greyscale lookup table
	table = vtk.vtkLookupTable()
	table.SetRange(0, 2000) # image intensity range
	table.SetValueRange(0.0, 1.0) # from black to white
	table.SetSaturationRange(0.0, 0.0) # no color saturation
	table.SetRampToLinear()
	table.Build()
	
	# Map the image through the lookup table
	color = vtk.vtkImageMapToColors()
	color.SetLookupTable(table)
	color.SetInputConnection(reslice.GetOutputPort())

	# Display the image
	actor = vtk.vtkImageActor()
	actor.GetMapper().SetInputConnection(color.GetOutputPort())

	renderer1 = vtk.vtkRenderer()
	renderer1.AddActor(actor)
	################
	    
	
	# set up cube actor with Orientation(R-L, A-P, S-O) using transform_cube
	# Set up to ALS (+X=A, +Y=S, +Z=L) source:
	cube = vtk.vtkAnnotatedCubeActor()
	cube.SetXPlusFaceText( "S" );
	cube.SetXMinusFaceText( "I" );
	cube.SetYPlusFaceText( "L" );
	cube.SetYMinusFaceText( "R" );
	cube.SetZPlusFaceText( "A" );
	cube.SetZMinusFaceText( "P" );
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
	# Set up the interaction
	interactorStyle = vtk.vtkInteractorStyleImage()
	interactor = vtk.vtkRenderWindowInteractor()
	interactor.SetInteractorStyle(interactorStyle)
	renWin1.SetInteractor(interactor)
	renWin1.Render()
	
	renderer1.Render()
	interactor.Start()
	renderer1.RemoveViewProp(axes)
						
	return transform_cube, zImagePlaneWidget.GetSliceIndex()
		 

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
		
		# *****************************************************************************
		# get_series and print which series to load
		print " "
		chdirname='Z:/Cristina/MassNonmass/'
		os.chdir(chdirname)  
		print os.getcwd()
		
		[abspath_ExamID, eID, SeriesID, studyFolder, dicomInfo] = get_series(StudyID, img_folder)
						
		# *****************************************************************************
		# Create a renderer, render window, and render window interactor to
		# display the results.
		renderer1 = vtk.vtkRenderer()
		renWin1 = vtk.vtkRenderWindow()
		iren1 = vtk.vtkRenderWindowInteractor()
		
		renWin1.SetSize(1000, 800);
		renWin1.AddRenderer(renderer1);
		iren1.SetRenderWindow(renWin1);
		
		
		# Set Up Camera view
		renderer1.SetBackground(0.0, 0.0, 0.0)
		camera = renderer1.GetActiveCamera()
					
 				
		visual_SerieID = raw_input('Enter Series# to load or x to exit \n\t input: ')
		chosen_folderID = abspath_ExamID+os.sep+str(SeriesID[int(visual_SerieID)])
		print chosen_folderID
		os.chdir(chosen_folderID)
		
		# Get total number of files
		list_files = get_only_filesindirectory(chosen_folderID)
		path_filenameID = chosen_folderID+os.sep+list_files[1]										
		dicomInfo_visual_SerieID = dicom.read_file(os.path.abspath(path_filenameID)) 
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
			dicomReader.SetDirectoryName( chosen_folderID ) # neccesarily, most VTK image readers flip the images from their original top-to-bottom rasterization to a bottom-to-top rasterization instead.
			dicomReader.FileLowerLeftOn()
			dicomReader.Update()
	
		# display the results.
		renderer1 = vtk.vtkRenderer()
		renWin1 = vtk.vtkRenderWindow()
		iren1 = vtk.vtkRenderWindowInteractor()
		
		renWin1.SetSize(1000, 800);
		renWin1.AddRenderer(renderer1);
		iren1.SetRenderWindow(renWin1);
		
							
		[transform_cube, z_visualSerie] = display(dicomReader.GetOutput(), image_pos_pat_visual_SerieID, image_ori_pat_visual_SerieID)
		
								
finally:
	file_ids.close()		
		
			
		
		
