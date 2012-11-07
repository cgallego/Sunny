#!/usr/bin/env python

import os
import os.path
import sys
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
from dicom.contrib.pydicom_PIL import show_PIL

import gdcm
import vtkgdcm
import vtkgdcmPython

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
		1) nth-time-points per Location folders
		2) bilateral post-contrast Volumes when phases are not separated. (i.e when Ph1/Ph2/Ph3.. is unavailable)
% Aug 22 -      Added function get_locs_for_all_phases() 
		3) Loopps through each one of the Ph1/Ph2/Ph3/Ph4 series and finds equal locations at 4 post-contrast timepoints in the T13D Gad Dynamic sequences by
		
-----------------------------------------
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
		path_arrangedF_ID = abspath_SeriesID+'/'+arrangedF
		#print path_SeriesID
		
		# Get total number of files
		listSeries_files = get_only_filesindirectory(path_arrangedF_ID)
					
		# Use only the first slice one file and get DICOM DICTIONARY
		if(listSeries_files != []):
			path_filenameID = abspath_SeriesID+'/'+arrangedF+'/'+listSeries_files[0]
						
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
	
def get_series(StudyID,img_folder):
	# obtain subdires in the StudyID directory, correspond to ExamsID 
	# check for one in subdirectory tree, (e.g Mass or NonMass) 
	
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
		path_ExamID = img_folder+StudyID+'/'+eID
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

			path_SeriesID = img_folder+StudyID+'/'+eID+'/'+sID
							
			abspath_SeriesID = os.path.abspath(path_SeriesID)
			
			# Get total number of files
			listSeries_files = get_only_filesindirectory(abspath_SeriesID)
						
			# Use only the first slice one file and get DICOM DICTIONARY
			path_filenameID = img_folder+StudyID+'/'+eID+'/'+sID+'/'+listSeries_files[0]
						
			dicomInfo = dicom.read_file(os.path.abspath(path_filenameID)) 
					
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
			path_ExamID = img_folder+StudyID+'/'+eID
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
	
				path_SeriesID = img_folder+StudyID+'/'+eID+'/'+sID
				
				abspath_SeriesID = os.path.abspath(path_SeriesID)
				#print abspath_SeriesID
				
				# Get total number of files
				listSeries_files = get_only_filesindirectory(abspath_SeriesID)
			
				# Use only the first slice one file and get DICOM DICTIONARY
				if(listSeries_files != []):
					path_filenameID = img_folder+StudyID+'/'+eID+'/'+sID+'/'+listSeries_files[0]
					dicomInfo = dicom.read_file(os.path.abspath(path_filenameID)) 
										
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
				
	return abspath_ExamID, eID, SeriesID, studyFolder
		
			
def get_slices_at_all_locs(img_folder, abspath_SeriesID, len_listSeries_files, StudyID, eID, SeriesID, choseSerie, listSeries_files, num_series, abspath_ExamID ):
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
								
	print "\n------------------------------------------------------------------"						
	print '''FINISH get_slices_at_all_locs \n'''	
	os.chdir(abspath_ExamID)
	return

def get_slices_for_volumes(img_folder, abspath_SeriesID, len_listSeries_files, StudyID, eID, SeriesID, choseSerie, listSeries_files, num_series, abspath_ExamID ):
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
	
	
	print '''\nJust get the number of diff SLICE LOCATIONS '''
	FileNms_slices_stack = reshape(FileNms_slices, [len_listSeries_files,2])
	FileNms_slices_sorted = sorted(FileNms_slices_stack, key=itemgetter(0))
	#print FileNms_slices_sorted
	#time.sleep(5) #will sleep for 5 seconds
	
	FileNms_slices_sorted_stack = reshape(FileNms_slices_sorted, [len_listSeries_files,2])
	#print FileNms_slices_sorted_stack
	#time.sleep(5) #will sleep for 5 seconds
	
	current_slice = FileNms_slices_sorted_stack[0,1]
	print current_slice
	num_locations = 0
	stack_num_locations = []
	stack_num_locations.append(current_slice)
	
	for sliceS in FileNms_slices_sorted_stack:
		# Get the num_series = 5 name_byLocations for a given Location
		while( sliceS[1] not in stack_num_locations ):
			print "Name: %s" % sliceS[0]
			print "Slice_loc: %s" % sliceS[1]
			num_locations = num_locations+1
			stack_num_locations.append(sliceS[1])
			
		# Finish getting all series for a given Location
		else:
			break;
				
	print '''\nPROCESS STACKS BY SLICE LOCATIONS '''
	print "Total number of locations found: %d" % num_locations
	FileNms_slices_stack = reshape(FileNms_slices, [len_listSeries_files,2])
	FileNms_slices_sorted = sorted(FileNms_slices_stack, key=itemgetter(1))
	FileNms_slices_sorted_stack = reshape(FileNms_slices_sorted, [len_listSeries_files,2])
	
	stack_byLocation = []
	name_byLocation = []
	scount = 0
	
	# Extracting the number of slices per bilateral 3D volumes based on num_series for 280/5 = 56
	numSlices_perBilat_vol = len_listSeries_files/num_series
	print "------\tNumber of Slices per Bilateral Volume: %d " % numSlices_perBilat_vol
	stack_byBilatVol = []
	
	# Get the folder names based on num_series
	if ( len_listSeries_files%num_series != 0):	
		stack_byBilatVol.append('pre-Contrast')
		
		for k in range(1,num_series):
			# Now link slices at location to folder
			stack_byBilatVol.append('post_Contrast-'+str(k))
	
	if ( len_listSeries_files%num_series == 0):
		for k in range(num_series):
			# Now link slices at location to folder
			stack_byBilatVol.append('post_Contrast-'+str(k+1))
		
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

def get_locs_for_all_phases(img_folder, SeriesID, choseSerie, StudyID, eID, abspath_ExamID, len_listSeries_files, listSeries_files ):
	# enter folder of Series selection /examID/S_XXX 
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
		i=1
	
	###### OBTAIN INFO ABOUT LOCATION
	# Make sure you're in root ExamId dir
	os.chdir(abspath_ExamID)
	scountPhase = 1
	
	# Obtain all locations in all phases_series
	for phaseS in phases_series:
		if 'S' in str(testSID):
			phaseS = 'S'+str(phaseS)
			
		abspath_PhaseID = abspath_ExamID+'/'+str(phaseS) # this will return last element on list (so phase1)
		print abspath_PhaseID
	
		# Get total number of files
		listPhase_files = sorted(list(get_only_filesindirectory(abspath_PhaseID)))
			
		print "Total images in phase: %d " % len(listPhase_files)
		len_listPhase_files = len(listPhase_files)
	
		slices_phase = []
		FileNms_phase_slices =  []
	
		for n in range(len_listPhase_files):
			# Use all DICOM slices on series
			''' EXTRACT DICOM SLICE LOCATION IN PHASEs '''
			abspath_PhaseLocID = abspath_PhaseID+'/'+listPhase_files[n]
								
			dInfo_phase = dicom.read_file(abspath_PhaseLocID)
			slices_phase.append(dInfo_phase.SliceLocation)
			FileNms_phase_slices.append(listPhase_files[n])
			FileNms_phase_slices.append(dInfo_phase.SliceLocation)
		
		print '''\nPROCESS STACKS BY SLICE LOCATIONS '''
		FileNms_phase_slices_stack = reshape(FileNms_phase_slices, [len_listPhase_files,2])
		FileNms_phase_slices_sorted = sorted(FileNms_phase_slices_stack, key=itemgetter(1))
		FileNms_slices_sorted_stack = reshape(FileNms_phase_slices_sorted, [len_listPhase_files,2])
			
		# change to next phase
		os.chdir(abspath_PhaseID)
		
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
				'''-----\t NOW HAVE ALL SLICES IN A stack_byLocation '''
				# Get the new Location folder
				# To verify "Above lengths should be = num_series + 1" len(stack_byLocation) len(name_byLocation)'''				
				current_loc = stack_byLocation[0]
										
				# Makedir of current loc
				os.chdir(abspath_ExamID)
				if not os.path.exists('DynPhases'):
					os.makedirs('DynPhases')
	
				os.chdir('DynPhases')	
				print os.getcwd()
				
				# if current loc folder doesn't exist create it		
				if not os.path.exists(str(current_loc)):
					os.makedirs(str(current_loc))
															
				# Get inside location directory
				os.chdir(str(current_loc))
				
				# Save the file list to read as series later
				filename_series = 'DIRCONTENTS.txt'
				file_series = open(filename_series, 'a')
				#print os.getcwd()
				
				# Now link slices at location to folder
				filename = str(name_byLocation[0])
				file_ending = '.MR.dcm'
									
				# list to read as series later
				for j in range(scount):
					link_to = '../../'+str(phaseS)+'/'+name_byLocation[j]
					name4link_to = '00'+str(scountPhase)+file_ending
					print "linking file: %s to: %s" % (link_to, name4link_to)
					ln_subp = subprocess.Popen(['ln', link_to, name4link_to], stdout=subprocess.PIPE)
					ln_subp.wait()
					file_series.write(str(name4link_to)+'\n')
													
				# Get back inside the  Series directory
				os.chdir(abspath_ExamID)
				
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
				
		scountPhase = scountPhase+1
		file_series.close()
		# Makedir of current loc
		os.chdir(abspath_ExamID)
		
	print "\n------------------------------------------------------------------"						
	print '''FINISH get_locs_for_all_phases \n'''	
	os.chdir('DynPhases')
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
                      Pathology.PathoResultType.PathoResultTypeID, Pathology.PathoResultType.Description AS Expr5, 
                      Radiology.NMassLesion.LesionFindingID, Lesion.NMassEnhancingType.NMassEnhanceTypeID, Lesion.NMassEnhancingType.Description AS Expr2, 
                      Lesion.NMassIntEnhType.NMassIntEnhanceTypeID, Lesion.NMassIntEnhType.Description AS Expr3, Lesion.NMassEnhSymmetryType.NMassEnhSymmTypeID, 
                      Lesion.NMassEnhSymmetryType.Description AS Expr4
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
                      Radiology.NMassLesion ON Radiology.Lesions.LesionFindingID = Radiology.NMassLesion.LesionFindingID AND 
                      Radiology.Lesions.LesionFindingID = Radiology.NMassLesion.LesionFindingID AND 
                      Radiology.Lesions.LesionFindingID = Radiology.NMassLesion.LesionFindingID AND 
                      Radiology.Lesions.LesionFindingID = Radiology.NMassLesion.LesionFindingID INNER JOIN
                      Lesion.NMassEnhancingType ON Radiology.NMassLesion.NMassEnhanceTypeID = Lesion.NMassEnhancingType.NMassEnhanceTypeID AND 
                      Radiology.NMassLesion.NMassEnhanceTypeID = Lesion.NMassEnhancingType.NMassEnhanceTypeID AND 
                      Radiology.NMassLesion.NMassEnhanceTypeID = Lesion.NMassEnhancingType.NMassEnhanceTypeID AND 
                      Radiology.NMassLesion.NMassEnhanceTypeID = Lesion.NMassEnhancingType.NMassEnhanceTypeID AND 
                      Radiology.NMassLesion.NMassEnhanceTypeID = Lesion.NMassEnhancingType.NMassEnhanceTypeID AND 
                      Radiology.NMassLesion.NMassEnhanceTypeID = Lesion.NMassEnhancingType.NMassEnhanceTypeID AND 
                      Radiology.NMassLesion.NMassEnhanceTypeID = Lesion.NMassEnhancingType.NMassEnhanceTypeID AND 
                      Radiology.NMassLesion.NMassEnhanceTypeID = Lesion.NMassEnhancingType.NMassEnhanceTypeID AND 
                      Radiology.NMassLesion.NMassEnhanceTypeID = Lesion.NMassEnhancingType.NMassEnhanceTypeID AND 
                      Radiology.NMassLesion.NMassEnhanceTypeID = Lesion.NMassEnhancingType.NMassEnhanceTypeID AND 
                      Radiology.NMassLesion.NMassEnhanceTypeID = Lesion.NMassEnhancingType.NMassEnhanceTypeID AND 
                      Radiology.NMassLesion.NMassEnhanceTypeID = Lesion.NMassEnhancingType.NMassEnhanceTypeID INNER JOIN
                      Lesion.NMassEnhSymmetryType ON Radiology.NMassLesion.NMassEnhSymmTypeID = Lesion.NMassEnhSymmetryType.NMassEnhSymmTypeID AND 
                      Radiology.NMassLesion.NMassEnhSymmTypeID = Lesion.NMassEnhSymmetryType.NMassEnhSymmTypeID AND 
                      Radiology.NMassLesion.NMassEnhSymmTypeID = Lesion.NMassEnhSymmetryType.NMassEnhSymmTypeID AND 
                      Radiology.NMassLesion.NMassEnhSymmTypeID = Lesion.NMassEnhSymmetryType.NMassEnhSymmTypeID AND 
                      Radiology.NMassLesion.NMassEnhSymmTypeID = Lesion.NMassEnhSymmetryType.NMassEnhSymmTypeID AND 
                      Radiology.NMassLesion.NMassEnhSymmTypeID = Lesion.NMassEnhSymmetryType.NMassEnhSymmTypeID AND 
                      Radiology.NMassLesion.NMassEnhSymmTypeID = Lesion.NMassEnhSymmetryType.NMassEnhSymmTypeID AND 
                      Radiology.NMassLesion.NMassEnhSymmTypeID = Lesion.NMassEnhSymmetryType.NMassEnhSymmTypeID AND 
                      Radiology.NMassLesion.NMassEnhSymmTypeID = Lesion.NMassEnhSymmetryType.NMassEnhSymmTypeID AND 
                      Radiology.NMassLesion.NMassEnhSymmTypeID = Lesion.NMassEnhSymmetryType.NMassEnhSymmTypeID AND 
                      Radiology.NMassLesion.NMassEnhSymmTypeID = Lesion.NMassEnhSymmetryType.NMassEnhSymmTypeID AND 
                      Radiology.NMassLesion.NMassEnhSymmTypeID = Lesion.NMassEnhSymmetryType.NMassEnhSymmTypeID INNER JOIN
                      Lesion.NMassIntEnhType ON Radiology.NMassLesion.NMassIntEnhanceTypeID = Lesion.NMassIntEnhType.NMassIntEnhanceTypeID AND 
                      Radiology.NMassLesion.NMassIntEnhanceTypeID = Lesion.NMassIntEnhType.NMassIntEnhanceTypeID AND 
                      Radiology.NMassLesion.NMassIntEnhanceTypeID = Lesion.NMassIntEnhType.NMassIntEnhanceTypeID AND 
                      Radiology.NMassLesion.NMassIntEnhanceTypeID = Lesion.NMassIntEnhType.NMassIntEnhanceTypeID AND 
                      Radiology.NMassLesion.NMassIntEnhanceTypeID = Lesion.NMassIntEnhType.NMassIntEnhanceTypeID AND 
                      Radiology.NMassLesion.NMassIntEnhanceTypeID = Lesion.NMassIntEnhType.NMassIntEnhanceTypeID AND 
                      Radiology.NMassLesion.NMassIntEnhanceTypeID = Lesion.NMassIntEnhType.NMassIntEnhanceTypeID AND 
                      Radiology.NMassLesion.NMassIntEnhanceTypeID = Lesion.NMassIntEnhType.NMassIntEnhanceTypeID AND 
                      Radiology.NMassLesion.NMassIntEnhanceTypeID = Lesion.NMassIntEnhType.NMassIntEnhanceTypeID AND 
                      Radiology.NMassLesion.NMassIntEnhanceTypeID = Lesion.NMassIntEnhType.NMassIntEnhanceTypeID AND 
                      Radiology.NMassLesion.NMassIntEnhanceTypeID = Lesion.NMassIntEnhType.NMassIntEnhanceTypeID AND 
                      Radiology.NMassLesion.NMassIntEnhanceTypeID = Lesion.NMassIntEnhType.NMassIntEnhanceTypeID
                WHERE Radiology.RadiologyReport.DicomExamNumber=%d     
			''', int(ExamID))
		for row in cur:
			print  row # "\tStudyNumber %d\t DicomExamNumber %d\t Description: %s\n AdditonalDescription: %s\n BreastSide: %s\t ClockStart: %s\t ClockEnd: %s\t MRISeries %s\t\n LesionCenterPixelX: %s\t LesionCenterPixelY: %s\t LesionFindingType.Description: %s\t\n\n" % (row[0], row[1], row[3], row[4], row[5], row[7], row[8], row[9], row[10], row[11], row[12])
			print len(row)
		if(row==[]):
			# initialize database vector
			row = [0 for x in range(26)]
			print  row[16]
			print len(row)
			
		conn.close()
		
		# *****************************************************************************
		# get_series and print which series to load
		print " "
		print "Current folder time"
		chdirname='Z:/Cristina/MassNonmass/'
		os.chdir(chdirname)  
		print os.getcwd()
		
		[abspath_ExamID, eID, SeriesID, studyFolder] = get_series(StudyID, img_folder)
		
		# *****************************************************************************
		# ask for which series to load
		print "\n----------------------------------------------------------"
		choseSerie = raw_input('Enter n Series to load (0-n), or x to exit: ')
		
		while(choseSerie!='x'):
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
			#time.sleep(3) #will sleep for 5 seconds
			
			# Use only the first slice one file and get DICOM DICTIONARY
			path_filenameID = img_folder+StudyID+'/'+eID+'/'+SeriesID[int(choseSerie)]+'/'+listSeries_files[0]
						
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
			
			if([0x0019,0x10a4] in dicomInfo):
				print dicomInfo[0x0019,0x10a4]	# SAT Fat/water/bone
				SAT_Fat_water_bone = dicomInfo[0x0019,0x10a4].value
			else:	SAT_Fat_water_bone=''
					
			print dicomInfo[0x0008,0x0060]	#print("Modality.........:", dicomInfo.Modality)
			print dicomInfo[0x0008,0x0030]	#print("Study Time.......:", dicomInfo.StudyTime)
			print dicomInfo[0x0008,0x0031]	#print("Series Time.......:", dicomInfo.SeriesTime)
			print dicomInfo[0x0008,0x0033]	# Imagetime
			print dicomInfo[0x0020,0x0011] 	# SeriesNumber
			print dicomInfo[0x0008,0x103e]	#print("SeriesDescription.......:", dicomInfo.SeriesDescription)
			print dicomInfo[0x0020,0x000E]	#Series Instance UID
							
			if([0x0018,0x1050] in dicomInfo):
				print dicomInfo[0x0018,0x1050]	#print("SpatialResolution.......:", dicomInfo.SpatialResolution)
				SpatialResolution=dicomInfo[0x0018,0x1050].value
			else:	SpatialResolution=''
			if([0x0018,0x1050] in dicomInfo):
				print dicomInfo[0x0020,0x0110]	#print("TemporalResolution.......:", dicomInfo.TemporalResolution)
				TemporalResolution = dicomInfo[0x0020,0x0110].value	
			else:	TemporalResolution=''	
			#print dicomInfo[0x0020,0x0100]	# Temporal Position Identifier
			#print dicomInfo[0x0020,0x0105]	# Number of Temporal Positions
			print dicomInfo[0x0018,0x0080]	# TR #print("Repetition Time.......:", dicomInfo.AcquisitionTime)
			print dicomInfo[0x0018,0x0081]	# TE #print("Echo Time.......:", dicomInfo.AcquisitionTime)
			print dicomInfo[0x0018,0x1314]	# Flip Angle
			#print dicomInfo[0x0018,0x0082]	# Inversion Time
			print dicomInfo[0x0018,0x0091]	# Echo Train Length
			
			print dicomInfo[0x0018,0x0050] 	#print("SliceThickness.......:", dicomInfo.SliceThickness)
			print dicomInfo[0x0018,0x0088] 	#print("SpacingBetweenSlices.......:", dicomInfo.SpacingBetweenSlices)
			if([0x0019,0x1017] in dicomInfo):
				print dicomInfo[0x0019,0x1017]	# Series if Exits
				## Get number of series, + 1 for precontrast
				num_series = dicomInfo[0x0019,0x1017].value
				if (is_number(num_series) ):
					print "\n------------------------------------------------------------------\n"
					if ( len_listSeries_files%num_series == 0): # means that there's no added pre-contrast phase (current series only concern to post contrast
						num_series = num_series
					else:
						num_series = num_series+1
					print "Number of series: %d" % num_series
			else:	
				num_series=0
				SeriesNum =''
			
			#///////////////////////////////////////////////////////////////////////////////////////////#
			# See what kind of Series are dealing with
			flag_dyn = 0
			
			print "choseto firs time"
			print os.getcwd()
			
			choseTodo = raw_input('Enter option to load: \n\t1) get_slices_at_all_locs or \n\t2) get_slices_for_volumes or \n\t3) get_locs_for_all_phases\n\tx) to exit \n\t input: ')
			while(choseTodo!='x'):
				if(choseTodo == '1'):
					print "\n------------------------------------------------------------------"						
					print '''GETTING SLICE AT ALL LOCS \n'''	
					get_slices_at_all_locs(img_folder, abspath_SeriesID, len_listSeries_files, StudyID, eID, SeriesID, choseSerie, listSeries_files, num_series, abspath_ExamID )
				if(choseTodo == '2'):	
					print "\n------------------------------------------------------------------"						
					print '''NOW GET SLICES FOR VOLUMES \n'''	
					get_slices_for_volumes(img_folder, abspath_SeriesID, len_listSeries_files, StudyID, eID, SeriesID, choseSerie, listSeries_files, num_series, abspath_ExamID )
				if(choseTodo == '3'):
					print "\n------------------------------------------------------------------"						
					print '''GETTING LOCS FOR ALL PHASES \n'''
					flag_dyn = 1
					get_locs_for_all_phases(img_folder, SeriesID, choseSerie, StudyID, eID, abspath_ExamID, len_listSeries_files, listSeries_files )
					
					
				choseTodo = raw_input('Enter option: \n\t1) get_slices_at_all_locs or \n\t2) get_slices_for_volumes or \n\t3) get_locs_for_all_phases\n\tx) to exit \n\t input: ')
						
			# *****************************************************************************
			print "Current folder time"
			chdirname='Z:/Cristina/MassNonmass/'+img_folder
			os.chdir(chdirname)  
			print os.getcwd()

			# Save the file list to read as series later
			summary_cases = 'SUMMARY_CASES.txt'
			file_summary_cases = open(summary_cases, 'a')
			
			# Already run header line - now just append infor
			#file_summary_cases.write(str("MRN\t")+' '+str("StudyNumber\t")+' '+str("DicomExamNumber \t")+' '+str("Description \t")+' '+str("AdditonalDescription \t")+' '+str("AllFindings \t")+' '+str("BreastSide \t")+' '+str("QuadrantIQ \t")+' '+str("ClockStart \t")+' '+str("ClockEnd \t")+' '+str("MRISeries \t")+' '+str("LesionCenterPixelX \t")+' '+str("LesionCenterPixelY \t")+' '+str("LesionFindingType.Description \t")+' '+str("WorstHistoType \t")+' '+str("PathoResultTypeID \t")+' '+str("PathoResultTypeID.Des \t")+' '+str("ACRScoreID \t")+' '+str("DCEInitialPhaseEnhTypeID \t")+' '+str("DCEInitialPhaseEnhTypeID.Desc \t")+' '+str("DCEDelayedPhaseEnhID \t")+' '+str("DCEDelayedPhaseEnhID.Desc \t")+' '+str("MassLesionIntEnhanceType \t")+' '+str("MassLesionIntEnhanceType.Desc \t")+' '+str("MassLesionShapeType \t")+' '+str("MassLesionShapeType.Desc \t")+' '+str("MassMarginTypeID \t")+' '+str("MassMarginTypeID.Desc \t")+' '+str("StudyDate\t")+' '+str("Patient id (MRN)\t")+' '+str("ClinicalTrialID\t")+' '+str("StudyID\t")+' '+str("PatientAge\t")+' '+str("Pulse sequence Name\t")+' '+str("Patient Position\t")+' '+str("Image Position (patient)\t")+' '+str("Image Orientation (patient)\t")+' '+str("First Scan location\t")+' '+str("Pulse sequence\t")+' '+str("SAT Fat/water/bone\t")+' '+str("Modality\t")+' '+str("SeriesNumber\t")+' '+str("SeriesDescription\t")+' '+str("SpatialResolution\t")+' '+str("TemporalResolution\t")+' '+str("TR\t")+' '+str("TE\t")+' '+str("Flip Angle\t")+' '+str("ETL\t")+' '+str("SliceThickness\t")+' '+str("SpacingBetweenSlices\t")+' '+str("Series\t")+'\n')
			if len(row) != 0:
				print "%s" % str("MRN: ")+str(row[0])+'\n'+str("StudyNumber: ")+str( row[1])+'\n'+str("DicomExamNumber:  ")+str(row[2])+'\n'
				print "%s" % str("Description: ")+str(row[4])+'\n'
				print "%s" % str("AdditonalDescription: ")+str( row[5])+'\n'
				print "%s" % str("BreastSide: ")+str(row[6])+'\n'+str("QuadrantIQ: ")+str(row[7])+'\n'+str("ClockStart: ")+str(row[8])+'\n'+str("ClockEnd: ")+str(row[9])+'\n'
				print "%s" % str("MRISeries: ")+str(row[10])+'\n'+str("LesionCenterPixelX: ")+str( row[11])+'\n'+str("LesionCenterPixelY: ")+str( row[12])+'\n'+str("LesionFindingType.Description ")+str( row[13])+'\n'
				print "%s" % str("WorstHistoType: ")+str( row[14])+'\n'+str("PathoResultTypeID: ")+str( row[15])+'\n'+str("PathoResultTypeID.Des: ")+str( row[16])+'\n'
				print "%s" % str("NMassLesionID: ")+str( row[17])+'\n'+str("NMassEnhanceTypeID.Desc ")+str( row[19])+'\n'+str("NMassIntEnhType.Desc: ")+str( row[21])
				print "%s" % str("NMassEnhSymmetryType.Desc: ")+str( row[23])+'\n'
			print "%s" % str("StudyDate: ")+str(dicomInfo[0x0008,0x0020].value)+'\n'+str("Patient id (MRN): ")+str(dicomInfo[0x0010,0x0020].value)+'\n'+str("ClinicalTrialID: ")+str(dicomInfo[0x0012,0x0040].value)+'\n'
			print "%s" % str("StudyID: ")+str(dicomInfo[0x0020,0x0010].value)+'\n'+str("PatientAge: ")+str(dicomInfo[0x0010,0x1010].value)+'\n'+str("Pulse sequence Name: ")+str(dicomInfo[0x0019,0x109e].value)+'\n'
			print "%s" % str("Patient Position: ")+str(dicomInfo[0x0018,0x5100].value)+'\n'+str("Image Position (patient): ")+str(dicomInfo[0x0020,0x0032].value)+'\n'
			print "%s" % str("Image Orientation (patient): ")+str(dicomInfo[0x0020,0x0037].value)+'\n'+str("First Scan location: ")+str(dicomInfo[0x0019,0x1019].value)+'\n'+str("Pulse sequence: ")+str(dicomInfo[0x0019,0x109c].value)+'\n'
			print "%s" % str("SAT Fat/water/bone: ")+str(SAT_Fat_water_bone)+'\n'+str("Modality: ")+str(dicomInfo[0x0008,0x0060].value)+'\n'+str("SeriesNumber: ")+str(dicomInfo[0x0020,0x0011].value)+'\n'+str("SeriesDescription: ")+str(dicomInfo[0x0008,0x103e].value)+'\n'
			print "%s" % str("SpatialResolution: ")+str(SpatialResolution)+'\n'+str("TemporalResolution: ")+str(TemporalResolution)+'\n'+str("TR: ")+str(dicomInfo[0x0018,0x0080].value)+'\n'+str("TE: ")+str(dicomInfo[0x0018,0x0081].value)+'\n'+str("Flip Angle: ")+str(dicomInfo[0x0018,0x1314].value)+'\n'+str("ETL: ")+str(dicomInfo[0x0018,0x0091].value)+'\n'
			print "%s" % str("SliceThickness: ")+str(dicomInfo[0x0018,0x0050].value)+'\n'+str("SpacingBetweenSlices: ")+str(dicomInfo[0x0018,0x0088].value)+'\n'+str("Series: ")+str(dicomInfo[0x0019,0x1017].value)+'\n'
			print "%s" % str("Total images in series: ")+str(len_listSeries_files)+'\n\n'
				
			if len(row) != 0:
				file_summary_cases.write(str("MRN: ")+str(row[0])+'\n'+str("StudyNumber: ")+str( row[1])+'\n'+str("DicomExamNumber:  ")+str(row[2])+'\n')
				file_summary_cases.write(str("Description: ")+str(row[4])+'\n')
				file_summary_cases.write(str("AdditonalDescription: ")+str( row[5])+'\n')
				file_summary_cases.write(str("BreastSide: ")+str(row[6])+'\n'+str("QuadrantIQ: ")+str(row[7])+'\n'+str("ClockStart: ")+str(row[8])+'\n'+str("ClockEnd: ")+str(row[9])+'\n')
				file_summary_cases.write(str("MRISeries: ")+str(row[10])+'\n'+str("LesionCenterPixelX: ")+str( row[11])+'\n'+str("LesionCenterPixelY: ")+str( row[12])+'\n'+str("LesionFindingType.Description ")+str( row[13])+'\n')
				file_summary_cases.write(str("WorstHistoType: ")+str( row[14])+'\n'+str("PathoResultTypeID: ")+str( row[15])+'\n'+str("PathoResultTypeID.Des: ")+str( row[16])+'\n')
				file_summary_cases.write(str("NMassLesionID: ")+str( row[17])+'\n'+str("NMassEnhanceTypeID.Desc ")+str( row[19])+'\n'+str("NMassIntEnhType.Desc: ")+str( row[21])+'\n')
				file_summary_cases.write(str("NMassEnhSymmetryType.Desc: ")+str( row[23])+'\n')
							
			file_summary_cases.write(str("StudyDate: ")+str(dicomInfo[0x0008,0x0020].value)+'\n'+str("Patient id (MRN): ")+str(dicomInfo[0x0010,0x0020].value)+'\n'+str("ClinicalTrialID: ")+str(dicomInfo[0x0012,0x0040].value)+'\n')
			file_summary_cases.write(str("StudyID: ")+str(dicomInfo[0x0020,0x0010].value)+'\n'+str("PatientAge: ")+str(dicomInfo[0x0010,0x1010].value)+'\n'+str("Pulse sequence Name: ")+str(dicomInfo[0x0019,0x109e].value)+'\n')
			file_summary_cases.write(str("Patient Position: ")+str(dicomInfo[0x0018,0x5100].value)+'\n'+str("Image Position (patient): ")+str(dicomInfo[0x0020,0x0032].value)+'\n')
			file_summary_cases.write(str("Image Orientation (patient): ")+str(dicomInfo[0x0020,0x0037].value)+'\n'+str("First Scan location: ")+str(dicomInfo[0x0019,0x1019].value)+'\n'+str("Pulse sequence: ")+str(dicomInfo[0x0019,0x109c].value)+'\n')
			file_summary_cases.write(str("SAT Fat/water/bone: ")+str(SAT_Fat_water_bone)+'\n'+str("Modality: ")+str(dicomInfo[0x0008,0x0060].value)+'\n'+str("SeriesNumber: ")+str(dicomInfo[0x0020,0x0011].value)+'\n'+str("SeriesDescription: ")+str(dicomInfo[0x0008,0x103e].value)+'\n')
			file_summary_cases.write(str("SpatialResolution: ")+str(SpatialResolution)+'\n'+str("TemporalResolution: ")+str(TemporalResolution)+'\n'+str("TR: ")+str(dicomInfo[0x0018,0x0080].value)+'\n'+str("TE: ")+str(dicomInfo[0x0018,0x0081].value)+'\n'+str("Flip Angle: ")+str(dicomInfo[0x0018,0x1314].value)+'\n'+str("ETL: ")+str(dicomInfo[0x0018,0x0091].value)+'\n')
			file_summary_cases.write(str("SliceThickness: ")+str(dicomInfo[0x0018,0x0050].value)+'\n'+str("SpacingBetweenSlices: ")+str(dicomInfo[0x0018,0x0088].value)+'\n'+str("Series: ")+str(dicomInfo[0x0019,0x1017].value)+'\n')
			file_summary_cases.write(str("Total images in series: ")+str(len_listSeries_files)+'\n\n')
			file_summary_cases.close()
				
			# Save the file list to read as series later
			brief_cases = 'BRIEF_CASES.txt'
			file_brief_cases = open(brief_cases, 'a')
			#file_brief_cases.write(str("MRN\t")+str("StudyDate\t")+str("Study#\t")+str("DicomExam#\t\t")+str("SeriesDescription\t")+str("#Images\t")+str("Pathology\t")+str("ACRScore\t")+str("TR\t")+str("TE\t")+str("FlipA\t")+str("ETL\t")+str("SliceThickness\t")+str("SpacingBetweenSlices\t")+str("Series\t")+'\n')
			file_brief_cases.write(str(dicomInfo[0x0010,0x0020].value)+'\t\t'+str(dicomInfo[0x0008,0x0020].value)+'\t'+str(dicomInfo[0x0012,0x0040].value)+'\t'+str(ExamID)+'\t\t'+str(dicomInfo[0x0008,0x103e].value)+'\t'+str(len_listSeries_files)+'\t\t'+str(row[14])+'\t'+str(row[16])+'\t\t\t'+str(dicomInfo[0x0018,0x0080].value)+'\t'+str(dicomInfo[0x0018,0x0081].value)+'\t'+str(dicomInfo[0x0018,0x1314].value)+'\t'+str(dicomInfo[0x0018,0x0091].value)+'\t'+str(dicomInfo[0x0018,0x0050].value)+'\t\t'+str(dicomInfo[0x0018,0x0088].value)+'\t\t\t'+str(dicomInfo[0x0019,0x1017].value)+'\n')
			file_brief_cases.close()
			
			# ask for which series to display
			print "\n----------------------------------------------------------"
			
			# Enter folder
			os.chdir(abspath_ExamID)
			if (flag_dyn == 1):				
				os.chdir('DynPhases')
				abspath_SeriesID='Z:/Cristina/MassNonmass/'+img_folder+StudyID+'/'+eID+'/'+'DynPhases/'
				
			arranged_folders = []
			arranged_folders_dyn = get_display(abspath_SeriesID)
			
			if(arranged_folders_dyn == [] and flag_dyn == 0):
				print "No Folders creater -- Change to higher lever folder"
				os.chdir(abspath_ExamID)
				
				print os.getcwd()
				arranged_folders = get_display_series(abspath_ExamID)
				
			if(arranged_folders_dyn == [] and flag_dyn == 1):
				print os.getcwd()
				arranged_folders = get_display_series(abspath_SeriesID)
				os.chdir(abspath_ExamID)
					
			choseDisplay = raw_input('Enter n Folder to load (0-n), or x to exit: ')
						
			while(choseDisplay!='x'):
				print "\n------------------------------------------------------------------"
				print ''' LOAD IMAGES '''
				if(arranged_folders_dyn != []):
					chosen_folder = arranged_folders_dyn[int(choseDisplay)]
					abspath_FolderID = abspath_SeriesID+'/'+chosen_folder
					os.chdir(abspath_FolderID)
					
				if(arranged_folders != []):
					chosen_folder = arranged_folders[int(choseDisplay)]
					abspath_FolderID = abspath_ExamID+'/'+chosen_folder
					os.chdir(abspath_FolderID)
				
				#proc = subprocess.Popen('ls', stdout=subprocess.PIPE)
				#listchosen_images = proc.stdout.read().split()
				#listchosen_images = list(listchosen_images)
				#print listchosen_images
				
				'''-----------------------------'''
				vtkStringArray = vtk.vtkStringArray()
				filename_series = 'DIRCONTENTS.txt'
				
				#check firs if it exists
				if(FileCheck(filename_series) ):
					file_series = open(filename_series,'r')
					files_toRead = file_series.read()
					i = 0
					try:
						for slicename in files_toRead.split():
							#print "adding : %s" %  slicename
							vtkStringArray.InsertNextValue( slicename )
					finally:
						file_series.close()
					
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
				else:
					os.chdir(abspath_ExamID)
					dicomReader  = vtk.vtkDICOMImageReader()
					dicomReader.SetDirectoryName( chosen_folder )
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
			print "Current folder time"
			chdirname='Z:/Cristina/MassNonmass/'
			os.chdir(chdirname)  
			print os.getcwd()
			
			[abspath_ExamID, eID, SeriesID, studyFolder] = get_series(StudyID, img_folder)
			
			print '''Ask for another series or continue \n'''
			   
			choseSerie = raw_input('Enter n Series to load (0-n), or x to exit: ')
			
	
						
						
finally:
	file_ids.close()		
		
		
		
		
