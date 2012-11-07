#!/usr/bin/env python

import os
import os.path
import sys
from sys import argv, stderr, exit
import shlex, subprocess
from numpy import *
from scipy.io import loadmat, savemat

# ******************************************************************************
print "usage: importDicom.py /folder-with-dicoms/.txt STUDYNo \t EXMAANo"
print '''
Run on MassNonmass root folder organizebyDicomSeries.py based on list of 
patient_list.txt(list of StudyIDs in folder mass or nonmass)

"usage:" 
python organizebyDicomSeries.py patient_list.txt
-- This script will recognize 'S' and 'I' tags in the filenames and create folders for 
each of the series Si and move the images 'Si'xxIj into Si:
	iterate over lines in patient_list.txt, each corresponding to a patient
	creates:
		abspath_studyID
		abspath_ExamID
		
	Enters each StudyID folder and finds subfolders (n Exams) with
		ExamsID = get_immediate_subdirectories(path_studyID);

	iterates ExamsID and moves corresponding files

% Copyright (C) Cristina Gallego, University of Toronto, 2012
% June 12/12 - Now works with added ExamsNo on second column of patient_list.txt
% June 27/12 - Now works with more added columns afer ExamsNo on TotalFromSharmila.txt
-
-----------------------------------------
'''

def get_immediate_subdirectories(dirname):
    return [name for name in os.listdir(dirname)
            if os.path.isdir(os.path.join(dirname, name))]
            
def find(strng, ch):
    index = 0
    while index < len(strng):
        if strng[index] == ch:
            return index
        index += 1
    return -1
    
# Get Root folder ( the directory of the script being run)
path_rootFolder = os.path.dirname(os.path.abspath(__file__))
print path_rootFolder
            
# Open filename list
file_ids = open(sys.argv[1],"r")
try:
	for line in file_ids:
		
		# Enter the studyID folder
		StudyID = line.lstrip()
		indStudyID = find(StudyID, '\t')
		StudyID = StudyID[0:indStudyID]
		print " "               
		print "StudyID: %s" % StudyID
		
		# obtain subdires in the StudyID directory, correspond to ExamsID 
		# check for one in subdirectory tree, (e.g Mass or NonMass) 
			
		path_studyID = 'mass/'+StudyID
		ExamsID = get_immediate_subdirectories(path_studyID);
		print ExamsID
		
		abspath_studyID = os.path.abspath(path_studyID)
		print abspath_studyID
		
		for eID in ExamsID:
			
			path_ExamID = 'mass/'+StudyID+'/'+eID
			print path_ExamID
			abspath_ExamID = os.path.abspath(path_ExamID)
			print abspath_ExamID
			
			# obtain total filenames in directory
			list_files = os.listdir(abspath_ExamID)
			os.chdir(abspath_ExamID)
			
			# Get total number of files
			print len(list_files)

			# Iterate and split based on image Series
			for filename in list_files:	
				# get image file names
				#print filename
				findSeries_tag = filename.find ( 'S' );
				findImages_tag = filename.find ( 'I' );
	
				# get only the series name
				series_name = filename[findSeries_tag:findImages_tag];
	
				if not os.path.exists(str(series_name)):
					os.makedirs(str(series_name))
		
				os.chdir(str(series_name))
				mv_examFolder = str(abspath_ExamID)+'/'+filename
				mv_examFolder2 = str(abspath_ExamID)+'/'+str(series_name)
				print mv_examFolder
				
				mv = subprocess.Popen(['mv', mv_examFolder, mv_examFolder2 ], stdout=subprocess.PIPE)
				os.chdir(abspath_ExamID)
			
			os.chdir(path_rootFolder)
			print "go back:"
			print path_rootFolder
		
		os.chdir(path_rootFolder)
   
finally:
	file_ids.close()    
	 	
		
	
	
	

