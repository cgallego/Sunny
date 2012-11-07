#!/usr/bin/env python

import os
import os.path
import sys
from sys import argv, stderr, exit
import shlex, subprocess


# ******************************************************************************
print "usage: importDicom.py /folder-with-dicoms/.txt"
print " "

def get_immediate_subdirectories(dir):
    return [name for name in os.listdir(dir)
            if os.path.isdir(os.path.join(dir, name))]
            
# Open filename list
file_ids = open(sys.argv[1],"r")
try:
	for line in file_ids:
		# Enter the studyID folder
		StudyID = line;
		ExamsID = get_immediate_subdirectories(StudyID);
		
		print ExamsID
		
		os.chdir(str(StudyID))
			
		# obtain total filenames in directory
		list_files = os.listdir(os.path.abspath(ExamsID))
		
		files_dir = os.path.abspath(StudyID)
		os.chdir(files_dir)

		print files_dir

		# Get total number of files
		print len(list_files)

		# Iterate and split based on image Series
		for filename in list_files:	
			# get image file names
			print filename
			findSeries_tag = filename.find ( 'S' );
			findImages_tag = filename.find ( 'I' );
	
			# get only the series name
			series_name = filename[findSeries_tag:findImages_tag];
	
			if not os.path.exists(str(series_name)):
				os.makedirs(str(series_name))
		
			os.chdir(str(series_name))
			mv_examFolder = str(files_dir)+'/'+filename
			mv_examFolder2 = str(files_dir)+'/'+str(series_name)
			print mv_examFolder
	
			mv = subprocess.Popen(['mv', mv_examFolder, mv_examFolder2 ], stdout=subprocess.PIPE)
			os.chdir(files_dir)	
   
	
		
	
	
	

