import sys 
import os
import re
import colorama
colorama.init()

import refactoring # Note: refactoring.py need to be in the current working directory

paths = [
	"/Users/petst/Work/Projects/Inviwo-Developent/Private/Inviwo-dev"
#	"/Users/petst/Work/Projects/Inviwo-Developent/Private/Inviwo-research"
]

#paths = [
#	"C:/Users/petst55/Work/Inviwo/Inviwo-dev",
#	"C:/Users/petst55/Work/Inviwo/Inviwo-research"
#]

excludespatterns = ["*/ext/*", "*moc_*", "*/proteindocking/*", "*/proteindocking2/*", 
					"*/genetree/*", "*.DS_Store", "*DS_mapp",
					"*.png", "*.ttf", "*.tif", "*.pyc", "*.raw", "*.bmp", "*.wav", "*.ico", "*.icns",
					"*.qch", "*.qhc", "*.exr", "*.pwm", "*.pvm", "*.pdf", "*.otf", "*.exe"]

files = refactoring.find_files(paths, ['*'], excludes=excludespatterns)

def replace(pattern, replacement) :
	print("Matches:")
	matches = refactoring.find_matches(files, pattern)
	
	print("\n")
	print("Replacing:")
	refactoring.replace_matches(matches, pattern, replacement)


copyright_replacements  = {
	r"(\s*[*#]\s+Copyright \(c\) 201\d-)201[123456]( Inviwo Foundation\s*)" : r"\g<1>2017\g<2>",
	r"(\s*[*#]\s+Copyright \(c\) )(201[123456])( Inviwo Foundation\s*)" : r"\g<1>\g<2>-2017\g<3>"
}

print("Looking in " + str(len(files)) + " files")

for k,v in copyright_replacements.items():
	replace(k, v)



