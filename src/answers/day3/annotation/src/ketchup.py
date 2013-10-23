import sys

outputfile = sys.argv[1]

readoutput = open(outputfile)

sbjcount = 0 # counter to make sure it only takes the first subject start coordinate

# runs through all the lines
# collects chromosome ID, start and end coordinate of first hit
for line in readoutput:
	if line[:4] == ">lcl":
		chr = line[5:].strip() # get chromosome ID from >lcl|chr... line
	elif line.startswith("Sbjct"): # the coordinates on the subject are part of this line
		linesave = line # current line is saved to use in the future (in case the line contains the stop coordinate)
		if sbjcount == 0: # only collect start coordinates for the first "Sbjct" line
			start = line.split()[1] # get start coordinate (2nd element in line)
			sbjcount = 1 # prevents going through this if loop again
	elif line.startswith(" Score") and sbjcount == 1: # only goes through this loop if start coordinate is already filled
		stop = linesave.split()[3] # get stop coordinate (4th element in line)
		break # exits loop entirely to make sure the information belongs to the first hit

print "%s:%s-%s" % (chr,start,stop)	

