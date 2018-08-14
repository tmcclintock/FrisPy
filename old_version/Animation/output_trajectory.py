import numpy as np

def write_trajectory_to_file(positions,n_times):
    outfile_name = "data/test_data.txt"
    outfile = open(outfile_name,"w")
    
    #We will only output the positions and time, not the angles
    #or the velocities
    headerline = "x\ty\tz\tt\n"
    outfile.write(headerline)
    
    for i in range(n_times):
        outline = "%.6E"%positions[i,1]+"\t%.6E"%positions[i,2]+"\t%.6E"%positions[i,3]+"\t%.6E"%positions[i,0]+"\n"
        outfile.write(outline)
        continue
    
    print "Output to file:\n\t"+outfile_name
    outfile.close()
    return
