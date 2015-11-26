import numpy as np
import cv2

#read in the video
cap = cv2.VideoCapture('/home/tom/Desktop/Disc_research/meanshifttest.mp4')
print cap.isOpened()

cap2 = cv2.VideoCapture(0)
print cap2.isOpened()

#take the first frame
ret,frame=cap.read()
print ret,frame

#set up the initial location of the window
r,c,h,w=550,630,150,150
track_window=(c,r,w,h)

#set up the region of interest (ROI) for tracking
roi=frame[r:r+h, c:c+w]
hsv_roi=cv2.cvtColor(roi, cv2.COLOR_BGR2HSV)
mask=cv2.inRange(hsv_roi, np.array((0,60,32)),np.array((180,255,255)))
roi_hist=cv2.calcHist([hsv_roi],[0],mask,[180],[0,180])
cv2.normalize(roi_hist,roi_hist,0,255,cv2.NORM_MINMAX)

#set the termination criteria for the algorithm: 10 iterations or 1 pt
term_crit=(cv2.TERM_CRITERIA_EPS | cv2.TERM_CRITERIA_COUNT, 10, 1)

output = open("positions.txt","w")

while(1):
	ret, frame=cap.read()
	
	if ret==True:
		hsv=cv2.cvtColor(frame,cv2.COLOR_BGR2HSV)
		dst=cv2.calcBackProject([hsv],[0],roi_hist,[0,180],1)
		
		#apply the meanshift algorithm to get the new location
		ret, track_window=cv2.meanShift(dst, track_window, term_crit)
		
		#draw the window on the image
		x,y,w,h=track_window
		img2=cv2.rectangle(frame, (x,y), (x+w, y+h), 255,2)
		cv2.imshow('img2',img2)

		print "X = ",x+w/2., "  Y = ",y+h/2.
		output.write(str(x+w/2.)+"\t"+str(y+h/2.)+"\n")
		k = cv2.waitKey(30) & 0xff
    	#manually terminate with cntrl+C - termination criteria not working

output.close()
