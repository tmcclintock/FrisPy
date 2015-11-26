import numpy as np
import cv2

#read in the video
cap = cv2.VideoCapture('camshifttest.mp4')

#take the first frame
ret,frame=cap.read()

#set the tracking window. The values are hand entered
r,c,h,w=150,600,200,200
#r and c are position and h and w are the dimensions of the box
track_window=(c,r,w,h)

#set up the region of interest (ROI) for tracking
roi=frame[r:r+h, c:c+w]
hsv_roi=cv2.cvtColor(roi, cv2.COLOR_BGR2HSV)
mask=cv2.inRange(hsv_roi, np.array((0,60,32)),np.array((180,255,255)))
roi_hist=cv2.calcHist([hsv_roi],[0],mask,[180],[0,180])
cv2.normalize(roi_hist,roi_hist,0,255,cv2.NORM_MINMAX)

#termination criteria for the camshift algorithm: 10 iterations or 1 point
term_crit=(cv2.TERM_CRITERIA_EPS | cv2.TERM_CRITERIA_COUNT, 10, 1)

while(1):
	ret, frame=cap.read()
	if ret==True:
		hsv=cv2.cvtColor(frame,cv2.COLOR_BGR2HSV)
		dst=cv2.calcBackProject([hsv],[0],roi_hist,[0,180],1)
		
		#apply camshift to find the new location
		ret, track_window=cv2.CamShift(dst, track_window, term_crit)

		#draw it on the image
		pts=cv2.boxPoints(ret)
		pts=np.int0(pts)
		img2=cv2.polylines(frame,[pts],True,255,2)
		cv2.imshow('img2',img2)

		k = cv2.waitKey(30) & 0xff
		#have to end the program manually with control+c - couldn't get the termination criteria to work
    
