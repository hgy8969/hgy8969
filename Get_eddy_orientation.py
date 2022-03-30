# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 14:55:41 2022
@author: hgy8969
Please acknowledge us in your future publication/presentation 
by citing our paper and mentioning the source code. Thank you so much.
https://doi.org/10.1016/j.rse.2021.112348
"""
import numpy as np
from math import cos,pi
import math 

def ellipse(X,Y):
    a=[pow(x,2) for x in X]
    b=[pow(y,2) for y in Y]
    c=[X[i]*Y[i] for i in range(len(X))]
    L=np.transpose(np.array([a,b,c,X,Y]))
    l=np.ones(len(X)).reshape([len(X),1])
    para=np.linalg.lstsq(L,l,rcond=-1)[0]     
    return para

def get_angle(eX_r,eY_r):
    para=ellipse(eX_r,eY_r)
    A=-para[0][0];B=-para[2][0];C=-para[1][0];D=-para[3][0];E=-para[4][0]
    
    delta=A*C-B**2
    if delta>0:
        if B==0 and abs(A)<abs(C):
            orientation=0 
        elif B==0 and abs(A)>abs(C):            
            orientation=np.pi/2          
        elif B!=0 and abs(A)<abs(C):           
            ex=0.5*np.arctan(B/(A-C))
            orientation=ex          
        elif B!=0 and abs(A)>abs(C):
            ex=0.5*np.arctan(B/(A-C))
            orientation=np.pi/2+ex
            
        if orientation>np.pi/2:
            orientation=orientation-np.pi           
    else:
        A=np.array([eX_r,eY_r])
        B=np.cov(A)  
        a,b= np.linalg.eig(B) 
         
        if a[0]>a[1]:
            orientation=(np.arctan(b[1][0]/b[0][0]))
        else:
            orientation=(np.arctan(b[1][1]/b[0][1]))
            
    orientation=orientation*180/np.pi
    return orientation

def cau_deg(list_contour_lon,list_contour_lat,list_center_lon,list_center_lat):
    ''''''
    cx=[];cy=[]
    for k in range(list_contour_lon):
        eX,eY=list_contour_lon[k],list_contour_lat[k]            
        center_lon,center_lat=list_center_lon[k],list_center_lat[k]
        
        eX_km=(np.array(eX)-center_lon)*111.7*np.cos(eY*pi/180)
        eY_km=(np.array(eY)-center_lat)*111.7
        
        angles=get_angle(eX_km,eY_km)
        orientation=angles*np.pi/180

        if center_lon<0:
            center_lon=center_lon+360
            
        para=ellipse(eX,eY)
        A=-para[0][0];B=-para[2][0];C=-para[1][0];D=-para[3][0];E=-para[4][0]
        
        delta=A*C-B**2
        if delta>0:          
            xinxx=(B*E-2*C*D)/(4*A*C-B*B);xinyy=(B*D-2*A*E)/(4*A*C-B*B)

            az=math.sqrt((2*(A*xinxx*xinxx+C*xinyy*xinyy+B*xinxx*xinyy-1))/(A+C+math.sqrt((A-C)**2+B**2)))
            bz=math.sqrt((2*(A*xinxx*xinxx+C*xinyy*xinyy+B*xinxx*xinyy-1))/(A+C-math.sqrt((A-C)**2+B**2)))              
            c=math.sqrt(abs(az**2-bz**2))

        else:
            A=np.array([eX,eY])
            B=np.cov(A)
            a,b= np.linalg.eig(B) 
            
            xr=np.zeros(len(eX));yr=np.zeros(len(eY))
        
            xr=eX*(np.cos(orientation))+eY*(np.sin(orientation))
            yr=-eX*(np.sin(orientation))+eY*(np.cos(orientation))
            xmin=min(xr);xmax=max(xr);ymin=min(yr);ymax=max(yr)
                     
            length=(xmax-xmin)/2
            height=(ymax-ymin)/2
            
            if length<height:
                t=length;length=height;height=t                    
            c=abs(math.sqrt(length**2-height**2))
            
        if k==0:
            cx.append(c*np.cos(orientation))
            cy.append(c*np.sin(orientation))
        else:
            cx1=c*np.cos(orientation)
            cy1=c*np.sin(orientation)
            cx2=-c*np.cos(orientation)
            cy2=-c*np.sin(orientation)
            
            dis1=math.sqrt(math.pow(float(cx1)-float(cx[len(cx)-1]),2)+math.pow(float(cy1)-float(cy[len(cy)-1]),2))
            dis2=math.sqrt(math.pow(float(cx2)-float(cx[len(cx)-1]),2)+math.pow(float(cy2)-float(cy[len(cy)-1]),2))
            
            if dis1 <= dis2:
                cx.append(cx1)
                cy.append(cy1)
            else:
                cx.append(cx2)
                cy.append(cy2)
###########################################################################################    
    angle=[];day=[]
    for u in range(0,len(cx)):  
        x=float(cx[u])
        y=float(cy[u])
        
        if x==0 and y>0: 
            angle_=0
        elif x>0 and y==0: 
            angle_=90
        elif x==0 and y<0: 
            angle_=180
        elif x<0 and y==0: 
            angle_=270
        elif x>0 and y>0: 
            angle_=np.pi/2-np.arctan(y/x)
        elif x>0 and y<0:
            angle_=np.pi/2-np.arctan(y/x)
        elif x<0 and y<0: 
            angle_=np.pi*1.5-np.arctan(y/x)
        elif x<0 and y>0: 
            angle_=np.pi*1.5-np.arctan(y/x)
        angle_=angle_*180/np.pi
        angle_-=90
        if angle_<0:
            angle_+=360
        angle.append(angle_)
        day.append(u+1)

    return angle

def roate_boudary(list_contour_lon,list_contour_lat,list_orientation):
    list_contour_lon_new=[];list_contour_lat_new=[]
    for i in range(len(list_contour_lon)):
        eX=np.array(list_contour_lon[i]);eY=np.array(list_contour_lat[i])
        orientation=-list_orientation[i]*pi/180
        xr=eX*(np.cos(orientation))+eY*(np.sin(orientation))
        yr=-eX*(np.sin(orientation))+eY*(np.cos(orientation))
        list_contour_lon_new.append(xr)
        list_contour_lat_new.append(yr)
    return list_contour_lon_new,list_contour_lat_new
    
    
if __name__ == "__main__":
    '''if your eddy is tracked 10 days, then lists should contain 10 eddy boundaries and centers'''
    list_contour_lon=[] #list of longitude of eddy boundaries
    list_contour_lat=[] #list of latitude of eddy boundaries
    list_center_lon=[] #list of lontitude of eddy centers
    list_center_lat=[] ##list of latitude of eddy centers
    '''get the origial eddy orientaion of eddy, which is 0 at the 
       positive-x-axis direction and Increase clockwise'''
    list_orientation=cau_deg(list_contour_lon,list_contour_lat,list_center_lon,list_center_lat)
    
    '''rotate eddy boundaries into horizontal direction'''
    list_contour_lon,list_contour_lat=roate_boudary(list_contour_lon,list_contour_lat,list_orientation)
    
    '''final step'''
    '''you should average the rotated eddy boundaries, and get its semimajor axis of positive x-axis (a1) 
    and negative x-axis (a2). If a1<a2, list_orientation should add 180 degree. Then, you can get the 
    final orientaions of each eddy snapshot'''
