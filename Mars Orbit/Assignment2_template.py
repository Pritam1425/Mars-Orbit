import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
from math import radians, cos, sin, tan, sqrt, fmod, atan, degrees
import datetime
from scipy import optimize as op

def get_times(data):
    times = [0.0]
    for i in range(1,12):
        d1 = datetime.datetime(data[i][0],data[i][1],data[i][2],data[i][3],data[i][4])
        d2 = datetime.datetime(data[i-1][0],data[i-1][1],data[i-1][2],data[i-1][3],data[i-1][4])
        days = d1-d2
        times.append(times[-1]+days.days+days.seconds/(3600*24))
    return times

def get_oppositions(data):
    oppositions = data[:,5]*30+data[:,6]+data[:,7]/60+data[:,8]/3600
    return oppositions

def dist(ex,ey,m,c,r):
    cx = np.cos(np.radians(c))
    cy = np.sin(np.radians(c))
    
    tangent = np.tan(np.radians(m))
    
    a = 1+tangent**2
    b = -2*cx+2*tangent*(ey-cy-ex*tangent)
    c = (ey-cy-ex*tangent)**2+cx**2-r**2
    
    discrm = np.sqrt(b**2-4*a*c)
    
    x1 = (-b-discrm)/(2*a)
    x2 = (-b+discrm)/(2*a)
    
    y1 = ey+(x1-ex)*tangent
    y2 = ey+(x2-ex)*tangent
    
    if 0<=m<=90 or 270<=m<=360:
        x = x1 if x1>=0 else x2
        y = y1 if x1>=0 else y2
    else :
        x = x1 if x1<=0 else x2
        y = y1 if x1<=0 else y2
    return x,y

def MarsEquantModel(c,r,e1,e2,z,s,times,oppositions):
    ol = list()
    el = list()
    
    for o in oppositions:
        x,y = dist(0,0,o,c,r)
        ol.append([x,y])
    for e in times:
        ex = e1*np.cos(np.radians(e2+z))
        ey = e1*np.sin(np.radians(e2+z))
        x,y = dist(ex,ey,fmod(e*s+z,360),c,r)
        el.append([x,y])
    de = list()
    do = list()
    
    for x in el:
        de.append(np.degrees(np.arctan2(x[1],x[0])))
    for o in ol:
        do.append(np.degrees(np.arctan2(o[1],o[0])))
    
    err = np.array(de)-np.array(do)
    errmax = np.max(np.abs(err))
    return err,errmax

def bestOrbitInnerParams(r,s,times,oppositions):
    cf = 120
    e1f = 1
    e2f = 20
    def obj(X0,oppositions):
        c,e1,e2,z = X0
        err,errmax = MarsEquantModel(c,r,e1,e2,z,s,times,oppositions)
        return errmax
    for i in range(3):
        z = np.linspace(0,360,360)
        zerr = np.array([obj((cf,e1f,e2f,i),oppositions) for i in z])
        zf = z[zerr.argmin()]
        
        e2 = np.linspace(zf,360,360)
        e2err = np.array([obj((cf,e1f,i,zf),oppositions) for i in e2])
        e2f = e2[e2err.argmin()]
        
        c = np.linspace(0,360,360)
        cerr = np.array([obj((i,e1f,e2f,zf),oppositions) for i in c])
        cf = c[cerr.argmin()]
        
        e1 = np.linspace(0,0.5*r,300)
        e1err = np.array([obj((cf,i,e2f,zf),oppositions) for i in e1])
        e1f = e1[e1err.argmin()]
    res = op.minimize(fun=obj,x0=[cf,e1f,e2f,zf],args=oppositions,method='Nelder-Mead', options={'xatol' : 1e-5 ,'disp':False, 'return_all' :False})
    cf,e1f,e2f,zf = res.x
    err,errmax = MarsEquantModel(cf,r,e1f,e2f,zf,s,times,oppositions)
    return cf,e1f,e2f,zf,err,errmax

def bestS(r,times,oppositions):
    maxError = 720
    sgood = 360/687
    iter=1
    for s in np.linspace(360/(687+0.1),360/(687-0.1),12):
        iter+=1
        _c,_e1,_e2,_z,_errors,_maxError = bestOrbitInnerParams(r,s,times,oppositions)
        if maxError > _maxError:
            c,e1,e2,z,errors,maxError=_c,_e1,_e2,_z,_errors,_maxError
            sgood = s
    c,e1,e2,z,errors,maxError = bestOrbitInnerParams(r,sgood,times,oppositions)
    return sgood,errors,maxError

def bestR(s,times,oppositions):
    errmax = 720
    best_rf = 5
    times = np.array(times)
    oppositions = np.array(oppositions)
    for o_rf in np.linspace(5,10,5):
        _c,_e1,_e2,_z,_err,_errmax = bestOrbitInnerParams(o_rf,s,times,oppositions)
        if errmax>_errmax:
            errmax = _errmax
            best_rf = o_rf
        cx = np.cos(np.radians(_c))
        cy = np.sin(np.radians(_c))
        ex = _e1 * np.cos(np.radians(_e2+_z))
        ey = _e1 * np.sin(np.radians(_e2+_z))
        
        eqLine = (times * s + _z) % 360
        
        x = (ey-ex*np.tan(np.radians(eqLine))) / (np.tan(np.radians(oppositions))-np.tan(np.radians(eqLine)))
        y = x * np.tan(np.radians(oppositions))
        
        dst = np.sqrt((x - cx) ** 2 + (y - cy) ** 2)
        nrf = np.mean(dst)
        
        for rf in np.linspace(nrf,min(10,nrf+0.5),5):
            _c,_e1,_e2,_z,_errors,_maxError = bestOrbitInnerParams(rf,s,times,oppositions)
            
            if errmax>_errmax:
                errmax = _errmax
                best_rf = rf
    c,e1,e2,z,errors,maxError = bestOrbitInnerParams(best_rf,s,times,oppositions)
    return best_rf,errors,maxError

def bestMarsOrbitParams(times,oppositions):
    r = 10
    s = 360/687
    emaxt = 1
    def obj(X0,oppositions):
        c,r,e1,e2,z,s = X0
        err,errmax = MarsEquantModel(c,r,e1,e2,z,s,times,oppositions)
        return errmax
    while(emaxt>4/60):
        s,err,errmax = bestS(r,times,oppositions)
        r,err,errmax = bestR(s,times,oppositions)
        
        c,e1,e2,z,err,errmax = bestOrbitInnerParams(r,s,times,oppositions)
        res = op.minimize(fun=obj,x0=[c,r,e1,e2,z,s],args=oppositions,method='Nelder-Mead',options={'xatol' : 1e-5 ,'disp':False, 'return_all' :False})
        c,r,e1,e2,z,s = res.x
        emaxt = errmax
        if emaxt-errmax<0.001:
            break
    c,e1,e2,z,err,errmax = bestOrbitInnerParams(r,s,times,oppositions)
    return r,s,c,e1,e2,z,err,errmax
if __name__ == "__main__":

    # Import oppositions data from the CSV file provided
    data = np.genfromtxt(
        "../data/01_data_mars_opposition_updated.csv",
        delimiter=",",
        skip_header=True,
        dtype="int",
    )

    # Extract times from the data in terms of number of days.
    # "times" is a numpy array of length 12. The first time is the reference
    # time and is taken to be "zero". That is times[0] = 0.0
    times = get_times(data)
    assert len(times) == 12, "times array is not of length 12"

    # Extract angles from the data in degrees. "oppositions" is
    # a numpy array of length 12.
    oppositions = get_oppositions(data)
    assert len(oppositions) == 12, "oppositions array is not of length 12"

    # Call the top level function for optimization
    # The angles are all in degrees
    
    r, s, c, e1, e2, z, errors, maxError = bestMarsOrbitParams(
        times, oppositions
    )
    assert max(list(map(abs, errors))) == maxError, "maxError is not computed properly!"
    print(
        "Fit parameters: r = {:.4f}, s = {:.4f}, c = {:.4f}, e1 = {:.4f}, e2 = {:.4f}, z = {:.4f}".format(
            r, s, c, e1, e2, z
        )
    )
    print("The maximum angular error = {:2.4f}".format(maxError))
