# -*- coding: utf-8 -*-
"""
Created on Sun Sep  4 18:46:09 2022

@author: Pritam Sarkar
"""

import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

# To remove Instruction matplotlib.pyplot commnet PlotOut Function And Line:578
import scipy.optimize as sp

# Note :Code written for fix no of 12 Opposition Data

# Assignment Q1 Return Value in Maxerror(int) Error (List)
def MarsEquantModel(cf, rf, e1f, e2f, zf, sf, oppositionsf):
    '''
    Parameter Info:
    (cx,cy)= Center of Circle calculated by Polar Coordinate C_dis=1 (distance from Sun) cf(angle from Aris)
    (ex,ey)= Equant position cal by Polar Coordinate e1f (Distance from Sun) e2f (angle from Aris)
    rf= Radius Of Orbit
    sf= Degree per Day
    Dotlinef= Calculation of Strokes angle from Aris and Rotate to Equant 0 by adding Angle Zf
    (X_Circle, Y_Circle)= point Intersection Point of Circle and 12 Strokes in the Direction of the respective Stroke
    Cal_Angle= aTan(Y_Circle/X_Circle) in (-pi,pi] range
    Act_Angle= Given Latitude Angle converted in (-pi,pi] range

    Logic:
    Solving Equation
    (x - cx) ** 2 + (y - cy)**2 =r **2
    (y- ey) = (x-ex) * Tan( Dotlinef )

    simplifying to 2nd Order Equation and solve it
    From (x1,y1), (x2,y2)
    discarding one point by checking Dotline (Stroke angle) lies in which XY quadrant
    ( First And Fourth Quadrant X1 >0 and opposite in Third and Second)

    Return Value:
    Error =Calculated Angle(Cal_Angle)- Actual Longitude Angle (Act_Angle)
    Max Error= Max of absolute Error
    '''
    oppositionsf=np.array(oppositionsf).T # Coverting to 2x12 Array (Daydiff, Actual Longitude)
    C_dis = 1
    cx = C_dis * np.cos(np.radians(cf))
    cy = C_dis * np.sin(np.radians(cf))
    ex = e1f * np.cos(np.radians(e2f+zf))
    ey = e1f * np.sin(np.radians(e2f+zf))
    # Doted Line Angle At Z=0;
    Dotlinef = (oppositionsf[0] * sf + zf) % 360
    # print('Doted line Angle=', Dotlinef)

    X_circle = np.array([np.nan for x in range(0,np.shape(Dotlinef)[0])])
    Y_circle = np.array([np.nan for x in range(0,np.shape(Dotlinef)[0])])

    for k in range(0, np.shape(Dotlinef)[0]): # For Intersect (X,Y) Coordinate
        tanValue = np.tan(np.radians(Dotlinef[k]))
        a = (1 + tanValue ** 2)
        b = -2 * cx + 2 * tanValue * (ey - cy - ex * tanValue)
        c = (ey - cy - ex * tanValue) ** 2 + cx ** 2 - rf ** 2
        delta = np.sqrt(b ** 2 - 4 * a * c)
        x1 = (-b - delta) / (2 * a)
        x2 = (-b + delta) / (2 * a)
        y1 = ey + (x1 - ex) * tanValue
        y2 = ey + (x2 - ex) * tanValue
        if 0 <= Dotlinef[k] <= 90 or 270 <= Dotlinef[k] <= 360:
            X_circle[k] = x1 if x1 >= 0 else x2
            Y_circle[k] = y1 if x1 >= 0 else y2
        else:
            X_circle[k] = x1 if x1 <= 0 else x2
            Y_circle[k] = y1 if x1 <= 0 else y2
        # plt.plot(X_circle[k],Y_circle[k] , 'bo')
    Cal_angle = np.degrees(np.arctan2(Y_circle,X_circle))
    Act_Angle1 = np.array([i if i <= 180  else i-360 for i in oppositionsf[1]])
    # print('Dotline',Dotlinef,'\nX_circle',X_circle,'\n Y_circle',Y_circle,'\nCal_angle=',Cal_angle,'\n Act_Angle',Act_Angle1)
    # errors = np.array([min(np.absolute(np.subtract(Act_Angle1,Cal_angle[i]))) for i in range(12)])
    errors = np.subtract(Cal_angle,Act_Angle1)
    return errors, np.max(np.absolute(errors))
# End Q1

def PlotOut(oppositionsf, zf, cf, rf, e1f, e2f,sf):
    oppositionsf = np.array(oppositionsf).T
    T = 360/sf
    C_dis = 1
    cx = C_dis * np.cos(np.radians(cf))
    cy = C_dis * np.sin(np.radians(cf))
    xp = np.linspace(0, rf+3, 200)
    xn = np.linspace(-(rf+3), 0, 200)
    ex = e1f * np.cos(np.radians(e2f +zf))
    ey = e1f * np.sin(np.radians(e2f +zf))
    Dotline_A = ((oppositionsf[0] * 360 / T) % 360 + zf) % 360
    Solidline_A= oppositionsf[1]
    # print('c',cf,'r',rf,'e1',e1f,'e2',e2f,'z',zf,'s',s)
    # print('\n Data _ opposition',oppositionsf)
    # print('\ncx',cx,'cy',cy,'ex',ex,'ey',ey)
    for i in range(0, np.shape(Dotline_A)[0]):
        if 0 <= Solidline_A[i] <= 90 or 270 <= Solidline_A[i] <= 360:
            plt.plot(xp, xp * np.tan(np.radians(Solidline_A[i])), label='Line' + str(i), linestyle="-",
                     linewidth=1 + i * .05)
            plt.text((rf-3)/2, (rf-3)/2 * np.tan(np.radians(Solidline_A[i])), str(i), fontsize=8,color='green')
        else:
            plt.plot(xn, xn * np.tan(np.radians(Solidline_A[i])), label='Line' + str(i), linestyle="-",
                     linewidth=1 + i * .05)
            plt.text(-(rf - 3) / 2, -(rf - 3) / 2 * np.tan(np.radians(Solidline_A[i])), str(i), fontsize=8, color='green')
    xp = np.linspace(ex, rf+3, 200)
    xn = np.linspace(-(rf+3), ex, 200)
    for i in range(0, np.shape(Dotline_A)[0]):
        if 0 <= Dotline_A[i] <= 90 or 270 <= Dotline_A[i] <= 360:
            plt.plot(xp, ((xp - ex) * np.tan(np.radians(Dotline_A[i])) + ey), label='Line' + str(i),
                     linestyle="--", linewidth=1 + i * .05)
            plt.text((rf - 1) / 2, (((rf - 1) / 2 - ex) * np.tan(np.radians(Dotline_A[i])) + ey), str(i), fontsize=8, color='black')
        else:
            plt.plot(xn, (xn - ex) * np.tan(np.radians(Dotline_A[i])) + ey, label='Line' + str(i),
                     linestyle="--", linewidth=1 + i * .05)
            plt.text(-(rf - 1) / 2, ((-(rf - 1) / 2 - ex) * np.tan(np.radians(Dotline_A[i])) + ey), str(i), fontsize=8, color='black')
    X_circle = np.array([np.nan for x in range(0,np.shape(Dotline_A)[0])])
    Y_circle = np.array([np.nan for x in range(0,np.shape(Dotline_A)[0])])
    for k in range(0, np.shape(Dotline_A)[0]):
        tanValue = np.tan(np.radians(Dotline_A[k]))
        a = (1 + tanValue ** 2)
        b = -2 * cx + 2 * tanValue * (ey - cy - ex * tanValue)
        c = (ey - cy - ex * tanValue) ** 2 + cx ** 2 - rf ** 2
        delta = np.sqrt(b ** 2 - 4 * a * c)
        x1 = (-b - delta) / (2 * a)
        x2 = (-b + delta) / (2 * a)
        y1 = ey + (x1 - ex) * tanValue
        y2 = ey + (x2 - ex) * tanValue
        plt.plot(x1, y1, 'go')
        plt.plot(x1, y1, 'ro')
        if 0 <= Dotline_A[k] <= 90 or 270 <= Dotline_A[k] <= 360:
            X_circle[k] = x1 if x1 >= 0 else x2
            Y_circle[k] = y1 if x1 >= 0 else y2
        else:
            X_circle[k] = x1 if x1 <= 0 else x2
            Y_circle[k] = y1 if x1 <= 0 else y2
        plt.plot(X_circle[k],Y_circle[k] , 'bo')
    plt.annotate('Rotating Axis Ref', xy=(.5, (-ex + .5) * np.tan(np.radians(zf)) + ey), xytext=(rf*.5, rf*.5),
                 arrowprops=dict(arrowstyle="->", connectionstyle="angle3,angleA=0,angleB=-90"))
    plt.annotate('Equant', xy=(ex, ey), xytext=(2, 2),
                 arrowprops=dict(arrowstyle="->", connectionstyle="angle3,angleA=0,angleB=-90"))
    theta = np.radians(np.linspace(0, 360 * 2, 100))
    plt.plot(cx + rf * np.cos(theta), cy + rf * np.sin(theta) )
    plt.plot(cx,cy,'ro')
# End Plot

# Assignment Q2 Return Value in cf, e1f, e2f, zf, Error (List), Max Error (int)
def bestOrbitInnerParams(rf,sf,oppositions):
    '''
    Parameter Info:
    c=  cf(Circle center angle from Aris)
    (e1,e2)= Polar Coordinate e1f (Distance from Sun) e2f (angle from Aris)
    rf= Radius Of Orbit
    sf= Degree per Day
    Oppositions = (Daydiff, Latitude Angle) (12 x 2 ) Array
    Dotlinef= Calculation of Strokes angle from Aris and Rotate to Equant 0 by adding Angle Zf

    Logic

    1)Approx Calculation:
    x0i (Intial Approx Functional Parameter Variable)
        Logic:(Grid Search Method by Keeping One Parameter Varied Over range fix remain at time)

    2)Fine Control Calculation:

    Optimized Method Used by Scipy Method : Nelder-Mead

    LoopValue= Control the Flow of Search
            Condition LoopValue Remains True
            If Max Error of Function is (Error>= 100,10<= Error<=100,0.3<= Error<=10)
            correspond incremental Initial parameter given to optimize fuction By Method 1,2,3 Respectively

    Method 1:
        If Max Error is >=100 than Initial Value Start From X0i with Inc In Parameter
            Increment is Define By Zf -(0,360) if its Exhust than Next (Cf, E2f ) -((0,360),(0,360))  Than E2 Varies (1,.5*rf)

    Method 2:
        If Max Error range (10,100)
            Increment(Step .3 of Last Value Optimize function Calculated) is Define by changeing Cf  & E2f to check the conversion is finite

    Method 3:
        If Max Error range (1,10)
            Passing Same Optimize Value to Optimize function And Checking
            Divergence Value(lfun(Second Last Max Error)-result.fun(Last Max Error))

    Total Fine Search Iteration Limit to 30

    Return Para:
    (cf, e1f, e2f, zf, Error, Max Error ) cal using Last Result

    '''

    # Approx Intial Value Deciding factor
    def intialValue(rf, sf, oppositionsf):
        '''

        :param rf: Fix orbit Radius
        :param sf: Fis Mars Deg/Day
        :param oppositionsf: Oppositon Data

        Logic:
            (C,E1,E2,Z) Max range ((0,360),(0,.5*Rf),(0,360),(0,360))
            Find Minerror by Varying each parameter Over a Bound Keeping Other Parameter Fixed in sequence of Z,E2,C,E1
            Loop:
                (Intial [C,E1,E2,Z] )--> [C,E1,E2,Z_New]-->[C,E1,E2_New,Z_New]-->[C_New,E1,E2_New,Z_New]-->[C_New,E1_New,E2_New,Z_New])
                New Value Corresponds to Min error at that parameter Varied
        :return:
            Approx Initial Guess [C,E1,E2,Z]
        '''

        def maxerror(cf, e1f, e2f, zf):
            Er, MEr = MarsEquantModel(cf, rf, e1f, e2f, zf, sf, oppositionsf)
            return MEr
        #Initial Guess
        cf = 10
        e1f = 1
        e2f = 10
        for i in range(0, 3):
            z1 = np.linspace(0, 360, 360)
            zerror = np.array([maxerror(cf, e1f, e2f, z1[i]) for i in range(0, z1.shape[0])])
            # print(zerror.argmin(),np.min(zerror),zerror)
            zf = z1[zerror.argmin()]
            e2f1 = np.linspace(zf, 360, 360)
            e2ferror = np.array([maxerror(cf, e1f, e2f1[i], zf) for i in range(0, e2f1.shape[0])])
            # print(e2ferror.argmin(), np.min(e2ferror), e2ferror)
            e2f = e2f1[e2ferror.argmin()]
            cf1 = np.linspace(0, 360, 360)
            cferror = np.array([maxerror(cf1[i], e1f, e2f, zf) for i in range(0, cf1.shape[0])])
            cf = cf1[cferror.argmin()]
            e1f1 = np.linspace(0, .5 * rf, 300)
            e1ferror = np.array([maxerror(cf, e1f1[i], e2f, zf) for i in range(0, e1f1.shape[0])])
            e1f = e1f1[e1ferror.argmin()]
            # print(e1ferror.argmin(), np.min(e1ferror), e1ferror)
        return cf, e1f, e2f, zf
    # End Intial Parameter code

    def maxerror(x):
        '''

        :param x: [c,e1,e2,z] Internal Function to avoid Iteration

        :return: Reuturn the Max Error to given c, e1,e2,z,s,Opposition Data
        '''
        c, e1, e2, z = x
        Er, MEr = MarsEquantModel(c, rf, e1, e2, z, sf, oppositions)
        return MEr

    x0i=[intialValue(rf,sf,oppositions)] # Intial Value
    x0=x0i
    LoopValue=True
    Iter = int(1)
    incf=[0,0,0,0]
    lfun=10000
    # print('Initial Value [c,e1,e2,z]=', x0)
    while(LoopValue==True and Iter <= 15):
        # print('Initial Value [c,e1,e2,z]=', x0)
        result = sp.minimize(maxerror,x0,method='Nelder-Mead', options={'xatol' : 1e-5 ,'disp':False, 'return_all' :False})
        if result.success or result.message=='Maximum number of function evaluations has been exceeded.':
            Iter += 1
            cf, e1f, e2f, zf = result.x
            lfun = result.fun
            if result.fun >= 100 or np.isnan(result.fun):
                if  x0i[3] == 359:
                    x0i[3] = 0
                    x0i[1]=x0i[1] if x0i[1]<= .5 *rf else 1
                    incf = [1,0,1,0]
                elif x0i[3] == 359:
                    x0i[0] = 0
                    x0i[2] = 0
                    incf=[0,.1,0,0]
                else:
                    incf= [0,0,0,.5]
                x0= [(a+b) % 360 for a, b in zip(x0i,incf)]
                x0i=x0
                LoopValue=True
            elif result.fun >=10:
                x0 = [(cf+.3) % 360, e1f if e1f >= 1 else 1, (e2f+.3) % 360, zf % 360]
                LoopValue = True
            elif result.fun >=1:
                if np.absolute(lfun-result.fun) <= 1e-3:
                    LoopValue = False
                else:
                    x0 = [cf+.1 % 360, e1f if e1f >= 1 else 1, e2f % 360, zf % 360]
                lfun=result.fun
                x0i = x0
        else:
            LoopValue=False
            print(result)
            raise ValueError(result.message)
    # print(result)
    # print('FInal Value [c,e1,e2,z]=', result.x)
    cf, e1f, e2f, zf = result.x
    Er, MEr = MarsEquantModel(cf, rf, e1f, e2f, zf, sf, oppositions)

    return cf, e1f, e2f, zf, Er, MEr
# End Q2 BestOebitInnerParams

# Assignment Q3 Return Value in s,errors,maxError Si=array(start,end) if Search Range Given
def bestS(rf, oppositions, Si=None):
    '''
    :param rf: Radius of Orbit Value Fixed
    :param oppositions:  Oppositions = (Daydiff, Latitude Angle) (12 x 2 ) Array
    :param Si: Range Of Si Search (Default:[360/(687-.5),360/(687+.5)])

    Logic:

    Searching minerror in range of Si in step of Precision_Control Factor(Default:30)
    deciding New Range Over

    New range =Si(Minerror Index to -1, Minerror Index to +1)

                            Mer Mer'' Mer'
    __________________________'_'_'____________________
            '       '       '      '      '       '
            S1      S1'     S1''  s2''   S2'     S2
    Loop 1 Range :(S1,S2) ==> Mer
    Loop 2 Range :(S1',S2') ==> Mer'
    Loop 3 Range :(S1'',S2'') ==?Mer''

    :return
        Si Value at Mer''
        Error =Error of Each Opposition
        Max Error= Max of Error

    '''
    Precision_Control=int(20)
    if Si == None:
        Timep=687
        Si = np.array([360/(Timep-1),360/(Timep+1)])
    for i in range(0,3):
        dis_s=np.linspace(Si[0],Si[1],Precision_Control)
        MErf=np.array([np.nan for i in range(0,Precision_Control)])
        # discretised search over si
        def maxerrors(sf):
            cf1, e1f1, e2f1, zf1, Erf1, MErf1 = bestOrbitInnerParams(rf, sf, oppositions)
            return MErf1
        for i in range(0,Precision_Control):
            # print(dis_s[i],rf)
            MErf[i] = maxerrors(dis_s[i])
        Si=[dis_s[MErf.argmin() if not(MErf.argmin()) else MErf.argmin()-1],dis_s[MErf.argmin() if (MErf.argmin() +1)==MErf.shape[0] else MErf.argmin() + 1]]
        OptS=dis_s[MErf.argmin()]
    cf, e1f, e2f, zf, Erf, MErf = bestOrbitInnerParams(rf, OptS , oppositions)
    print('S', OptS, 'MErf',MErf)
    return OptS,Erf,MErf
#End Q3 BestS

# Assignment Q4 Return r,errors,maxError =
def bestR(sf,oppositionsf,int_rf=None,Error_Enable=None):
    '''

    :param sf:Fixed Orbit Rotational Speed (Deg/Day)
    :param oppositionsf: Opposition Data
    :param int_rf: Initial Guess Of Orbit Radius (Default :5 AU)
    :param Error_Enable: Return Additional Parameter Of functional Optimized Successfully
                    (Return Error=False)  or Diverge (Return Error=True)

    Calculation Rf:
        1) Solving Equation Y=X *tan(Act_Line)  and (Y-ey)=(X-ex) *tan(Stroke)
            Simplified Equation X=(ey - ex * tan(Stroke)) / (tan(Act_Line) - tan(Stroke))
        2) Calculate Each (X,Y) distance from (Cx,Cy)
        3) New_Rf= is mean of Each distance between (X,Y) and (Cx,Cy)

    Search Control:
        calculate New_rf until the Error Start increasing continuously

        Implement flow:

            1)Appending Each iteration Value
                i) R_New (In r_Vari)
                ii) Max Error (In Mer)
                iii) Error_Change ( Append Value True if error_tol became positive else false)
                    Error_tol: it describe error change from last value
            2)Checking Divergence Condition:(Success than Error=True)
                i) If R_New decreasing ==> Break Internal Loop & Increment Factor += & Set Error True
                ii) If R_New Exponentially Increasing ==>Break Internal Loop & Increment Factor -= Set Error True
                iii) Above Condition Fail Than Set Error False to deactivate Outer Loop
            3) Check if more than  (LoopValue/2) (Default 25) no of times error tol become positive Inner While loop breaks
                    (Note : error_tol positive mean --error curve(Intial Decresing ,Makes Min, Started increasing) in last
                    part of started raising)
            4) Final Rf correspond to min Max Error (Find by Mer,r_Vari)
            5) Outer Loop Control By Error Value + (Max Search range of Rf(Defaul 20)) + (LoopValue= Max No. time innerloop
                Diverge)

    :return:
        Return 3 Parameter Rf (Optimized),Error,Max Error

        Additional Parameter Enable if Status of Convergence if Error_Change=True
    '''

    C_dis = 1
    oppositions_2x12 = np.array(oppositionsf).T
    if int_rf==None:
        int_rf = 5
    if Error_Enable==None:
        Error_Enable=False

    Loop_Value = 20 # Loop control Value if Diverge the Convergence
    error_change=np.array([False for i in range(0, Loop_Value)])
    Mer=np.array([360 for i in range(0,Loop_Value)])
    r_vari=np.array([int_rf for i in range(0,Loop_Value)])
    last_MEr1=0
    inc_fact=.2
    Error=True
    New_rf = int_rf
    Count_Error=1
    while (Error==True and 0<New_rf<20 and Count_Error<=Loop_Value):
        error_change = np.array([False for i in range(0, Loop_Value)])
        while(np.unique(error_change, return_counts=True)[1][0] > (int((Loop_Value)/2) +3) ):
            cf, e1f, e2f, zf, Er1f, MEr1=bestOrbitInnerParams(New_rf,sf,oppositionsf)
            error_tol= MEr1 - last_MEr1
            # print('Cal_rf', New_rf, 'MaxEr', MEr1, 'error_tol', error_tol, '\nerror_Change', error_change)
            last_MEr1=MEr1
            r_vari=np.append(r_vari,New_rf)
            r_vari=np.delete(r_vari,0)
            Mer = np.append(Mer, MEr1)
            Mer = np.delete(Mer, 0)
            error_change= np.append(error_change, True if error_tol >= 0 else False)
            error_change=np.delete(error_change, 0)
            # Setting Error Message
            R_Check = r_vari[-4:] - r_vari[-5:-1]
            if (np.any(R_Check >= 5) or np.all(R_Check < 0)):
                Count_Error +=1
                Error = True
                inc_fact += .2
                # print('BestR could Find R is Div. to given s Value, Try with diff. Parameter,\n'
                #       'Return initial value\n'
                #       'Error', Error)
                r_vari = np.delete(r_vari, -1)
                Mer = np.delete(Mer, -1)
                r_vari = np.insert(r_vari,0,int_rf)
                Mer = np.insert(Mer,0,360)
                break
            else:
                Error = False
            cx = C_dis * np.cos(np.radians(cf))
            cy = C_dis * np.sin(np.radians(cf))
            ex = e1f * np.cos(np.radians(e2f+zf))
            ey = e1f * np.sin(np.radians(e2f +zf))
            # Doted Line Angle At Z=0;
            Dotlinef = (oppositions_2x12[0] * sf + zf) % 360
            # print('Doted line Angle=', Dotlinef)
            X_Line = (ey-ex*np.tan(np.radians(Dotlinef))) / (np.tan(np.radians(oppositions_2x12[1]))-np.tan(np.radians(Dotlinef)))
            Y_Line = X_Line * np.tan(np.radians(oppositions_2x12[1]))
            dis_C = np.sqrt((X_Line - cx) ** 2 + (Y_Line - cy) ** 2)
            New_rf = np.mean(dis_C)
            # print(New_rf)
            # print(X_Line, '\n', Y_Line, '\n', dis_C, '\n', New_rf, MEr1)
        # print(r_vari,'\n',Mer)
        New_rf = int_rf + inc_fact
        # print(New_rf)
    # plt.plot(X_Line,Y_Line,'go')
    rf=r_vari[Mer.argmin()]
    # print(rf,Mer.argmin(),Mer[Mer.argmin()])

    cf, e1f, e2f, zf, Erf, MEr = bestOrbitInnerParams(rf, sf, oppositionsf)
    # print('Rf',rf, 'Mer',MEr)
    if Error_Enable:
        return rf, Erf, MEr,Error
    else:
        return rf, Erf, MEr
# End Code Q4 BestR

# Assignment Q5 Best Parameter
def bestMarsOrbitParams(oppositions):
    '''

    :param oppositions: 12 Opposition Input Data (Daydiff , Latitude)

    Logic:
    Simultaneously  running the BestR and BestS function with intial Guess rf(Default : 3) and sf(Default: 360/687) until constrain gets true
    Constrain:
        1) Max Error < 4 Min. (Implemneted using while condition)

        2) if Max Error change diverge in last 2 complete iteration in way
            i) Remain Constant
            ii) After each BestR and BestS use Max error get increased from its last iteration
                (Implemented by traking Error_Change last 4 Parameter)


    :return:
        Success :Return All parameter
        Fail:Return All parameter with Error Message Printed
    '''
    rf=3
    sf=360/687
    MEr1f=1
    error_change=np.array([0,0,0,0])
    last_Mer=1
    while(MEr1f>(4/60)):
        rf,Er1f,MEr1f,Error=bestR(sf,oppositions,int_rf=rf,Error_Enable=True)
        print('Using Best R:\n','s=',sf,'r=',rf,'Max Error',MEr1f,'\n')
        error_change=np.append(error_change,MEr1f-last_Mer)
        last_Mer=MEr1f
        print(error_change)
        sf, Er, MEr2f = bestS(rf, oppositions)
        print('Using Best S:\n', 's=', sf, 'r=', rf, 'Max Error', MEr2f, '\n')
        error_change = np.append(error_change, MEr2f - last_Mer)
        last_Mer = MEr2f
        print(error_change)
        if (np.all(error_change[-4:-1]>=0)):
            print('Optimization Can not converged')
            break
    cf, e1f, e2f, zf, Erf, MErf = bestOrbitInnerParams(rf,sf,oppositions)
    return rf,sf,cf,e1f,e2f,zf,Erf,MErf
# End Code of Q5 BestMarsOrbitParas


def main():
    # Data Import .CSV is in same Path
    Year, Month, Day, Hour, Minute, ZodiacIndex, Degree, Minute1, Second, LatDegree, LatMinute, ZodiacIndexAverageSun, DegreeMean, MinuteMean, SecondMean = \
        np.genfromtxt("data_mars_opposition_updated.csv", dtype='int', delimiter=",", skip_header=1, unpack=True, )
    # Please Change Path Accroding to your Execution Directory

    # Opposition Data & # Cal Daydiff
    Ref_Date = [Year[0], Month[0], Day[0], Hour[0], Minute[0]]
    Daydiff = np.array([0 for x in range(np.shape(Year)[0])], dtype='float')
    for i in range(1, np.shape(Year)[0]):
        To_Date = [Year[i], Month[i], Day[i], Hour[i], Minute[i]]
        delta = datetime(Year[i], Month[i], Day[i], Hour[i], Minute[i])-datetime(Year[0], Month[0], Day[0], Hour[0], Minute[0])
        Daydiff[i] = float(delta.total_seconds()/(24*3600))
    # print('Daydiff=',Daydiff)
    Solidline = ZodiacIndex * 30 + Degree + Minute1 / 60 + Second / (60 * 60)
    print('Actual Longitude (Solidline Angle)=', Solidline)
    print('Day Difference #Taking Reference as 1st longitude Datatime=',Daydiff)
    Opposition_Data = tuple((Daydiff[i],Solidline[i])for i in range(0,np.shape(Solidline)[0]))
    # End Code
    Opposition_Data
    #Intial Guess for Q1
    '''c = 120
    e1 = 1
    e2 = 93
    z = 57
    r = 10
    s = 360 / 687

    # Q1 Test Case
    print('Q1 Error(deg.) & absolute Max error (deg.) Test Case Result')
    Er1, MEr1 = MarsEquantModel(c, r, e1, e2, z, s, Opposition_Data)
    print('\nc=', c, '\ne1=', e1, '\ne2=', e2, '\nz=', z, '\nr=', r, '\ns=', s, '\nT=', 360 / s, '\nError=', Er1,
          '\nMaxError=', MEr1)
    # Q2 Test Case
    print('Q2 BestInnerOrbitPara  Test Case using Q1 Initial Para (r=10, S=360/687) ')
    c,e1,e2,z,Er1,MEr1=bestOrbitInnerParams(r,s,Opposition_Data)
    print('\nc=', c, '\ne1=', e1, '\ne2=', e2, '\nz=', z, '\nr=', r, '\ns=', s, '\nT=', 360 / s, '\nError=', Er1,
          '\nMaxError=', MEr1)

    # Q3 Test Case
    r=10
    print('Q3 Best S (Deg/Day) Test Case using Q1 Initial Para (r=10)')
    s,Er1,MEr1 = bestS(r,Opposition_Data)
    print('\nc=', c, '\ne1=', e1, '\ne2=', e2, '\nz=', z, '\nr=', r, '\ns=', s, '\nT=', 360 / s, '\nError=', Er1,
          '\nMaxError=', MEr1)

    # Q4 Test Case
    s=360/687
    print('Q4 Best R (respect to 1 AU B/W Sun and Orbit Center) Test Case using Q1 Initial Para (s=360/687)')
    r,Er1,MEr1 = bestR(s,Opposition_Data)
    print('\nc=', c, '\ne1=', e1, '\ne2=', e2, '\nz=', z, '\nr=', r, '\ns=', s, '\nT=', 360 / s, '\nError=', Er1,
          '\nMaxError=', MEr1)

    # Q5 Test Case
    print('Q5 Best All Para Test Case using initial Guess r=3, s=360/687 ')
    r, s, c, e1, e2, z, Er1, MEr1 = bestMarsOrbitParams(Opposition_Data)
    print('\nc=', c, '\ne1=', e1, '\ne2=', e2, '\nz=', z, '\nr=', r, '\ns=', s, '\nT=', 360 / s, '\nError=', Er1,
          '\nMaxError=', MEr1)
    PlotOut(Opposition_Data, z, c, r, e1, e2, s)

    plt.xlim(-(r + 3), r + 3)
    plt.ylim(-(r + 3), r + 3)
    plt.grid()
    plt.show()'''

if __name__== "__main__":
    main()



