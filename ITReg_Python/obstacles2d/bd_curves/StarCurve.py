import numpy as np

def StarCurve(name):
    curve["name"] = name;
    curve["bd_eval"] = StarCurve_bd_eval
    curve["radial"] = StarCurve_radial
    return curve

def StarCurve_bd_eval(curve,n,der):
    t=2*np.pi*np.array(list(range(0:n)))/n
    cost = np.cos(t)
    sint = np.sin(t)
    curve["q"] = np.zeros((der+1)*n).reshape(der+1,n)
    # remember that matlab indexes from 1 and not from 0
    #TODO: get the curve.name function!
    def FindMe(t,j):
        print("FixMe")
        return 0
    cnf = FindMe
    for j in range(0,der+1):
        curve["q"][j,] = cnf(t,j)
    q = curve["q"]
    curve["z"] = np.vstack((q[0,]*cost,q[0,]*sint))
    if der>1:
        curve["zp"] = np.vstack((q[1,]*cost - q[0,]*sint, q[1,]*sint + q[0,]*cost))
        curve["zpabs"] = np.sqrt(curve["zp"][0,]**2 + curve["zp"][1,]**2)
        #outer normal vector
        curve["normal"] = np.vstack((curve["zp"][1,], -curve["zp"][0,]]))

    if der>=1:
        curve["zp"] = np.vstack((q[1,]*cost - q[0,]*sint, q[1,]*sint + q[0,]*cost))
        curve["zpabs"] = np.sqrt(curve["zp"][0,]**2 + curve["zp"][1,]**2)
        #outer normal vector
        curve["normal"] = np.vstack((curve["zip"][1,],-curve["zp"][0,]))

    if der>=2:
        curve["zpp"] = np.vstack((q[2,]*cost - 2*q[1,]*sint - q[0,]*cost,
                                q[2,]*sint + 2*q[1,]*cost - q[0,]*sint))

    if der>=3:
        curve["zppp"] = np.vstack((q[3,]*cost - 3*q[2,]*sint - 3*q[1,]*cost + q[0,]*sint,
                                q[3,]*sint + 3*q[2,]*cost - 3*q[1,]*sint - q[0,]*cost))

    if der>3:
        #TODO: Error
        print('only derivatives up to order 3 implemented')

    return curve

def StarCurve_radial(curve,n):
    #TODO: Try/Catch with function usage
    t=2*np.pi*np.array(list(range(0:n)))/n
    rad = curve["name"](t,0)
    return rad

def peanut(t,der):
    if der == 0:
        res = 1/2*(3*np.cos(t)**2+1)**(1/2)
    elif der == 1:
        res = -3/2/(4*np.cos(t)**2+sin(t)**2)**(1/2)*np.cos(t)*np.sin(t)
    elif der == 2:
        res = -3/2*(3*np.cos(t)**4+2*np.cos(t)**2-1)/(3*np.cos(t)**2+1)**(3/2)
    elif der == 3:
        res = 3/2*np.cos(t)*np.sin(t)*(9*np.cos(t)**4+6*np.cos(t)**2+13)/(3*np.cos(t)**2+1)**(5/2)
    else:
        #TODO: Error
        print('derivative not implemented')
    return res

def round_rect(t,der):
    #TODO: Format this nicer
    co = 2/3
    if der == 0:
        res = (np.sin(t)**10 + (co*np.cos(t))**10)**(-0.1)
    elif der == 1:
        res = -1/10/(np.sin(t)**10+co**10*np.cos(t)**10)**(11/10)*(10*np.sin(t)**9*np.cos(t)-10*co**10*np.cos(t)**9*np.sin(t))

    elif der == 2:
        res = 11/100/(np.sin(t)**10+co**10*np.cos(t)**10)**(21/10)*(10*np.sin(t)**9*np.cos(t)-10*co**10*np.cos(t)**9*np.sin(t))**2-1/10/(np.sin(t)**10+co**10*np.cos(t)**10)**(11/10)*(90*np.sin(t)**8*np.cos(t)**2-10*np.sin(t)**10+90*co**10*np.cos(t)**8*np.sin(t)**2-10*co**10*np.cos(t)**10)
    elif der == 3:
        res = -231/1000/(np.sin(t)**10+co**10*np.cos(t)**10)**(31/10)*(10*np.sin(t)**9*np.cos(t)-10*co**10*np.cos(t)**9*np.sin(t))**3+33/100/(np.sin(t)**10+co**10*np.cos(t)**10)**(21/10)**(10*np.sin(t)**9*np.cos(t)-10*co**10*np.cos(t)**9*np.sin(t))*(90*np.sin(t)**8*np.cos(t)**2-10*np.sin(t)**10+90*co**10*np.cos(t)**8*np.sin(t)**2-10*co**10*np.cos(t)**10)-1/10/(np.sin(t)**10+co**10*np.cos(t)**10)**(11/10)*(720*np.sin(t)**7*np.cos(t)**3-280*np.sin(t)**9*np.cos(t)-720*co**10*np.cos(t)**7*np.sin(t)**3+280*co**10*np.cos(t)**9*np.sin(t))
    else:
        #TODO Error
        print('derivative not implemented');
    return res

def apple(t,der):

    if der == 0:
        res = (0.5+0.4*np.cos(t)+0.1*np.sin(2*t))/(1+0.7*np.cos(t))
    elif der == 1:
        res = (-2/5*np.sin(t)+1/5*np.cos(2*t))/(1+7/10*np.cos(t))+7/10*(1/2+2/5*np.cos(t)+1/10*np.sin(2*t))/(1+7/10*np.cos(t))**2*np.sin(t)
    elif der == 2:
        res = (-2/5*np.cos(t)-2/5*np.sin(2*t))/(1+7/10*np.cos(t))+7/5*(-2/5*np.sin(t)+1/5*np.cos(2*t))/(1+7/10*np.cos(t))**2*np.sin(t)+49/50*(1/2+2/5*np.cos(t)+1/10*np.sin(2*t))/(1+7/10*np.cos(t))**3*np.sin(t)**2+7/10*(1/2+2/5*np.cos(t)+1/10*np.sin(2*t))/(1+7/10*np.cos(t))**2*np.cos(t)
    elif der == 3:
        res = (2/5*np.sin(t)-4/5*np.cos(2*t))/(1+7/10*np.cos(t))+21/10*(-2/5*np.cos(t)-2/5*np.sin(2*t))/(1+7/10*np.cos(t))**2*np.sin(t)+147/50*(-2/5*np.sin(t)+1/5*np.cos(2*t))/(1+7/10*np.cos(t))**3*np.sin(t)**2+21/10*(-2/5*np.sin(t)+1/5*np.cos(2*t))/(1+7/10*np.cos(t))**2*np.cos(t)+1029/500*(1/2+2/5*np.cos(t)+1/10*np.sin(2*t))/(1+7/10*np.cos(t))**4*np.sin(t)^3+147/50*(1/2+2/5*np.cos(t)+1/10*np.sin(2*t))/(1+7/10*np.cos(t))**3*np.sin(t)*np.cos(t)-7/10*(1/2+2/5*np.cos(t)+1/10*np.sin(2*t))/(1+7/10*np.cos(t))**2*np.sin(t)
    else:
        #TODO Error
        print('derivative not implemented')
    return res

def three_lobes(t,der):

    if der == 0:
        res = 0.5 + 0.25*np.exp(-np.sin(3*t)) - 0.1*np.sin(t)
    elif der == 1:
        res = -3/4*np.cos(3*t)*np.exp(-np.sin(3*t))-1/10*np.cos(t)
    elif der == 2:
        res = 9/4*np.sin(3*t)*np.exp(-np.sin(3*t))+9/4*np.cos(3*t)**2*np.exp(-np.sin(3*t))+1/10*np.sin(t)
    elif der == 3:
        res = 27/4*np.cos(3*t)*np.exp(-np.sin(3*t))-81/4*np.sin(3*t)*np.cos(3*t)*np.exp(-np.sin(3*t))-27/4*np.cos(3*t)**3*np.exp(-np.sin(3*t))+1/10*np.cos(t)
    else:
        #TODO: Error
        print('derivative not implemented')
    return res

def pinched_ellipse(t,der):

    if der == 0:
        res = 3/2*np.sqrt(1/4*np.cos(t)**2 + np.sin(t)**2)
    elif der == 1:
        res = 9/4/(-3*np.cos(t)**2+4)**(1/2)*np.cos(t)*np.sin(t)
    elif der == 2:
        res = 9/4*(3*np.cos(t)**4-8*np.cos(t)**2+4)/(3*np.cos(t)**2-4)/(-3*np.cos(t)**2+4)**(1/2)
    elif der == 3:
        res = -9/4*np.cos(t)*np.sin(t)*(9*np.cos(t)**4-24*np.cos(t)**2+28)/(3*np.cos(t)**2-4)**2/(-3*np.cos(t)**2+4)**(1/2)
    else:
        #TODO Error
        print('derivative not implemented')
    return res

def smoothed_rectangle(t,der):
function res = smoothed_rectangle(t,der)
    if der == 0:
        res = (np.cos(t)**10 +2/3*np.sin(t)**10)**(-1/10)
    elif der == 1:
        res = -1/10/(np.cos(t)**10+2/3*np.sin(t)**10)**(11/10)*(-10*cos(t)**9*np.sin(t)+20/3*np.sin(t)**9*np.cos(t))
    elif der == 2:
        res = 11/100/(np.cos(t)**10+2/3*np.sin(t)**10)**(21/10)*(-10*np.cos(t)**9*np.sin(t)+20/3*np.sin(t)**9*np.cos(t))**2-1/10/(np.cos(t)**10+2/3*np.sin(t)**10)**(11/10)*(90*np.cos(t)**8*np.sin(t)**2-10*np.cos(t)**10+60*np.sin(t)**8*np.cos(t)**2-20/3*np.sin(t)**10)
    elif der == 3:
        res = -231/1000/(np.cos(t)**10+2/3*np.sin(t)**10)**(31/10)*(-10*np.cos(t)**9*np.sin(t)+20/3*np.sin(t)**9*np.cos(t))**3+33/100/(np.cos(t)**10+2/3*np.sin(t)**10)**(21/10)*(-10*np.cos(t)**9*np.sin(t)+20/3*np.sin(t)**9*np.cos(t))*(90*np.cos(t)**8*np.sin(t)**2-10*np.cos(t)**10+60*np.sin(t)**8*np.cos(t)**2-20/3*np.sin(t)**10)-1/10/(np.cos(t)**10+2/3*np.sin(t)**10)**(11/10)*(-720*np.cos(t)**7*np.sin(t)**3+280*np.cos(t)**9*np.sin(t)+480*np.sin(t)**7*np.cos(t)**3-560/3*np.sin(t)**9*np.cos(t));
    else:
        #TODO: Error
        print('derivative not implemented')
return res

def nonsym_shape(t,der):

    if der == 0:
        res =(1 + 0.9*np.cos(t) + 0.1*np.sin(2*t))/(1 + 0.75*np.cos(t))
    elif der == 1:
        res = 4/5*(-3*np.sin(t)+8*np.cos(t)**2-4+3*np.cos(t)**3)/(16+24*np.cos(t)+9*np.cos(t)**2)
    elif der == 2:
        res = -4/5*(12*np.cos(t)-9*np.cos(t)**2+64*np.sin(t)*np.cos(t)+36*np.sin(t)*np.cos(t)**2+9*np.sin(t)*np.cos(t)**3+24*np.sin(t)+18)/(64+144*np.cos(t)+108*np.cos(t)**2+27*np.cos(t)**3)
    elif der == 3:
        res = -4/5*(144*np.sin(t)*np.cos(t)+114*np.sin(t)-40+240*np.cos(t)**3+192*np.cos(t)-27*np.sin(t)*np.cos(t)**2+368*np.cos(t)**2+144*np.cos(t)**4+27*np.cos(t)**5)/(256+768*np.cos(t)+864*np.cos(t)**2+432*np.cos(t)**3+81*np.cos(t)**4)
    else:
        #TODO: Error
        print('derivative not implemented')
    return res

def circle(t,der):
   if der==0:
       res = np.ones(len(t))
   else:
       res= zeros(len(t))
   return res
