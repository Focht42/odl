import StarCurve
import numpy as np
def GenCurve(name):

    if name == "kite":
        curve["name"] = name
        curve["bd_eval"] = GenCurve_bd_eval
    else:
        curve = StarCurve.StarCurve(name)
    return curve


def GenCurve_bd_eval(curve,n,der):
    t = 2*np.pi*np.array(list(range(0,n)))/n
    #TODO: set the correct curve name function here! cnf = curve.name
    cnf = kite
    curve["z"] = cnf(t,0)

    if der>=1:
        curve["zp"] = cnf(t,1)
        curve["zpabs"] = np.abs(curve["zp"][0,]) + np.abs(curve["zp"][1,])
        # outer normal vector
        curve["normal"] = np.vstack((curve["zp"][1,],-curve["zp"][0,]))

    if der>=2:
        curve["zpp"] = cnf(t,2)

    if der>=3:
        curve["zppp"] = cnf(t,3)
    if der>3
        print('only derivatives up to order 3 implemented')
    return curve

def kite(t,der):
    if der == 0:
        res = np.vstack((np.cos(t)+0.65*np.cos(2*t)-0.65, 1.5*np.sin(t)))
    elif der == 1:
        res = np.vstack((-np.sin(t)-1.3*np.sin(2*t), 1.5*np.cos(t)))
    elif der == 2:
        res = np.vstack((-np.cos(t)+2.6*np.cos(2*t), -1.5*np.sin(t)))
    elif der == 3:
        res = np.vstack((np.sin(t)+5.2*np.sin(2*t), -1.5*np.cos(t)))
    else:
        #TODO: Error!
        print('derivative not implemented');
    return res
