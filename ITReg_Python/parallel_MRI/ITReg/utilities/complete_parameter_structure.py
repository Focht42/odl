def complete_parameter_structure(p,p_ref):

    """ compares parameters specified in a dict structure p to a reference structure
    p_ref and returns a parameter structure p_new in which missing entries in p
    are taken from p_ref"""
    #print("DBG p: ", p,"\np_ref: ", p_ref)
    p_new = p
    if not isinstance(p_ref,dict):
        return p_new
    for k in p_ref.keys():
        v = p_ref[k]
        if not isinstance(p_new,dict):
            p_new[k] = v
        elif not k in p_new.keys():
            p_new[k] = v
        else:
            #if the value is a dict we need recursive application
            p_new[k] = complete_parameter_structure(p_new[k],v)

    #TODO: make this real warnings
    for k in p_new.keys():
        if not k in p_ref.keys():
            print("WARNING: No reference value specified for parameter ",k," !")

    return p_new
