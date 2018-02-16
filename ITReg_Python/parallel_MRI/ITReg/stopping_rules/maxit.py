'''
% Parameter choice rule. Chooses last index

function stoprule = maxit(par)

stoprule.rec = [];
stoprule.stop = @stop;
stoprule.select_index = @select_index;
end

function [bool,stoprule] = stop(stoprule, F, xn,data,par)
    stoprule.rec = [stoprule.rec xn];
    if (par.it_step < par.N_max_it)
        bool = 0;
    else
        bool = 1;
    end;
end

function [maxit,res] = select_index(stoprule,F)
    n = size(stoprule.rec);
    maxit = n(2);
    res = stoprule.rec(:,n(2));
end
'''

#% Parameter choice rule. Chooses last index

def maxit(par):

    stoprule = dict()
    stoprule["rec"] = []
    stoprule["stop"] = stop
    stoprule["select_index"] = select_index
    return stoprule

def stop(stoprule, F, xn,data,par):
    stoprule["rec"].append(xn)
    if (par["it_step"] < par["N_max_it"]):
        rval = 0;
    else:
        rval = 1;

    return rval,stoprule

def select_index(stoprule,F):
    res = stoprule["rec"][-1]

    return maxit,res
