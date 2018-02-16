"""
function [stat,F] = error_stat_output(F,pcond,stat,x_k,y_obs,y_k,x_0,par)
if par.it_step==0
    stat.Xerr=[];
    stat.Yerr=[];
    if ~isempty(par.store_rec)
        stat.x = {};
    end
end
stat.Yerr = [stat.Yerr,sqrt(real((y_obs'-y_k') * F.applyGramY(F,(y_obs-y_k))))];
condPrint(pcond,'----------------------------------------------\n');
condPrint(pcond,'It.%i: resi=%1.3e',par.it_step,stat.Yerr(end));
if F.syntheticdata_flag
    stat.Xerr = [stat.Xerr,sqrt(real((F.xdag-x_k)' * F.applyGramX(F,(F.xdag-x_k))))];
    condPrint(pcond,', err=%1.3e',stat.Xerr(end));
    if isfield(F,'other_X_err')
        if (par.it_step==0) stat.Xerr2 = []; end
        for l = 1:size(F.other_X_err,2)
           err = F.other_X_err{1,l}(F,x_k);
           stat.Xerr2 = [stat.Xerr2,err];
           condPrint(pcond,[', ',F.other_X_err{2,l}]);
           condPrint(pcond,'=%1.3e',err);
        end
    end
    if ~isempty(find(par.plot_steps==par.it_step))
       F = F.plot(F,x_k,x_0,y_k,y_obs,par.it_step);
    end
else
    if ~isempty(find(par.plot_steps==par.it_step))
       F = F.plot(F,x_k,x_0,y_k,y_obs,par.it_step);
    end
end;
if ~isempty(intersect(par.store_rec,par.it_step))
         stat.x = [stat.x,x_k];
end
end
"""
import ITReg.utilities.condPrint as cndP
import numpy as np

def error_stat_output(F,pcond,stat,x_k,y_obs,y_k,x_0,par):
    if par["it_step"]==0:
        stat["Xerr"]=[]
        stat["Yerr"]=[]
        if len(par["store_rec"]) < 1:
            stat["x"] = dict()


    stat["Yerr"].append(np.sqrt(np.real((y_obs.conj().transpose()-y_k.conj().transpose()) * F["applyGramY"](F,(y_obs-y_k)))))
    cndP.condPrint(pcond,'----------------------------------------------\n')
    cndP.condPrint(pcond,'It.%i: resi=%1.3e',par["it_step"],stat["Yerr"][-1])
    if F["syntheticdata_flag"]:
        stat["Xerr"].append(np.sqrt(np.real((F["xdag"]-x_k).conj().transpose().dot(F["applyGramX"](F,(F["xdag"]-x_k))))))
        cndP.condPrint(pcond,', err=%1.3e',stat["Xerr"][-1])
        if 'other_X_err' in F.keys():
            if (par["it_step"]==0):
                stat["Xerr2"] = []
            for l in range(F["other_X_err"].shape[1]):
                #TODO:we don't have this yet
                #err = F.other_X_err{1,l}(F,x_k);
                #stat.Xerr2 = [stat.Xerr2,err];
                #cndP.condPrint(pcond,[', ',F.other_X_err{2,l}]);
                #cndP.condPrint(pcond,'=%1.3e',err);
                pass


        if len(np.nonzero(par["plot_steps"]==par["it_step"])) > 0:
            F = F["plot"](F,x_k,x_0,y_k,y_obs,par["it_step"])

    else:
        if len(np.nonzero(par["plot_steps"]==par["it_step"])) > 0:
            F = F["plot"](F,x_k,x_0,y_k,y_obs,par["it_step"])


    if len(np.intersect1d(par["store_rec"],par["it_step"]))>0:
         stat["x"].append(x_k)

    return stat,F
