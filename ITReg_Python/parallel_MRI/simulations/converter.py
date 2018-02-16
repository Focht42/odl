import scipy.io
import numpy as np
fn = "ksp2x2.mat"
data = scipy.io.loadmat(fn)
print(data)
for i in data:
    if '__' not in i and 'readme' not in i:
        np.savetxt(("file.csv"),data[i],delimiter=',')
