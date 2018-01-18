#!/usr/bin/env python2

import math
import matplotlib.pyplot as plt
import sys
import json

def read_file(filename):
    raw = file(filename).read()
    try: 
        return json.loads(raw)
    except:
        raw += '"empty":[]}'
    return json.loads(raw)
def sort(l1, l2):
    return (list(t) for t in zip(*sorted(zip(l1, l2))))
def get_xplot(data, index, setname, vecidx):
    x = []
    y = []
    for dataset in data["%.3f"%index][setname]:
        if dataset["x"] >= 100:
            continue
        x.append(dataset["x"])
        y.append(dataset["v"][vecidx])
    return sort(x, y)

def get_iyplot(data, volidx, setname, vecidx):
    x = []
    y = []
    global pos
    for dataset in data:
        if dataset == "empty" or dataset == "meta":
            continue
        pos = data[dataset][setname][volidx]["x"]
        x.append(float(dataset))
        d = data[dataset][setname][volidx]
        y.append(-1*d["v"][vecidx])
    return sort(x, y)

def model_i_lim(D, z, L, c0):
    kappa = 0
    D_mean = 0
    for i, v in enumerate(c0):
        kappa += z[i]**2*D[i]*c0[i]
        D_mean += D[i]
    D_mean /= float(len(c0))
    return kappa/L*D_mean/D[1]


fname = sys.argv[1]
V = float(sys.argv[2])
data = read_file(fname)
nv = data["meta"]["nv"]
D0 = 1.0e-9
D = [Di/D0 for Di in data["meta"]["diffcoeff"]]
z = data["meta"]["charge"]
c = [[]]*(nv-1)
for i in xrange(nv-1):
    [x,c[i]] = get_xplot(data, V, "stt", i)
[x3,cd] = get_xplot(data, V, "flx", nv-1)
i_lim_DNS = cd[0]
c0 = [c[i][0] for i in xrange(nv-1)]
c_fix = 0.001
print "c_fix = %g" % c_fix
print "alpha = %g" % (c0[0]/c_fix)
print c0

i_lim_ANA = model_i_lim(D, z, 100, c0)
print "Current density:"
print "DNS: %g" % i_lim_DNS
print "Model: %g" % i_lim_ANA
print "Rel. Err: %g" % (abs(i_lim_ANA-i_lim_DNS)/i_lim_DNS)

[x, f] = get_iyplot(data, 1055, "flx", nv-1)
for i, v in enumerate(f):
    f[i] = -v
cd_ana = [i_lim_ANA]*len(x)
fig = plt.figure()
plt.title("Simulated and analytical limiting current density")
plt.plot(x, f)
plt.plot(x, cd_ana)
#plt.show()
