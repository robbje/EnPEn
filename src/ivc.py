#!/usr/bin/env python2

import matplotlib.pyplot as plt
import json
import sys
import math

pos = 0.0

def read_file(filename):
    raw = file(filename).read()
    data = dict()
    try: 
        data = json.loads(raw)
    except: 
        raw += '"empty":[]}'
        data = json.loads(raw)
    return data

def diagnosis(data):
    indices = []
    x = []
    y = []
    z = []
    for idx in data:
        if idx == "empty":
            continue
        indices.append(float(idx))
        states = data[idx]["stt"]
        for s in states:
            x.append(s["x"])
            y.append(s["y"])
            z.append(s["z"])
    print "idx %f - %f" % (min(indices), max(indices))
    print "x %f - %f" % (min(x), max(x))
    print "y %f - %f" % (min(y), max(y))
    print "y %f - %f" % (min(z), max(z))

def sort(l1, l2):
    return (list(t) for t in zip(*sorted(zip(l1, l2))))

def get_xplot(data, index, setname, vecidx):
    x = []
    y = []
    for dataset in data["%.3f"%index][setname]:
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

def get_phi_for_idx(data, index, volidx):
    return data["%.3f" % index]["stt"][volidx]["v"][3]

def get_index_for_x(xvec, xval):
    for i,v in enumerate(xvec):
        if xval < v:
            return i

filename = "outdata/sweep.json"
if len(sys.argv) > 1:
    filename = sys.argv[1]
data = read_file(filename)
nv = data["meta"]["nv"]
charge = data["meta"]["charge"]
names = data["meta"]["species"]
D = data["meta"]["diffcoeff"]
c0 = data["meta"]["left_boundary"]

y = [[0]]*nv
for i in xrange(nv):
    [x, y[i]] = get_iyplot(data, 923, "flx", i)
for i,v in enumerate(y[nv-1]):
    y[nv-1][i] = math.fabs(y[nv-1][i])

#[x_x, c] = get_xplot(data, 0, "stt", 0)
#for j, w in enumerate(x_x):
#    if c[j] >= 2*c0[0]:
#        i1 = j
#
#v_real = []
#for i, v in enumerate(x):
#    [x_x, voltage] = get_xplot(data, v, "stt", 2)
#    v_real.append(voltage[i1] - voltage[0])
fig = plt.figure()
plot = plt.plot(x, y[nv-1], 'g-')
#plot = plt.plot(x, v_real, 'r-')
plt.title("Current at x=%g" % pos)
plt.show()
for i,v in enumerate(x):
     print "%g, %g" % (v, y[nv-1][i])


### CALCULATE TRANSPORTNUMBERS ###
t = []
legnames = []
for i, z in enumerate(charge):
    if z == 0:
        continue
    empty = []
    for j in xrange(len(x)):
        try:
            empty.append(abs(z*y[i][j]/y[nv-1][j]))
        except:
            empty.append(0)
    t.append(empty)
    legnames.append(names[i])

#print t[0][0], t[1][0], t[2][0]
#print c0[0], c0[1]

#V = int(sys.argv[2])
#print "%g" % (t[0][V]*c0[1]/(t[1][V]*c0[0]))

##### PLOT TRANSPORT NUMBERS ###
fig = plt.figure()
plots = []
for tn in t:
    plot, = plt.plot(x[1:], tn[1:])
    plots.append(plot)
plt.title("Transport numbers for species at x=%g" % pos)
plt.legend(plots,legnames)
plt.axes().set_ylim((0, 1))
plt.show()

#print "V, H+, OH-, Na+, Ca2+, Cl-"
#for i, v in enumerate(x):
#    print "%g, %g, %g, %g, %g" % (v, t[0][i], t[1][i], t[2][i], t[3][i])
### PLOT the location of 1-D(-)/D(+)Kw/[H+]^2 = 0 ###
#debye = 1.57781e-07
#potential = []
#positions = []
#for dataset in data:
#    if dataset == "meta" or dataset == "empty":
#        continue
#    [x,H] = get_xplot(data, float(dataset), "stt", 1)
#    pos = 0
#    for i, pos in enumerate(x):
#        val = 1-D[2]/D[1]*1e-14/H[i]**2
#        if val < 0:
#            positions.append((100-pos))
#            potential.append(float(dataset))
#            break
#[potential, positions] = sort(potential, positions)
#fig = plt.figure()
#plt.title("Distance of 1-a*Kw/[H+]^2==0 from interface in micron")
#plt.text(60,12,"Debye length = %gm" % debye)
#plt.plot(potential, positions, 'b-')
#plt.show()

