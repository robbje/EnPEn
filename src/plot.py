#!/usr/bin/env python2

import matplotlib.pyplot as plt
import json
import sys
import math

def read_file(filename):
    raw = file(filename).read()
    try: 
        return json.loads(raw)
    except:
        raw += '"empty":[]}'
    return json.loads(raw)

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

def get_emp(data, index):
    x = []
    y = []
    for dataset in data["%.3f"%index]["stt"]:
        x.append(dataset["x"])
        y.append(dataset["emp"])
    return sort(x, y)

def get_iyplot(data, volidx, setname, vecidx):
    x = []
    y = []
    for dataset in data:
        if dataset == "empty":
            continue
        x.append(float(dataset))
        d = data[dataset][setname][volidx]
        y.append(d["v"][vecidx])
    return sort(x, y)
def get_idx_at(x, x0):
    for i, v in enumerate(x):
        if x0 < v:
            return i

filename = "outdata/single.json"
if len(sys.argv) >= 2:
    filename = sys.argv[1]

data = read_file(filename)
nv = data["meta"]["nv"]
charge = data["meta"]["charge"]
names = data["meta"]["species"]
D = data["meta"]["diffcoeff"]
F = 96485.3
R = 8.31446
T = 298.15

V = 0
if len(sys.argv) == 3:
    V = float(sys.argv[2])

### GET ALL STATES ###
y = [[0]]*nv
for i in xrange(nv):
    [x,y[i]] = get_xplot(data, V, "stt", i)

#### PLOT EMP ###
#[asdf,emp] = get_emp(data, V)
#fig = plt.figure()
#ans, = plt.plot(x, emp)
#plt.title("EMP at %f" % V)
#plt.show()
#points = [15.4+0.5*15.4]
#indices = map(lambda p: get_idx_at(x,p), points)
#mid = get_idx_at(x, 150.0)
#print "%g, %g" % (V, (emp[mid]-emp[0])*1000.0*R*T*1e-5)
#for i,v in enumerate(x):
#    print "%g, %g" % (v, R*T*1000.0*emp[i]*1e-5)

### Print pH at a certain position ###
#print "%g, %g" % (V, -math.log10(y[1][get_idx_at(x, 90)]))
#print x[2234]
#print get_idx_at(x, 250)

### PLOT pH ###
#pH = []
#for v in y[1]:
#    try:
#        pH.append(-math.log10(v))
#    except:
#        pH.append(-1)
#fig = plt.figure()
#ans, = plt.plot(x, pH)
#plt.title("pH at %f" % V)
#plt.show()
#for i, v in enumerate(x):
#    print "%g, %g" % (v, pH[i])

### PLOT WATER IONS ###
fig = plt.figure()
plots = []
for i in range(nv-5,nv-3):
    plot, = plt.plot(x, y[i])
    plots.append(plot)
plt.title("Concentrations at %f" % V)
plt.legend(plots, names[nv-5:nv-3])
plt.show()

#for i, v in enumerate(x):
#    if v > 199.9 and v < 200.5:
#        print "%.12f, %.12f, %.12f, %.12f" % (v, y[0][i], y[1][i], y[2][i])

#def phi(x, phi0=-1):
#    
#    k = (x*1e-6/3.03232e-08)
#    h = (1-math.exp(0.5*phi0))/(1+math.exp(0.5*phi0))
#    phi = 2 * math.log((1-h*math.exp(0.5*phi0)*math.exp(-k))/(1+math.exp(0.5*phi0)*math.exp(-k)))
#    return phi
#
#y_anal = [[],[],[]]
#for i,v in enumerate(x):
#    p = phi(100.0-v)
#    y_anal[0].append(1e-3*math.exp(-1.0*p))
#    y_anal[1].append(1e-3*math.exp(1.0*p))
#    y_anal[2].append(p)

### PLOT ALL NON-WATER SPECIES ###
#fig = plt.figure()
#plots = []
#for i in range(0,nv-1):
#    plot, = plt.plot(x, y[i])
#    plots.append(plot)
#plt.title("Concentrations at %f" % V)
#plt.legend(plots, names[0:nv-1])
#plt.show()

rho = []
new_x = []
for i,v in enumerate(x):
    if v < 20.0: continue
    new_x.append(v)
    new_rho = y[0][i]-y[1][i]
    rho.append(new_rho)
    #print '%.10f, %.10f' % (v, new_rho)

#fig = plt.figure()
#plt.plot(new_x, rho)
#plt.show()


#for i,v in enumerate(x):
#    print '%.8f, %.8f, %.8f' % (v, y[0][i], y[1][i])
##print "x, %s, %s, %s" % (names[1],names[3],names[4])
##print "mu, M"
#for i,v in enumerate(x):
#    if v > 101.0:
#        continue
#    print "%.8f, %.8f, %.8f" % (v, y[0][i], y[1][i])


### PLOT POTENTIAL ###
#fig = plt.figure()
#phi, = plt.plot(x, y[nv-1])
#plt.title("Potential at %f" % V)
#plt.legend([phi], [names[nv-1]])
#plt.show()
#
#for i, xval in enumerate(x):
#    if xval > 98 and xval < 101:
#        print "%.8f, %.8f" % (xval, y[nv-1][i])

### PLOT CONDUCTIVITY ###
#sigma = []
#for i, pos in enumerate(x):
#    sig = 0
#    for k, z in enumerate(charge):
#        sig += F**2/(R*T)*z**2*D[k]*y[k][i]*1000
#    if pos >= 100.0 and pos <= 200.0:
#        sig *= 0.2
#    sigma.append(sig*10)
#fig = plt.figure()
#sigplot, = plt.plot(x, sigma)
#plt.title("Conductivity [mS/cm]")
#plt.show()

### PLOT SOURCE TERM
#R = []
#k_f = 1e2
#k_b = 1e16
#N = 20 
#for i, pos in enumerate(x):
#    src = 0
#    total_x = 0
#    for k in range(i-N,i+N+1):
#        if k < N or k >= len(x)-N:
#            continue
#        dx = 0.5*abs(x[k-1]-x[k]) + 0.5 * abs(x[k+1]-x[k])
#        src += dx*(k_f - k_b*(y[1][k]*y[2][k]))
#        total_x += dx
#    src /= total_x
#    R.append(src)
#fig = plt.figure()
#hplt, = plt.plot(x, R)
#plt.title("Source term");
#plt.show()

#### CALCULATE SPACE CHARGE ###
#ix = []
#rho = []
#for i, pos in enumerate(x):
#    r = 0
#    s = 0
#    for j, z in enumerate(charge):
#        r += z*y[j][i]
#        s += 0.5*z**2*y[j][i]
#    if pos >= 100.0 and pos <= 200.0:
#        r -= 1.0
#    rho.append(r)
#    ix.append(s)
#fig = plt.figure()
#strength, = plt.plot(x, ix, 'b-')
#rhoplt, = plt.plot(x, rho, 'r-')
#plt.legend([strength,rhoplt],["Ionic Strength", "Charge density"])
#plt.title("rho, Ix at %f" % V)
#plt.show()
#
#devw = []
#for i, v in enumerate(x):
#    kw = (1e-14)-y[1][i]*y[2][i]
#    devw.append(kw)
#print "Maximum rel. deviation from Water equilibrium: %g" % max(devw)

### GET ALL FLUXES ### 
#f = [[0]]*nv
#for i in xrange(nv):
#    [x2,f[i]] = get_xplot(data, V, "flx", i)

### KIDBRANE VALIDATION
#def flow_analytic(x, j, u=0.01, c0=1e-3):
#    return (c0 - 0.5 * j / u) * math.exp(u * x) + 0.5 * j / u
#
#j = f[nv-1][20]
#c = []
#for i, v in enumerate(x):
#    c.append(flow_analytic(v, j))
#
#fig = plt.figure()
#plots = []
#for i in range(0,nv-1):
#    plot, = plt.plot(x, y[i])
#    plots.append(plot)
#plot, = plt.plot(x, c)
#plots.append(plot)
#plt.title("Concentrations at %f" % V)
#plt.legend(plots, names[0:nv-1]+["analytical"])
#plt.show()
#for i,v in enumerate(x):
#    if v <= 100:
#        print "%g, %g, %g" % (v, y[0][i], y[1][i])

#### CALCULATE TRANSPORTNUMBERS ###
#t = []
#legnames = []
#for i, z in enumerate(charge):
#    if z == 0:
#        continue
#    empty = []
#    for j in xrange(len(x2)):
#        try:
#            empty.append(z*f[i][j]/f[nv-1][j])
#        except:
#            empty.append(0)
#    t.append(empty)
#    legnames.append(names[i])

### PLOT TRANSPORT NUMBERS ###
#fig = plt.figure()
#plots = []
#for tn in t:
#    plot, = plt.plot(x2, tn)
#    plots.append(plot)
#plt.title("Transport numbers for species over domain at %f" % V)
#plt.legend(plots,legnames)
#plt.axes().set_ylim((0, 1))
#plt.show()

#### PLOT SOME FLUXES ###
#fig = plt.figure()
#plots = []
#legnames = []
#for i in range(nv-1,nv):
#    for j, v in enumerate(f[i]):
#        f[i][j] = f[i][j]
#    plot, = plt.plot(x2, f[i])
#    plots.append(plot)
#    legnames.append(names[i])
#plt.title("Fluxes for species over domain at %f" % V)
#plt.legend(plots, legnames)
#plt.show()

#fig = plt.figure()
#plot,= plt.plot(x2, f[nv-1])
#plt.title("current density")
#plt.show()
