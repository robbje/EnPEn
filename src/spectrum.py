#!/usr/bin/env python2

import matplotlib.pyplot as plt
import numpy as np
import scipy
import scipy.signal
from scipy.optimize import leastsq
import json
import sys
import math
import argparse
import os

class IVC:
    def __init__(self, fname, pos=0.0):
        self.filename = fname
        self.data = self.parse(fname) 
        self.meta = self.data['meta']

    def parse(self, fname):
        rd = file(fname).read()
        data = {}
        try:
            data = json.loads(rd)
        except:
            rd += '"empty":[]}'
            data = json.loads(rd)
        return data

    def sort(self, l1, l2):
        return (list(t) for t in zip(*sorted(zip(l1, l2))))

    def xplot(self, index, setname, vecidx):
        x = []
        y = []
        for d in self.data['%.3f' % index][setname]:
            x.append(d['x'])
            y.append(d['v'][vecidx])
        return self.sort(x, y)

    def yplot(self, volidx, setname, vecidx):
        x = []
        y = []
        for d in self.data:
            if d in ['empty','meta']: continue
            x.append(float(d))
            ds = self.data[d][setname][volidx]
            y.append(ds['v'][vecidx])
        return self.sort(x, y)

    def plot_ivc(self, show=True):
        x, f = self.yplot(1055, 'flx', self.meta['nv']-1)
        f = np.array(f)*-1.0
        if show:
            fig = plt.figure()
            plt.plot(x,f)
            plt.show()
        return x, f

    def plot_conc(self, V, show=True, indices=[]):
        y = [[0]]*(self.meta['nv'])
        for i in xrange(self.meta['nv']):
            x, y[i] = self.xplot(V, 'stt', i)
        if indices == []: 
            indices = xrange(self.meta['nv'])
        names = [self.meta['species'][i] for i in indices]
        if show:
            fig = plt.figure()
            plots = []
            for i in indices:
                plot, = plt.plot(x, y[i])
                plots.append(plot)
            plt.legend(plots, names)
            plt.show()
        return x, [y[i] for i in indices], names


class Spectrum:
    def __init__(self, fnames):
        self.spectrum = [EIS(*self.parse(f)) for f in fnames]
        self.spectrum.sort(key = lambda x: x.frequency)
        if len(self.spectrum) < 2:
            self.f, Z = self.spectrum[0].impedance()
        else:
            self.f, Z = map(list,map(None,*[z.impedance() for z in self.spectrum]))
        self.omega = 2 * math.pi * np.array(self.f)
        self.Re = np.real(Z)
        self.Im = np.imag(Z)
        self.eqc = lambda omega: 0
        self.p = []

    def parse(self, fname):
        raw = file(fname).readlines()
        try:
            meta = json.loads("".join(raw[:2]))
        except ValueError as e:
            raise ValueError("%s: %s" % (fname, e.message))
        data = []
        for d in raw[3:-1]:
            data.append([float(x.strip()) for x in d.split(',')])
        if not data:
            raise ValueError("%s: Empty file" % fname)
        return map(list,map(None,*data)) + [meta]

    def interpolate(self, edf):
        self.Re = np.interp(edf, self.f, self.Re)
        self.Im = np.interp(edf, self.f, self.Im)
        self.f = np.array(edf)
        self.omega = 2 * math.pi * self.f

    def csv(self):
        output = ""
        output += "%i\r\n" % len(self.f)
        for i,v in enumerate(self.f):
            output += "%.14G %.14G %.14G\r\n" % (self.Re[i], -1.0*self.Im[i], v)
        return output

    def plot_nyquist(self, show=True, eqc=False, data=True):
        extrema = self.get_extrema()
        if show:
            plt.figure()
            plt.gca().invert_yaxis()
            plt.gca().set_aspect('equal')
            if data:
                plt.plot(self.Re, self.Im, 'b-o')
                for x in extrema:
                    plt.annotate('%g'%x[0],
                                 xy=(x[1],x[2]), xycoords='data')
                    plt.plot(x[1],x[2],'ko')
            if eqc:
                z = [self.eqc(w) for w in self.omega]
                plt.plot(np.real(z), np.imag(z), 'rx-')
            plt.title("Nyquist diagram")
            plt.show()
        return self.Re, self.Im

    def plot_impedance(self, eqc=False, data=True):
        plt.figure()
        if data:
            re_plot, = plt.semilogx(self.omega, self.Re,'b+-')
            im_plot, = plt.semilogx(self.omega, self.Im,'g+-')
        if eqc:
            z = [self.eqc(w) for w in self.omega]
            plt.semilogx(self.omega, np.real(z),'rx-')
            plt.semilogx(self.omega, np.imag(z),'rx-')
        plt.legend([re_plot,im_plot],['Real','Imag'])
        plt.show()

    def get_extrema(self):
        im = self.Im
        indices = scipy.signal.argrelextrema(im,np.greater)
        extrema = [(self.omega[i],self.Re[i], self.Im[i]) for i in indices[0]]
        im = self.Im * -1.0
        indices = scipy.signal.argrelextrema(im,np.greater)
        extrema += [(self.omega[i],self.Re[i], self.Im[i]) for i in indices[0]]
        return extrema

    def add_rc(self, R, C, parallel=False):
        self.p += [R, C]
        n = len(self.p)
        i = [n-2, n-1]
        prefix = "|| " if parallel else "+ "
        RC = lambda omega: self.p[i[0]] / \
                           (1 + 1j*omega*self.p[i[0]]*self.p[i[1]])
        tmp = self.eqc
        if parallel:
            self.eqc = lambda omega: 1/(1/tmp(omega) + 1/RC(omega))
        else:
            self.eqc = lambda omega: tmp(omega) + RC(omega)
        print "%sRC(%g,%g) at [%i, %i]" % (prefix, R, C, i[0], i[1])
        return i

    def add_r(self, R, parallel=False):
        self.p += [R]
        n = len(self.p)
        i = [n-1]
        tmp = self.eqc
        prefix = "|| " if parallel else "+ "
        if parallel:
            self.eqc = lambda omega: 1/(1/tmp(omega)+1/self.p[i[0]])
        else:
            self.eqc = lambda omega: tmp(omega) + self.p[i[0]]
        print "%sR(%g) at [%i]" % (prefix, R, i[0])
        return i

    def add_cpe(self, Q, N, parallel=False):
        self.p += [Q,N]
        n = len(self.p)
        i = [n-2,n-1]
        tmp = self.eqc
        prefix = "|| " if parallel else "+ "
        CPE = lambda omega: 1/(self.p[i[0]]*np.power(1j*omega,self.p[i[1]]))
        if parallel:
            self.eqc = lambda omega: 1.0/(1.0/tmp(omega)+1.0/CPE(omega))
        else:
            self.eqc = lambda omega: tmp(omega) + CPE(omega)
        print "%sCPE(%g,%g) at [%i, %i]" % \
                (prefix, Q, N, i[0], i[1])
        return i

    def add_rcpe(self, R, Q, N, parallel=False):
        self.p += [R,Q,N]
        n = len(self.p)
        i = [n-3, n-2, n-1]
        tmp = self.eqc
        prefix = "|| " if parallel else "+ "
        CPE = lambda omega: 1/(self.p[i[1]]*np.power(1j*omega,self.p[i[2]]))
        RCPE = lambda omega: 1/(1/self.p[i[0]]+1/CPE(omega))
        if parallel:
            self.eqc = lambda omega: 1.0/(1.0/tmp(omega)+1.0/RCPE(omega))
        else:
            self.eqc = lambda omega: tmp(omega) + RCPE(omega)
        print "%sRCPE(%g,%g,%g) at [%i, %i, %i]" % \
                (prefix, R, Q, N, i[0], i[1], i[2])
        return i
    def add_rl(self, R, L, parallel=False):
        self.p += [R,L]
        n = len(self.p)
        i = [n-2,n-1]
        tmp = self.eqc
        prefix = "|| " if parallel else "+ "
        RL = lambda omega: self.p[i[0]]/ \
                           (1+self.p[i[0]]/(1j*omega*self.p[i[1]]))
        if parallel:
            self.eqc = lambda omega: 1.0/(1.0/tmp(omega)+1.0/RL(omega))
        else:
            self.eqc = lambda omega: tmp(omega) + RL(omega)
        print "%sRL(%g,%g) at [%i, %i]" % \
                (prefix, R, L, i[0], i[1])
        return i

    def add_warb(self, R, T, parallel=False,N=0.5):
        self.p += [R, T]
        n = len(self.p)
        i = [n-2,n-1]
        tmp = self.eqc
        prefix = "|| " if parallel else "+ "
        def W2(omega):
            tan = 0
            if self.p[i[1]]*omega < 1e4:
                tan = np.tanh(np.power(1j*self.p[i[1]]*omega,N))
            else:
                tan = 1
            return self.p[i[0]] * tan/np.power(1j*self.p[i[1]]*omega,N)
        W = lambda omega: W2(omega)
                
        if parallel:
            self.eqc = lambda omega: 1.0/(1.0/tmp(omega)+1.0/W(omega))
        else:
            self.eqc = lambda omega: tmp(omega) + W(omega)
        print "%sW(%g, %g) at [%i, %i]" % \
                (prefix, R, T, i[0], i[1])
        return i

    def add_zarc(self, R, T, N, parallel=False):
        self.p += [R, T, N]
        n = len(self.p)
        i = [n-3,n-2,n-1]
        tmp = self.eqc
        prefix = "|| " if parallel else "+ "
        ZARC = lambda omega: self.p[i[0]]/ \
                            (1+np.power(self.p[i[1]]*omega*1j,self.p[i[2]]))
        if parallel:
            self.eqc = lambda omega: 1.0/(1.0/tmp(omega)+1.0/ZARC(omega))
        else:
            self.eqc = lambda omega: tmp(omega) + ZARC(omega)
        print "%sZARC(%g,%g,%g) at [%i, %i, %i]" % \
                (prefix, R, T, N, i[0], i[1], i[2])
        return i


    
    def reset_eqc(self):
        self.eqc = lambda omega: 0
        self.p = []

    def fit(self):
        p0 = self.p 
        def residuals(p,s):
            s.p = np.abs(p)
            res = []
            for i,w in enumerate(s.omega):
                z = s.eqc(w)
                res.append(np.real(z)-s.Re[i])
                res.append(np.imag(z)-s.Im[i])
            return res
        plsq = leastsq(residuals, p0,args=(self,), \
                       full_output=True, xtol=1e-20, ftol=1e-12)
        self.p = np.abs(plsq[0])
        print self.p
        print plsq[3]


class EIS:
    def __init__(self, t, v, i, meta):
        self.t = np.array(t)
        self.v = np.array(v)
        self.i = -1*np.array(i)
        self.meta = meta
        self.frequency = meta['imp']['f']

    def complex(self, signal):
        # Discard the last entry of signal, to get periodic:
        # signal = signal[:-1]
        # Do Fourier transform
        nfft = signal.size
        n = np.floor(nfft/2)
        s_ft = np.fft.fft(signal)
        # Discard the first entry (omega = 0)
        S = np.abs(s_ft[1:n])
        S_ft = s_ft[1:]
        # Get complex value
        idx = np.argmax(S)
        z = S_ft[idx]/n
        return z

    def impedance(self):
        if not hasattr(self, 'Z'):
            z_v = self.complex(self.v)
            z_i = self.complex(self.i)
            self.Z = z_v/z_i
        return (self.frequency, self.Z)

    def __str__(self):
        f, Z = self.impedance()
        return "%g, %g, %g" % (f, np.real(Z), np.imag(Z))

    def plot_signal(self):
        fig, ax1 = plt.subplots()
        voltage,= ax1.plot(self.t, self.v,'r-')
        ax1.set_xlabel('time t [s]')
        ax1.set_ylabel('voltage V [-]')
        ax2 = ax1.twinx()
        ax2.set_ylabel('current density [-]')
        current,= ax2.plot(self.t, self.i,'b-')
        plt.legend([voltage,current],["V","i"])
        plt.title("f = %g [Hz]" % self.meta['imp']['f'])
        plt.show()

parser = argparse.ArgumentParser(description="EnPEn - EIS plot and analysis")
parser.add_argument('-s', '--single_plot', nargs=1, type=str, action='store',
                    help='Plot signal of a single simulation')
parser.add_argument('-S', '--spectrum', nargs='*', type=str, action='store',
                    help='Plot entire spectrum')
parser.add_argument('-m', '--movify', nargs=1, type=str, action='store')
parser.add_argument('-e', '--equivalent', action='store_true')
args = parser.parse_args()

if args.single_plot:
    s = Spectrum(args.single_plot)
    s.spectrum[0].plot_signal()
    sys.exit(0)
if args.spectrum:
    s = Spectrum(args.spectrum)
    edf = np.power(10.0, np.arange(-5,10,0.1))
    s.interpolate(edf)
    s.add_r(1e5, False)
    s.add_zarc(5e4,2e-6,0.99, False)
    #s.add_rc(2e4,1e-10, False)
    s.add_warb(1e4,1, False)
    s.add_rl(5e4,5e-5, False)
    s.add_rl(1e5,1e-5, False)
    s.fit()
    s.plot_impedance(eqc=args.equivalent)
    #print s.csv()
    s.plot_nyquist(data=True, eqc=args.equivalent)
    sys.exit(0)

if args.movify:
    fname = args.movify[0]
    ivc = IVC(fname)
    v, i = ivc.plot_ivc(False) 
    spectra = {}
    c = {}
    fig = plt.figure(figsize=(12.8,9.6))
    index = 0
    for bias in v:
        modulus = 0.001
        directory = "%s_%.2f-%.3f/" % (fname[:-5], bias, modulus)
        fnames = ['%s/%s' % (directory, f) for f in os.listdir(directory)]
        spectra[bias] = Spectrum(fnames)
        x, c[bias], names = ivc.plot_conc(bias, False, indices=[0,1])
        sys.stdout.write('.')
        sys.stdout.flush()

        ### Depletion zone plot ###
        ax1 = plt.subplot(221)
        ax1.set_xlim([0,105.0])
        ax1.set_ylim([0,0.0013])
        plt.title("Depletion boundary layer")
        plots = []
        for conc in c[bias]:
            p, = plt.plot(x, conc)
            plots.append(p)
        plt.legend(plots, names)
        plt.xlabel("x [micron]")
        plt.ylabel("c [M]")

        ### Enrichment zone plot ###
        ax2 = plt.subplot(223)
        ax2.set_xlim([195.0,300.0])
        ax2.set_ylim([0,0.003])
        plt.title("Enrichment boundary layer")
        plots = []
        for conc in c[bias]:
            p, = plt.plot(x, conc)
            plots.append(p)
        plt.legend(plots, names)
        plt.xlabel("x [micron]")
        plt.ylabel("c [M]")

        ### IVC
        ax3 = plt.subplot(222)
        plt.xlabel("V [-]")
        plt.ylabel("i [-]")
        plt.plot(v, i, 'k-')
        plt.plot(bias, i[v.index(bias)], 'ro')

        ### Nyquist plot
        ax4 = plt.subplot(224)
        plt.gca().invert_yaxis()
        plt.gca().set_aspect('equal')
        plt.title("Nyquist diagram")
        plt.xlabel("Real part")
        plt.ylabel("Imaginary part")
        #ax4.set_xlim([0,200000])
        #ax4.set_ylim([-200000,100000])
        re, im = spectra[bias].plot_nyquist(False)
        plt.plot(re, im)
        extrema = spectra[bias].get_extrema()
        for x in extrema:
            plt.annotate('%g'%x[0],
                         xy = (x[1],x[2]),
                         xycoords = 'data')
            plt.plot(x[1],x[2],'bo')

        plt.savefig("dump/conc-%05i.png" % index, dpi=125)
        plt.clf()
        index += 1
    print
    sys.exit(0)
