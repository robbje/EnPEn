#!/usr/bin/env python2

import sys
import subprocess
import os
import threading
import Queue
import numpy as np
import time

eis_bin = './eis'
MAX_THREADS = 6

def run_eis(frequency, charge=0.1, bias=0.0, modulus=0.001, name='impedance'):
    args = {}
    args['-eis-bias'] = bias
    args['-eis-modulus'] = modulus
    args['-eis-frequency'] = frequency
    args['-charge'] = charge
    if charge == 0.0:
        dirname = 'eis/%s_unmod_%.2f-%g' % (name, bias, modulus)
    else:
        dirname = 'eis/%s_mod%.2f_%.2f-%g' % (name, charge, bias, modulus)
    try:
        os.mkdir('outdata/%s' % dirname)
    except:
        pass
    filename = '%s/%gHz' % (dirname, frequency)
    args['-eis-name'] = filename
    args['-snes_atol'] = '1e-8'
    args['-snes_rtol'] = '1e-30'
    args['-snes_stol'] = '1e-16'
    argv = [eis_bin]
    for arg in [[k, str(args[k])] for k in args]:
        argv += arg
    try:
        subprocess.check_call(argv)
    except Exception as e:
        print "[!] Problem occurred simulating frequency: %g - %s" % (frequency, e)
        os.unlink('outdata/%s.json' % filename)
        return

def worker(q):
    while True:
        args = q.get()
        run_eis(**args)
        q.task_done()

if __name__ == '__main__':
    frequencies = np.power(10.0, np.arange(-6,10,0.1))
    for charge in [0.1, 0.5]:
        print "\n[+] Charge: %g\n" % charge
        for bias in np.arange(0.0, 0.9, 1.0):
            q = Queue.Queue()
            for f in frequencies:
                q.put({'frequency':f,'bias':bias,'modulus':0.001,'charge':charge,
                       'name':'asymetric/review/trimix_3'})
            for i in xrange(MAX_THREADS):
                t = threading.Thread(target=worker, args=(q,))
                t.daemon = True
                t.start()
            q.join()
