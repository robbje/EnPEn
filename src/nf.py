#!/usr/bin/env python2

import sys
import subprocess
import os
import threading
import Queue
import numpy as np
import time

eis_bin = './flow'
MAX_THREADS = 6
TOLERANCES = {'-snes_atol':'1e-8', '-snes_stol':'1e-30', '-snes_stol':'1e-16'}

def run_simul(args):
    filename = 'nf/ITA_NAOH_u{}_naoh{}_c{}'.format(args['-velocity'], args['-naoh_conc'], args['-charge'])
    args['-name'] = filename
    args.update(TOLERANCES)
    argv = [eis_bin]
    for arg in [[k, str(args[k])] for k in args]:
        argv += arg
    subprocess.check_call(argv)

def worker(q):
    while True:
        args = q.get()
        run_simul(args)
        q.task_done()

if __name__ == '__main__':
    for charge in np.arange(0.001, 0.050, 0.001):
        print "\n[+] Charge: %g\n" % charge
        for velocity in np.arange(0.001, 0.05, 0.001):
            q = Queue.Queue()
            for conc in np.arange(0.0001, 0.003, 0.0001):
                q.put({'-velocity':velocity,
                       '-naoh_conc':conc,
                       '-charge':charge})
            for i in xrange(MAX_THREADS):
                t = threading.Thread(target=worker, args=(q,))
                t.daemon = True
                t.start()
            q.join()
