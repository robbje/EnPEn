#!/usr/bin/env python2

import sys
from plex import *
from plex.traditional import re
import re as rex

group_z = {"H+":1, "COOH":0,"COO-":-1, "OH-":-1, "H2O":0, "a-CH2":0, "HC=CH2":0, "CH2":0, "NH2":0,"NH3+":1}

class Reaction():
    def __init__(self, kf, num_species):
        self.kf = kf
        self.nu = [0]*num_species
        self.kb = 0
    def set_kb(self, kb):
        self.kb = kb
    def __str__(self):
        if self.kb == 0:
            print "Error: No backward reaction constant defined"
        s = ""
        for coeff in self.nu:
            s += "%i " % coeff
        s += "%g %g" % (self.kf, self.kb)
        return s

class Species():
    def __init__(self, identifier):
        self.identifier = identifier
        self.charge = 0
        self.groups = []
        self.val_next = self.set_mw
    def set_mw(self, mw):
        self.mw = mw
        self.val_next = self.set_D
    def set_D(self, D):
        self.D = D
        self.val_next = self.set_N
    def set_N(self, N):
        self.N = N
    def add_group(self, num, group):
        global group_z
        for i in xrange(num):
            self.groups.append(group)
        self.charge += num * group_z[group]
    def __str__(self):
        s = "%s %i %g %g %g"%(self.identifier,self.charge,self.mw,self.D,self.N)
        return s

class MyParser():
    species = []
    reactions = []
    cur_species = None
    rere = rex.compile("([-+]?[0-9]*)([a-zA-Z][-+=a-zA-Z0-9]*)")
    def new_species(self, scanner, text):
        self.species.append(Species(text))
    def add_species_info(self, scanner, text):
        self.species[-1].val_next(float(text))
    def new_group(self, scanner, text):
        for s in self.species:
            if s.identifier == text:
                self.cur_species = s
                return
        print "Species not found: %s" % text
    def add_gdef(self, scanner, text):
        m = self.rere.findall(text)[0]
        self.cur_species.add_group(int(m[0]), m[1])
    def new_reac(self, scanner, text):
        kf = float(text[2:-1])
        self.reactions.append(Reaction(kf, len(self.species)))
    def add_back(self, scanner, text):
        kb = float(text[2:-1])
        if text[0] == 'k':
            kb = self.reactions[-1].kf/(10**-kb)
        self.reactions[-1].set_kb(kb)
    def add_reac(self, scanner, text):
        m = self.rere.findall(text)[0]
        for index, s in enumerate(self.species):
            if s.identifier == m[1]:
                break
        coeff = int(m[0])
        self.reactions[-1].nu[index] = coeff

class MyScanner(Scanner):
    coeff = re("[-+]?[0-9]*")
    fp = re("[-+]?[0-9]+(\.[0-9]*)?([eE]?[-+]?[0-9]+)?")
    specname = re("[a-zA-Z][-+a-zA-Z0-9=]*")
    funcgroup = re("[0-9]+")+specname
    reactand = coeff+specname
    comment = re("#[^\n]*\n")
    forw = Str("f=")+fp+Str(";")
    back = Str("b=")+fp+Str(";")
    logk = Str("k=")+fp+Str(";")
    w = Any(" \t\n")
    endl = Str("\n")
    p = MyParser()

    id_spec = re("S")
    id_group = re("G")
    id_reac = re("REAC")
    lexicon = Lexicon([
        (Bol+id_spec, Begin('species')),
        (Bol+id_group, Begin('group')),
        (Bol+id_reac, Begin('reaction')),
        (comment, IGNORE),
        (w, IGNORE),
        (State('species', [
            (comment, Begin('')),
            (endl, Begin('')),
            (specname, p.new_species),
            (w, IGNORE),
            (fp, p.add_species_info)
            ])),

        (State('group', [
            (comment, Begin('')),
            (endl, Begin('')),
            (w, IGNORE),
            (specname, p.new_group),
            (funcgroup, p.add_gdef),
            ])),

        (State('reaction', [
            (comment, Begin('')),
            (endl, Begin('')),
            (w, IGNORE),
            (forw, p.new_reac),
            (back, p.add_back),
            (logk, p.add_back),
            (reactand, p.add_reac),
            ]))
        ])

    def __init__(self, f, name):
        Scanner.__init__(self, self.lexicon, f, name)
    def __str__(self):
        global group_z
        res = "%i %i %i\n" % \
                (len(self.p.species), len(self.p.reactions), len(group_z))
        for s in self.p.species:
            res += "%s\n" % s
        for r in self.p.reactions:
            res += "%s\n" % r
        for z in group_z:
            res += "%i\n" % group_z[z]
        for s in self.p.species:
            group = ""
            for g in group_z:
                group += "%i " % s.groups.count(g)
            res += "%s\n" % group
        res += "298.15\n"  # Temperature
        res += "0.001\n"   # Viscosity
        res += "997\n"     # Density
        res += "14.9\n"    # PDH rho
        res += "96485.3\n" # Faraday
        res += "8.31446\n" # Gas constant
        res += "6.90612e-10" # epsilon_r epsilon_0
        return res

test = MyScanner(open(sys.argv[1], "r"), "test")
test.read()
print test
