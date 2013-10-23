#!/usr/bin/env python

import sys
import os
from subprocess import Popen, PIPE

if len(sys.argv) != 5:
    print "Usage: %s ref_genome coord_start coord_end chromosome" % sys.argv[0]
    sys.exit(1)

ref = sys.argv[1]
coordA = int(sys.argv[2])
coordB = int(sys.argv[3])
chromo = sys.argv[4]

def get_name(line):
    ids = line[8].split(";")
    for id in ids:
        if id.startswith("Name="):
            name = id.split("=", 1)[-1]
            sys.stderr.write("Getting name for gene %s\n" % name)
            sys.stderr.flush()
            p = Popen(["perl", "src/id2name.pl", name], stdout=PIPE)
            stdout, stderr = p.communicate()
            if stdout:
                print "GeneName %s" % stdout.rstrip()
            else:
                print "GeneID %s" % name
       
            return

    sys.stderr.write("Gene found but coordinates don't match\n")
    sys.exit(1)

with open(ref) as fh:
    for line in fh:
        line = line.rstrip().split("\t")
        if len(line) != 9:
            continue

        if line[2] != "gene":
            continue

        else:
            if line[0] != chromo:
                continue

            ref_coordA = int(line[3])
            ref_coordB = int(line[4])

            if (coordA >= ref_coordA and coordA <= ref_coordB) or (
               coordB >= ref_coordA and coordB <= ref_coordB) or (
               coordA >= ref_coordA and coordB <= ref_coordB):
                get_name(line)

            if coordB < ref_coordA:
                break
