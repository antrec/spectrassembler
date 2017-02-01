#!/usr/bin/env python

import os
import argparse


def cmd_printer(exec_line, *args, **kwargs):
    """
    Creates string to be printed to bash in order to run executable.

    Parameters
    ----------
    exec_line : str (execution line including path/to/exec)
    extra args : str (option to be passed, e.g. "-v 1")
    extra kwargs : may include "output" kwarg if redirection to "output" file
    is demanded.

    Returns
    ----------
    cmd : str (block of text to be printed to a bash script to run exec_line
    [with options and with output redirected if specified in *args and *kwargs]
    when the bash script is ran.)
    """
    cmd = "time %s " % exec_line
    for arg in args:
        cmd += "%s " % arg
    if kwargs is not None:
        for key, value in kwargs.iteritems():
            if key in ["out", "output", "output_file"]:
                cmd += "> %s" % value
    cmd += "\n"

    return cmd

parser = argparse.ArgumentParser(description="creates bash script"
                                "to test assemblers on a dataset")
parser.add_argument("reads", help="path to reads file (fasta or fastq)")
parser.add_argument("-p", "--prefix", required=True, help="prefix for output")
parser.add_argument("-d", "--dir",required=True, help="directory for output")
parser.add_argument("--minimap", required=True,
                    help="path/to/minimap (executable)."
                    "Required for spectrassembler.")
parser.add_argument("--minimap_opts", default="-Sw5 -L100 -m0",
                    help="options for minimap")
parser.add_argument("--miniasm", help="path/to/miniasm (executable)")
parser.add_argument("--miniasm_opts", default="",
                    help="options for miniasm")
parser.add_argument("--racon", help="path/to/racon (executable)")
parser.add_argument("--canu", help="path/to/canu (executable)")
parser.add_argument("--canu_opts", default="",
                    help="options for canu. Read file needs to be in it through"
                    "e.g. -nanopore-raw reads.fasta")
parser.add_argument("--spectrassembler_opts", default="",
                    help="options for spectrassembler")
parser.add_argument("--spoa", help="path/to/spoa/")

args = parser.parse_args()
tools_path = os.path.realpath(__file__)
spectrassembler = tools_path.split('/tools')[0]

# Open script file
script_fn = "%s_script.sh" % (args.prefix)
print "Create bash script %s" %	(script_fn)
fh = open(script_fn, 'wb')

# Prelude
fh.write("#!/bin/bash\n")
fh.write("reads=\"%s\"\n" % args.reads)
fh.write("if [ ! -d \"%s\" ]; then mkdir %s; fi\n" % (args.dir, args.dir))

# Check read ids are unique
fh.write("# Check that read ids are unique\n")
fh.write("python %s/tools/check_reads.py %s\n" % (spectrassembler, args.reads))
readspp = '.'.join(args.reads.split('.')[:-1]) + '.pp.' + args.reads.split('.')[-1]
fh.write("if [ -e \"%s\" ]; then reads=\"%s\"; fi\n" % (readspp, readspp))

# Run minimap
fh.write("# Run minimap\n")
minimap_out = "%s/%s.mini.paf" % (args.dir, args.prefix)
mycmd = cmd_printer(args.minimap, args.minimap_opts, "$reads",
"$reads", output=minimap_out)
fh.write(mycmd)

# Run miniasm if provided
if args.miniasm:
    fh.write("# Run miniasm\n")
    miniasm_out = "%s/%s.miniasm.gfa" % (args.dir, args.prefix)
    mycmd = cmd_printer(args.miniasm, args.miniasm_opts, "-f", args.reads,
    minimap_out, output=miniasm_out)
    fh.write(mycmd)
    miniasm_fasta = "%s/%s.miniasm.fasta" % (args.dir, args.prefix)
    mycmd = "awk '/^S/{print \">\"$2\"\\n\"$3}' %s | fold > %s\n" % (miniasm_out, miniasm_fasta)
    fh.write(mycmd)

# Run spectrassembler
fh.write("# Run spectrassembler\n")
spectral_out = "%s/%s.spectral.fasta" % (args.dir, args.prefix)
spectral_dir = "%s/spectral" % (args.dir)
mycmd = cmd_printer("python %s/spectrassembler.py" % (spectrassembler), "-f $reads",
            "-m", minimap_out, "-r", args.dir, "--spoapath", args.spoa,
            args.spectrassembler_opts, "-r", spectral_dir, output=spectral_out)
fh.write(mycmd)

# Run Racon if provided
if args.racon and args.minimap:
    # Racon + spectral
    fh.write("# Check that reads are .fastq for racon\n")
    mycmd = "if [ \"${reads:(-1)}\" == \"q\" ]; then\
    \n  readsq=\"$reads\"\
    \nelse\
    \n  echo \"Reads are in fasta format. Convert to fastq before calling racon...\"\
    \n  readsq=\"${reads/.f*a/.fastq}\"\
    \n  python %s/tools/fasta2fastq.py \"$reads\" \"$readsq\"\
    \nfi\n" % (spectrassembler)
    fh.write(mycmd)

    fh.write("# Run racon + spectral\n")
    # Map reads to miniasm assembly with minimap
    spectral_mappings_out = "%s/%s.spectral_mappings.paf" % (args.dir, args.prefix)
    mycmd = cmd_printer(args.minimap, spectral_out, "$readsq",
    output=spectral_mappings_out)
    fh.write(mycmd)

    # Run Racon
    spectral_racon_out = "%s/%s.spectral_racon.fasta" % (args.dir, args.prefix)
    mycmd = cmd_printer(args.racon, "$readsq", spectral_mappings_out,
    spectral_out, spectral_racon_out)
    fh.write(mycmd)

    # Racon + miniasm
    if args.miniasm:
        fh.write("# Run racon + miniasm\n")
        # Map reads to miniasm assembly with minimap
        mappings_out = "%s/%s.miniasm_mappings.paf" % (args.dir, args.prefix)
        mycmd = cmd_printer(args.minimap, miniasm_fasta, "$readsq", output=mappings_out)
        fh.write(mycmd)

        # Run Racon
        racon_out = "%s/%s.miniasm_racon.fasta" % (args.dir, args.prefix)
        mycmd = cmd_printer(args.racon, "$readsq", mappings_out, miniasm_out, racon_out)
        fh.write(mycmd)

# Run Canu if provided
if args.canu:
    fh.write("# Run Canu\n")
    canu_dir = "%s/canu" % (args.dir)
    mycmd = cmd_printer(args.canu, "-p", args.prefix, "-d", canu_dir,
    args.canu_opts)
    fh.write(mycmd)

fh.close()
