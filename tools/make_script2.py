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
    comment = False
    cmd = "%s " % exec_line
    for arg in args:
        cmd += "%s " % arg
    if kwargs is not None:
        for key, value in kwargs.iteritems():
            if key in ["out", "output", "output_file"]:
                cmd += "> %s" % value
            if key == 'comment' and value == True:
                comment=True
    cmd += "\n"

    if comment:
        cmd = '#' + cmd

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
parser.add_argument("--python", default="python", help="path/to/python")

args = parser.parse_args()
tools_path = os.path.realpath(__file__)
spectrassembler = tools_path.split('/tools')[0]

# In order not to recompute old things
do_comment = False

# Create main script
runall_script = "%s_runall.sh" % (args.prefix)
runall_fh = open(runall_script, 'wb')
runall_fh.write("#!/bin/bash\n")

# Create preprocess script file
pp_script = "%s_preprocess.sh" % (args.prefix)
pp_fh = open(pp_script, 'wb')

# Prelude
pp_fh.write("#!/bin/bash\n")
pp_fh.write("reads=\"%s\"\n" % args.reads)
pp_fh.write("if [ ! -d \"%s\" ]; then mkdir %s; fi\n" % (args.dir, args.dir))

# Check read ids are unique
pp_fh.write("# Check that read ids are unique\n")
pp_fh.write("%s %s/tools/check_reads.py %s\n" % (args.python, spectrassembler, args.reads))
readspp = '.'.join(args.reads.split('.')[:-1]) + '.pp.' + args.reads.split('.')[-1]
pp_fh.write("if [ -e \"%s\" ]; then reads=\"%s\"; fi\n" % (readspp, readspp))

# make fastq file for racon
if args.racon and args.minimap:
    # Racon + spectral
    pp_fh.write("# Check that reads are .fastq for racon\n")
    mycmd = "if [ \"${reads:(-1)}\" == \"q\" ]; then\
    \n  readsq=\"$reads\"\
    \nelse\
    \n  echo \"Reads are in fasta format. Convert to fastq before calling racon...\"\
    \n  readsq=\"${reads/.f*a/.fastq}\"\
    \n  %s %s/tools/fasta2fastq.py \"$reads\" \"$readsq\"\
    \nfi\n" % (args.python, spectrassembler)
    pp_fh.write(mycmd)

mycmd = \
"/usr/bin/time -v /bin/bash %s 2> %s.log 1> %s.out \n" % \
(pp_script, pp_script.split('.sh')[0], pp_script.split('.sh')[0])
runall_fh.write(mycmd)


# Create minimap script
minimap_script = "%s_minimap.sh" % (args.prefix)
minimap_fh = open(minimap_script, 'wb')
minimap_fh.write("#!/bin/bash\n")
minimap_fh.write("reads=\"%s\"\n" % args.reads)
minimap_fh.write("if [ -e \"%s\" ]; then reads=\"%s\"; fi\n" % (readspp, readspp))
minimap_fh.write("# Run minimap\n")
minimap_out = "%s/%s.mini.paf" % (args.dir, args.prefix)
mycmd = cmd_printer(args.minimap, args.minimap_opts, "$reads",
"$reads", output=minimap_out, comment=do_comment)
minimap_fh.write(mycmd)

mycmd = \
"/usr/bin/time -v /bin/bash %s 2> %s.log 1> %s.out \n" % \
(minimap_script, minimap_script.split('.sh')[0],
 minimap_script.split('.sh')[0])
runall_fh.write(mycmd)

# Create miniasm script
if args.miniasm:
    miniasm_script = "%s_miniasm.sh" % (args.prefix)
    miniasm_dir = "%s/miniasm" % (args.dir)
    miniasm_fh = open(miniasm_script, 'wb')
    miniasm_fh.write("#!/bin/bash\n")
    miniasm_fh.write("if [ ! -d \"%s\" ]; then mkdir %s; fi\n" % (miniasm_dir, miniasm_dir))
    miniasm_fh.write("reads=\"%s\"\n" % args.reads)
    miniasm_fh.write("if [ -e \"%s\" ]; then reads=\"%s\"; fi\n" % (readspp, readspp))
    miniasm_fh.write("# Run miniasm\n")
    miniasm_out = "%s/%s.miniasm.gfa" % (miniasm_dir, args.prefix)
    mycmd = cmd_printer(args.miniasm, args.miniasm_opts, "-f $reads",
    minimap_out, output=miniasm_out, comment=do_comment)
    miniasm_fh.write(mycmd)
    miniasm_fasta = "%s/%s.miniasm.fasta" % (miniasm_dir, args.prefix)
    mycmd = "awk '/^S/{print \">\"$2\"\\n\"$3}' %s | fold > %s\n" % (miniasm_out, miniasm_fasta)
    miniasm_fh.write(mycmd)

mycmd = \
"/usr/bin/time -v /bin/bash %s 2> %s.log 1> %s.out \n" % \
(miniasm_script, miniasm_script.split('.sh')[0], miniasm_script.split('.sh')[0])
runall_fh.write(mycmd)

# Run spectrassembler
spectral_script = "%s_spectral.sh" % (args.prefix)
spectral_dir = "%s/spectral" % (args.dir)
spectral_fh = open(spectral_script, 'wb')
spectral_fh.write("#!/bin/bash\n")
spectral_fh.write("if [ ! -d \"%s\" ]; then mkdir %s; fi\n" % (spectral_dir, spectral_dir))
spectral_fh.write("reads=\"%s\"\n" % args.reads)
spectral_fh.write("if [ -e \"%s\" ]; then reads=\"%s\"; fi\n" % (readspp, readspp))
spectral_fh.write("# Run spectral\n")
spectral_out = "%s/%s.spectral.fasta" % (spectral_dir, args.prefix)
mycmd = cmd_printer("%s %s/spectrassembler.py" % (args.python, spectrassembler), "-f $reads",
            "-m", minimap_out, "-r", spectral_dir,
            args.spectrassembler_opts, output=spectral_out, comment=do_comment)
spectral_fh.write(mycmd)

mycmd = \
"/usr/bin/time -v /bin/bash %s 2> %s.log 1> %s.out \n" % \
(spectral_script, spectral_script.split('.sh')[0], spectral_script.split('.sh')[0])
runall_fh.write(mycmd)

# Run Racon if provided
if args.racon and args.minimap:
    raconspectral_script = "%s_raconspectral.sh" % (args.prefix)
    raconspectral_dir = "%s/raconspectral" % (args.dir)
    raconspectral_fh = open(raconspectral_script, 'wb')
    raconspectral_fh.write("#!/bin/bash\n")
    raconspectral_fh.write("if [ ! -d \"%s\" ]; then mkdir %s; fi\n" % (raconspectral_dir, raconspectral_dir))
    raconspectral_fh.write("reads=\"%s\"\n" % args.reads)
    raconspectral_fh.write("# Check that reads are .fastq for racon\n")
    mycmd = "if [ \"${reads:(-1)}\" == \"q\" ]; then\
    \n  readsq=\"$reads\"\
    \nelse\
    \n  readsq=\"${reads/.f*a/.fastq}\"\
    \nfi\n"
    raconspectral_fh.write(mycmd)
    raconspectral_fh.write("# Run racon + spectral\n")
    # Map reads to miniasm assembly with minimap
    spectral_mappings_out = "%s/%s.spectral_mappings.paf" % (raconspectral_dir, args.prefix)
    mycmd = cmd_printer(args.minimap, args.minimap_opts, spectral_out, "$readsq",
    output=spectral_mappings_out, comment=do_comment)
    raconspectral_fh.write(mycmd)

    # Run Racon
    spectral_racon_out = "%s/%s.spectral_racon.fasta" % (raconspectral_dir, args.prefix)
    mycmd = cmd_printer(args.racon, "$readsq", spectral_mappings_out,
    spectral_out, spectral_racon_out, comment=do_comment)
    raconspectral_fh.write(mycmd)

    mycmd = \
    "/usr/bin/time -v /bin/bash %s 2> %s.log 1> %s.out \n" % \
    (raconspectral_script, raconspectral_script.split('.sh')[0],
    raconspectral_script.split('.sh')[0])
    runall_fh.write(mycmd)


    # Racon + miniasm
    if args.miniasm:
        raconminiasm_script = "%s_raconminiasm.sh" % (args.prefix)
        raconminiasm_dir = "%s/raconminiasm" % (args.dir)
        raconminiasm_fh = open(raconminiasm_script, 'wb')
        raconminiasm_fh.write("#!/bin/bash\n")
        raconminiasm_fh.write("if [ ! -d \"%s\" ]; then mkdir %s; fi\n" % (raconminiasm_dir, raconminiasm_dir))
        raconminiasm_fh.write("reads=\"%s\"\n" % args.reads)
        raconminiasm_fh.write("# Check that reads are .fastq for racon\n")
        mycmd = "if [ \"${reads:(-1)}\" == \"q\" ]; then\
        \n  readsq=\"$reads\"\
        \nelse\
        \n  readsq=\"${reads/.f*a/.fastq}\"\
        \nfi\n"
        raconminiasm_fh.write(mycmd)
        raconminiasm_fh.write("# Run racon + miniasm\n")
        # Map reads to miniasm assembly with minimap
        mappings_out = "%s/%s.miniasm_mappings.paf" % (raconminiasm_dir, args.prefix)
        mycmd = cmd_printer(args.minimap, args.minimap_opts, miniasm_fasta, "$readsq", output=mappings_out, comment=do_comment)
        raconminiasm_fh.write(mycmd)

        # Run Racon
        racon_out = "%s/%s.miniasm_racon.fasta" % (raconminiasm_dir, args.prefix)
        mycmd = cmd_printer(args.racon, "$readsq", mappings_out, miniasm_out, racon_out, comment=do_comment)
        raconminiasm_fh.write(mycmd)

        # Round 2 of Racon + miniasm
        raconminiasm_fh.write("# Run racon + miniasm 2nd iteration\n")
        # Map reads to miniasm assembly with minimap
        mappings2_out = "%s/%s.racon2_mappings.paf" % (raconminiasm_dir, args.prefix)
        mycmd = cmd_printer(args.minimap, args.minimap_opts, racon_out, "$readsq", output=mappings2_out)
        raconminiasm_fh.write(mycmd)

        # Run Racon
        racon2_out = "%s/%s.miniasm_racon2.fasta" % (raconminiasm_dir, args.prefix)
        mycmd = cmd_printer(args.racon, "$readsq", mappings2_out, racon_out, racon2_out)
        raconminiasm_fh.write(mycmd)

        mycmd = \
        "/usr/bin/time -v /bin/bash %s 2> %s.log 1> %s.out \n" % \
        (raconminiasm_script, raconminiasm_script.split('.sh')[0],
        raconminiasm_script.split('.sh')[0])
        runall_fh.write(mycmd)



# Run Canu if provided
if args.canu:
    canu_script = "%s_canu.sh" % (args.prefix)
    canu_dir = "%s/canu" % (args.dir)
    canu_fh = open(canu_script, 'wb')
    canu_fh.write("#!/bin/bash\n")
    canu_fh.write("if [ ! -d \"%s\" ]; then mkdir %s; fi\n" % (canu_dir, canu_dir))
    canu_fh.write("reads=\"%s\"\n" % args.reads)
    canu_fh.write("if [ -e \"%s\" ]; then reads=\"%s\"; fi\n" % (readspp, readspp))
    canu_fh.write("# Run canu\n")
    mycmd = cmd_printer(args.canu, "-p", args.prefix, "-d", canu_dir,
    args.canu_opts, comment=do_comment)
    canu_fh.write(mycmd)

    mycmd = \
    "/usr/bin/time -v /bin/bash %s 2> %s.log 1> %s.out \n" % \
    (canu_script, canu_script.split('.sh')[0],
    canu_script.split('.sh')[0])
    runall_fh.write(mycmd)
