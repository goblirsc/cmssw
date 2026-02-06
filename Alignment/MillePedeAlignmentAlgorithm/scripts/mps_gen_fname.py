#!/usr/bin/env python3

from Alignment.MillePedeAlignmentAlgorithm.mpslib.mps_formats import mps_formats 
from argparse import ArgumentParser

if __name__ == "__main__":
    parser = ArgumentParser("mps_gen_fname")
    parser.add_argument("-m", dest="milleFormat", help="print the CMSSW mille output file name", action="store_true") 
    parser.add_argument("-o", dest="milleOutFile", help="print the overall mille output file name", action="store_true") 
    parser.add_argument("-p", dest="pedeInFile", help="print the pede input file name", action="store_true") 
    parser.add_argument("-c", dest="copyCommand", help="print the copy command", action="store_true") 
    parser.add_argument("-g", dest="compressCommand", help="print the compression command", action="store_true") 
    parser.add_argument("-r", dest="cleanupCommand", help="print the cleanup command", action="store_true") 
    parser.add_argument("--eos", dest="eosInstance", help = "EOS instance to use", default="eoscms.cern.ch", type=str)
    parser.add_argument("fname",type=str, help="file name stub without extension")
    parser.add_argument("directory",type=str, help="target directory for file")
    parser.add_argument("format",type=str, help="file extension identifier string")


    args = parser.parse_args() 

    formats = mps_formats(args.eosInstance)
    theFormat = formats.identifyFormat(args.format)
    if args.milleFormat:
        print( formats.milleJobTarget(args.fname, theFormat))
    if args.milleOutFile:
        print( formats.milleJobTarget(args.fname, theFormat)) 
    if args.pedeInFile:
        print( formats.pathForPede(args.directory, args.fname, theFormat)) 
    if args.copyCommand:
        print( formats.copyCommand(args.directory, args.fname, theFormat)) 
    if args.compressCommand:
        print( formats.compressCommand(args.fname, theFormat)) 
    if args.cleanupCommand:
        print( formats.cleanupCommand(args.fname, theFormat)) 
