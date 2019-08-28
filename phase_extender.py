import os 
import sys
import time
import argparse

from phaser import phase_converter
from arg_builders import get_args
from val_extractor import args_to_val

def main():
    parser = argparse.ArgumentParser()
    args_namespace = get_args(parser)
    print(args_namespace)
    parsed_args = args_to_val(args_namespace)
    # soi, outputdir, nt, input_file, lods_cut_off, snp_threshold, num_of_hets, max_is, maxed_as, use_bed, bed_file, use_refhap, refhap, use_sample ,hapstats, writelod, addmissingsites

    phase_converter(*parsed_args)


if __name__ == "__main__":
    main()