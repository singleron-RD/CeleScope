#!/bin/env python
#coding=utf8

import argparse

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='fastq to matrix')
    subparsers = parser.add_subparsers()

    from barcode import barcode, get_opts1
    parser1 = subparsers.add_parser('barcode', description='extract barcode and umi')
    get_opts1(parser1,True)
    parser1.set_defaults(func=barcode)

    from cutadapt import cutadapt, get_opts2
    parser2 = subparsers.add_parser('cutadapt', description='cutadapt')
    get_opts2(parser2,True)
    parser2.set_defaults(func=cutadapt)

    from STAR import STAR, get_opts3
    parser3 = subparsers.add_parser('STAR')
    get_opts3(parser3,True)
    parser3.set_defaults(func=STAR)

    from featureCounts2 import featureCounts, get_opts4
    parser4 = subparsers.add_parser('featureCounts2')
    get_opts4(parser4,True)
    parser4.set_defaults(func=featureCounts)

    from count import count, get_opts5
    parser5 = subparsers.add_parser('count')
    get_opts5(parser5,True)
    parser5.set_defaults(func=count)

    from run import run, get_opts6
    parser6 = subparsers.add_parser('run',conflict_handler='resolve')
    get_opts1(parser6,False)
    get_opts2(parser6,False)
    get_opts3(parser6,False)
    get_opts4(parser6,False)
    get_opts5(parser6,False)
    get_opts6(parser6,False)
    parser6.set_defaults(func=run)

    args = parser.parse_args()
    args.func(args)
