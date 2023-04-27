import pysam
import pandas as pd
import re
import math
import gzip

def rev_comp(seq):
    """ Returns the reverse complement of a DNA sequence,
        retaining the case of each letter"""
    complement = ""

    for base in seq:
        if base == "A": complement += "T"
        elif base == "T": complement += "A"
        elif base == "G": complement += "C"
        elif base == "C": complement += "G"
        elif base == "N": complement += "N"
        elif base == "a": complement += "t"
        elif base == "t": complement += "a"
        elif base == "g": complement += "c"
        elif base == "c": complement += "g"
        elif base == "n": complement += "n"
        elif base == "*": complement += "*"
        else:
            complement += base
            print("Warning: reverse complement function encountered unknown base " + "'" + base + "'")

    reverseComplement = complement[::-1]

    return reverseComplement

def flip_fastq(infile, outfile):

    reverse = ''

    if infile.endswith('.gz'):
        f = gzip.open(infile, 'rt')
    else:
        f = open(infile, 'r')
    o = open(outfile, 'w')

    for line in f:
        line = line.strip()
        # if line starts with @, write identifier to ouput file
        if re.match('^@.{19,49}', line):
            o.write(line + "\n")
        elif re.match('^[ATCG].*', line):
            seq = line
            # if there is a TTTTTT repeat at begining of read, reverse sequence then write to output file
            if re.match('TTTTTTTTTTTTTTTTTT*', seq):
                line = rev_comp(seq)
                o.write(line + '\n')
                reverse = True
            # if there is no TTTTT repeat at beginng of read, write to output file
            else:
                o.write(line + '\n')
                reverse = False
        else:
            if not re.match('^@ATCG([+]\n)', line):
                qual = line
                # if the sequence line was reversed, reverse the quality score line then write to ouput
                if reverse == True:
                    line = qual[::-1]
                    o.write(line + '\n')
                # if the sequence line was not reversed, write quality score line to output
                else:
                    o.write(line + '\n')
    o.close()
    f.close()

def reverse_alignment(infile, outfile, threads=1):

    reverse_strand = {0: 16, 16: 0}

    if infile.endswith('.bam'):
        in_mode = 'rb'
    else:
        in_mode = 'r'
    if outfile.endswith('bam'):
        out_mode = 'wb'
    else:
        out_mode = 'w'
    input =  pysam.AlignmentFile(infile, in_mode, threads=threads)
    output = pysam.AlignmentFile(outfile, out_mode, template=input, threads=threads)

    for read in input:
        if read.has_tag('ts') and read.flag in reverse_strand:
            if read.get_tag('ts') == '-':
                read.flag = reverse_strand[read.flag]
                read.set_tag('ts', '+')

        output.write(read)
    input.close()
    output.close()
