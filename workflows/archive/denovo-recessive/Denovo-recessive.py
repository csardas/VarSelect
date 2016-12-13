#!/usr/bin/python

import os, sys, getopt

def main(argv):
    usage = "usage: python Denovo-recessive.py -v [vcf file] -p [ped file]"

    try:
        opts, args = getopt.getopt(argv,"hv:p:",["vcf=","ped="])
    except getopt.GetoptError:
        print (usage)
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            print (usage)
            sys.exit()
        elif opt in ("-v", "--vcf"):
            vcffile = arg
            outputfile = os.path.basename(os.path.splitext(arg)[0]) + '_Denovo-recessive.txt'
        elif opt in ("-p", "--ped"):
            pedfile = arg

    if len(opts) != 2:
        print ("try 'python Denovo-recessive.py -h' for help")
        sys.exit(2)


    f = open(pedfile, 'r')
    for line in f:
        if line.startswith('#'):
            continue

        tmp = line.strip().split('\t')
        if tmp[2] != '0' and tmp[3] != '0' and tmp[5] == '2':
            c_ = tmp[1]
            f_ = tmp[2]
            m_ = tmp[3]
            # only take first affected child
            break
    f.close()

    f2 = open(vcffile, 'r')

    tmp = []
    for line in f2:
        if '#CHROM' in line:
            cfm = line.strip().split('\t')
            c_index = cfm.index(c_)
            f_index = cfm.index(f_)
            m_index = cfm.index(m_)
        elif line.startswith('#'):
            pass
        else:
            tmp.append(line)
    f2.close()

    vlist = [x.strip().split('\t') for x in tmp]


    f3 = open(outputfile, 'w')
    f3.write('#chrom\tpos\tref\talt\tis_DR\t'+ c_ + '\t' + f_ + '\t' + m_ + '\n')

    HOM_ALT = ['1/1', '2/2', '3/3', '1/2', '2/1', '1/3', '3/1', '2/3', '3/2']
    for i in range(len(vlist)):
        tmp = '\t'.join(vlist[i])
        tmp = tmp.replace('.', '0').split('\t')

        child = tmp[c_index].split(':')[0]
        father = tmp[f_index].split(':')[0]
        mother = tmp[m_index].split(':')[0]

        if child == '0':
            child = '0/0'
        if father == '0':
            father = '0/0'
        if mother == '0':
            mother = '0/0'

        true = vlist[i][0] + '\t' + vlist[i][1] + '\t' + vlist[i][3] + '\t' + vlist[i][4] + '\t' + '1' + '\t' + vlist[i][c_index].split(':')[0] + '\t' + vlist[i][f_index].split(':')[0] + '\t' + vlist[i][m_index].split(':')[0] + '\n'
        false = vlist[i][0] + '\t' + vlist[i][1] + '\t' + vlist[i][3] + '\t' + vlist[i][4] + '\t' + '0' + '\t' + vlist[i][c_index].split(':')[0] + '\t' + vlist[i][f_index].split(':')[0] + '\t' + vlist[i][m_index].split(':')[0] + '\n'

        # pass if chrom = MT or in Non PAR region
        if vlist[i][0] == 'MT' or (vlist[i][0] == 'X' and ((int(vlist[i][1]) < 60001) or (2699520 <= int(vlist[i][1]) <= 154931044) or (int(vlist[i][1]) > 155270560))):
            f3.write(false)
            continue

        if (child == '0/0' and father in HOM_ALT and mother in HOM_ALT) or (child in HOM_ALT and father == '0/0' and mother == '0/0'):
            f3.write(true)
        else:
            f3.write(false)

    f3.close()

if __name__ == "__main__":
    main(sys.argv[1:])
