#!/usr/bin/python

import os, sys, getopt

def main(argv):
    usage = "usage: python Autosomal-recessive.py -v [vcf file] -p [ped file]"

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
            outputfile = os.path.basename(os.path.splitext(arg)[0]) + '_Autosomal-recessive.txt'
        elif opt in ("-p", "--ped"):
            pedfile = arg

    if len(opts) != 2:
        print ("try 'python Autosomal-recessive.py -h' for help")
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
            break # only take one affected child
    f.close()

    # make vlist from VCF
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


    # Filtering
    f3 = open(outputfile, 'w')
    f3.write('#chrom\tpos\tref\talt\tis_AR\t'+ c_ + '\t' + f_ + '\t' + m_ + '\n')
    
    HET = ['0/1', '1/0', '0/2', '2/0', '0/3', '3/0']
    for i in range(len(vlist)):
        tmp = '\t'.join(vlist[i])
        # consider '.' as '0'
        tmp = tmp.replace('.', '0').split('\t')

        child = tmp[c_index].split(':')[0]
        father = tmp[f_index].split(':')[0]
        mother = tmp[m_index].split(':')[0]

        # if genotype in VCF is only one '.' or '0', consider as '0/0'
        if child == '0':
            child = '0/0'
            tmp[c_index] = '0/0'
        if father == '0':
            father = '0/0'
            tmp[f_index] = '0/0'
        if mother == '0':
            mother ='0/0'
            tmp[m_index] = '0/0'

        c1 = tmp[c_index][0]
        c2 = tmp[c_index][2]

        true = vlist[i][0] + '\t' + vlist[i][1] + '\t' + vlist[i][3] + '\t' + vlist[i][4] + '\t' + '1' + '\t' + vlist[i][c_index].split(':')[0] + '\t' + vlist[i][f_index].split(':')[0] + '\t' + vlist[i][m_index].split(':')[0] + '\n'
        false = vlist[i][0] + '\t' + vlist[i][1] + '\t' + vlist[i][3] + '\t' + vlist[i][4] + '\t' + '0' + '\t' + vlist[i][c_index].split(':')[0] + '\t' + vlist[i][f_index].split(':')[0] + '\t' + vlist[i][m_index].split(':')[0] + '\n'


        # pass if chrom = MT or in Non PAR region
        if vlist[i][0] == 'MT' or (vlist[i][0] == 'X' and ((int(vlist[i][1]) < 60001) or (2699520 <= int(vlist[i][1]) <= 154931044) or (int(vlist[i][1]) > 155270560))):
            f3.write(false)
            continue

        if c1 == c2 != '0' and c1 in father and c1 in mother and (father in HET) and (mother in HET):
            f3.write(true)
        else:
            f3.write(false)

    f3.close()

if __name__ == "__main__":
    main(sys.argv[1:])
