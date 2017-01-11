#!/usr/bin/python

"""
    gender: True -> Male
            False -> Female

    is_F: Father is affected
    is_M: Mother is affected
"""

import os, sys, getopt

def main(argv):
    usage = "usage: python X-linked.py -v [vcf file] -p [ped file]"

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
        elif opt in ("-p", "--ped"):
            pedfile = arg

    if len(opts) != 2:
        print ("try 'python X-linked.py -h' for help")
        sys.exit(2)

    # get all parents from PED
    f = open(pedfile, 'r')
    parents = []
    for line in f:
        if line.startswith('#'):
            continue

        tmp = line.strip().split('\t')
        if tmp[2] != '0':
            if tmp[2] not in parents:
                parents.append(tmp[2])
        if tmp[3] != '0':
            if tmp[3] not in parents:
                parents.append(tmp[3])
    # Re-read file
    f.seek(0)

    is_female = False
    # get first affected sample
    for line in f:
        if line.startswith('#'):
            continue

        tmp = line.strip().split('\t')
        # affected must be Male
        if tmp[1] not in parents and tmp[5] == '2':
            if tmp[4] == '1':
                c_ = tmp[1]
                f_ = tmp[2]
                m_ = tmp[3]
                # only take first affected child
                break
            else:
                c_ = tmp[1]
                f_ = tmp[2]
                m_ = tmp[3]
                is_female = True
##        if tmp[1] not in parents and tmp[4] == '1' and tmp[5] == '2':
#            c_ = tmp[1]
#            f_ = tmp[2]
#            m_ = tmp[3]
#            # only take first affected child
#            break

#    # if no Male affected sample found, exit program
#    try:
#        c_
#    except NameError:
#        print ("[Notice] None affected sample found, or it's a Female.\n[Notice] Use other genetic model instead.")
#        sys.exit()

    # Re-read file
    f.seek(0)
    # check if father or mother affected
    is_F = False
    is_M = False
    for line in f:
        if line.startswith('#'):
            continue

        tmp = line.strip().split('\t')
        if tmp[1] == f_:
            if tmp[5] == '2':
                is_F = True
        elif tmp[1] == m_:
            if tmp[5] == '2':
                is_M = True
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
    print ('#chrom\tpos\tref\talt\tis_XL\t'+ c_ + '\t' + f_ + '\t' + m_)

    for i in range(len(vlist)):
        tmp = '\t'.join(vlist[i])
        tmp = tmp.replace('.', '0').split('\t')

        if '/' not in tmp[c_index]:
            child = tmp[c_index][0]
            father = tmp[f_index][0]
            mother = tmp[m_index].split(':')[0]

            if mother == '0':
                mother = '0/0'
                tmp[m_index] == '0/0'

            m1 = tmp[m_index][0]
            m2 = tmp[m_index][2]

            true = vlist[i][0] + '\t' + vlist[i][1] + '\t' + vlist[i][3] + '\t' + vlist[i][4] + '\t' + '1' + '\t' + vlist[i][c_index][0] + '\t' + vlist[i][f_index][0] + '\t' + vlist[i][m_index].split(':')[0]
            false = vlist[i][0] + '\t' + vlist[i][1] + '\t' + vlist[i][3] + '\t' + vlist[i][4] + '\t' + '0' + '\t' + vlist[i][c_index][0] + '\t' + vlist[i][f_index][0] + '\t' + vlist[i][m_index].split(':')[0]

            # if affected sample is female, all variants are false
            if is_female:
                print (false)

            # pass if chrom = MT or Y
            if vlist[i][0] == 'MT' or vlist[i][0] == 'Y':
                print (false)
                continue

            # 1. both parents are unaffected
            if not(is_F) and not(is_M):
                if child in mother and child != '0' and father == '0' and (m1 == '0' or m2 == '0'):
                    print (true)
                else:
                    print (false)

            # 2. Father is affected, mother is unaffected
            elif is_F and not(is_M):
                if child in mother and child != '0' and father != '0' and (m1 == '0' or m2 == '0'):
                    print (true)
                else:
                    print (false)

            # 3. Father is unaffected, mother is affected
            elif not(is_F) and is_M:
                if child in mother and child != '0' and father == '0' and (m1 != '0' and m2 != '0'):
                    print (true)
                else:
                    print (false)
        else:
            child = tmp[c_index].split(':')[0]
            father = tmp[f_index].split(':')[0]
            mother = tmp[m_index].split(':')[0]
            print (vlist[i][0] + '\t' + vlist[i][1] + '\t' + vlist[i][3] + '\t' + vlist[i][4] + '\t' + '0' + '\t' + child + '\t' + father + '\t' + mother)


if __name__ == "__main__":
    main(sys.argv[1:])
