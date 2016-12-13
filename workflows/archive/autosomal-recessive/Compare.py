#!/usr/bin/python

import os, sys, getopt

def main(argv):
    usage = "usage: python Compare.py -v [vcf file] -p [ped file] -c [compare file]"

    try:
        opts, args = getopt.getopt(argv,"hv:p:c:",["vcf=","ped=","cfile="])
    except getopt.GetoptError:
        print (usage)
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            print (usage)
            sys.exit()
        elif opt in ("-v", "--vcf"):
            vcffile = arg
            outputfile = os.path.basename(os.path.splitext(arg)[0]) + '_Compared.txt'
        elif opt in ("-p", "--ped"):
            pedfile = arg
        elif opt in ("-c", "--cfile"):
            comparefile = arg

    if len(opts) != 3:
        print ("try 'python Compare.py -h' for help")
        sys.exit(2)


    # get first child whom has been filterd and make clist
    f = open(comparefile, 'r')
    FL = f.readline().strip().split('\t') # First Line [list]
    FC = FL[5] # First Child
    FATHER = FL[6] # Father genotype
    MOTHER = FL[7] # Mother genotype
    f.close()

    # get all other child, no matter it's affected or unaffected
    f2 = open(pedfile, 'r')

    OC = {} # Other Child [dict]
    for line in f2:
        if line.startswith('#'):
            continue

        tmp = line.strip().split('\t')
        # it could be aunt or uncle
#        if tmp[1] != FC and tmp[2] != '0' and tmp[3] != '0':
        if tmp[1] != FC and tmp[1] != FATHER and tmp[1] != MOTHER:
            OC[tmp[1]] = tmp[5]

    if not OC:
        print ("[Notice] Only one child was founded.")
        sys.exit()
    f2.close()


    f = open(comparefile, 'r')
    tmp = []
    for line in f:
        if line.startswith('#'):
            pass
        else:
            tmp.append(line)
    f.close()

    clist = [x.strip().split('\t') for x in tmp]


    # seperate other child into two categories
    affected_child = []
    unaffected_child = []
    for k, v in OC.items():
        if v == '2':
            affected_child.append(k)
        else:
            unaffected_child.append(k)

    # parse is_model=1 data, else store to other
    data = {}
    other = {}
    for i in range(len(clist)):
        if clist[i][4] == '1':
            data[clist[i][0]+'-'+clist[i][1]] = clist[i][5]
        else:
            other[clist[i][0]+'-'+clist[i][1]] = clist[i][5]

    # get new child indexes
    f4 = open(vcffile, 'r')

    FL = '\t'.join(FL)
    affected_indexes = []
    unaffected_indexes = []
    for line in f4:
        if '#CHROM' in line:
            tmp = line.strip().split('\t')

    if affected_child:
        for i in affected_child:
            FL = FL + '\t' + i
            affected_indexes.append(tmp.index(i))

    if unaffected_child:
        for j in unaffected_child:
            FL = FL + '\t' + j
            unaffected_indexes.append(tmp.index(j))

    f4.close()

    # comparison
    removed_count = 0
    removed_variants = {}
    ALT = ['1', '2', '3', '1/1', '2/2', '3/3', '1/2', '2/1', '1/3', '3/1', '2/3', '3/2']
    for a in affected_indexes:
        f4 = open(vcffile, 'r')
        for line in f4:
            if line.startswith('#'):
                pass
            else:
                tmp = line.strip().split('\t')
                chrom = tmp[0]
                pos = tmp[1]
                genotype = tmp[a].split(':')[0]

#                if genotype == '.':
#                    genotype = './.'

                if chrom+'-'+pos in data:
                    if genotype not in ALT:
#                        del data[chrom+'-'+pos]
                        removed_count += 1
                        removed_variants[chrom+'-'+pos] = genotype
#                        print ('[Removing] Chrom: ' + chrom + '\tPos: ' + pos)
                    else:
                        data[chrom+'-'+pos] = genotype
                # add is_model=0 data back
                if chrom+'-'+pos in other:
                    data[chrom+'-'+pos] = genotype

        f4.close()
        for i in range(len(clist)):
            if clist[i][0] + '-' + clist[i][1] in data:
                clist[i].append(data[clist[i][0]+'-'+clist[i][1]])

    for u in unaffected_indexes:
        f4 = open(vcffile, 'r')
        for line in f4:
            if line.startswith('#'):
                pass
            else:
                tmp = line.strip().split('\t')
                chrom = tmp[0]
                pos = tmp[1]
                genotype = tmp[u].split(':')[0]

#                if genotype == '.':
#                    genotype = './.'

                if chrom+'-'+pos in data:
                    if genotype in ALT:
#                        del data[chrom+'-'+pos]
                        removed_count += 1
                        removed_variants[chrom+'-'+pos] = genotype
#                        print ('[Removing] Chrom: ' + chrom + '\tPos: ' + pos)
                    else:
                        data[chrom+'-'+pos] = genotype
                # add is_model=0 data back
                if chrom+'-'+pos in other:
                    data[chrom+'-'+pos] = genotype

        for i in range(len(clist)):
            if clist[i][0] + '-' + clist[i][1] in data:
                clist[i].append(data[clist[i][0]+'-'+clist[i][1]])
    f4.close()

    # output compared data
    f5 = open(outputfile, 'w')
    f5.write(FL + '\n')
    for i in range(len(clist)):
        if clist[i][0] + '-' + clist[i][1] in data:
            if clist[i][0] + '-' + clist[i][1] in removed_variants:
                clist[i][4] = '0'
            f5.write('\t'.join(clist[i]) + '\n')
#        elif clist[i][0] + '-' + clist[i][1] in removed_variants:
#            clist[i][4] = '0'
#            f5.write('\t'.join(clist[i]) + '\t' + removed_variants[clist[i][0]+'-'+clist[i][1]] + '\n')
#        if clist[i][0] + '-' + clist[i][1] in data:
#            f5.write('\t'.join(clist[i]) + '\t' + data[clist[i][0]+'-'+clist[i][1]] + '\n')
    f5.close()
    print ('Removed variants: ' + str(removed_count))

if __name__ == "__main__":
    main(sys.argv[1:])
