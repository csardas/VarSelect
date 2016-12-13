#!/usr/bin/python

import os, sys, getopt, subprocess

def main(argv):
    usage = "usage: python Compare.py -v [vcf file] -p [ped file] -d [db file] -c [compare file]"

    try:
        opts, args = getopt.getopt(argv,"hv:p:d:c:",["vcf=","ped=","db=","cfile="])
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
        elif opt in ("-d", "--db"):
            dbfile = arg
        elif opt in ("-c", "--cfile"):
            comparefile = arg

    if len(opts) != 4:
        print ("try 'python Compare.py -h' for help")
        sys.exit(2)


    # get first child whom has been filterd and make clist
    f = open(comparefile, 'r')
    FL = f.readline().strip().split('\t') # First Line [list]
    FC = FL[6] # First Child
    FATHER = FL[7]
    MOTHER = FL[8]
    f.close()

    # get all other child, no matter it's affected or unaffected, also father and mother
    f2 = open(pedfile, 'r')

    OC = {} # Other Child [dict]
    for line in f2:
        if line.startswith('#'):
            continue

        tmp = line.strip().split('\t')
#        if tmp[1] != FC and tmp[2] != '0' and tmp[3] != '0':
        if tmp[1] != FC and tmp[1] != FATHER and tmp[1] != MOTHER:
            OC[tmp[1]] = tmp[5]
        elif tmp[1] == FATHER:
            f_ = tmp[1]
        elif tmp[1] == MOTHER:
            m_ = tmp[1]

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

    # get new child indexes, also father and mother
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

    f_index = tmp.index(f_)
    m_index = tmp.index(m_)

    f4.close()

    # comparison
    gene_dic = {}
    removed_count = 0
    removed_variants = {}
    REF = ['.', './.', '0/0']
    HET = ['0/1', '1/0', '0/2', '2/0', '0/3', '3/0']
    for a in affected_indexes:
        f4 = open(vcffile, 'r')
        for line in f4:
            if line.startswith('#'):
                pass
            else:
                tmp = line.strip().split('\t')
                chrom = tmp[0]
                pos = tmp[1]
                father = tmp[f_index].split(':')[0]
                mother = tmp[m_index].split(':')[0]

#                CSQ = tmp[7].split('CSQ')[1]
#                gene = CSQ.split('|')[3]

                genotype = tmp[a].split(':')[0]

                if chrom+'-'+pos in data:
                    gene = data[chrom+'-'+pos]
                    data[chrom+'-'+pos] = genotype
                    if genotype not in HET:
#                        del data[chrom+'-'+pos]
                        removed_count += 1
                        removed_variants[chrom+'-'+pos] = genotype
#                        print ('[Removing] Chrom: ' + chrom + '\tPos: ' + pos)
                    else:
                        if gene in gene_dic:
                            if (gene_dic[gene] == 'F' and father in REF and mother in HET) or (gene_dic[gene] == 'M' and father in HET and mother in REF):
                                gene_dic[gene] = 1
                        else:
                            if father in HET and mother in REF:
                                gene_dic[gene] = 'F'
                            elif father in REF and mother in HET:
                                gene_dic[gene] ='M'

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
                father = tmp[f_index].split(':')[0]
                mother = tmp[m_index].split(':')[0]

#                CSQ = tmp[7].split('CSQ')[1]
#                gene = CSQ.split('|')[3]

                genotype = tmp[u].split(':')[0]

                if chrom+'-'+pos in data:
                    gene = data[chrom+'-'+pos]
                    data[chrom+'-'+pos] = genotype
                    if genotype not in REF:
#                        del data[chrom+'-'+pos]
                        removed_count += 1
                        removed_variants[chrom+'-'+pos] = genotype
#                        print ('[Removing] Chrom: ' + chrom + '\tPos: ' + pos)
                    else:
                        if gene in gene_dic:
                            if (gene_dic[gene] == 'F' and father in REF and mother in HET) or (gene_dic[gene] == 'M' and father in HET and mother in REF):
                                gene_dic[gene] = 1
                        else:
                            if father in HET and mother in REF:
                                gene_dic[gene] = 'F'
                            elif father in REF and mother in HET:
                                gene_dic[gene] = 'M'

                # add is_model=0 data back
                if chrom+'-'+pos in other:
                    data[chrom+'-'+pos] = genotype

        f4.close()
        for i in range(len(clist)):
            if clist[i][0] + '-' + clist[i][1] in data:
                clist[i].append(data[clist[i][0]+'-'+clist[i][1]])

    # output compared data
    f5 = open(outputfile, 'w')
    f5.write(FL + '\n')
    for i in range(len(clist)):
        if clist[i][0] + '-' + clist[i][1] in data:
            # is_model=0
            if clist[i][4] == '0':
#                f5.write('\t'.join(clist[i]) + '\t' + data[clist[i][0]+'-'+clist[i][1]] + '\n')
                f5.write('\t'.join(clist[i]) + '\n')
            # after checking there are both father and mother in the same gene
            elif clist[i][5] in gene_dic:
                if gene_dic[clist[i][5]] == 1:
#                    f5.write('\t'.join(clist[i]) + '\t' + data[clist[i][0]+'-'+clist[i][1]] + '\n')
                    f5.write('\t'.join(clist[i]) + '\n')
                else:
                    clist[i][4] = '0'
                    f5.write('\t'.join(clist[i]) + '\n')
            elif clist[i][0] + '-' + clist[i][1] in removed_variants:
                clist[i][4] = '0'
                f5.write('\t'.join(clist[i]) + '\n') 
#        else:
#            clist[i][4] = '0'
#            f5.write('\t'.join(clist[i]) + '\t' + removed_variants[clist[i][0]+'-'+clist[i][1]] + '\n')
    f5.close()
    print ('Removed variants: ' + str(removed_count))

if __name__ == "__main__":
    main(sys.argv[1:])
