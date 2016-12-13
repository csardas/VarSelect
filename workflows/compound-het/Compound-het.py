#!/usr/bin/python

import os, sys, getopt, subprocess

def main(argv):
    usage = "usage: python Compound-het.py -v [vcf file] -p [ped file] -d [db file]"

    try:
        opts, args = getopt.getopt(argv,"hv:p:d:",["vcf=","ped=","db="])
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
        elif opt in ("-d", "--db"):
            dbfile = arg

    if len(opts) != 3:
        print ("try 'python Compound-het.py -h' for help")
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
    f.close()

    # get fist affected sample
    f = open(pedfile, 'r')
    gender = False
    for line in f:
        if line.startswith('#'):
            continue

        tmp = line.strip().split('\t')
        if tmp[1] not in parents and tmp[5] == '2':
            c_ = tmp[1]
            f_ = tmp[2]
            m_ = tmp[3]
            if tmp[4] == '1':
                gender = True
            break # only take one affected child
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

    # store gene to list from gemini db
    output = subprocess.Popen(['gemini', 'query', '-q', 'SELECT gene FROM variants', dbfile], stdout=subprocess.PIPE).communicate()[0].strip()
    gene_list = output.split('\n')

    gene_dic = {}
    def is_HET(genotype):
        tmp = genotype.split('/')
        if len(tmp) != 2:
            return False

        g1 = tmp[0]
        g2 = tmp[1]

        if g1 == '.':
            g1 = '0'
        if g2 == '.':
            g2 = '0'

        if (g1 =='0' and g2 != '0') or (g2 == '0' and g1 != '0'):
            return True
        else:
            return False


    for i in range(len(vlist)):
        gene = gene_list[i]

        if gene == 'None':
            continue

        tmp = '\t'.join(vlist[i])
        tmp = tmp.replace('.', '0').split('\t')

        child = tmp[c_index].split(':')[0]
        father = tmp[f_index].split(':')[0]
        mother = tmp[m_index].split(':')[0]

        if child == '0':
            child = '0/0'
            tmp[c_index] = '0/0'
        if father == '0':
            father = '0/0'
            tmp[f_index] = '0/0'
        if mother == '0':
            mother = '0/0'
            tmp[m_index] = '0/0'

        c1 = tmp[c_index][0]
        c2 = tmp[c_index][2]

        # pass if chrom = Y or MT or in Non PAR region
        if vlist[i][0] == 'Y' or vlist[i][0] == 'MT':
            vlist[i].append('0')
            continue

        if (vlist[i][0] == 'X' and ((int(vlist[i][1]) < 60001) or (2699520 <= int(vlist[i][1]) <= 154931044) or (int(vlist[i][1]) > 155270560))):
            vlist[i].append('0')
            continue

        if (is_HET(child)) and ((child == father and mother == '0/0') or (child == mother and father == '0/0')):
            if gene in gene_dic:
                if (gene_dic[gene] == 'F' and ((c1 == '0' and c2 in mother) or (c2 == '0' and c1 in mother))) or (gene_dic[gene] == 'M' and ((c1 == '0' and c2 in father) or (c2 == '0' and c1 in father))):
                    gene_dic[gene] = 1
            else:
                if (c1 == '0' and c2 in father or c2 == '0' and c1 in father):
                    gene_dic[gene] = 'F'
                elif (c1 == '0' and c2 in mother or c2 == '0' and c1 in mother):
                    gene_dic[gene] = 'M'

            vlist[i].append('1')
        else:
            vlist[i].append('0')


    print ('#chrom\tpos\tref\talt\tis_CH\tgene\t'+ c_ + '\t' + f_ + '\t' + m_)

    for j in range(len(vlist)):
        gene = gene_list[j]

        tmp = '\t'.join(vlist[j])
        tmp = tmp.replace('.', '0').split('\t')

        child = tmp[c_index].split(':')[0]
        father = tmp[f_index].split(':')[0]
        mother = tmp[m_index].split(':')[0]

        true = vlist[j][0] + '\t' + vlist[j][1] + '\t' + vlist[j][3] + '\t' + vlist[j][4] + '\t' + '1' + '\t' + gene + '\t' + vlist[j][c_index].split(':')[0] + '\t' + vlist[j][f_index].split(':')[0] + '\t' + vlist[j][m_index].split(':')[0]
        false = vlist[j][0] + '\t' + vlist[j][1] + '\t' + vlist[j][3] + '\t' + vlist[j][4] + '\t' + '0' + '\t' + gene + '\t' + vlist[j][c_index].split(':')[0] + '\t' + vlist[j][f_index].split(':')[0] + '\t' + vlist[j][m_index].split(':')[0]

        if vlist[j][-1] == '1':
            if gene in gene_dic:
                if gene_dic[gene] != 1:
                    print (false)
                else:
                    print (true)
        else:
            print (false)

   
if __name__ == "__main__":
    main(sys.argv[1:])
