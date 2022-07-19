#!/usr/bin/env python3

# example call:
# /home/ak/anaconda3/bin/python -u windowtrees.py inputfile_test.tab -o testrun --binary --outgroup outroup-name -w 100000 -N 0.1 --cpu 14 -lw 80000
 
from pyfaidx import Fasta
import re
import os
import subprocess
from multiprocessing.pool import ThreadPool
from ete3 import Tree
import argparse
import time

def parsargs():
    parser = argparse.ArgumentParser(
        description='WindowTrees: Calculate window trees for whole genomes')
    # compulsory parameters 
    req = parser.add_argument_group(title='required arguments')
    req.add_argument('inputfile', type=str,
                        help='Colon seperated input file with short names and full path per line',default=False)
    req.add_argument('-o', dest='outdir', type=str,required=True,
                        help='Output directory')
    req.add_argument('--outgroup', dest='outgroup',required=True,
                        type=str, help='Specify outgroup with short name as given in input file')
    req.add_argument('-w', dest='windowsize', type=int,required=True,
                        help='Specify window size')  
    
    # optional parameters
    parser.add_argument('-lw',dest='gapsize',type=int,
                        help='Positive (gap) or negative (overlap) integer [0]', default=0)    
    parser.add_argument('--binary', dest='mode', action="store_true",
                        help='Activate binary mode with BINGAMMA model to only use transversions')
    parser.add_argument('-N', dest='nthreshold', type=float,
                        help='Ratio of allowed Ns per sequence in each window [0.1]',default=0.1)
    parser.add_argument('--cpu', dest='ncpus', type=int,
                        help='Number of CPUs to use [2]',default=2)
 
    args = parser.parse_args()
    return(args)

def read_input(inputfile):
    input_list = list()
    infile = open(inputfile, 'r')
    for line in infile:
        input_list.append(line.strip().split(':'))
    print(len(input_list))
    return(input_list)

def split_to_msa(inputfile, outdir):
    # read information from inputfile
    inputlist = read_input(inputfile)
    for i in range(len(inputlist)):
        print("Splitting "+inputlist[i][1]+"...")
        scaffolds = Fasta(inputlist[i][1])
        getkeys = scaffolds.keys()
        # split scaffolds in scaffoldwise MSA fasta files
        # according to same name (key)
        for key in getkeys:
            # create filename and write header and sequence
            filename = os.path.join(outdir, key+'.fa')
            with open(filename, 'a+') as output_msa:
                output_msa.write('>' + inputlist[i][0] + '\n')
                output_msa.write(str(scaffolds[key])+'\n')

def window_checker(indir, outdir, windowsize, nthreshold):
    # get all files from MSA folder
    files = os.listdir(indir)
    # write window statistics header
    with open(indir+'_windowstats.out', 'w') as out:
        print('## Output directory: '+outdir[:-4], file=out)
        print('## Window size: '+str(windowsize) +
              '\tmax(N): ' + str(nthreshold*100)+'%', file=out)
        print('## Sequence;CountN;Filter;CountA;CountT;CountG;CountC;', file=out)

    for data in files:
        if data.endswith('.fa') or data.endswith('.fasta'):
            # index fasta file for fast access
            specimen = Fasta(indir+'/'+data)
            # set start parameter
            start = 0
            stop = windowsize
            getkeys = specimen.keys()
            algnln = len(specimen[1])
            totsize = windowsize + gapsize
            if algnln > windowsize:
                # single window scaffold
                numwindows = 1
                if algnln > totsize:
                    numwindows = algnln//totsize
                for window in range(numwindows):
                    badwindow = False
                    # empty list
                    keylist = list()
                    msalist = list()
                    statslist = list()
                    filterlist = list()
                    for key in getkeys:
                        curr_seq = str(specimen[key][start:stop])
                        # get nucleotide composition and included N
                        nts = nt_comp(curr_seq)
                        # check for binary mode
                        if binary_mode == True:
                            curr_seq = binary_seq(curr_seq)
                        # append key msa and stats to lists
                        keylist.append(key)
                        msalist.append(curr_seq)
                        statslist.append(nts)
                        # check if N threshold exceeded
                        if nts[0]/windowsize > nthreshold:
                            badwindow = True
                            filterlist.append('Failed')
                        else:
                            filterlist.append('Passed')
            
                    if badwindow == True:
                        # don't check varsites for windows that did not pass threshold
                        varsites = "NA"
                        print('Window discarded, nthreshold exceeded: ' +
                              data+':'+str(start+1)+'-'+str(stop))
                    else:
                        # determine variable sites for window msa
                        varsites = check_varsites(msalist)
                        print('Using window: '+data+':' +
                              str(start+1)+'-'+str(stop))
                        # print msa to temp file
                        with open(outdir+'/'+data+'_VS-'+str(varsites)+'_'+str(start+1)+'-'+str(stop)+'.fa', 'w') as f:
                            for i in range(len(keylist)):
                                print('>'+keylist[i]+'\n'+msalist[i], file=f)
                    # print window stat output
                    with open(indir+'_windowstats.out', 'a') as out:
                        print('#'+data+'_'+str(start+1) +'-'+str(stop)+'\tPassNtreshold: '+str(not badwindow)+'\tVarSites: '+str(varsites), file=out)
                        for i in range(len(keylist)):
                            print(keylist[i] + ';' + str(statslist[i][0]) + ';'+filterlist[i]+';' + str(statslist[i][1]) +
                                ';' + str(statslist[i][2]) + ';' + str(statslist[i][3]) + ';' + str(statslist[i][4]) + ';', file=out)
                    # assign new start and stop
                    start = stop+gapsize
                    stop = start+windowsize
            else:
                print(data + ' Window discarded, scaffold smaller than specified window size')

def check_varsites(msalist):
    counter = 0
    for i in range(len(msalist[0])):
        usedlist = ['n','?']
        for j in range(len(msalist)):
            if msalist[j][i].casefold() not in usedlist:
                usedlist.append(msalist[j][i].casefold())
                if len(usedlist)>3:
                    counter = counter +1
                    break
    return(counter)

def nt_comp(sequence):
    ntlist = ['n','a','t','g','c']
    clist = list()
    for n in ntlist:
        clist.append(sequence.casefold().count(n))
    return(clist)

def binary_seq(sequence):
    sequence = sequence.casefold().replace('n', '?')
    sequence = sequence.casefold().replace('a', '0')
    sequence = sequence.casefold().replace('t', '1')
    sequence = sequence.casefold().replace('g', '0')
    sequence = sequence.casefold().replace('c', '1')
    return(sequence)


def purge(dir, pattern):
    for f in os.listdir(dir):
        if re.search(pattern, f):
            os.remove(os.path.join(dir, f))

def raxml_worker(inpath, outpath, outgroup, file):
    # create window MSA in size
    name = os.path.splitext(file)[0]
    print("Start RAxML: "+str(name))
    # select model
    if binary_mode == True:
        model = 'BINGAMMA'
    else:
        model = 'GTRGAMMA'

    process = subprocess.Popen(['raxmlHPC-PTHREADS', '-T', '2', '-m', model, '-p', '12345', '-s',
                                inpath+'/'+file, '-o', outgroup, '-n', name, '-w', outpath], stdout=subprocess.PIPE)
    # print stdout to terminal
    output = process.stdout.read()
    # grab output and write to file
    with open(outpath+'/'+file+'_outfile.txt', 'w') as f:
        f.write(output.decode())

def generate_trees(inpath, outpath, outgroup, threads):
    # get all temp msa files
    dirlist = os.listdir(inpath)
    que = list()
    for file in dirlist:
        if file.endswith('.fa') or file.endswith('.fasta'):
            que.append((str(inpath), str(outpath), str(outgroup), str(file)))
    # setup ThreadPool for parallel RAxML
    thpool = ThreadPool(int(threads))
    for sample in que:
        thpool.apply_async(raxml_worker, sample)
    thpool.close()
    thpool.join()

def get_similar_trees(path):
    treelist = list()
    dirlist = os.listdir(path)
    for file in dirlist:
        if file.startswith('RAxML_bestTree'):
            treelist.append([Tree(path+'/'+file), file])
    counter = 0
    while len(treelist) > 0:
        counter = counter+1
        tree_num = 1
        currtreedir = outdir+'/tree_'+str(counter)
        os.mkdir(currtreedir)
        curr_tree = treelist.pop(0)
        curr_tree[0].write(format=0, outfile=currtreedir+'/tree_' +
                           str(counter)+'_'+str(tree_num)+'_'+curr_tree[1]+'.nw')
        keeplist = list()
        for i in range(0, len(treelist)):
            comp_tree = treelist[i][0]
            rf_res = curr_tree[0].compare(comp_tree)
            if rf_res['rf'] == 0:
                tree_num = tree_num+1
                comp_tree.write(format=0, outfile=currtreedir+'/tree_' +
                                str(counter)+'_'+str(tree_num)+'_'+treelist[i][1]+'.nw')
            else:
                keeplist.append(treelist[i])

        treelist = keeplist
        treetopo = curr_tree[0].write(format=9)
        with open(outdir+'/foundTrees.txt', 'a') as f:
            f.write('tree_'+str(counter)+'\t'+treetopo+'\t'+str(tree_num)+'\n')
        print('tree_'+str(counter)+'\t'+treetopo+'\t'+str(tree_num))

def run_windowtrees(inputfile, outdir, outgroup, windowsize, nthreshold, ncpus):
    
    # split different input files to scaffoldwise msa
    msa_path = os.path.join(outdir, 'MSA')
    os.mkdir(msa_path)
    split_to_msa(inputfile, msa_path)
    
    # check windows
    tmp_win_msa_path = os.path.join(outdir, 'tmp')
    os.mkdir(tmp_win_msa_path)
    window_checker(msa_path, tmp_win_msa_path, windowsize, nthreshold)

    # remove scaffoldwise msa folder
    purge(msa_path, '.')
    os.rmdir(msa_path)
    
    # calculate trees
    treepath = os.path.join(outdir, 'trees')
    os.mkdir(treepath)
    generate_trees(tmp_win_msa_path, treepath, outgroup, ncpus/2)
    
    # remove window msa folder
    purge(tmp_win_msa_path, '.')
    os.rmdir(tmp_win_msa_path)
    
    # get sorted trees
    get_similar_trees(treepath)
    
    # remove trees temp folder
    #purge(treepath, '.')
    #os.rmdir(treepath)

# add check for at least 5 input files
# %%
if __name__ == '__main__':
    start_time = time.time()
    # parse args provided by call
    args = parsargs()

    inputfile = args.inputfile
    outdir = os.path.abspath(args.outdir)
    binary_mode = args.mode
    outgroup = args.outgroup
    windowsize = args.windowsize
    gapsize = args.gapsize
    nthres = args.nthreshold
    ncpus = args.ncpus

    run_windowtrees(inputfile, outdir, outgroup, windowsize, nthres, ncpus)
    print("Total: --- %s seconds ---" % (time.time() - start_time))
