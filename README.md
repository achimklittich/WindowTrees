# WindowTrees README
## Intro
WindowTrees is a Python RAxML wrapper. WindowTrees allows the generation and evaluation of a large number of sliding window trees (trees built from only a certain window of the sequence) in short time. As input it uses consensus sequences in fasta format. These consensus sequences must have been created by mapping reads of different individuals or species to the same reference. All scaffolds or contigs that are included need to have the same name in every input file.

## What does WindowTrees do?
WindowTrees uses the consensus fasta files as input and creates windows of a previously specified size (-w) by "sliding" the window along the sequence. It is irrelevant whether short scaffolds or complete genomes are provided as input.

Afterwards the base content is determined. It is possible to specify a N-threshold to exclude windows with high numbers of missing data ("N"). Furthermore, WindowTrees determines the number of variable positions in the used windows.

The generated and checked window will then be used as input for RAxML to calculate the respective tree. It is necessary to provide an outgroup. In order to increase the speed, the tree generation step is parallelized.

In the final output WindowTrees provides folders that contain trees with the same topology. Additionally, it also counts the number of trees found for each topology.

WindowTrees provides several optional arguments to be adaptable to the specific purpose of different analyses. Beside the already mentioned N-threshold those are the possibility to specify gaps or overlaps between the windows (-lw) and a binary mode in RAxML (BINGAMMA) to only use transversions (--binary). For the latter the fasta files must be non-binary. WindowTrees will do the transformation. This can be especially useful when working with ancient DNA data.

The output trees can be used for various analyses like determining the direction of gene flow. Examples for its usage are Hempel et al. (2022), where the tool in its current form was first applied or Barlow et al. (2018), where the idea has first been applied. It also gives a good estimate on the possible topologies found in the sequences.

---

## Dependencies
### Python 3.7:
- ETE3 (Tree evaluation) [[ETEToolkit website]](http://etetoolkit.org/)
- pyfaidx (Sequence slicing) [[pyfaidx GitHub]](https://github.com/mdshw5/pyfaidx)

### Other Software
- RAxML [[RAxML website)]](https://cme.h-its.org/exelixis/web/software/raxml/)

---

## Installation
A simple and fast method to provide the required dependencies is to use Anaconda [[Anaconda install instructions]](https://docs.anaconda.com/anaconda/install/index.html) or miniconda [[Miniconda install instructions]](https://docs.conda.io/en/latest/miniconda.html#).

After the respective conda version has been installed, it is possible to install all required dependencies with the following command in the CLI:
```sh
# create conda environment called windowtrees and install dependencies from bioconda and conda-forge channel

conda create -n windowtrees -c bioconda -c conda-forge ete3 pyfaidx raxml python=3.7

# activate environment
source activate windwotrees

```

---

## Analysis Workflow
1. Provide several consensus sequences in fasta format that have been mapped to the same reference. The names of the scaffolds in the different consensus files must have the same name.
2. Create a colon seperated list of short_name:/path/to/fasta. Each row of the file can contain information for one sample with the according path. (see example input_test.tab)
3. Start WindowTrees with appropriate parameters for your anaylsis. 

Example WindowTrees call:
```sh
python3 -u windowtrees.py -o testrun --outgroup outgroup-name -w 200000 inputfile_test.tab
```
This will start the window tree analysis using WindowTrees. We define the window size with -w to be 200kB, the outgroup is defined by using the short_name from the input file. The generated output and the temporary files will be written into the folder testrun (-o).

4. Inspect and further process the generated and topology sorted trees.

## Usage
Calling windowtrees:
```sh
# example call
python3 -u windowtrees.py -o testrun --outgroup outgroup-name -w 200000 inputfile_test.tab
```

```sh
# help
python3 windowtrees.py -h
```

```
usage: python3 windowtrees.py [-h] -o OUTDIR --outgroup OUTGROUP -w WINDOWSIZE [-lw GAPSIZE] [--binary] [-N NTHRESHOLD] [--cpu NCPUS] inputfile

windowtrees.py: error: the following arguments are required: inputfile, -o, --outgroup, -w

WindowTrees: Calculate window trees for whole genomes to determine gene flow

required arguments:
  inputfile            Input file with short name and full path per line, colon separated
  -o OUTDIR            Output directory
  --outgroup OUTGROUP  Specify outgroup with short name as given in input file
  -w WINDOWSIZE        Specify window size

optional arguments:
  -h, --help           show this help message and exit
  -lw GAPSIZE          Positive (gap) or negative (overlap) integer [0]
  --binary             Activate binary mode with BINGAMMA model [GTRGAMMA] to only use transversions, ! input for binary mode must be non-binary fasta file !
  -N NTHRESHOLD        Ratio of allowed Ns per sequence in each window [0.1]
  --cpu NCPUS          Number of CPUs [2]


```

## Output
The program will generate the following output:
- as many folders as trees were found by the program numbered consecutively (e.g., tree_1, tree_2...). These folders contain the trees that were found in newick format. The name of each tree contains the following: the number of that tree topology, a consecutive number of this tree, RaxML_bestTree, the scaffold number.fa, the number of variable sites and the position of the window (example: tree_1_1_RAxML_bestTree.Scaffold_5.fa_VS-137_37400001-37420000.nw".
- a folder called "trees" which contains the raxml outfile for each of the found trees for each window. The name consists of the scaffold name, number of variable sites and window position (example: Scaffold_1.fa_VS-70_60300001-60320000.fa_outfile.txt)
- a file called "MSA_windowstats.out"
  This file gives the path to the output directory, the applied window size and the allowed N threshold at the top. Next it provides the statistics for each window giving the scaffold name, the window position, the number of Ns, if the N threshold was passed (failed/passed), amount of A, T, G and C, if the N threshold was passed (True/False) and the number of variable sites in case the N threshold was passed. 
- a file called "foundTrees.txt" gives a summary of the topologies found and the amount of each one.
 
Here is an example output of MSA_windowstats.out
```
## Output directory: /path/to/output/dir
## Windowsize: 20000	max(N): 50.0%
## Sequence;CountN;Filter;CountA;CountT;CountG;CountC;
#scaffold_1.fa_1-20000	PassNtreshold: False	VarSites: NA
Indiv1;20000;Failed;0;0;0;0;
Indiv2;20000;Failed;0;0;0;0;
Indiv3;20000;Failed;0;0;0;0;
Indiv4;20000;Failed;0;0;0;0;
.
.
.
#scaffold_25.fa_400001-420000	PassNtreshold: True	VarSites: 167
Indiv1;4795;Passed;4538;4089;3238;3340;
Indiv2;1444;Passed;5410;5050;3974;4122;
Indiv3;1722;Passed;5295;4987;3918;4078;
Indiv4;21;Passed;5773;5485;4317;4404;

```
Here is an example output of foundTrees.txt
```
tree_1	(((Indiv1,Indiv2),Indiv3),Indiv4);	4000
tree_2	(((Indiv1,Indiv3),Indiv2),Indiv4);	7000
tree_3	((Indiv1,(Indiv3,Indiv2)),Indiv4);	2000
```

## How to cite
When you use this tool please cite this github site and Hempel, E., Bibi, F., Faith, J.T., Koepfli, K.-P., Klittich, A.M., Duchêne, D.A., Brink, J.S., Kalthoff, D.C., Dalén, L., Hofreiter, M. & Westbury, M.V. 2022 Blue turns to gray - Paleogenomic insights into the evolutionary history and extinction of the blue antelope (Hippotragus leucophaeus). Molecular Biology and Evolution 39 (12): msac241; doi: [https://doi.org/10.1093/molbev/msac241](https://doi.org/10.1093/molbev/msac241).

## References
  - Barlow, A. Cahill, J.A., Hartmann, S., Theunert, Chr., Xenikoudakis, G., Fortes, G.G., Paijmans, J.L.A., Rabeder, G., Frischauf, Chr., Grandal-d'Anglade, A., Garcia-Vàzquez, A., Murtskhvaladze, M., Saarma, U., Anijalg, P., Skrbinšek, T., Bertorelle, G., Gasparian, B., Bar-Oz, G., Pinhasi, R., Slatkin, M., Dalén, L., Shapiro, B. & Hofreiter, M. 2018, Partial genomic survival of cave bears in living brown bears. Nature Ecology & Evolution 2: 1563-1570. doi: [https://doi.org/10.1038/s41559-018-0654-8](https://doi.org/10.1038/s41559-018-0654-8). 
  - Stamatakis, A. 214. RAxML Version 8: A tool for Phylogenetic Analysis and Post-Analysis of Large Phylogenies. Bioinformatics 30: 1312-1313. doi: [https://doi.org/10.1093/bioinformatics/btu033](https://doi.org/10.1093/bioinformatics/btu033).
  - Huerta-Cepas, J., Serra, F. & Bork, P. 2016. ETE 3: Reconstruction, analysis and visualization of phylogenomic data. Molecular Biology and Evolution 33: 1635-1638. doi: [https://doi.org/10.1093/molbev/msw046](https://academic.oup.com/mbe/article/33/6/1635/2579822).
  - Hempel, E., Bibi, F., Faith, J.T., Koepfli, K.-P., Klittich, A.M., Duchêne, D.A., Brink, J.S., Kalthoff, D.C., Dalén, L., Hofreiter, M. & Westbury, M.V. 2022 Blue turns to gray - Paleogenomic insights into the evolutionary history and extinction of the blue antelope (Hippotragus leucophaeus). Molecular Biology and Evolution 39 (12): msac241; doi: [https://doi.org/10.1093/molbev/msac241](https://doi.org/10.1093/molbev/msac241).
  - Shirley, M.D., Ma, Z., Pedersen, B.S. & Wheelan, S.J.. 2015. Efficient "pythonic" access to FASTA files using pyfaidx. PeerJ PrePrints 3: e970v1. doi: [https://doi.org/10.7287/peerj.preprints.970v1](https://doi.org/10.7287/peerj.preprints.970v1).
