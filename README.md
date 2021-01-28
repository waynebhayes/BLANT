# BLANT: Basic Local Alignment for Networks Tool
## Analogy with **BLAST**
If you are in the bioinformatics field, you have probably heard of the tool named *BLAST*, the *Basic Local Alignment Search Tool*. BLAST is an algorithm for quickly finding local alignments in genomic (or proteomic) sequences. BLAST works by first creating a comprehensive database of all *k*-letter sequences (called "*k*-mers") that appear in the corpus of the sequence to be searched and/or aligned. Such *k*-mers can be used to "seed" a local alignment between two distant regions of sequence. Below we show a hypothetical alignment between two distant regions of sequence, both of which contain the **boldfaced** *k*-mer "**GAGACCGT**":

  ACTAGAT*C*CAC*C*TCTAGGCAGGTA **GAGACCGT** GTTCTTCA*G*AGGTGA*A*GGAGACGCACAAACGGGCCC
 
  ACTAGAT*A*CAC*G*TCTAGGCAGGTA **GAGACCGT** GTTCTTCA*T*AGGTGA*C*GGAGACGCACAAACGGGCCC

By storing every *k*-mer and its location, BLAST can "line up" the regions around two identical *k*-mers, and then check to see if this local alignment "extends" further beyond the *k*-mers. In the case above, even though the sequences contain minor differences (highlighted with *italics*), the fact that they contain the depicted *k*-mer seed means we can still find the near-perfect match. BLAST is extremely fast at performing the above operations, which is the reason BLAST has become the near-ubiquitous tool for comparing and aligning sequences that contain billions of letters. BLAST automatically chooses the appropriate value of *k* to create *k*-mer seeds in a particular search and alignment task, and uses a sophisticated extend algorithm to create full seed-and-extend local alignments.

Our new tool, called BLANT (*Basic local Aligment Network Tool*), is intended to form the basis of a seed-and-extend local alignment algorithm, but for networks: given an undirected network *G*, and a value of *k*, it samples connected *k*-node subgraphs called *k-graphlets*. Since the number of *k*-graphlets in a graph of *n* nodes is exponential in both *k* and *n*, BLANT does not exhaustively enumerate all *k*-graphlets, but instead randomly samples as many as the user specifies. Furthermore, since uniform random sampling of *k*-graphlets is difficult, there are several choices among sampling methods, each with different trade-offs. Finally, BLANT allows for several different methods of output: it can produce *orbit-degree vectors* (ODVs) like [ORCA](https://academic.oup.com/bioinformatics/article/30/4/559/205331), or graphlet frequencies, or an explicit list of *k*-graphlets that can be used as seeds for later extension. At present, BLANT does not provide an "extend" functionality; there are *many* seed-and-extend local alignment algorithms in the literature, each with its own method of seeding and extending. Although BLANT currently is by far the fastest method of producing a large number of *seeds*, we have not yet tested how the various extend algorithms perform using our seeds; this is a clear area of future work, and suggestions are welcome. Despite the lack of an "extend" feature, BLANT is still capable of useful bioinformatics, as described in our first tool paper in the journal *Bioinformatics* (citation to be filled later).

## USAGE
### Quick Start guide
#### Stack Size
Before starting *anything* below, you need to ensure your OS allows your programs enough stack space. Some machines today still ship with a default limit of 8MB---a *ridiculously* small limit on any machine built after about 1999.  You do not require sudo privileges to change this. If you're running Linux or MacOS, type "ulimit -s unlimited" to your Bash shell; if you're running any other system, you're on your own.

### Building BLANT for the first time
To make and then test *everything* just type

    ./regression-test-all.sh -make
  
Be warned it may take up to an hour, especially creating the k=8 lookup tables.
Note that doing the above performs a "make pristine", which is even cleaner than "make clean". In particular, "clean" cleans all executables but doesn't touch the lookup tables; use "make pristine" to remove even the lookup tables (which shouldn't be necessary since they never change, though they take up disk space).

Once the above is done, you should have an executable called "blant" in the main repo directory, as well as many files in the directory *canon_maps*; these are the lookup tables and associated files that allow BLANT to sample graphlets so fast. Except for major BLANT upgrades, you normally shouldn't *ever* need to touch anything in the canon_maps directory; make it once and forget about it. Note that even "make clean" doesn't nuke the contents of canon_maps; to do that, type "make pristine".

### Required Command-line arguments
SYNOPSIS (things inside {} are mandatory; those inside [] are optional; if no source network is given, it is read from the standard input):

    blant {-k graphletSize} {-n numSamples} {-m outputMode} {-s samplingMethod} [-d displayMode] [-t THREADS] [-r randomSeed] [sourceNetwork]
#### -k *k*: graphlet size (from *k*=3 up to *k*=8) of graphlets to sample
*k* is an integer that can be 3 through 8 inclusive, and represents the number of nodes in the graphlets that will be sampled. For now BLANT almost always samples only connected graphlets, returning a disconnected "graphette" only if there is a connected component in the input graph that has fewer than *k* nodes; such disconnected graphlets are *never* returned by the MCMC method.

#### -n {integer}: the number of random graphlet samples to take
*n* is a non-negative integer representing the total number of graphlets that should be sampled. 0 is allowed; values into the billions are feasible (more if you use parallelism).

#### -m{outputMode}: how to use the sampled graphlets on output. All output goes to the standard output (aka "terminal")
*outputMode* is a single character, one of the following:
##### -m*f*D: frequency mode (with display D option)
BLANT's output will consist of two columns: the first column is a frequency, the 2nd column is the ID of the canonical graphlet. Lines are output in order of the canonical graphlets, and the number of lines of output is exactly equal to the number of canonical graphlets that exist for the given value of *k*; even graphlets with a zero count are output. Furthermore, another character **D** can be appended to -mf, which can be either an *i* (integer) representing that the frequency should be displayed as a raw integer count, or *d* (density) in which case the frequencies are normalized to a "density" (aka "concentration"); the latter is useful in comparing frequencies that used different numbers of samples (cf. -n option). If D is omitted, it defaults to *i* for all sampling methods except MCMC, which defaults to *d*.

##### -m*o*: ORCA/ODV (*O*rbit *D*egree *V*ector) format
The output here is similar to ORCA's output: the number of lines is exactly equal to the number of nodes in the input graph, and the number of columns is equal to the number of orbits that exist for the given value of *k*. The integer in each location is the number of times that node "touched" the specified orbit. Unlike ORCA, where this value is a deterministic constant since ORCA performs exhaustive enumeration, the value output by BLANT will be stochastic and depend on the sampling method, the number of samples, and any bias inherent in the random sampling method. However, these values should be roughly proportional to ORCA's output (modulo sample size), and the similarity should improve with increasing sample size.

##### -m*g*: GDV (*G*raphlet *D*egree *V*ector) format
Similar to ODV format above, except for each node, the columns are the number of times that node touches a particular *graphlet, **independent of which orbit is touched**.* **NOTE**: beware that many authors use the term GDV when they are in fact referring to an ODV. BLANT is more precise in clearly distinguishing between the two.

##### -m*i*: graphlet indexing mode
This is the mode that is most unique to BLANT: it outputs as many lines as there are samples (see -n above)---beware that if *n* is large, this can produce *huge* output files. Each line consists of exactly (*k*+1) columns: first, the canonical ID of the graphlet that was sampled, and then exactly *k* node identifiers (using whatever node naming scheme that was in the network input file). **Most importantly**, BLANT guarantees that for any two lines that have the same ID in the first column (ie., they are the same graphlet), the remaining columns are in an order that imposes an exact local alignment between the two. (There may be more than one local alignment depending on the orbits in said graphlet, but BLANT only outputs one such local alignment.) This is the mode that can be used to create a database of *k*-graphlets that can act as seeds for future implementations of an "extend" algorithm to produce local alignments with more than *k* nodes.  (Note, however, that while BLAST can produce an *exhaustive* list of all *k*-mers in a sequence corpus since their number is linear in the amount of sequence, BLANT only produces a random sample of *k*-graphlets; this means an extend algorithm cannot assume that the list of identical *k*-graphlets produced by BLANT is exhaustive. However, with enough samples, every node in a graph corpus can be "covered" by multiple *k*-graphlet seeds, and we believe this multiple coverage can be leveraged to produce larger local alignments.)

##### -m*j*: orbit indexing mode
Similar to -m*i* above, except nodes that occupy the same orbit are output separated by colons. Thus, the number of space-separated columns (excluding the first column, which is still the canonical ID of the sampled graphlet) is the number of orbits in the graphlet. Unlike -m*i* mode, the order of the nodes is *not* guaranteed to impose a local alignment, but instead all that is guaranteed is that if the same set of nodes is sampled more than once, they will be output in a constant order; this allows duplicates to be detected.

#### -s {samplingMethod}: the method used to sample graphlets
BLANT can sample graphlets in many different ways, each with advantages and disadvantages. The allowed values are:

##### NBE (Node Based Expansion)
BLANT picks an edge from the input network uniformly at random. The two endpoint nodes initialize the set *S* of nodes. The remaining (*k*-2) nodes are added to *S* by finding all nodes that are adjacent to the current set of nodes (excluding those already in *S*), and then picking one such node uniformly at random. We return the sampled graphlet when |*S*|=*k*. If at any point the set of nodes adjacent to *S* (excluding nodes inside *S*) is null, then *S* is retained by selecting a new starting location chosen uniformly at random. This can only happen if the initial edge is picked from a connected component *C* that has fewer than *k* nodes; in that case, *S* will be a disconnected graphlet (aka *graphette*) that consists of *C*, plus (*k*-|*C*|) nodes from elsewhere in the input network. The NBE method is not guaranteed to produce graphlet samples that are unbiased, although empirically we have found it is not terribly biased. Each sample produced by NBE starts at a new random edge, and so NBE tends to "see" most regions of the network even if the number of samples is small.

##### EBE (Edge Based Expansion)
Node based expansion (NBE) can be slow if the mean degree of the network is large. For networks with high mean degree, the EBE sampling method is asymptotically faster than NBE, though it may be more biased. Like NBE, it starts each sample with a freshly chosen edge from the network, so that widely separated regions are sampled with ease. Then, one edge is chosen uniformly at random from all edges emanating from *S*. Note that, while NBE choses among nodes adjacent to *S* with equal probability, EBE follows *edges* emanating from *S* with equal probability. This means that (a) the next node added to *S* is more likely to be a neighbor of a node in *S* with many edges emanating from *S* than one with few edges emanating from *S*---in other words, the choice of next node added to *S* is dominated by "hubs" already in *S*; and (b) a node adjacent to *S* is chosen with probability proportional to the number of edges that reach it from inside *S*.

##### RES (Reservior sampling)
Starting with an NBE-created sample of *S* nodes, we take a small number of "steps" where one node is deleted from *S* and simultaneously a randomly chosen node adjacent to *S* is added. This random walk erases the memory of where NBE started and tends to reduce the bias of NBE.

##### MCMC (Markov Chain Monte Carlo)
This method is the fastest *per sample*, and in the long run is also guaranteed to produce unbiased samples (ie., those whose relative frequencies approach those of ORCA's exhaustive enumeration). It starts with a randomly chosen edge, then builds *S* similar to NBE, and outputs one sample. However, unlike the other methods, it takes a *long* random walk by deleting one node from *S* and randomly adding one node adjacent to *S*, and then outputting that new set as a new sample. This has the effect that, on short timescales, each sample differs from the previous sample by only one (or a few) nodes. Thus, for example, in indexing mode (-m*i*), adjacent lines in the output will have many nodes in common. This has the disadvantage that it may take a long time for the MCMC method to "visit" all regions of the network. On the other hand, MCMC is very fast *per sample*, because only one (or a few) nodes change between samples. Furthermore, the MCMC method computes the bias on-the-fly, and in the long run the -m*f* mode will output concentrations that are asymptotically correct---although indexing mode should not be used to estimate graphlet frequencies because indexing mode does not take these computed biases into account.

##### AR (accept/reject)
AR is a horrendously slow but asymptotically correct method of picking *k* nodes uniformly at random, and accepting the sample only if the *k* nodes form a connected graphlet. For large values of *k* and for most typical input networks that are sparse, this results in the vast majority of *k*-node sets being rejected. This method is not recommended since it is exponentially slow with large values of *k*.

### Optional command-line arguments
#### -t*N*: run BLANT in parallel mode with *N* threads.
For all modes except indexing (-mi) mode, speedup should be linear (this has been tested on a machine with 64 cores, and speedup is linear.) However, with indexing mode, the processes are I/O bound and speedup is only linear up to about 8 threads, and little speedup is observed beyond 8 threads.
#### -d*m*: displayMode
Display mode determines how the graphlet ID of the sampled graphlet is displayed: default is *i*, BLANT's internal integer ordinal of the canonical (order similar to that described in [Hasan, Chung, Hayes 2017](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0181570) except using the lower rather than upper triangle, to be more compatible with [Jesse](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0147078). Other formats include *d*, the decimal value of the canonical; *b* the same integer dispalyed as binary; *j* use Jesse's ID; and *o* use ORCA's ID.

#### -m*w*: Windowing mode
Experimental, may be discontinued.

#### -w*l*: Window sampling method
Currently unused and experimental; may be removed.

# Brief description of internals and source code structure
## Directories
### Directory libwayne
An extensive C-language library by Wayne Hayes, similar to a "class" library in C++ , containing old and well-tested code implementing data structures and algorithms for graphs, sets, multisets, combinatorics, heaps, stacks, compressed flies, searching + sorting, linked lists, event-driven simulation and numerical analysis (the latter few not used much in BLANT but included for completeness). See [Wayne's Little Data Structures and Algorithms library](http://www.cs.toronto.edu/~wayne/libwayne). This library is necessary for BLANT and is created when "make" or "make all" is performed; libwayne was written mostly before C++ was invented, and most of its functionality is significantly faster and more memory efficient than equivalent C++ classes.

### Directory networks
A directory containing many example input networks in edge list (2-columns, space or tab separated) format. While BLANT also accepts other input formats such as GML and LEDA, edge lists are the simplest and most recommended.

### Directory orca_jesse_blant_table
Many methods of automatically (and manually) identifying all the various graphlets exist. There is the original manually-generated numbering of graphlets of size 3, 4, and 5 nodes ([Przulj *et al* 2004](https://academic.oup.com/bioinformatics/article/20/18/3508/202438)) which is also used by the popular program [ORCA](https://academic.oup.com/bioinformatics/article/30/4/559/205331), and the automatic methods from [Jesse](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0147078), and our own ([Hasan, Chung, Hayes](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0181570)) which has been modified from using an upper- to lower-triangle representation for BLANT, following Jesse's lead. This directory contains code to translate between all these different naming conventions, including both the graphlet and orbit numbering schemes.

### Directory regression_tests
As the name applies, this directory contains regression tests, in the following format: each subdirectory contains code and data for one regression test. Any filename in such a directory whose name ends in ".sh" will be automatically run by the program regression-test-all.sh, which is run nightly by the senior author on a Jenkins machine. If you add functionality to BLANT, you are encouraged to create a directory here with your own regression test, and then do a pull request to me to add your code to the main BLANT repo (it's best to send me an email first at whayes@uci.edu to discuss it first).

### Directory canon_maps.correct
Contains correct output for comparison to files in canon_maps to ensure your installation of BLANT is working correctly.

### Directory "canon_maps"
Though not in the repo, this crucial directory is built when you install BLANT using "make" or "make all". It contains the lookup tables from all possible graphlets to the "canonical" ones; these lookup tables are the secret to BLANT's speed. See [Hasan, Chung, Hayes (2017)](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0181570).

### Directory "Draw"
Contains code to create PDFs that depict any graphlet BLANT uses. Type "make Draw" to create the program Draw/graphette2dot. Then run ./Draw/graphette2dot for help and options. All graphette2dot does is create an input file and a command line for use by the "neato" command, which performs the actual creation of a PDF. So for example, to draw the 3-graphlet which is a triangle, run

./Draw/graphette2dot -k 3 -i 3 -o triangle | sh

That will create a file called "triangle.dot" and the command line for neato to convert triangle.dot to triangle.pdf; piping the output to sh simply runs the neato command, which actually generates the triangle.pdf file.

## Brief manifest and description of source files
### Makefile
Highlights: *k*=8 files may take quite a while to create (up to an hour); if you are only interested in sampling graphlets up to *k*=7, comment out the "EIGHT" variable near the top of the Makefile, and the "make" will finish in a few minutes rather than an hour. ("make all" also performs some extensive testing and will take longer even without EIGHT.)

The first time you install BLANT, it's best to run "make all". This may take up to an hour (or a few minutes if you exclude *k*=8). This will create all the files necessary to run BLANT, as well as run some tests. Thereafter, any changes you make should require only to type "make". You can run sanity and regression tests using "make test_blant".

### blant-sanity.c
Run during "make test_blant" from the Makefile, it takes the output of BLANT's index mode and verifies that each graphlet with the same ID is actually identical.

### blant.c and blant.h
Main source code and header files that takes user input that specifies the value of *k*, the desired sampling method and number of samples, output mode, and input network. It then creates the internal graph, starts parallel invokations of BLANT if requested by the user (for increased speed), then performs graphlet sampling of the requested number of *k*-graphlets using the requested sampling method. Along the way it reads or *mmap()*'s several large files from the **canon_maps** directory.

### compare_canon_maps.cpp
Small code changes can sometimes result in different (but still correct) permutations between non-canonical and canonical graphlets (See [Hasan, Chung, Hayes (2017)](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0181570).) This program checks two different canonical permutation mappings to ensure that they are equivalent. That is, one file is correct if and only if the other is correct.

### compute-alphas-{MCMC|NBE}
As a part of the unbiased MCMC graphlet sampling method, and as a preliminary implementatino of an unbiased NBE method, the *alpha* values determine the expected over/under representation of each type of graphlet sampled. These programs compute those *alpha* values.

### convert.cpp
Code that allows BLANT to take input of various graph representations such as edge list, GML, and LEDA formats.

### create-bin-data.c
This program reads the text-format output of fast-canon-map.c (see below), creating the internal lookup and permutation tables. Then, it simply dumps those tables into binary files that can be quickly read or *mmap()*'d; *mmap()*'ing these binary files is *hundreds* of times faster than reading the text files, which makes BLANT's startup time virtually instantaneous even for the 1.25GB canon_map and permutation files required for *k*=8.

### fast-canon-map.c
This file creates the text version of the non-canonical to canonical lookup table and permutation map for any value of *k* from 3 to 8. It's only run once in the Makefile to create files in the **canon_maps** directory; these text output files will be converted to binary files by create-bin-data.c (see above). fast-canon-map runs virtually instantaneously for all values of *k* up to 6; *k*=7 takes less than a minute, while *k*=8 can take anywhere from 5 minutes to an hour depending on the speed of your computer. These files normally only need to be created once. If you are not interested in *k*=8, you can comment out the variable **EIGHT** in the Makefile.

### libblant.c
A collection of often-used routines needed by most of the other C files including blant.c, fast-canon-map.c, create-bin-data.c, etc. Prototypes for the functions contained herein are all in blant.h

### magictable.cpp
Code to convert between all the various naming/numbering schemes of graphlets (see **Directory orca_jesse_blant_table** above).

### make-orbit-maps.c
After creating the list of all non-canonical and canonical graphlets for a given value of *k*, this program enumerates all the orbits of each canonical, putting the data in to the file canon_maps/orbit_map*k*.txt

### make-subcanon_maps.c
Each canonical graphlet of size *k* contains exactly *k* subgraphs of size (*k*-1) (not all of which may be connected). This program computes these subgraphs and their canonical IDs. Although not yet used, these sub-canon maps will be used later for more efficient search and alignment.

### makeEHD.c
*EHD* stands for *Edge Hamming Distance*: it is the minimum number of edges required to convert one canonical graphlet to another, and although not yet used by BLANT, they will be used later to implement partial / approximate graphlet matching (ie., matching graphlets with missing or exta edges). This file creates the lookup table representing the EHD between any two canonical graphlets.

### regression-test-all.sh
Script to iterate through each directory in **Directory regression-tests**, running any script inside whose name ends in ".sh"; such scripts should perform one or more regression tests using other files in the given directory; a regression test that fails should return a non-zero exit status, which will cause the parent regression to fail as well.

### slow-canon-maps.c
The old (deprected) algorithm described in [Hasan, Chung, Hayes](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0181570). Replaced by fast-canon-maps.c above, which is exponentially faster.

### test-bin-data.c
Testing code for files created by create-bin-data.c above.

# References
#### Our first paper on how BLANT performs so quickly on up to k=8 node graphlets: Hasan, Adib, Po-Chien Chung, and Wayne Hayes. ["Graphettes: Constant-time determination of graphlet and orbit identity including (possibly disconnected) graphlets up to size 8." PloS one 12, no. 8 (2017): e0181570.](https://doi.org/10.1371/journal.pone.0181570)

#### First conference presentation about BLANT: Hayes, W. and Maharaj, S., 2018. [BLANT: sampling graphlets in a flash. q-bio, Rice University, Houston, Texas, USA.](http://q-bio.org/wp/wp-content/uploads/2018/04/Maharaj_Sridevi.pdf)

#### Journal "tool" announcement in BioInformatics: Maharaj, Sridevi, Brennan Tracy, and Wayne B. Hayes. ["BLANT—fast graphlet sampling tool." Bioinformatics 35, no. 24 (2019): 5363-5364.](https://doi.org/10.1093/bioinformatics/btz603)

