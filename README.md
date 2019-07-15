# BLANT: Basic Local Alignment for Networks Tool
## Analogy with **BLAST**
If you are in the bioinformatics field, you have probably heard of the tool named *BLAST*, the *Basic Local Alignment Search Tool*. BLAST is an algorithm for quickly finding local alignments in genomic (or proteomic) sequences. BLAST works by first creating a comprehensive database of all *k*-letter sequences (called "*k*-mers") that appear in the corpus of sequence to be searched and/or aligned. Such *k*-mers can be used to "seed" a local alignment between two distant regions of sequence. Below we show a hypothetical alignment between two distance regions of sequence, both of which contain the *boldfaced** *k*-mer "**GAGACCGT**":

  ACTAGATCCACCTCTAGGCAGGTAA **GAGACCGT** TGTTCTTCAGAGGTGAAGGAGACGCACAAACGGGCCC
  ACTAGATACACGTCTAGGCAGGTAG **GAGACCGT** AGTTCTTCATAGGTGACGGAGACGCACAAACGGGCCC

By storing every *k*-mer and its location, BLAST can "line up" the regions around two identical *k*-mers, and then check to see if this local alignment "extends" further beyond the *k*-mers. In the case above, even though the sequences contain minor differences, the fact that they contain the depicted *k*-mer seed means we can still find the near-perfect match. BLAST is extremely fast at performing the above operations, which is the reason BLAST has become the near-ubitquitous tool for comparing and aligning sequences that contain billions of letters. BLAST automaticall chooses the appropriate value of *k* to use in a particular search and alignment task. The general name for methods like this is *seed-and-extend* methods.

BLANT (*Basic local Aligment Network Tool*) is intended for form the basis of a tool similar to BLAST, but for networks: given an undirected network *G* as a list of edges, and a value of *k*, it samples *k*-node subgraphs called *k-graphlets*. Since the number of *k*-graphlets in a graph of *n* nodes is exponential in both *k* and *n*, BLANT does not exhaustively enumerate all *k*-graphlets, but instead samples as many as the user specifies. Furthermore, since uniform random sampling of *k*-graphlets is difficult, there are several choices among sampling methods, each with different trade-offs. Finally, BLANT allows for several different methods of output: it can produce *orbit-degree vectors* (ODVs) like *ORCA* (citation), or it can produce an explicit list of *k*-graphlets that can be used as seeds for later extension. At present, BLANT does not provide an "extend" functionality; there are *many* seed-and-extend local alignment algorithms in the literature, and although BLANT currently is by far the fastest method of producing a large number of seeds, we have not yet implemented or tested any extend algorithms; this is an area of future work. Despite the lack of an "extend" feature, BLANT is still capable of useful bioinformatics, as described in our first tool paper in the journal *Bioinformatics* (citation).

## USAGE
### Command-line arguments
bla bla bla

# Brief description of internals and source code structure
## Directories
### Directory libwayne
An extensive C-language library by Wayne Hayes, similar to a "class" library in C++ (but written mostly before C++ was invented), containing old and well-tested code implementing data structures and algorithms for graphs, sets, multisets, combinatorics, heaps, stacks, compressed flies, searching + sorting, linked lists, event-driven simulation and numerical analysis (the latter few not used much in BLANT but included for completeness). See [Wayne's Little Data Structures and Algorithms library](http://www.cs.toronto.edu/~wayne/libwayne). This library is necessary for BLANT and is created when "make" or "make all" is performed.

### Directory networks
A directory containing many example input networks in edge list (2-colunmns, space or tab separated) format. While BLANT also accepts other input formats such as GML and LEDA, edge lists are the simplest and most recommended.

### Directory orca_jesse_blant_table
Many methods of automatically (and manually) identifying all the various graphlets exist. There is the original manually-generated numbering of graphlets of size 3, 4, and 5 nodes ([Przulj *et al* 2004](https://academic.oup.com/bioinformatics/article/20/18/3508/202438)) which is also used by the popular program [ORCA](https://academic.oup.com/bioinformatics/article/30/4/559/205331), and the automatic methods from [Jesse](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0147078), and our own ([Hasan, Chung, Hayes](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0181570)) which has been modified from using an upper- to lower-triangle representation for BLANT, following Jesse's lead. This directory contains code to translate between all these different naming conventions.

### Directory regression_tests
As the name applies, this directory contains regression tests, in the following format: each subdirectory contains code and data for one regression test. Any filename in such a directory whose name ends in ".sh" will be automatically run by the program regression-test-all.sh, which is run nightly by the senior author on a Jenkins machine. If you add functinoality to BLANT, you are encouraged to create a directory here with your own regression test, and then do a pull request to me to add your code to the main BLANT repo (it's best to send me an email first at whayes@uci.edu to discuss it first).

### Directory canon_maps.correct
Contains correct output for comparison to files in canon_maps to ensure your installation of BLANT is working correctly.

### Directory "canon_maps"
Though not in the repo, this crucial directory is built when you install BLANT using "make" or "make all". It contains the lookup tables from all possible graphlets to the "canonical" ones; these lookup tables are the secret to BLANT's speed.  See [Hasan, Chung, Hayes (2017)](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0181570).

### Directory "Draw"
Contains code to create PDFs that depict any graphlet BLANT uses. Type "make Draw" to create the program Draw/graphette2dot. Then run ./Draw/graphette2dot for help and options. All graphette2dot does is create an input file and a command line for use by the "neato" command, which performs the actual creation of a PDF. So for example, to draw the 3-graphlet which is a triangle, run

./Draw/graphette2dot -k 3 -i 3 -o triangle | sh

That will create a file called "triangle.dot" and the command line for neato to convert triangle.dot to triangle.pdf; piping the output to sh simply runs the neato command, which actually generates the triangle.pdf file.

## Brief manifest and description of source files
### Makefile
Highlights: *k*=8 files may take quite a while to create (up to an hour); if you are only interested in sampling graphlets up to *k*=7, comment out the "EIGHT" variable near the top of the Makefile.

The first time you install BLANT, it's best to run "make all". This may take up to an hour (or a few minutes if you exclude *k*=8). This will create all the files necessary to run BLANT, as well as run some tests. Thereafter, any changes you make should require only to type "make". You can run sanity and regression tests using "make test_blant".

### blant-sanity.c
Run during "make test_blant" from the Makefile, it takes the output of BLANT's index mode and verifies that each graphlet with the same ID is actually identical.

### blant.c and blant.h
Main source code and header files that takes user input that specifies the value of *k*, the desired sampling method and number of samples, output mode, and input network. It then creates the internal graph, starts parallel invokations of BLANT if requested by the user (for increased speed), then performs graphlet sampling of the requested number of *k*-graphlets using the requested sampling method. Along the way it reads or mmap's several large files from the **canon_maps** directory.

### compare_canon_maps.cpp
Small code changes can sometimes result in different (but still correct) permutations between non-canonical and canonical graphlets (See [Hasan, Chung, Hayes (2017)](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0181570).) This program checks two different canonical permutation maps to ensure that they are equivalent module permutations. That is, one file is correct if and only if and the other is correct.

### compute-alphas-*
As a part of the unbiased MCMC graphlet sampling method, the *alpha* values determine the expected over/under representation of each type of graphlet sampled. These files compute those *alphas*.

### convert.cpp
Code that allows BLANT to take input of various graph representations such as edge list, GML, and LEDA formats.

### create-bin-data.c
This program reads the text-format output of fast-canon-map.c (see below), creating the internal lookup and permutation tables. Then, it simple dumps those tables into binary files that can be quickly read or mmap'd; mmap'ing these binary files is *hundreds* of times faster than reading the text files, which makes BLANT's startup time virtuall instantaneous even for the 1.25GB canon_map and permutation files required for *k*=8.

### fast-canon-map.c
This file creates the text version of the non-canonical to canonical lookup table and permutation map for any value of *k* up to 8. It's only run once in the Makefile to create files in the **canon_maps** directory; these text output files will be converted to binary files by create-bin-data.c (see above). fast-canon-map runs virtually instantaneously for all values of *k* up to 6; *k*=7 takes less than a minute, while *k*=8 can take anywhile from 5 minutes to an hour depending on the speed of your computer. These files normally only need to be created once. If you are not interested in *k*=8, you can comment out the variable **EIGHT** in the Makefile.

### libblant.c
A collection of often-used routines needed by most of the other C files including blant.c, fast-canon-map.c, create-bin-data.c, etc. Prototypes for the functions contain herein are all in blant.h

### magictable.cpp
Code to convert between all the various naming/numbering schemes of graphlets (see **Directory orca_jesse_blant_table** above).

### make-orbit-maps.c
After creating the list of all non-canonical and canonical graphlets for a given value of *k*, this program enumerates all the orbits of each canonical, putting the data in to the file canon_maps/orbit_map*k*.txt

### make-subcanon_maps.c
Each canonical graphlet of size *k* contains exactly *k* subgraphs of size (*k*-1) (not all of which may be connected). This program computes these subgraphs and their canonical IDs. Although not yet used, these sub-canon maps will be used later for more efficient search and alignment.

### makeEHD.c
*EHD* stands for *Edge Hamming Distance*: it is the minimum number of edges required to convert one canonical graphlet to another, and although not yet used by BLANT, will be used later to implement partial / approximate graphlet matching (ie., matching graphlets with missing or exta edges). This file creates the lookup table representing the EHD between any two canonical graphlets.

### regression-test-all.sh
Script to iterate through each directory in **Directory regression-tests**, running any script inside whose name ends in ".sh"; such scripts should perform one or more regression tests using other files in the given directory.

### slow-canon-maps.c
The old (deprected) algorithm described in [Hasan, Chung, Hayes](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0181570). Replaced by fast-canon-maps.c above, which is exponentially faster.

### test-bin-data.c
Testing code for files created by create-bin-data.c above.
