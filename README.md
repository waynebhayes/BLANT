#BLANT - Short Documentation
## Analogy with **BLAST**
If you are in the bioinformatics field, you have probably heard of the tool named *BLAST*, the *Basic Local Alignment Search Tool*. BLAST is an algorithm for quickly finding local alignments in genomic (or proteomic) sequences. BLAST works by first creating a comprehensive database of all *k*-letter sequences (called "*k*-mers") that appear in the corpus of sequence to be searched and/or aligned. Such *k*-mers can be used to "seed" a local alignment between two distant regions of sequence. Below we show a hypothetical alignment between two distance regions of sequence, both of which contain the *boldfaced** *k*-mer "**GAGACCGT**":

  ACTAGATCCACCTCTAGGCAGGTAA **GAGACCGT** TGTTCTTCAGAGGTGAAGGAGACGCACAAACGGGCCC
  ACTAGATACACGTCTAGGCAGGTAG **GAGACCGT** AGTTCTTCATAGGTGACGGAGACGCACAAACGGGCCC

By storing every *k*-mer and its location, BLAST can "line up" the regions around two identical *k*-mers, and then check to see if this local alignment "extends" further beyond the *k*-mers. In the case above, even though the sequences contain minor differences, the fact that they contain the depicted *k*-mer seed means we can still find the near-perfect match. BLAST is extremely fast at performing the above operations, which is the reason BLAST has become the near-ubitquitous tool for comparing and aligning sequences that contain billions of letters. BLAST automaticall chooses the appropriate value of *k* to use in a particular search and alignment task. The general name for methods like this is *seed-and-extend* methods.

BLANT (*Basic local Aligment Network Tool*) is intended for form the basis of a tool similar to BLAST, but for networks: given an undirected network *G* as a list of edges, and a value of *k*, it samples *k*-node subgraphs called *k-graphlets*. Since the number of *k*-graphlets in a graph of *n* nodes is exponential in both *k* and *n*, BLANT does not exhaustively enumerate all *k*-graphlets, but instead samples as many as the user specifies. Furthermore, since uniform random sampling of *k*-graphlets is difficult, there are several choices among sampling methods, each with different trade-offs. Finally, BLANT allows for several different methods of output: it can produce *orbit-degree vectors* (ODVs) like *ORCA* (citation), or it can produce an explicit list of *k*-graphlets that can be used as seeds for later extension. At present, BLANT does not provide an "extend" functionality; there are *many* seed-and-extend local alignment algorithms in the literature, and although BLANT currently is by far the fastest method of producing a large number of seeds, we have not yet implemented or tested any extend algorithms; this is an area of future work. Despite the lack of an "extend" feature, BLANT is still capable of useful bioinformatics, as described in our first tool paper in the journal *Bioinformatics* (citation).

## BLANT USAGE
### Command-line arguments
bla bla bla

# BLANT: brief description of internals and source code structure
## Directory libwayne
An extensive C-language library by Wayne Hayes, similar to a "class" library in C++ (but written mostly before C++ was invented), containing old and well-tested code implementing data structures and algorithms for graphs, sets, multisets, combinatorics, heaps, stacks, compressed flies, searching + sorting, linked lists, event-driven simulation and numerical analysis (the latter few not used much in BLANT but included for completeness). See [Wayne's Little Data Structures and Algorithms library](http://www.cs.toronto.edu/~wayne/libwayne).

## Directory networks
A directory containing many example input networks in edge list (2-colunmns, space or tab separated) format. While BLANT also accepts other input formats such as GML and LEDA, edge lists are the simplest and most recommended.

## Directory orca_jesse_blant_table
Many methods of automatically (and manually) identifying all the various graphlets exist. There is the original manually-generated numbering of graphlets of size 3, 4, and 5 nodes ([Przulj *et al* 2004](https://academic.oup.com/bioinformatics/article/20/18/3508/202438)) which is also used by the popular program [ORCA](https://academic.oup.com/bioinformatics/article/30/4/559/205331), and the automatic methods from [Jesse](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0147078), and our own ([Hasan, Chung, Hayes](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0181570)) which has been modified from using an upper- to lower-triangle representation for BLANT, following Jesse's lead. This directory contains code to translate between all these different naming conventions.

## Directory regression_tests
As the name applies, this directory contains regression tests, in the following format: each subdirectory contains code and data for one regression test. Any filename in such a directory whose name ends in ".sh" will be automatically run by the program regression-test-all.sh, which is run nightly by the senior author on a Jenkins machine. If you add functinoality to BLANT, you are encouraged to create a directory here with your own regression test, and then do a pull request to me to add your code to the main BLANT repo (it's best to send me an email first at whayes@uci.edu to discuss it first).

## Directory canon_maps.correct
Contains correct output for comparison to files in canon_maps to ensure your installation of BLANT is working correctly.

## Directory "Draw"
Contains code to create PDFs that depict any graphlet BLANT uses. Type "make Draw" to create the program Draw/graphette2dot. Then run ./Draw/graphette2dot for help and options. All graphette2dot does is create an input file and a command line for use by the "neato" command, which performs the actual creation of a PDF. So for example, to draw the 3-graphlet which is a triangle, run

./Draw/graphette2dot -k 3 -i 3 -o triangle | sh

That will create a file called "triangle.dot" and the command line for neato to convert triangle.dot to triangle.pdf; piping the output to sh simply runs the neato command, which actually generates the triangle.pdf file.
