The original network is HI-union.el, in this directory. You probably don't need the details, but if interested see http://www.interactome-atlas.org/

To predict edges, we basically sample small subgraphs called "graphlets". If we observe that a particular graphlet usually appears appears complete, but occasionally appears with a certain edge missing, then we can predict that that edge should actually exist. The URL above describes the case of when an L3 (path of length three) can be "closed" to form a square--the missing edge tha connects the ends of the L3 path; graphlets are a generalization of this idea.

There are several version of the data, depending on how big the graphlets are that we sample; the variable "k" is used to denote the size of sampled subgraphs, so that k=4 means 4-node subgraph samples, k=8 is 8-node, etc. For example the k45 dataset contains graphlets of size 4 and 5; k456 contains graphlets up to size 6, etc. I have found that larger graphlets seem to produce better predictions, but the catch is that it takes exponentially longer to sample larger graphlets, and also the number of features grows exponentially with k. In the largest dataset there are about 200,000 possible features (ie., 200,000 possible columns where one column is one feature).

Each row of the data contains one *pair of nodes*, say (u,v). Since the original network has about 9,000 nodes, there are (9000 choose 2) possible pairs--about 40 million pairs. Most of the folds have only 10-20 million lines, since the rest of pairs have counts of zero in all columns. Also, to avoid requiring 200,000 columns on every line, the "columns" are separated by tabs, and within one tab-delimited column are two space-separated strings: a feature ID of the form k:o:p where k is as above, and o and p are "orbit" identifiers (again details probably not relevant); and a value which is the feature value. Some lines (ie., node pairs (u,v)) have have only one nonzero feature; other node pairs have dozens of nonzero features.

The very first column contains the node pair as u:v, followed by a space, followed by a boolean (0 or 1) defining whether that pair of nodes has a known edge. So for example one line looks like this:

ENSG00000150783:ENSG00000068976 0       4:11:11 0.173909

so ENSG00000150783 is one node, ENSG00000068976 is the other, there is no edge (0); then there's a tab; there is one nonzero feature named 4:11:11, with value 0.173909.  A longer line might look like this:

ENSG00000182173:ENSG00000154743 1       5:40:40 0.192355        4:11:11 0.511344        5:60:59 0.021941

so this pair has an edge, and three non-zero features. Note that the feature IDs are *not* sorted, so feature ID 4:11:11 appears in the middle of these three nonzero features.

There are 10 folds, numbered 0 through 9.  Within each fold, 10% of the edges have been removed from HI-union.el and placed into a file called test HI-union-testN.el (N from 0 through 9); the remaining edges are the "training" edge list, called HI-union-trainN.el. Note that YOU CANNOT CHANGE THE FOLDS OR CREATE YOUR OWN FOLDS. This is because the feature values must be computed from scratch for every fold independently. Furthermore, the edges that you want to predict *are already in the feature file*, with a (false) zero in the boolean (edge) column. Thus, technically, this problem can be interpreted as a data "cleaning" problem rather than prediction. This is because *every pair of nodes*, whether they have an edge or not, must be considered when creating the feature vectors. Thus, all edges can be considered reliable, but every non-edge is considered suspect; your job is to "detect" (ie., predict) which non-edges are falsely labeled.

Note that the conglomerate of all the 10 test sets equal the original network.

There are two versions of the data; a non-normalized version, and a normalized one. In my own tests, I have found that the normalized set gives slightly better predictions, but I've provided the raw (*.out) files as well since I can't be sure that my normalization is the best possible.

The file L3.predictions.top1000.el contains the top 1,000 edge predictions (in order of L3 count, not normalized in any way) from the raw L3 method discussed at the above URL of the original data.

More formally, the columns of the training file are of the following form: tab-separated "major" columns, each with two sub-columns space-separated. That is,

    Node1:Node2 [space] B [tab] MotifPairID [space] count [tab] MotifPairID [space] count [tab] .......

where:
    Node1:Node2 is the two protein names, separated by a colon; they are listed in REVERSE "alphabetical" order, ie., Node1 is always lexicographically larger than Node2.
    B is a boolean digit: 0 meaning the nodes are not connected, 1 if they form an edge

    MotifPairID: the name of the orbit pair. For historical reasons there are two possible naming conventions. The files in an older version of the data used the following naming scheme: k:c:q:r, where k is the number of nodes in the motif (4 through 8 inclusive); c is the "canonical ordinal ID" of the graphlet/motif (which is essentially the "line number" of the orbit_map[k].txt file where line 0 is the line that's all zeros in that file); and q and r are the *columns* on that line, numbered from zero, and q>r.
    
    The other convention, used in the current files, is "k:o:p", where k is the number of nodes in the motif (4 through 8 inclusive), and o and p are orbits using the orbit IDs from the file canon_maps/orbit_list[k].txt from a build of BLANT. 

    count is the count of how frequently orbits o and p connect Node1 with Node2 (like the L3 path count).

There are as many (MotifPair count) columns as were observed by BLANT while processing that fold's training network.

Across all 10 folds, there are about 200,000 observed orbit pairs, so your ML algorithm will need to handle that many features.

BLANT is available (but not needed for the ML portion of this project) at https://github.com/waynebhayes/BLANT

