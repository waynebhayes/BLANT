*** 22 October 2023 (Wayne Hayes)
*** This code is from https://bitbucket.org/sjain12/cliquecounting/src/master. I modified it slightly to more easily compile
*** on multiple systems, then added it to the BLANT repository (github.com/waynebhayes/BLANT).

************************************************************************************************
*** original README.md from https://bitbucket.org/sjain12/cliquecounting/src/master is below ***
************************************************************************************************
UPDATE 1/22/2018
This code has been updated with a fix for a minor bug that made the algorithm cliques_turan_shadow give biased, erroneous estimates in certain cases. For questions, please contact the author.

=================================

# CliqueCounting
Code for the algorithms and experiments in "A Fast and Provable Method for Estimating Clique Counts Using Turán's Theorem" by S. Jain and C. Seshadhri

Steps to run the code:

* Store the graph to be analyzed as a list of edges in a .txt file in the graphs folder. For example, you could download and uncompress one of the graphs from the SNAP library, say amazon0601. cd into the python folder and run:

  python sanitize.py ../graphs amazon0601.txt
  
  This will create a file 'amazon0601.edges' in the graphs folder.
  
* To run the TuránShadow algorithm, run:
  
  ./test_cliques_turan_shadow -i amazon0601.edges -m 5 -M 5 -r 10 -n 50000
  
  This will save the tree stats in '/amazon0601_5_50000_params' csv file in results/cliques/ and the estimates of the 10 runs in 'amazon0601_5_50000_10_data'
  
* To run the Edge Sampling algorithm, run:

  ./test_cliques_edge_sampling -i amazon0601.edges -m 5 -M 5 -r 10 -p 0.75
  
  This will save the estimates from 10 runs of the algorithm in amazon0601_5_0.75_10_es in results/cliques/
  
* To run the GRAFT algorithm, run:

  ./test_cliques_graft -i amazon0601.edges -m 5 -M 5 -r 10 -n 100000
  
  This will save the estimates from 10 runs of the algorithm in amazon0601_5_100000_10_graft in results/cliques/
  
* To run the Color Coding algorithm, run:

  ./test_cliques_color_coding -i amazon0601.edges -m 5 -M 5 -r 10 
  
  This will save the estimates from 10 runs of the algorithm in amazon0601_5_10_cc in results/cliques/
  
* To run the brute force clique counting algorithm, run:

  ./test_cliques_brute_force -i amazon0601.edges -m 5 -M 5 -r 10 
  
  This will save the estimates from 10 runs of the algorithm in amazon0601_5_10_bf in results/cliques/

For questions, contact Shweta Jain at sjain12@ucsc.edu.
