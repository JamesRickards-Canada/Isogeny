# Isogeny

Code to efficiently compute the directed supersingular isogeny graph G(p, L), written in PARI. This allows for L = integer or vector of integers, currently restricted to at most 17. For larger l, it requires downloading and processing the relevant modular polynomials (see "modpol_processraw" for more details).

There is currently code built-in to Sage to do this, however this project allows for non-prime l and is significantly faster than the Sage implementation:
* l = 2, about 67 times as fast
* l = 3, about 107 times as fast
* l = 5, about 211 times as fast
* l = 7, about 217 times as fast

## Main methods
gp: 
* ssl_graph: computes the supersingular isogeny graph. Does both the l and L isogeny graphs.
* ssl_graphadjmat: returns the adjacency matrix of said graph
* ssl_graph_scipy: does not return the graph but saves it to a file. This file can be used to import this data into python, and efficiently compute its eigenvalues.

python:
* sparseadjmtx: load the saved file (from ssl_graph_scipy) as a sparse array in scipy
* secondeval: returns the second largest eigenvalue
* secondabseval: returns the second largest eigenvalue in absolute value
* allevals: returns all eigenvalues

## Installation
1. git clone this repository

2. You need to know where the version of PARI/GP you want to use is installed. The default location is inside /usr/local, but this may change based on your Linux distro, or if you want to use it through SageMath. If you think it is installed in the default location, you can simply call "make".

3. Otherwise, call "make setup" to search for the correct files. They are likely somewhere in /usr, but you can change this if not. If the program finds potential matches, it will ask you to confirm which files are correct, and saves them to "pari_loc.txt". Once this step completes, a call to "make" will compile the project! Modifying the program (e.g. via git pull) won't require redoing this setup, unless the version of PARI/GP or Sage you use changes.

4. Call "gp isogeny" to start gp and load the methods. ?supersingular accesses the help

5. Call "make clean" to clean up the files created.

## Sage integration
You need to make sure that when you call "make setup", you find the version of PARI/GP installed with the version of Sage you are using. For example, on CoCalc, it should be in /ext/sage/VERSION/local. Once the project is built:

1. Open Sage/Jupyter, and call gp.read('isogeny.gp') to load the methods

2. To use the methods, the syntax is "g = gp.ssl_graph(101, [2, 3, 11])" (which computes the graph for p=103 and l=[2, 3, 11]). Note that this returns a PARI/GP object, which you may have to modify a bit to use with other sage methods

3. ?gp.supersingular accesses the help
