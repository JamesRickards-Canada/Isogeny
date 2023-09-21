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

2. Enter the folder and call "make". If your version of PARI/GP is not installed in the "default location", or you want to use this from Sage, see bullet point 1 of the "Sage incorporation" below.

3. Call "gp isogeny" to start gp and load the methods. ?supersingular accesses the help

## Sage incorporation
These are extra steps / changes required to use this inside Sage.

1. Find where PARI/GP is installed in Sage (it's important to find this version, as it may differ in version number if you have it installed separately), i.e. the folder X where:  
   X/lib contains the file libpari.so;  
   X/lib/pari contains the file pari.cfg;  
   X/include/pari contains about 18 .h header files for PARI/GP.  
   On CoCalc, this is /ext/sage/default/local

2. In the repository, make a file named "pariloc.txt" containing the path to X

3. Call "make" from the terminal to build the project

4. Open Sage/Jupyter, and call gp.read('isogeny.gp') to load the methods

5. To use the methods, the syntax is "g = gp.ssl_graph(101, [2, 3, 11])" (which computes the graph for p=103 and l=[2, 3, 11]). Note that this returns a PARI/GP object, which you may have to modify a bit to use with other sage methods

7. ?gp.supersingular accesses the help
