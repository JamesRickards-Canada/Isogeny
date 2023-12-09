# Isogeny

Code to efficiently compute the directed supersingular isogeny graph G(p, L), written in PARI. This allows for L = integer or vector of integers, currently restricted to at most 17. For larger l, it requires downloading and processing the relevant modular polynomials (see "modpol_processraw" for more details).

There is currently code built-in to Sage to do this, however this project allows for non-prime l as well as L-isogeny graphs, and it is significantly faster than the Sage implementation:
* l = 2, about 67 times as fast
* l = 3, about 107 times as fast
* l = 5, about 211 times as fast
* l = 7, about 217 times as fast

## Main methods
gp: 
* ```ssl_graph```: computes the supersingular isogeny graph. Does both the l and L isogeny graphs.
* ```ssl_graphadjmat```: returns the adjacency matrix of said graph
* ```ssl_graph_scipy```: does not return the graph but saves it to a file. This file can be used to import this data into python, and efficiently compute its eigenvalues.

Python:
* ```sparseadjmtx```: load the saved file (from ```ssl_graph_scipy```) as a sparse array in scipy
* ```secondeval```: returns the second largest eigenvalue
* ```secondabseval```: returns the second largest eigenvalue in absolute value
* ```allevals```: returns all eigenvalues

Sage:
* ```ssl_adjmat```: computes the supersingular isogeny graph using the PARI/GP code, and makes the adjacency matrix as a sparse matrix in Sage. This is about twice as fast as ```ssl_graph```, as we don't need to make a DiGraph and add the edges.
* ```ssl_graph```: computes the supersingular isogeny graph using the PARI/GP code, and converts it to a Sage object. The output is a DiGraph with vertices labelled by supersingular j-invariants. This is the replacement for ```E.isogeny_ell_graph(l, directed = True, label_by_j = True)```: it will give an isomorphic graph, where the F_p^2 vertices may be presented differently (the choice of generator for F_p^2 may be different). Note that the vertices are labelled by the actual j-invariants in GF(p^2), and not the strings that give the j-invariants (which is what isogeny_ell_graph does).

## Installation
1. ```git clone``` this repository, and enter the folder created.

2. You need to know where the version of PARI/GP you want to use is installed. The default location is inside /usr/local, but this may change based on your Linux distro, or if you want to use it through SageMath. If you think it is installed in the default location, you can simply call ```make```.

3. Otherwise, call ```make setup``` to search for the correct files. By default the program searches in ```/usr```, but there is a chance it is not installed there (this sometimes happens on a server). If this is the case, you can supply an alternate location.

4. If the program finds potential matches, it will ask you to confirm which files are correct, and saves them to "pari_loc.txt". Once this step completes, a call to ```make``` will compile the project! Modifying the program (e.g. via ```git pull```) won't require redoing this setup, unless the version of PARI/GP or Sage you use changes.

5. Call ```gp isogeny``` to start gp and load the methods. ```?supersingular``` accesses the help

6. Call ```make clean``` to clean up the object files (.o) created.

## Troubleshooting installation

* **No library found**: the files were not found in the search location. Try asking to search in "/", which searches everywhere. This will be slow, but is guaranteed to find the correct files, if they exist.
* **Wrong version**: Maybe you found the libraries, but there were warnings with ```make```. The likely cause is the version of PARI/GP you found was too old. If there are multiple copies of PARI/GP on your computer, then perhaps you chose the wrong one! Check the shared object file created: it will be called "libisogeny-X-Y.so", where "X.Y" is the version of PARI/GP in the libraries. If this does not match the version you are using, then you found the wrong one!
* **Miscellaneous**: when you compile with ```make```, object files (.o) are created. However, if the underlying code did not change, then nothing will happen. If you change the version of PARI/GP you are working with (or perhaps, trying different installations), then you should call ```make clean``` in between to clear out these object files. Otherwise, recompiling will do nothing!
* If you are still having issues with installation, please get in touch, and I will try to help sort it out!
 
## Sage integration
You need to make sure that when you call ```make setup``` you find the version of PARI/GP installed with the version of Sage you are using. For example, on CoCalc, it should be in /ext/sage/VERSION/local. Once the project is built, you can open Sage and call:
```
load("ssl_pari.sage")
G = ssl_graph(101, 8)
M = ssl_adjmat(15013, [2, 3])
```
to create the supersingular 8-isogeny graph for p = 101 as a DiGraph, and the adjacency matrix for the [2, 3]-isogeny graph for p=15013. Note that this is a bit slower than the native PARI/GP methods, as there is overhead in converting them to Sage, as well as making the DiGraph object. If you want to play with the native PARI/GP methods in Sage, then:

1. Open Sage/Jupyter, and call ```gp.read('isogeny.gp')``` to load the methods

2. To use the methods, the syntax is ```g = gp.ssl_graph(101, [2, 3, 11])``` (which computes the graph for p=103 and l=[2, 3, 11]). Note that this returns a PARI/GP object, which you may have to modify a bit to use with other sage methods.

3. ```?gp.supersingular``` accesses the help

However, I suggest simply using PARI/GP instead if this is the case. Every call to ```gp.METHOD``` has a significant overhead (perhaps as much as 5ms).

## Working in a different directory
The easiest way to use this package in a different directory (to the installation location) is to create soft links to the relevant files. Let D be the directory containing the methods, and navigate to the folder where you want to run them from. Call:
```
ln -s D/isogeny.gp
ln -s D/libisogeny-X-Y.so
ln -s D/modpolys/
```
where X.Y is the version of PARI/GP you are using; it's in the file name created with ```make```. If you want to use it with sage, also link ```ln -s D/ssl_pari.sage```.

Now the program will work exactly the same as if you were running it from the installation directory.
