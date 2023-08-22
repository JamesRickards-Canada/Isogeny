# Isogeny

Code to efficiently compute the supersingular isogeny graph G(p, l), written in PARI.

To use in Sage or Jupyter, follow these steps:
1. git clone this repository

2. Find where PARI/GP is installed in Sage (it's important to find this version, as it may differ in version number if you have it installed separately), i.e. the folder X where:  
   X/lib contains the file libpari.so;  
   X/lib/pari contains the file pari.cfg;  
   X/include/pari contains about 18 .h header files for PARI/GP.  
   On CoCalc, this is /ext/sage/default/local

3. In the repository, make a file named "pariloc.txt" containing X

4. Call "make" from the terminal to build the project

5. Open Sage/Jupyter, and call gp.read('isogeny.gp') to load the methods

6. g=gp.ssl_graph(101, 3) computes the graph for p=103 and l=2

7. ?gp.ssl_graph accesses the help for this method
