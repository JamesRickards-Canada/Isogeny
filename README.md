# Isogeny

Code to efficiently compute the supersingular isogeny graph G(p, l), written in PARI.

To use in Sage or Jupyter, follow these steps:
1. Find where PARI/GP is installed, i.e. the folder X where:  
   X/lib contains the file libpari.so;  
   X/lib/pari contains the file pari.cfg;  
   X/include/pari contains about 18 .h header files for PARI/GP.  
   On CoCalc, this is /ext/sage/default/local

3. Make a file named "pariloc.txt" containing X

4. Call make to build the project

5. From Sage/Jupyter, call gp.read('isogeny.gp')

6. g=gp.ssl_graph(101, 3) makes a graph

7. ?gp.ssl_graph accesses the help for this method
