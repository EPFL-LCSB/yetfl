yETFL

- cobra (https://opencobra.github.io/cobrapy/)
- pytfa (https://github.com/EPFL-LCSB/pytfa)

The code can be run on different types of OS including UNIX and Windows.

You will need Git LFS to properly download some binary files:

git lfs install
git lfs pull

We recommend the use of Docker to set up a container instead of running the codes directly. To solve the genome-scale models, a licensed commercial solver, either Gurobi or CPLEX, should be installed. The free academic versions can be obtained through Gurobi for Acedemics and Researchers (https://www.gurobi.com/academia/academic-program-and-licenses/) or IBM Academic initiative (https://developer.ibm.com/academic/). Alternatively, to test the general performance of the code, the scripts in tests/ subdirectory can be run with the free solver GLPK. For further details, see README.rst.

To install the package the setup.py script should be run. The installation process will be fast (It takes <10 seconds).
