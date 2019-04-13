## This project implemented the simple random walk algorithms for barrier options pricing problems. The reference could be found here : https://arxiv.org/pdf/1211.5726.pdf

### We have implemented efficient pricing algorithms for Barrier cap/floor and BarrierSwaption

### Three types of algorithms can be chose :
    * mode = 1, the random walk algorithm 1 introduced in the paper, with weal error rate 1
    * mode = 2, the random walk algorithm 2 introduced in the paper, with weak error rate 1/2 
    * mode = 3, traditional monte calro method

### To run the simulation for barrier cap/floor, with 10000 times and algo type 1:
```
    make barrierCap
    ./barrierCap 10000 1
```

### To run the simulation for barrier swaption, with 10000 times and algo type 1:
```
    make barrierSwap
    ./barrierSwap 10000 1
```

### Directories/files description
    * Makefile : make rules, users should adapt the include path
    * Eigen/ : one C++ linear algebra library, used here for cholesky decomposition
    * optim/ : one light optimization library, used here for solving optimization problem   
                 when pricing barrier swaption. Notice armadillo should be installed to 
                 use this tool.
    * util.hpp/.cpp : implemented several utility functions like create files, cholesky     
                 decomposition, generate correlation matrix for multi libor case, etc.
    * process.hpp : implemented template classes / methods related to the stochastic
                 processes which are defined by euler schema. 
    * liborrates.hpp/.cpp : implemented Libor rate class which contains the update rule,
                 the reset rule, etc.

    * barrierCap.hpp/.cpp : implemented Barrier cap/floor class, which contains the 
                 implementation of three algos (including simple monte carlo), intrinsic value calculation, the true value calculation, etc.
    * barrierCapSimu.cpp : main file for running the simulations of Barrier cap/floor. We 
                 could use it for pricing barrier call or put, caplet or floor. Notice that 
                 we store the mean value and variance every 5000 simulations.

    * barrierSwaption.hpp/.cpp : implemented Barrier swaption class, which contains the 
                 simulation of different algos, intrinsic value calculation, apporimate va-
                 lue calculation, etc.
    * barrierSwapSimu.cpp : main file for running the simulations of Barrier cap/floor 
                 swaption, of call type or put type. Notice that we store the monte carlo result every 1000 simulations. 

    * output/ : contains the result of experiments.
                    