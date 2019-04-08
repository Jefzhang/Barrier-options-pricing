## This project implemented the simple random walk algorithms for barrier options pricing problems. The reference could be found here : https://arxiv.org/pdf/1211.5726.pdf

### We have implemented efficient pricing algorithms for Barrier cap/floor and BarrierSwaption

### Three types of algorithms can be chose :
    * mode = 1, the random walk algorithm 1 introduced in the paper, with weal error rate 1
    * mode = 2, the random walk algorithm 2 introduced in the paper, with weak error rate 1/2 
    * mode = 3, traditional monte calro method

### To run the simulation for barrier cap/floor :
'''cpp
    //run the simulation for 10000 times with algo type 1
    make barrierCap N=10000 mode=1
'''

### To run the simulation for barrier swaption :
'''cpp
    //run the simulation for 10000 times with algo type 1
    make barrierSwap N=10000 mode=1

