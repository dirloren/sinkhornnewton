# sinkhornnewton

This repository contains the MATLAB code that generates the figures in the paper "A Sinkhorn-Newton method for entropic optimal transport" available at <https://arxiv.org/abs/1710.06635>.

##### Authors:
- Christoph Brauer    (<ch.brauer@tu-braunschweig.de>)
- Dirk Lorenz    (<d.lorenz@tu-braunschweig.de>)
- Christian Clason    (<christian.clason@uni-due.de>)
- Benedikt Wirth (<benedikt.wirth@uni-muenster.de>

Contents
--------

##### Drivers (run these to generate figures):
    driver_sinkhorn.m        test script to generate figure 1
    driver_parameters.m      test script to generate figure 2
    driver_discretization.m  test script to generate figure 3

##### Routines called by the drivers:
    sinkhorn.m                implements the Sinkhorn-Knopp method
    sinkhorn_newton.m         implements the proposed Sinkhorn-Newton method
    sinkhorn_newton_primal.m  the same using only primal variables
    MNISTGroundMetric.m       generates a cost-metric used in driver_parameters
    loadMNISTImages.m         loads the images from the MNIST database

##### Dependencies
`driver_parameters.m` needs images from the MNIST database which can be obtained from <http://yann.lecun.com/exdb/mnist/>
