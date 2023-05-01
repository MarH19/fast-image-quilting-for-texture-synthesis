# Advanced Systems Lab -- Group 24

## Build Information

To build the tests run in a linux shell the following command.
```shell
$ make tests
$ ./bin/tests/test_l2norm.out  # a small test script
$ ./bin/tests/test.out         # play-function
$ make cleanall
```
To build the benchmark run in a linux shell the following command.
```shell
$ make benchmark
$ ./bin/benchmark/test.out  # runs benchmark tests
$ make cleanall
```


This compiles the code only if changes occur and does basically all the linking stuff on its own.
The manual compilation is of course still possible.

## Useful links

- project description: [link](https://acl.inf.ethz.ch/teaching/fastcode/2023/project/project-ideas/Image-Quilting.pdf)
- matlab implementation: [link](https://jmecom.github.io/projects/computational-photography/texture-synthesis/)
- original paper: [link](http://graphics.cs.cmu.edu/people/efros/research/quilting/quilting.pdf)
- lecture website: [link](https://acl.inf.ethz.ch/teaching/fastcode/2023/)