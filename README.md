## INFO
This is the repository of the project from the Advanced Systems Lab course offered at ETH Zurich.

## Build Information

To build the tests run in a linux shell the following command.
```shell
$ make run_tests               # runs + builds all tests
$ make tests                   # builds all tests
$ ./bin/tests/test_l2norm.out  # runs a small test script
$ ./bin/tests/test.out         # runs play-function
$ make cleanall                # deletes all created files
```
To build the benchmark run in a linux shell the following command.
```shell
$ make benchmarks
$ ./bin/benchmarks/test.out  # runs benchmark tests
$ make cleanall
```
`make benchmarks` and `make tests` use different compiler flags to compile their binaries.
But in Makefiles object files are **not** recompiled when some flags change. Therefore it is important to *always* run `make clean` or `make cleanall` between `make tests` and `make benchmarks`.



This compiles the code only if changes occur and does basically all the linking stuff on its own.
The manual compilation is of course still possible.

## Useful links

- project description: [link](https://acl.inf.ethz.ch/teaching/fastcode/2023/project/project-ideas/Image-Quilting.pdf)
- matlab implementation: [link](https://jmecom.github.io/projects/computational-photography/texture-synthesis/)
- original paper: [link](http://graphics.cs.cmu.edu/people/efros/research/quilting/quilting.pdf)
- lecture website: [link](https://acl.inf.ethz.ch/teaching/fastcode/2023/)
