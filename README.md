# Advanced Systems Lab -- Group 24

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
$ ./bin/benchmark/test.out  # runs benchmark tests
$ make cleanall
```

It is important to note that benchmarks and tests use different flags to compile.
But if you use `make tests` then the object files get **not** recompiled if you use `make benchmarks`.
Therefore if you switch between those to modes *always* run `make clean` or `make cleanall` beforehand. 


This compiles the code only if changes occur and does basically all the linking stuff on its own.
The manual compilation is of course still possible.

## Useful links

- project description: [link](https://acl.inf.ethz.ch/teaching/fastcode/2023/project/project-ideas/Image-Quilting.pdf)
- matlab implementation: [link](https://jmecom.github.io/projects/computational-photography/texture-synthesis/)
- original paper: [link](http://graphics.cs.cmu.edu/people/efros/research/quilting/quilting.pdf)
- lecture website: [link](https://acl.inf.ethz.ch/teaching/fastcode/2023/)