![passing](https://img.shields.io/badge/build-passing-green)
![version](https://img.shields.io/badge/cpp-8.3.0-brightgreen.svg)
![make](https://img.shields.io/badge/make-4.2.1-brightgreen.svg)
![cmake](https://img.shields.io/badge/cmake-3.18.4-brightgreen.svg)
# temporal bipartite graph

## How to compile the code:
Run `Makefile`:
```shell script
make
```
Run `cmake`: 
```shell
cmake --build ./cmake-build-debug --target abcore -- -j 6

./cmake-build-debug/abcore -dataset /dir/example -query /dir/querystream.txt
```