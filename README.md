# TriFMatch

## Compile
Under the root directory of the project, execute the following commands to compile the source code.

```zsh
mkdir build
cd build
cmake ..
make
```


## Execute
A running example. The time limit is 300 seconds in default;

```zsh
cd build/src
./main.o --data ../../dataset/hprd.graph --query ../../dataset/hprd_test_query.graph --num 100000
```

where:
- `--data`, the path of the data graph file
- `--query`, the path of the query graph file
- `--num`, number of results intended to find

## Input
Both the input query graph and data graph are vertex-labeled and in the same format.
The first line of each graph is formatted as 't N M' where N is the number of vertices and M is the number of edges
Then each vertex and edge is represented as 'v VertexID LabelID Degree' and 'e VertexID VertexID', respectively. The vertex id is started
from 0. Here is a input sample.

Example:

```zsh
t 5 6
v 0 0 2
v 1 1 3
v 2 2 3
v 3 1 2
v 4 2 2
e 0 1
e 0 2
e 1 2
e 1 3
e 2 4
e 3 4
```

## Configuration
You can configure the enumeration process by adjusting the following macros in 'configuration/config.h'

| Macro | Description |
| :-----------------------------------------------: | :-------------: |
|TIME_LIMIT| set the time limit for the enumeration process|
|FAILING_SET_PRUNING| enable the filtering technique FAILING SET PRUNING  |
|PROACTIVE_CANDIDATE_COMPUTING| enable the proactive candidate computing if disabled the whole algorithm is |
|CD_FILTERING | enable the containment-driven filtering |
|FD_FILTERING | enable the failure-driven filtering |
|PRINT_RESULT| print the results in details, the i-th element printed is the data vertex matching the i-th vertex of the order | 
|PRINT_MEM_INFO| print the peak memory |
|PRINT_LEAF_STATE | print the number of states being filtered |
