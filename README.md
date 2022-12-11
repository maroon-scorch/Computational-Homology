# Computational-Homology

## How to Run

To compute the homology of 3-dimensional Lens Spaces:
```
python main.py lens <number of tetrahedrons>
```
To compute the homology on a custom Delta Complex:
```
python main.py <path to file>
```

## Input Specification

We will explain how the input files should be formatted here, note that you can find some examples under the [example folder](example/).<br/>

The program takes a file of the following format:
- The first line of the file should be an non-negative integer indicating the number of vertices in the Delta Complex:
```
n
```
- The proceeding lines are as follows - The $i$-th line of file should contain a sequence of $i$-dimensional simplices, each separated by a semi-colon **;**: 
```
<simplex 1>; <simplex 2>; ...; <simplex k>
```
- For each simplex on line $i$, it should be represented by a list of indices, where the $j$-th index $v$ of the simplex represents the $j$-th face of the simplex. The index value $v$ means this face is the $v$-th simplex of line $i-1$:
```
<index_1> ... <index i>
```
For example, the following is the standard Delta Complex specification of the 2-torus:
```
1
0 0; 0 0; 0 0
0 1 2; 2 1 0
```
This structure only has 1 vertex, so all of its 3 edges start and end at the same point. Its two faces go in opposite orientations, hence their edges go in opposite directions.
