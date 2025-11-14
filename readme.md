# Revised Simplex Algorithm in C
STILL WORK IN PROGRESS!!!

TODO: 




<del>1- Add error handling.</del>


2- <del>Add infeasibility check while solving.</del>

3- <del>Implement proper memory de-allocation.</del>

4- <del>Might as well substitute .txt data input with .csv data instead. Maybe even work on better data parsing since it does more computational effort, in my opinion. </del>

5- <del> Add debug mode </del>.

6- Using Gauss' pivoting technique for matrix inversion is an effort of  $O(n^3)$, can I do better? 

# How to contribute 

Test this with as many LP models please. So far I gave it nothing but LPs from linear algebra/Operation research text books and internet forums. Be sure to report any issues if necessary.

# About
In this code I attempted to write a mathematical solver that utilizes the Revised Simplex Algorithm by George B. Dantzig (1953) using the Big-M method. The model input should be written (In general form) as such:

```
OBJECTIVE
Coeffs
Constraints
```
Where `OBJECTIVE` is the objective function's direction which takes either `MINIMIZE` and `MAXIMIZE` as keywords. `Coeffs` refers to the objective function coefficients for the variables all in one line, the next lines will be strictly for the constraints where each constraint will have its left hand sign contain nothing but variable coefficients preceding the constraint's symbol `<=`, `=` or `>=`. Make sure each entry is separated with a `,`. 
One last thing: Make sure the right-hand side is positive and that there is no empty line beneath the last constraint. Use 0s in the LHS of each constraints if a variable is irrelevant for the latter. This is because the program checks if each LHS has exactly $n$ number of variables, else the solver won't start.

Example 1: 

```
MAXIMIZE
9,7
10, 5, <=, 50
6, 6, <=, 36
4.5, 18, <=, 81
```
Which corresponds to:

$$
\begin{aligned}
\text{MAXIMIZE } & z = 9x_1 + 7x_2 \\
\text{subject to } 
& 10x_1 + 5x_2  \leq 50 \\
& 6x_2 + 6x_3 \leq 36 \\
& 4.5x_1 + 18x_2 \leq 81 \\
& x_1, x_2 \ge 0
\end{aligned}
$$


```
MINIMIZE
-3,1,1
1,-2,1,<=,11
-4,1,2,>=,3
2,0,1,=,1
```
Which corresponds to:

$$
\begin{aligned}
\text{MINIMIZE } \quad 
& z = -3x_1 + x_2 + x_3 \\
\text{subject to } \quad
& x_1 - 2x_2 + x_3 \le 11 \\
& -4x_1 + x_2 + 2x_3 \ge 3 \\
& 2x_1 + x_3 = 1 \\
& x_1, x_2, x_3 \ge 0
\end{aligned}
$$

Implicitly, all variables are non-negative (of course) you won't need to consider this. To run this program, simply pass the text file and objective arguments:

```./RSA "filepath" "-Debug"```


Where `filepath`is the `.csv` file path, `-Debug` is an optional flag that can be added as an argument, allowing you to see the solver operations step by step, the displayed indices are 0 indexed. The program will convert the problem to its canonical form then iteratively execute the algorithm until it terminates. The reason I am using the .csv file format is because of how portable it is + you can easily view whether the entries are valid using a Graphical reader to easily display whether data is missing.



