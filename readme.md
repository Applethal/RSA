# Revised Simplex Algorithm in C
STILL WORK IN PROGRESS!!!

TODO: 

1- Add error handling.


2- Add infeasibility check while solving.

3- <del>Implement proper memory de-allocation.</del>

4- Might as well substitute .txt data input with .csv data instead. Maybe even work on better data parsing since it does more computational effort, in my opinion.
# How to contribute 

Test this with as many LP models please. So far I gave it nothing but LPs from linear algebra text books. Be sure to report any issues if necessary.

# About
In this code I attempted to write a mathematical solver that utilizes the Revised Simplex Algorithm by George B. Dantzig (1953). Model inputs should be written (In general form) as such:

```
Coeffs
Constraints
```
Where `Coeffs` refers to the objective function coefficients for the variables (There are no sign restrictions here) all in one line, the next lines will be strictly for the constraints where each constraint will have its left hand sign contain nothing but variable coefficients preceding the constraint's symbol `<=`, `=` or `>=`. Make sure each line ends with a `;`. 

Example: 

```
9 7;
10 5 <= 50;
6 6 <= 36;
4.5 18 <= 81;
```
$$
\begin{aligned}
\text{Minimize } & z = 9x_1 + 7x_2 \\
\text{subject to } 
& 10x_1 + 5x_2  \leq 50 \\
& 6x_2 + 6x_3 \leq 36 \\
& 4.5x_1 + 18x_2 \leq 81 \\
& x_1, x_2 \ge 0
\end{aligned}
$$

Implicitly, all variables are non-negative, of course, you won't need to consider this. To run this program, simply pass the text file and objective arguments:

```./RSA "model.txt" "MINIMIZE"```


# Update 27-10-2025
As of today, the solver works but assumes that your text input follows the rules above. It utilizes' Gauss-Jordan elimination technique to invert the basics' matrix, if degeneracy is detected, the program terminates. So far it takes models that have proven optimality and has a check for unbounded LPs, though no infeasibility verification has been implemented as of yet.  
