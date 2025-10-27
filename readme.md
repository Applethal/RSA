# Revised Simplex Algorithm in C
STILL WORK IN PROGRESS!!!

TODO: 
1- Add error handling
2- Add infeasibility check while solving
3- Implement proper memory de-allocation

# How to contribute 

Test this with as many LP models please. So far I gave it nothing but LPs from linear algebra text books. Be sure to report any issues if necessary.

In this code I attempted to write a mathematical solver that utilizes the Revised Simplex Algorithm by George B. Dantzig (1953). Model inputs should be written as such:

```
OBJECTIVE Coeffs
Constraints
```
Where OBJECTIVE refers to the objective function and takes either the `M̀INIMIZE` or `M̀AXIMIZE` tokens and next to it, put each variable coefficient together with a space in between each coefficient (There are no sign restrictions here) all in one line, the next lines will be strictly for the constraints where each constraint will have its left hand sign contain nothing but variable coefficients preceding the constraint's symbol `<=`, `=` or `>=`.

Example: 

```
MINIMIZE 24 29 10 38;
1 4 5 0 = 60;
0 2 1 0 <= 12;
2 1 -1 4 >= 10;
```
$$
\begin{aligned}
\text{Minimize } & z = 24x_1 + 29x_2 + 10x_3 + 38x_4 \\
\text{subject to } 
& x_1 + 4x_2 + 5x_3 = 60 \\
& 2x_2 + x_3 \le 12 \\
& 2x_1 + x_2 - x_3 + 4x_4 \ge 10 \\
& x_1, x_2, x_3, x_4 \ge 0
\end{aligned}
$$

# Update 27-10-2025
As of today, the solver works but assumes that your text input follows the rules above. It utilizes' Gauss-Jordan elimination technique to invert the basics' matrix, if degeneracy is detected, the program terminates. So far it takes models that have proven optimality and has a check for unbounded LPs, though no infeasibility verification has been implemented as of yet.  
