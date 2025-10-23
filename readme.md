# Revised Simplex Algorithm in C
STILL WORK IN PROGRESS!!!
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


