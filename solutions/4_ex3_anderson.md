### Solution for the exercise "Anderson acceleration" in 4_Solving_the_SCF_problem.ipynb

In contrast to the damped iterations of exercise 2,
the iteration with Anderson acceleration converges for all repeats (2, 4, 6, 8)
with the damping value of 0.8 (which failed even for a single repeat in exercise 2!).

As the system size increases also the number of SCF iterations increases,
roughly as follows:

repeat | iterations
------ | -----------
   1   |  6
   2   |  8 
   4   | 11
   6   | 16
   8   | 29

The goal of a size-independent number of iterations has thus not yet been achieved.