# CR2 (Simple)

**Description**

* This is a simple and basic version of a CR2 robust variance estimator based on Pustejovsky & Tipton's implementation in R/Stata.
* Highly recommended for users to use the CR2 library instead of this specific implementation unless a very specific requirement.

**Using CR2 (Simple)**

* This version of the CR2 estimator does not use a modified adjustment matrice (i.e. Theorem 2; Pustejovsky & Tipton, 2018), for fringe cases where the theorem might not hold.
* For this reason, this is highly unoptimized and only used for testing purposes.

**Prerequisites**

* Python
* Standard python libraries (e.g. Pandas)
