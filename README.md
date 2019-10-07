# OPOT
Data files for comparing optimal transport algorithms

At the moment it's a minimalistic data set for comparing
dense optimal transport in single processor settings.
We plan to support Octave (open source MATLAB)
and GPUs in the near future.

We also plan to eventually expand the data set to include,
in order of priority:
+ high dimensional point data
+ low dimensional point data
+ sparse graphs

## Data format:
+ Line 1: n and m, the # of rows and columns respectively
+ Line 2: n numbers, demands (in integers) of the left side
+ Line 3: m numbers, demands (in integers) of the right side (total of line 2 equals total of line 3)
+ Line 4 - n + 3: m numbers per line, the costs of transporting something from i to j

## Summary of Results

In progress.
Latest version in OTML_19_paper.pdf.

## More Methods Compared
OTML_19_paper.pdf also compared again
+ https://github.com/JasonAltschuler/OptimalTransportNIPS17 
+ https://github.com/nathaniellahn/CombinatorialOptimalTransport
