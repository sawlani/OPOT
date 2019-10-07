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
Data is in plaintext format.
+ Line 1: n and m, the # of rows and columns respectively
+ Line 2: n numbers, demands (in integers) of the left side
+ Line 3: m numbers, demands (in integers) of the right side (total of line 2 equals total of line 3)
+ Line 4 - n + 3: m numbers per line, the costs of transporting something from i to j
A format checker for this format is in Data/Verify.cpp

We are working on separate converters
from these data into other formats (e.g. .mat).
However, we believe it's best to require the rows/columns
of data to be scrambled randomly.

## Summary of Results

In progress.
Latest version in OTML_19_paper.pdf.

At the moment the only evaluation criteria is
ratio to optimum ground truth cost value, vs. time.
We plan to also include k-nearest-neighbor classification
accuracy as an alternate measure.

## More Methods Compared
OTML_19_paper.pdf also compared with
+ https://github.com/JasonAltschuler/OptimalTransportNIPS17 
+ https://github.com/nathaniellahn/CombinatorialOptimalTransport

Due to MATLAB licensing issues, we are in the process of figuring out how run these on the same machine.
