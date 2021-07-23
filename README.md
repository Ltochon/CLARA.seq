[![Licence: MIT](https://img.shields.io/badge/Licence-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
# CLARA.seq
*R package for sequences clustering*
<br>
This package works with the *TraMineR* library which provides the *LCS* method for distance calculations.

## Installation
To use this package on R, this following lines are required :
```R
install.packages("devtools")
library(devtools)
install_github("Ltochon/CLARA.seq")
```

If the vignette are not refreshed and you can't see the CLARA.seq's one, use this line :

```R
devtools::install(build_vignettes = TRUE)
```

## Algorithms
This package contains 3 differents algorithms
- CLARA
- CLARANS
- CLARA-FUZZY

And a quality index 
- Davies-Bouldin Index

## Purpose of this package
*TraMineR's* package has a limitation of the number of sequences used (~46'300). To counter this limit, different algorithm has been implemented to use subsets of the entire dataset to extract the best clustering for the big dataset.
At this moment, a dataset with 227'000 sequences has been tested and perfectly clustered with the CLARA algorithm.

**WARNINGS** CLARA algorithm is the most efficient one. CLARANS and CLARA-FUZZY are still in test phase.   

## Documentation
The complete documentation is available in the folder *Documentation*
- [User guide](https://github.com/Ltochon/CLARA.seq/blob/master/Documentation/user_guide.pdf)
- [Package Documentation](https://github.com/Ltochon/CLARA.seq/blob/master/Documentation/CLARA.seq_1.1.1.pdf)
