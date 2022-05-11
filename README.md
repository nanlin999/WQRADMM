# WQRADMM for Windows users
R package for distributed quantile regression in longitudinal big data based on multi-block ADMM.

It can be used to reproduce the simulation studies in the following paper:

**Ye Fan, Nan Lin and Liqun Yu**. *Distributed Quantile Regression for Longitudinal Big Data.*

Two main functions are included: **WQRADMM( )** and **paraWQRADMM( )**, repectively designed for non-distributed and distributed environments.

**Note:** please install RTools and further put its location on the environment variable PATH before installing this package. 

**A sample code**

(```)

e = rnorm(100)

(```)
