# C3


This is the readme file for C3 and its usage and parameters

## Usage

C3 is a python script, and it can be run as a python script using the format `python c3Method.py --cancer --filter --percentile --up --down --J_percentile --alpha --date --a --bound --cond --address --output --c --network`


## Arguments

`--cancer` this is the cancer argument, or more or less your output tag

`--filter` this argument is a filters the minimum number of times a gene can be altered to count towards C3. Default value is `20` 

`--percentile` a similar argument to filter, it filters only the top percentile of genes by alteration frequency. It acts as a secondary filter paramter in lieu or combined with filter to select altered genes. Default value is `90`.

`--up` an argument that determines the upper threshold for a GISTIC copy number signal to count as significant. Default is `3`

`--down` an argument that determines the lower threshold for a GISTIC copy number signal to count as significant. Default is `-1`
 
`--J_percentile`  the `J_percentile` 4 Default value is `90`.

`--alpha` the `alpha` argument. The default argument is `0.29`. 

`--date` the date argument (for labelling purposes). The default date is `April 29` (the date our project started :) )

`--a` the `a` argument. The default argument is 4

`--bound` the argument that indicates the upward bound for the cluster size. Default is `5`

`--cond` which  positive and negative weights do you address. ME is always a negative weight, so this argument determines the positive weights. We currently support

1. `EX_CO_ME` positive weight being EX (expression) and CO (coverage)
2. `NI_CO_ME` positive weight being NI (expression) and CO (coverage)
3. `triple` using NI, EX and CO as separate weights

`--address` the address you use as your working folder where your files are

`--output` the address of your destination folder

`--c` This determines the weight of the first parameter in the case of a `EX_CO_ME` or `NI_CO_ME` weighting scheme

`--network` the destination of your network folder (if applicable)

`--cov` the weight of your coverage parameter. Only used if `--cond` is `triple`. Default is `0.33`

`--exp` the weight of your expression parameter. Only used if `--cond` is `triple`. Default is `0.33`

`--net` the weight of your network parameter. Only used if `--cond` is `triple`. Default is `0.33`


## Credits

Authors:

Jack Hou
Greg Puleo
Amin Emad
Jian Ma
Olgica Milenkovic
