# kCollectiveFP
The source code for our work published in MDM2019, namely k-Collective Influential Facility Placement over Moving Object.

## Environment

1.  These experiments are implementes in C++.

2.  IDE is VS 2013.

##Dataset

1.  There are two datasets which are recorded by two text files(`checkins-10162.txt` and `checkins-2321.txt`).

    *   checkins-10162.txt: 10162 users from Gowalla located in California.
    *   checkins-2321.txt: 2321 users from Foursquare located in Singapore.
2.  The trajectory of user consist of username and a series of check-in points.

3.  Candidate sets are randomly generated from check-in points.
4. Both datasets can be found in this page(https://github.com/lihuixidian/PINOCCHIO/tree/master/datasets)
).


##Algorithm

1. `PINOCCHIO`: We evaluate the inf(·) for all candidates, and select the top k candidates with the maximum inf(·) as the results.
2. `GreedyP`: The GreedyP algorithm in Algorithm 1 in paper.
3. `GreedyPS`: The GreedyPS algorithm in Algorithm 2 in paper.

##Supplement

1.  `PLS.cpp `: Main function. 

2.  `pino.cpp`: This file contains the implements of all algorithoms (eg., Baseline, GreedyP, GreedyPS) in our paper.

## Usage

1.  All data files should be placed in a local folder named as 'Release', e.g., '`D:\Experiment\PLS\Release`'.

2.  We should load `boost library` in this project which provides corresponding utilities with respect to R-tree.

3.  There are some precompiles in program.

    `CANDIDATES_GENERATION` : Generate candidates.
    `GEN_FROM_10162` : From 10162 (1) or 2321 (0) dataset.
    `PICK_FROM_UNIQUE` : In CANDIDATES_GENERATION, pick from unique coordinates (set) but not all check-in logs (vector).
    `CHECKINS_EXCLUDING`: Exclude check-ins that appear in candidates. And generate checkins.txt.
    `DATA_LOADING` : Load the datas about candidates and users.
    **They are all in PLS.cpp.**
4. `HUSH_NUM`: the number of bitmaps in `pino.h`.
