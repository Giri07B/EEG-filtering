# EEG-filtering
Removal of EOG artifacts from EEG using Adaptive Filtering.

Algorithm reference: Removal of Ocular Artifacts from Electro-encephalogram by adaptive filtering by P.He - Medical & Biological Engineering & Computing 2004, Vol. 42

Data used: Published by Manousos Klados and Panagiotis Bamidis - https://data.mendeley.com/datasets/wb6yvr725d/3

File Description
1) raw_data.m - Display data used.
2) claim1.m - This method is easy to implement and stable, converges fast and is suitable for on- line removal of EOG artifacts.
3) claim2.m - When the filtering algorithm is applied to an EEG recorded at a remote site, e.g. O2, that contains very few EOG artifacts,               the filters are basically shut down and the original EEG simply passes the system without visible changes.
4) claim3.m - The exact value of λ – “the forgetting factor” is not critical to the performance of the algorithm.
5) claim4.m - The performance of the adaptive filter is not sensitive to the choice of M: “the filter order”.
6) claim5.m - As the filter order increases, Mean square error (MSE) decreases.
