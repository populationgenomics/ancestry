################################################
#### scRNA-seq data

# this is assuming it is all QC'd, and we are only considering good-quality cells
# also any filtering on genes? E.g. highly variable / highly expressed genes?

# open scRNA-seq count matrix: GxN

# this should be a genes x cells matrix containing normalised, logged counts
# we do not expect missing values, but probably lots of zeroes
# best way to store this? Currently using "pickle"'d version of tables (e.g., csv, tsv)

# Cellular contexts (C): NxC

# option 1: generate expression count PCs
# from the full count matrix, extract cell-level PCs, at least 100 (C=100)
# can this be done especially fast if N is very large, e.g. using hail?
# this will only work if PCs (or other latent variable representations)
# capture axes of variation of interest (e.g. cell (sub)type)

# option 2: use known annotations as contexts
# cell (sub)type annotations - these are available, will probably not capture hyerarchy
# and continuums like latent variables, but probably more interpretable

# option 3: other/combination 

# Phenotype (y): Nx1

# for each test, only one gene at a time will be extracted and used as outcome variable (y)
# this means one row of the count matrix gets considered
# some further manipulation can be applied here, typically standardisation,
# or quantile normalisation

# Covariates (W): NxK

# these will be modeled as fixed effects and we will not test for interactions with these factors
# things like batch, age, sex - (some of) these metadata may actually come as part of the WGS
# data (below)

################################################
#### SNP data

# Genotypes (G): NxS

# this is typically actually well defined as donors x SNPs (DxS),
# but needs to be expanded to cell-level (NxS), such that all cells 
# from a given donor get assigned the same genotype values
# I am used to plink files and usually use  pandas_plink.read_plink1_bin
# but open to different file formats

# GRM (K): NxN

# genetic relatedness / kinship matrix
# How was this calculated? Was this SNP-based?
# Again this is well defined as DxD, but should be expanded out to NxN
# In the model, we actually only provide hK such that hK @ hK.T = K
# to avoid building / loading large matrices
# we use cholesky for this, but there may be different ways/implementations

################################################
#### WGS data

# Am not familiar with this at all at the moment, use hail?