# diversity / dissimilarity index calculation and plots:

using JLD2, DataFrames, CSV, Dates, Query
using LinearAlgebra, Clustering, StatsPlots

include("/home/kristen/Documents/coda_utilities/coda_ops_transforms.jl")

#load in the data:
@load "/home/kristen/Documents/V6V8_analysis/analysis_products/oligo_df_zero_plus.jld2"