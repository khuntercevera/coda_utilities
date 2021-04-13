#quickly export dataframe columns:

using JLD2, DataFrames, CSV, Dates
using RCall
include("/home/kristen/Documents/coda_utilities/coda_ops_transforms.jl")

#using Distributions


## Need to have oligo_df_plus loaded! This has zeros replaced with small amount:
@load "/home/kristen/Documents/V6V8_analysis/analysis_products/oligo_df_zero_plus.jld2"

##
#y=closure(convert(Array,oligo_df_plus[!,Between(:csO1,:csO6)]))
y=closure(convert(Array,oligo_df_plus[!,Between(:csO1,:csO6)]))
sbp0=[1 -1 1 1 1 1; 1 0 1 -1 -1 1; 0 0 0 1 -1 0; 1 0 1 0 0 -1; 1 0 -1 0 0 0]
ilr_y, bal=ilr(y,sbp0)

wisp=findall(dayofyear.(oligo_df_plus[:,:date]) .< 166)
su=findall((dayofyear.(oligo_df_plus[:,:date]) .>= 166) .& (dayofyear.(oligo_df_plus[:,:date]) .< 258))
fa=findall(dayofyear.(oligo_df_plus[:,:date]) .>= 258)

test_fa = y[fa,:]
test_wisp = y[wisp,:]
test_su = y[su,:]
## Ternary Plots, courtesy of the lovely Ternary package in R:

R"require(compositions)"
#R"require(png)"

R"""
W_tr=t($sbp0)
V=gsi.buildilrBase(W_tr)
y=acomp($test_fa)
colnames(y) <- c('O1','O2','O3','O4','O5','O6')
png(file="/home/kristen/Documents/V6V8_analysis/paper_V6V8/EM_submission/test_dendro_fa.png")
CoDaDendrogram(y,V)
dev.off()
"""
##

R"""
W_tr=t($sbp0)
V=gsi.buildilrBase(W_tr)
y=acomp($test_wisp)
colnames(y) <- c('O1','O2','O3','O4','O5','O6')
png(file="/home/kristen/Documents/V6V8_analysis/paper_V6V8/EM_submission/test_dendro_wisp.png")
CoDaDendrogram(y,V)
dev.off()
"""
##


##
R"""
W_tr=t($sbp0)
V=gsi.buildilrBase(W_tr)
y=acomp($test_su)
colnames(y) <- c('O1','O2','O3','O4','O5','O6')
png(file="/home/kristen/Documents/V6V8_analysis/paper_V6V8/EM_submission/test_dendro_su.png")
CoDaDendrogram(y,V)
dev.off()
"""




#break up year into seasons:

# X_wisp = convert(Array{Union{Float64,Missing}},oligo_df_plus[wisp,vars])
# Y_wisp = ilr_y[wisp,:]
# #jj=findall(ismissing.(X_wisp)) #-> this should be zero

# X_su = convert(Array{Union{Float64,Missing}},oligo_df_plus[su,vars])
# Y_su = ilr_y[su,:]
# jj=findall(ismissing.(X_su))
# qq=unique(map(x->x[1], jj)) #no nutrient values for one of the samples, exclude
# ii=setdiff(collect(1:size(X_su,1)),qq)
# X_su=X_su[ii,:]
# Y_su=Y_su[ii,:]