# Distance dendrogram

using JLD2, DataFrames, CSV, Dates, Query
using LinearAlgebra, Clustering, StatsPlots

include("/home/kristen/Documents/coda_utilities/coda_ops_transforms.jl")

#load in the data:
@load "/home/kristen/Documents/V6V8_analysis/analysis_products/oligo_df_zero_plus.jld2"
#should return dataframe oligo_df_plus

## closure on the dataset:
y=closure(convert(Array,oligo_df_plus[!,Between(:csO1,:csO6)])) #seasonally imputed zeros
# findall(y .== 0) #this should be 0 as a double check!
vyp=ait_var(y) #compute aitchison variation

#other metrics to consider:
# test=exp.(-vyp.^2 / 2)
# test=exp.(-sqrt.(vyp))

#visualize! but must be with Plots.jl!
hc = hclust(Symmetric(ait_var(y)), linkage=:single, uplo=:U)
#-> a few other linkage options, check hclust help

##
w=plot(hc,xticks=false,linewidth=2)
node_labels=["O1-I","O2-CB5","O3-II","O4-III/IV","O5-I","O6-I*"] #oligotypes 1-6
#node_labels[hc.order]
xticks!(w,([1,2,3,4,5,6],node_labels[hc.order]),fontsize=16) #maps from hc and matching up oligos to clades
ylabel!(w,"Distance, Aitchision variation")
##
savefig(w, "/home/kristen/Documents/V6V8_analysis/paper_V6V8/V6V8_Syn_CoDa_paper/figures/dendrogram_single.pdf")
