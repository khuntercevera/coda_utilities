# Zero Imputation for compositional data:

using JLD2, DataFrames, CSV, Dates

include("/home/kristen/Documents/coda_utilities/coda_ops_transforms.jl")

## load in data:
@load "/home/kristen/Documents/V6V8_analysis/analysis_products/mvco_syn_oligo_counts_plus_env_data.jld2"
@load "/home/kristen/Documents/V6V8_analysis/analysis_products/mvco_syndata.jld2" time_syn daily_syn syn_avg

## replace zeros with an appropriate estimate:

# A few different replacement options:
# r1= bayes_multi_zero_replace(test,"Geometric","modified")
# r2=bayes_multi_zero_replace(test,"Square-root","modified")
# r3=bayes_multi_zero_replace(test,"Square-root","uniform")
# r4=bayes_multi_zero_replace(test,"Bayes-Laplace","uniform")
#
# j=2
# ii=findall(testc[:,j] .== 0)
#show(stdout, "text/plain", [testc[ii,j] r[ii,j] r2[ii,j] r3[ii,j] r4[ii,j]])

#whole year:
temp=convert(Array,oligo_counts_df[!,Between(:Oligo1,:Oligo6)])
r_allyear= bayes_multi_zero_replace(temp,"Geometric","modified")

#by season -
# wisp=findall(dayofyear.(oligo_counts_df[!,:date]) .<=182)
# su=findall((dayofyear.(oligo_counts_df[!,:date]) .> 182) .& (dayofyear.(oligo_counts_df[!,:date]) .< 274))
# fa=findall(dayofyear.(oligo_counts_df[!,:date]) .>= 274)

#dayofyear.(DateTime(2003,June,15))
wisp=findall(dayofyear.(oligo_counts_df[!,:date]) .< 166)
su=findall((dayofyear.(oligo_counts_df[!,:date]) .>= 166) .& (dayofyear.(oligo_counts_df[!,:date]) .< 258))
fa=findall(dayofyear.(oligo_counts_df[!,:date]) .>= 258)

r_wisp= bayes_multi_zero_replace(temp[wisp,:],"Geometric","modified")
r_su= bayes_multi_zero_replace(temp[su,:],"Geometric","modified")
r_fa= bayes_multi_zero_replace(temp[fa,:],"Geometric","modified")

temp2=ones(size(temp))
temp2[wisp,:]=r_wisp
temp2[su,:]=r_su
temp2[fa,:]=r_fa


# add to database:
oligo_df_plus=deepcopy(oligo_counts_df[!,Not(Between(:Oligo7,:Oligo14))])
oligo_df_plus.cO1=r_allyear[:,1]
oligo_df_plus.cO2=r_allyear[:,2]
oligo_df_plus.cO3=r_allyear[:,3]
oligo_df_plus.cO4=r_allyear[:,4]
oligo_df_plus.cO5=r_allyear[:,5]
oligo_df_plus.cO6=r_allyear[:,6]

oligo_df_plus.csO1=temp2[:,1]
oligo_df_plus.csO2=temp2[:,2]
oligo_df_plus.csO3=temp2[:,3]
oligo_df_plus.csO4=temp2[:,4]
oligo_df_plus.csO5=temp2[:,5]
oligo_df_plus.csO6=temp2[:,6]

## if want to exclude too low abundance samples:
# ind=findall(oligo_counts_df[!,:syncount] .> 10)
# oligo_df_plus=oligo_df_plus[ind,:]

## if want to export any:
#CSV.write("/home/kristen/Documents/V6V8_analysis/scripts/rel_abnd_analysis_scripts/matlab_scripts/oligo_df_plus.csv", oligo_df_plus, delim = ',')
#CSV.write("/home/kristen/Documents/V6V8_analysis/scripts/rel_abnd_analysis_scripts/matlab_scripts/oligos_7_14.csv", oligo_counts_df[!,Between(:Oligo7,:Oligo14)], delim = ',')

#and save as .jld2 for future:
@save "/home/kristen/Documents/V6V8_analysis/analysis_products/oligo_df_zero_plus.jld2" oligo_df_plus


## Example in paper:

C=[6 3 9 15 3 24;
3  2  4  7  0  12;
24  12  36 60 12 96;
1  0  2  3  0  4]

X=closure(C)

r = bayes_multi_zero_replace(c,"Geometric","modified")
