# diversity / dissimilarity index calculation and plots:

using JLD2, DataFrames, CSV, Dates, Query
using LinearAlgebra, Clustering
using MAT
#using StatsPlots
using Gadfly

include("/home/kristen/Documents/coda_utilities/coda_ops_transforms.jl")

#load in the data:
@load "/home/kristen/Documents/V6V8_analysis/analysis_products/oligo_df_zero_plus.jld2"

## compute distance of each sample to all others, and plot by season:

y=closure(convert(Array,oligo_df_plus[!,Between(:csO1,:csO6)])) #rows are samples, columns are components
m=size(y,1)
oligo_dist = zeros(Float64,m,m)#Array{Float64,2}(undef,m,m)

for i=1:m
    for j=i:m
       oligo_dist[i,j] = ait_dist(y[i,:],y[j,:])
       oligo_dist[j,i] = ait_dist(y[i,:],y[j,:])
    end
end

## hmmm...perhaps, easier to organize by yearday?

ind=sortperm(dayofyear.(oligo_df_plus[!,:date])) #rev=true) #returns indicies...
# dayofyear.(oligo_df_plus[ind,:date]) #to check....
y_yrdy = y[ind,:]
oligo_byday = zeros(Float64,m,m)#Array{Float64,2}(undef,m,m)

for i=1:m
    for j=i:m
       oligo_byday[i,j] = ait_dist(y_yrdy[i,:],y_yrdy[j,:])
       oligo_byday[j,i] = ait_dist(y_yrdy[i,:],y_yrdy[j,:])
    end
end





## NEED TO WORK ON COLOR CODING, ETC...
file = matopen("/home/kristen/Documents/V6V8_analysis/paper_V6V8/matlab_figure_scripts/exports_from_julia/dissimilarity_matrix.mat", "w")
write(file, "oligo_dissim", oligo_dist)
write(file, "oligo_dissim_yrdy", oligo_byday)
write(file,"ind", ind)
write(file,"sample_yearday", dayofyear.(oligo_df_plus[!,:date]))
close(file)






## easier in a dataframe for plotting?
y=closure(convert(Array,oligo_df_plus[!,Between(:csO1,:csO6)])) #rows are samples, columns are components
m=size(y,1)
df_dist = DataFrame(x=Int64[], y=Int64[], dist=Float64[], date=DateTime[])

for i=1:m
    for j=i:m
        push!(df_dist,[i,j,ait_dist(y[i,:],y[j,:]),oligo_df_plus[i,:date]])
    end
end

## #yes - definitely easier for plotting! Very cool!

df_dist.yearday=dayofyear.(df_dist[!,:date])
plot(sort(df_dist,:yearday),x="x",y="y",color="dist", Geom.rectbin) #oops, no - this gives the same thing!



#CSV.write("/home/kristen/Documents/V6V8_analysis/paper_V6V8/matlab_figure_scripts/exports_from_julia/etc...", clade_array, delim = ',')














#what about a dendrogram? 
#visualize! but must be with Plots.jl!
hc = hclust(Symmetric(oligo_dist), linkage=:single, uplo=:U)
#-> a few other linkage options, check hclust help

w=plot(hc,xticks=false,linewidth=2)
#node_labels=["O1-I","O2-CB5","O3-II","O4-III/IV","O5-I","O6-I*"] #oligotypes 1-6
#node_labels[hc.order]
#xticks!(w,([1,2,3,4,5,6],node_labels[hc.order]),fontsize=16) #maps from hc and matching up oligos to clades
#ylabel!(w,"Aitchision distance")
#savefig(w, "/home/kristen/Documents/V6V8_analysis/paper_V6V8/V6V8_Syn_CoDa_paper/figures/dendrogram_single.pdf")










## very cool! now, let's check out the aitchsion square norm for inqueality:


y=closure(convert(Array,oligo_df_plus[!,Between(:csO1,:csO6)])) #rows are samples, columns are components
m=size(y,1)
oligo_ineq = zeros(Float64,m,1)#Array{Float64,2}(undef,m,m)

for i=1:m
    oligo_ineq[i] = square_ait_norm(y[i,:])
end

## export it!
file = matopen("/home/kristen/Documents/V6V8_analysis/paper_V6V8/matlab_figure_scripts/exports_from_julia/ineq_metric.mat", "w")
write(file, "oligo_ineq", oligo_ineq)
rr=1 .- exp.(-oligo_ineq)
write(file,"scaled_ineq", rr) #should already ahve the smaples dates somewhere else...
close(file)

## plot it up!

Gadfly.set_default_plot_size(30cm,10cm)

xticks=DateTime("2010-01-01"):Year(1):DateTime("2019-01-01")
xyearticks = Dates.value.(xticks)
xyearticks_label=map(x-> Dates.format(x,"yyyy"),xticks)
labels = Dict(zip(xyearticks, xyearticks_label))

fontsize=Theme(minor_label_font_size=16pt,major_label_font_size=16pt,plot_padding=[5mm, 5mm, 2mm, 2mm])

l1=layer(x=Dates.value.(oligo_df_plus[!,:date]),y=oligo_ineq,Geom.line,Geom.point)
p2=plot(l1,Guide.xticks(ticks=xyearticks),Guide.xlabel(""),Scale.x_continuous(labels = x -> labels[x]),fontsize)


## can also scale it from 0-1:

rr=1 .- exp.(-oligo_ineq)
l2=layer(x=Dates.value.(oligo_df_plus[!,:date]),y=rr,Geom.line,Geom.point)
p2=plot(l2,Guide.xticks(ticks=xyearticks),Guide.xlabel(""),Scale.x_continuous(labels = x -> labels[x]),fontsize)




        
        
        #Guide.yticks(ticks=[0 0.25 0.5 0.75 1])
        # Coord.Cartesian(ymin=0, ymax=1),
        # Guide.ylabel("Relative abundance"),fontsize,
        # Guide.annotation(compose(context(), text(Dates.value.(DateTime(2010,1,10)), 0.9, "B"))))
