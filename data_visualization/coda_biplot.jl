# Script to construct a biplot of compositional data:

using JLD2, DataFrames, CSV, Dates, Query
using LinearAlgebra
using Gadfly, Colors, Compose, ColorSchemes

include("/home/kristen/Documents/coda_utilities/coda_ops_transforms.jl")

#load in the data:
@load "/home/kristen/Documents/V6V8_analysis/analysis_products/oligo_df_zero_plus.jld2"
#should return dataframe oligo_df_plus

## center the data:
data=convert(Array,oligo_df_plus[!,Between(:csO1,:csO6)]) #convert columns of dataframe to just an array
y=closure(data)
cen=center(y) #calcuate center of dataset
centered=perturb(y,inverse(cen)) #perturb by inverse of this center

Z=clr(centered) #take clr transform
#sum(centered, dims=2) # should all be 1!

#Perform SVD:
bp = svd(clr(centered))
U, S, V = bp
#double checks:
#ztest = U * Diagonal(S) * V'
#Z ≈ U * Diagonal(S) * V' #this should be true!

#Construct the biplots:
#Z=FG^T
alpha=0 #α = 1, form biplot, α = 0; covariance biplot
F=U[:,1:2] * diagm(S[1:2].^alpha)
G=diagm(S[1:2].^(1-alpha)) * V'[1:2,:]
# simply scaled or unscaled singular values

#From CoDa book - doesn't say why scaling of sqrt(n) is there?
# nu=0
# a = sqrt(size(Z,1)) * U[:,1:2] * diagm(S[1:2].^(1-nu))
# b = (1/sqrt(size(Z,1))) * diagm(S[1:2].^nu) * V'[1:2,:]

#percentage of variance explained:
#values of S are the eigenvalues:
S[1].^2 / sum(S.^2) #PC 1
S[2].^2 / sum(S.^2) #PC 2
(S[2].^2 + S[1].^2) / sum(S.^2) #total var explained by first two PCs

rows=F
cols=G # / sqrt(128) #for covariance biplot, scale back

## Plot with scaled row and columns:
#  Rows are scaled by sqrt(n-1), column vectors are scaled by 1/sqrt(n-1):

Gadfly.set_default_plot_size(20cm,17cm)
s=7
coord = Coord.cartesian(xmin=-s, xmax=s, ymin=-s, ymax=s) # Geom.vector also needs scales for both axes
xsc  = Scale.x_continuous(minvalue=-s, maxvalue=s)
ysc  = Scale.y_continuous(minvalue=-s, maxvalue=s)

layer1 = layer(x=zeros(size(cols,2),1), y=zeros(size(cols,2),1), xend=(1/sqrt(128)).*cols[1,:], yend=(1/sqrt(128)).*cols[2,:], Geom.segment,style(default_color=RGB(0.4,0.4,0.4)))
layer2=layer(x=(1/sqrt(128)).*cols[1,:], y=(1/sqrt(128)).*cols[2,:],label=["O1","O2","O3","O4","O5","O6"],Geom.point,Geom.label,style(default_color=RGB(0.4,0.4,0.4)))
layer3=layer(x=sqrt(128).*rows[:,1], y=sqrt(128).*rows[:,2], color=dayofyear.(oligo_df_plus[!,:date]), Geom.point)

plt=plot(layer1,layer2,layer3,xsc, ysc,coord,
        Scale.color_continuous(minvalue=1,maxvalue=365,colormap=p -> get(ColorSchemes.jet1, p)),Guide.colorkey(title="Year Day"),
        Guide.xlabel("PCA1, 62.3%"),Guide.ylabel("PCA2, 29.5%"))

## gadfly magic:
#export:
plt |> PNG("/home/kristen/Documents/V6V8_analysis/paper_V6V8/V6V8_Syn_CoDa_paper/figures/Fig3_clrPCA.png"); #@async run(`eog /tmp/test.svg`)


##  export if needed somewhere else:
#CSV.write("/home/kristen/Documents/V6V8_analysis/scripts/rel_abnd_analysis_scripts/matlab_scripts/olgio_df_plus.csv", oligo_df_plus, delim = ',')
