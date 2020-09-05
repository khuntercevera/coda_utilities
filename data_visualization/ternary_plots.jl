# If want to explore ternary plots and line fitting:

using JLD2, DataFrames, CSV, Dates, Query
using LinearAlgebra, GLM
using RCall
using Gadfly, Colors

include("/home/kristen/Documents/coda_utilities/coda_ops_transforms.jl")

#load in the data:
@load "/home/kristen/Documents/V6V8_analysis/analysis_products/oligo_df_zero_plus.jld2"
#should return dataframe oligo_df_plus

## Ternary Plots, courtesy of the lovely Ternary package in R:

R"require(Ternary)"
R"require(png)"

x136=closure(convert(Array,oligo_df_plus[!,[:Oligo1, :Oligo3,:Oligo6]])) #convert Dataframe columns to array
R"""
TernaryPlot()
TernaryPlot(atip = "O1-IC", btip = "O3-II/XV", ctip = "O6-I?", lab.col=c("#0000FF","#FF7F00","#000099"),tip.col=c("#0000FF","#FF7F00","#000099"))
TernaryPoints($x136, pch=16, cex=0.8, col='black')
"""

x245=closure(convert(Array,oligo_df_plus[!,[:csO2, :csO4,:csO5]]))
R"""
TernaryPlot()
TernaryPlot(atip = "O2-CB5",btip = "O4-III/IV", ctip = "O5-IE", lab.col=c("#CC334C","#009900","#00CCFF"),tip.col=c("#CC334C","#009900","#00CCFF"),grid.lines = 10, grid.minor.lines = 0, grid.col = "lightgrey")
TernaryPoints($x245, pch=16, cex=0.8, col='black')
"""


## Fancier year day colors:

x123=closure(convert(Array,oligo_df_plus[!,[:csO1,:csO2,:csO3]]))
yearday=dayofyear.(oligo_df_plus[!,:date])

R"""
TernaryPlot()
jet <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                     "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
jmap <- jet(365)[$yearday] #first argument is max....I think...
TernaryPlot(atip = "O1-I", btip = "O2-CB5", ctip = "O3-II/XV", lab.col=c("#0000FF","#CC334C","#FF7F00"),tip.col=c("#0000FF","#CC334C","#FF7F00"),grid.lines = 10, grid.minor.lines = 0, grid.col = "lightgrey")
TernaryPoints($x123, pch=16, cex=0.8, col=jmap) #whoo hoo!!!!!
"""

## Fit a line through this data in the simplex

#Note, this is easiest to do with an ilr tranform, fitting the line with standard linear regression and then back transfrom to simplex:

#Line throuhg oligotype 1,2,3 data:
x123=closure(convert(Array,oligo_df_plus[!,[:csO1,:csO2,:csO3]]))
sbp=[-1 -1 1; 1 -1 -0] #to form basis
ilr_x, bal=ilr(x123,sbp) #ilr transform

#see what this transformed data looks like:
Gadfly.set_default_plot_size(35cm,25cm)
layer1=layer(x=ilr_x[:,2], y=ilr_x[:,1],Geom.point, style(default_color=RGB(0,0,1))) #okay, that seems to be working...
plt=plot(layer1)

#fit line:
T1=DataFrame(x = ilr_x[:,2], y = ilr_x[:,1])
test1 = lm(@formula(y ~ x), T1)
bf=zeros(2,1)
bf[2] = GLM.coef(test1)[2]
bf[1] = GLM.coef(test1)[1]

xf=collect(-2:0.1:8) #and see what this fit looks like:
yf=bf[1] .+ bf[2] .* xf
layer2=layer(x=xf, y=yf, Geom.line)
push!(plt,layer2) #check that line goes through points!


## Now plot this line back in the simplex and save:

ll=invilr([yf xf],bal) #back transform fitted line

R"""
TernaryPlot()
jet <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                     "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
jmap <- jet(365)[$yearday] #first argument is max....I think...
TernaryPlot(atip = "O1-I", btip = "O2-CB5", ctip = "O3-II/XV", lab.col=c("#0000FF","#CC334C","#FF7F00"),tip.col=c("#0000FF","#CC334C","#FF7F00"),grid.lines = 10, grid.minor.lines = 0, grid.col = "lightgrey")
TernaryPoints($x123, pch=16, cex=0.8, col=jmap) #whoo hoo!!!!!
TernaryLines($ll, col='black')
"""

#and if wanted to save:
R"""
png(file="/home/kristen/Documents/V6V8_analysis/paper_V6V8/V6V8_Syn_CoDa_paper/figures/tern1.png")
TernaryPlot()
jet <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                     "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
jmap <- jet(365)[$yearday] #first argument is max....I think...
TernaryPlot(atip = "O1-I", btip = "O2-CB5", ctip = "O3-II/XV", lab.col=c("#0000FF","#CC334C","#FF7F00"),tip.col=c("#0000FF","#CC334C","#FF7F00"),grid.lines = 10, grid.minor.lines = 0, grid.col = "lightgrey")
TernaryPoints($x123, pch=16, cex=0.8, col=jmap) #whoo hoo!!!!!
TernaryLines($ll, col='black')
dev.off()
"""




## Again for the remaining set:

x456=closure(convert(Array,oligo_df_plus[!,[:csO4,:csO5,:csO6]]))
sbp=[-1 1 -1; 1 0 -1]
ilr_x, bal=ilr(x456,sbp)
layer1=layer(x=ilr_x[:,2], y=ilr_x[:,1],Geom.point, style(default_color=RGB(0,0,1))) #okay, that seems to be working...
plt=plot(layer1)

#fit line just within date:
T1=DataFrame(x = ilr_x[:,2], y = ilr_x[:,1])
test1 = lm(@formula(y ~ x), T1)
bf=zeros(2,1)
bf[2] = GLM.coef(test1)[2]
bf[1] = GLM.coef(test1)[1]

xf=collect(-3:0.1:5)
yf=bf[1] .+ bf[2] .* xf
layer2=layer(x=xf, y=yf, Geom.line)
push!(plt,layer2) #check that line goes through points!

ll=invilr([yf xf],bal) #back transform fitted line...

R"""
png(file="/home/kristen/Documents/V6V8_analysis/paper_V6V8/V6V8_Syn_CoDa_paper/figures/tern2.png")
TernaryPlot()
jet <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                     "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
jmap <- jet(365)[$yearday] #first argument is max....I think...
TernaryPlot(atip = "O4-III/IV", btip = "O5-IE", ctip = "O6-I",lab.col=c("#009900","#00CCFF","#000099"),tip.col=c("#009900","#00CCFF","#000099"),grid.lines = 10, grid.minor.lines = 0, grid.col = "lightgrey")
TernaryPoints($x456, pch=16, cex=0.8, col=jmap) #whoo hoo!!!!!
TernaryLines($ll, col='black')
dev.off()
"""



## to get the simplified relationship:

a=ll[:,1]./ll[:,2]
b=ll[:,2]./ll[:,3]

# a=ll[:,1]./ll[:,3]
# b=ll[:,3]./ll[:,2]

T1=DataFrame(x = log.(a), y = log.(b))

test2 = lm(@formula(y ~ x), T1)
bf=zeros(2,1)
bf[2] = GLM.coef(test2)[2]
bf[1] = GLM.coef(test2)[1]

plot(x=log.(a),y=log.(b),Geom.point)
