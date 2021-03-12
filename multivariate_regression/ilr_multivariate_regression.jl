#Perform mulivariate regression on ilr transformed compositions

#Kristen R. Hunter-Cevera, Marine Biological Laboratory, Spring 2020

using JLD2, DataFrames, CSV, Dates, Query
#using StatsPlots, Colors
using LinearAlgebra, Clustering, GLM
using Gadfly, Colors, Compose
using RCall

using Distributions

include("/home/kristen/Documents/coda_utilities/coda_ops_transforms.jl")

## Need to have oligo_df_plus loaded! This has zeros replaced with small amount:
@load "/home/kristen/Documents/V6V8_analysis/analysis_products/oligo_df_zero_plus.jld2"


##############################################################################################

#F and Chi-squared approximations:
function approx_F(Λ,p,vE,vH)
    if (p * vH) == 2
        t = 1
    else
        t = sqrt((p^2 * vH^2 - 4) / (p^2 + vH^2 -5))
    end
    w = vE + vH - 0.5 * (p + vH +1)
    df1 = p * vH
    df2 = w * t - 0.5 * (p * vH -2)

    F = ((1 - Λ.^(1/t)) / (Λ.^(1/t))) .* (df2 / df1)
    return F, df1, df2
end


function approx_Chi(Λ,p,vE,vH)
    df = p * vH
    chi = -(vE - 0.5 *(p - vH +1)) * log(Λ)
    return chi, df
end

##############################################################################################

#transform compositions to multivariate real numbers with ilr:
y=closure(convert(Array,oligo_df_plus[!,Between(:csO1,:csO6)]))
sbp0=[1 -1 1 1 1 1; 1 0 1 -1 -1 1; 0 0 0 1 -1 0; 1 0 1 0 0 -1; 1 0 -1 0 0 0]
lb0=["O1,O3,O4,O5,O6 | O2", "O1,O3,O6 | O4,O5", "O4 | O5", "O1,03 | O6", "O1 | O3"] #labels

ilr_y, bal=ilr(y,sbp0) #ilr transform
lb=lb0 #set labels

## Full model regression test - manually choose the variable combination:

vars=[:temp,:lightweek,:PO4,:NH4,:SiOH,:NO3]

#break up year into seasons:
wisp=findall(dayofyear.(oligo_df_plus[:,:date]) .< 166)
su=findall((dayofyear.(oligo_df_plus[:,:date]) .>= 166) .& (dayofyear.(oligo_df_plus[:,:date]) .< 258))
fa=findall(dayofyear.(oligo_df_plus[:,:date]) .>= 258)

X_wisp = convert(Array{Union{Float64,Missing}},oligo_df_plus[wisp,vars])
Y_wisp = ilr_y[wisp,:]
#jj=findall(ismissing.(X_wisp)) #-> this should be zero

X_su = convert(Array{Union{Float64,Missing}},oligo_df_plus[su,vars])
Y_su = ilr_y[su,:]
jj=findall(ismissing.(X_su))
qq=unique(map(x->x[1], jj)) #no nutrient values for one of the samples, exclude
ii=setdiff(collect(1:size(X_su,1)),qq)
X_su=X_su[ii,:]
Y_su=Y_su[ii,:]

X_fa = convert(Array{Union{Float64,Missing}},oligo_df_plus[fa,vars])
Y_fa = ilr_y[fa,:]
#jj=findall(ismissing.(X_fa)) #-> this also should be zero


## Step wise progression to test significance of added variables, one at a time:

#find the minimum number of variables that explains the most variability for each season:


#record in a dataframe
stepwise_df=DataFrame(season=String[], step=Int64[], base=String[], var=String[], Lambda=Float64[], F=Float64[], df1=Float64[], df2=Float64[], F_pval=Float64[], B_hat=Array{Union{Missing, Float64},2}[])

# first, find out which variable results in lowest lambda value:
seasons=["winter/spring", "summer", "fall"]

for j=1:3

    global season=seasons[j]
    #season = "winter/spring"
    if season == "winter/spring"
        X0=X_wisp
        Y=Y_wisp
    elseif season == "summer"
        X0=X_su
        Y=Y_su
    elseif season == "fall"
        X0=X_fa
        Y=Y_fa
    end

    intercept=ones(size(X0,1),1) #setup design matrix
    n=size(X0,1)
    y_bar=mean(Y, dims =1)

    @show season

    for q=1:6 # number of predictor variables to run through

        #@show q
        X=hcat(intercept,X0[:,q]) #build design matrix with each variable

        B_hat= inv(X' * X) * (X' * Y) #retrieve parameters

        H_full = (Y' * Y) .- (B_hat' * X' * Y) #calculate lambda
        H_null = Y' * Y .- (n .* (y_bar' * y_bar)) #E+H in book parlance
        Λf = det(H_full) ./ det(H_null)

        vH = length(q) #number of X vars, intercept excluded
        vE = n-(length(q))-1 #n - x vars -1
        p = size(Y,2) #num y vars

        f, df1 ,df2=approx_F(Λf,p,vE,vH)
        pvalF=1-cdf(FDist(df1,df2),f)

        push!(stepwise_df,[season,0,string(vars[q]), "-", Λf, f, df1, df2, pvalF, B_hat])

end

#once we have a baseline, find which var explains the most variance and add from there!
r=findall(stepwise_df[!,:season] .== season)
global qbase=argmin(stepwise_df[r,:Lambda])
global qtest=setdiff(1:6,qbase)
global stepp=1

while stepp <= 6
    global qtest
    global qbase
    global stepp

    for q in qtest

        # @show q
        # @show qbase
        X=hcat(intercept,X0[:,cat(qbase,q,dims=1)]) #full model with variable q as candidate variable
        B_hat= inv(X' * X) * (X' * Y) #retrieve parameters
        H_full = (Y' * Y) .- (B_hat' * X' * Y) #calculate lambda
        H_null = Y' * Y .- (n .* (y_bar' * y_bar)) #E+H in book parlance
        Λf = det(H_full) ./ det(H_null)

        Xr=hcat(intercept,X0[:,qbase]) #Reduced model
        B_hat_r = inv(Xr' * Xr) * (Xr' * Y)
        H_r = (Y' * Y) .- (B_hat_r' * Xr' * Y)
        Λr = det(H_r) ./ det(H_null) #overall test of signifance, 0.86-0.89

        Λ_test = Λf ./ Λr #Lamba to compare addition of extra variable

        # approximate F and χ statistics:
        vH = length(q) #additional varsglobal step=1
        vE = n-(length(qbase))-1 #n - x vars in reduced model -1
        p = size(Y,2) #num y vars

        f, df1, df2 = approx_F(Λ_test,p,vE,vH)
        #chi, df = approx_Chi(Λ_test,p,vE,vH)

        #p-values:
        pvalF=1-cdf(FDist(df1,df2),f)
        #pvalchi = 1-cdf(Chisq(df),chi)
        if length(qbase) == 1
            push!(stepwise_df,[season, stepp, string(vars[qbase]), string(vars[q]), Λ_test, f, df1, df2, pvalF, B_hat])
        else
            push!(stepwise_df,[season, stepp, join(string.(vars[qbase]),"-"), string(vars[q]), Λ_test,f, df1, df2, pvalF, B_hat])
        end

    end #end qtest loop

    r=findall((stepwise_df[!,:step] .== stepp) .& (stepwise_df[!,:season] .== season))

    if any(stepwise_df[r,:F_pval] .< 0.05) #alright- another variable does add significant explanation!

        qbase=vcat(qbase,qtest[argmin(stepwise_df[r,:Lambda])])
        qtest=setdiff(1:6,qbase)
        stepp=stepp .+ 1

    else
        break
    end

end

end #seaons
## Sweet!  Now, backtransform for plotting if needed:

#if needed (see below for saving)
@load  "/home/kristen/Documents/V6V8_analysis/analysis_products/oligo_ilr_regression_output2.jld2"

##Parameter values back to simplex:

#sumwinter/spring: temp-lightweek
rr=findall(stepwise_df[!,:season] .== "winter/spring")
final_step=maximum(stepwise_df[rr,:step])-1
rr2=findall((stepwise_df[!,:season] .== "winter/spring") .& (stepwise_df[!,:step] .== final_step))
best_model_wisp=rr2[argmin(stepwise_df[rr2,:Lambda])]

bhat_wisp_simplex=invilr(stepwise_df[best_model_wisp,:B_hat],bal)


#summer: temp-PO4-lightweek
rr=findall(stepwise_df[!,:season] .== "summer")
final_step=maximum(stepwise_df[rr,:step])-1
rr2=findall((stepwise_df[!,:season] .== "summer") .& (stepwise_df[!,:step] .== final_step))
best_model_su=rr2[argmin(stepwise_df[rr2,:Lambda])]

bhat_summer_simplex=invilr(stepwise_df[best_model_su,:B_hat],bal)

#fall: temp-lightweek
rr=findall(stepwise_df[!,:season] .== "fall")
final_step=maximum(stepwise_df[rr,:step])-1
rr2=findall((stepwise_df[!,:season] .== "fall") .& (stepwise_df[!,:step] .== final_step))
best_model_fa=rr2[argmin(stepwise_df[rr2,:Lambda])]

bhat_fall_simplex=invilr(stepwise_df[best_model_fa,:B_hat],bal) #bhat_wisp = Bhat


## Fitted data back to simplex:

Y_hat = Array{Union{Missing, Float64}}(missing, size(ilr_y))

#populate by seasons:
x =hcat(ones(size(wisp,1),1), X_wisp[:,[1,2]])
b = stepwise_df[best_model_wisp,:B_hat]
temp = b' * x' #B_hat is full model
Y_hat[wisp,:] .= temp'

x =hcat(ones(size(ii,1),1), X_su[:,[1,3,2]]) #initial inexing into summer as one sample does not have all the nutrients!
b = stepwise_df[best_model_su,:B_hat]
temp = b' * x' #B_hat is full model
Y_hat[su[ii],:] .= temp'

x =hcat(ones(size(fa,1),1), X_fa[:,[1,2]])
b = stepwise_df[best_model_fa,:B_hat]
temp = b' * x' #B_hat is full model
Y_hat[fa,:] .= temp'

#and back transform to the simplex:
yhat_simplex = invilr(Y_hat,bal)



#Only temperature as parameter:
Y_hat_temp = Array{Union{Missing, Float64}}(missing, size(ilr_y))

#populate by seasons:
x =hcat(ones(size(wisp,1),1), X_wisp[:,[1]])
b = stepwise_df[1,:B_hat] #temperature only model
temp = b' * x' #B_hat is full model
Y_hat_temp[wisp,:] .= temp'

x =hcat(ones(size(ii,1),1), X_su[:,1]) #initial inexing into summer as one sample does not have all the nutrients!
b = stepwise_df[16,:B_hat]
temp = b' * x' #B_hat is full model
Y_hat_temp[su[ii],:] .= temp'

x =hcat(ones(size(fa,1),1), X_fa[:,1])
b = stepwise_df[34,:B_hat]
temp = b' * x' #B_hat is full model
Y_hat_temp[fa,:] .= temp'

#and back transform to the simplex:
yhat_simplex = invilr(Y_hat,bal)
yhat_simplex_temp = invilr(Y_hat_temp,bal)
## and save!

@save "/home/kristen/Documents/V6V8_analysis/analysis_products/oligo_ilr_regression_output2.jld2" stepwise_df fa wisp su ii X_wisp X_fa X_su Y_hat yhat_simplex bhat_fall_simplex bhat_summer_simplex bhat_wisp_simplex

## if need to add to oligo_df_zero_plus dataframe:

# oligo_df_plus.ilr1 = ilr_y[:,1]
# oligo_df_plus.ilr2 = ilr_y[:,2]
# oligo_df_plus.ilr3 = ilr_y[:,3]
# oligo_df_plus.ilr4 = ilr_y[:,4]
# oligo_df_plus.ilr5 = ilr_y[:,5]

#oligo_df_plus.yearday = dayofyear.(oligo_df_plus[!,:date])
#CSV.write("/home/kristen/Documents/V6V8_analysis/scripts/rel_abnd_analysis_scripts/matlab_scripts/olgio_df_plus_with_ilr.csv", oligo_df_plus, delim = ',')

b=convert(DataFrame,[X_su y[su[ii],:]]) #Y_fa
bnames=names(b)
for j=1:6
    rename!(b,bnames[j] => vars[j])
end

# rename!(b,:x7 => :ilr1)
# rename!(b,:x8 => :ilr2)
# rename!(b,:x9 => :ilr3)
# rename!(b,:x10 => :ilr4)
# rename!(b,:x11 => :ilr5)

rename!(b,:x7 => :O1)
rename!(b,:x8 => :O2)
rename!(b,:x9 => :O3)
rename!(b,:x10 => :O4)
rename!(b,:x11 => :O5)
rename!(b,:x12 => :O6)

CSV.write("/home/kristen/Documents/V6V8_analysis/scripts/R_scripts/su_data.csv", b, delim = ',')

# a=DataFrame(yhat_1=yhat_simplex[:,1],
# yhat_2=yhat_simplex[:,2],
# yhat_3=yhat_simplex[:,3],
# yhat_4=yhat_simplex[:,4],
# yhat_5=yhat_simplex[:,5],
# yhat_6=yhat_simplex[:,6],
# yhat_temp1=yhat_simplex_temp[:,1],
# yhat_temp2=yhat_simplex_temp[:,2],
# yhat_temp3=yhat_simplex_temp[:,3],
# yhat_temp4=yhat_simplex_temp[:,4],
# yhat_temp5=yhat_simplex_temp[:,5],
# yhat_temp6=yhat_simplex_temp[:,6],
# yhat_ilr1=Y_hat[:,1],
# yhat_ilr2=Y_hat[:,2],
# yhat_ilr3=Y_hat[:,3],
# yhat_ilr4=Y_hat[:,4],
# yhat_ilr5=Y_hat[:,5])
#
# CSV.write("/home/kristen/Documents/V6V8_analysis/scripts/rel_abnd_analysis_scripts/matlab_scripts/yhat_simplex.csv",a, delim = ',')


## If want to just models individually:

X0=X_fa, #X_wisp, X_su
Y=Y_fa # Y_su, Y_wisp
n=size(X0,1)
y_bar=mean(Y, dims =1)

intercept=ones(size(X0,1),1) #setup design matrix
X=hcat(intercept,X0[:,1]) #vars here

B_hat= inv(X' * X) * (X' * Y) #retrieve parameters

H_full = (Y' * Y) .- (B_hat' * X' * Y) #calculate lambda
H_null = Y' * Y .- (n .* (y_bar' * y_bar)) #E+H in book parlance
Λf = det(H_full) ./ det(H_null)

# Λ is distributed as q, p,n− p−1 as well as p,q,n−q−1, where:
# p is number of response variables
# q is the number of independent predictor variables
# n is the number of observations

# If Λf < Λ(p,q,n-q-1), then model is significant

#Can also estimate p-values with F or Chi-squared approximation:
## for full model:

vH = size(X,2)-1 #number of X vars, intercept excluded
vE = n-(size(X,2)-1)-1 #n - x vars -1
p = size(Y,2) #num y vars

f, df1 ,df2=approx_F(Λf,p,vE,vH)
pvalF=1-cdf(FDist(df1,df2),f)

#chi, df = approx_Chi(Λ_test,p,vE,vH)
#pvalchi = 1-cdf(Chisq(df),chi)

## parameter fits and projected model:

#If model is significant, see model estimates: populate once have run through each season:
seas_yhat=zeros(size(ilr_y))

Y_hat = B_hat' * X' #B_hat is full model
seas_yhat[ii[seas],:] = Y_hat'

# can also back transfrom fitted estimates to simplex:
# comp_yhat=invilr(seas_yhat,bal)
# comp_yhat_temp=invilr(seas_yhat_temp,bal)

#Can also back transform parameter fits into the simplex:

comp_bhat_wisp=invilr(bhat_wisp,bal) #bhat_wisp = Bhat


## correlations with phosphate:

using Statistics
using StatsBase

jj=findall(ismissing.(oligo_df_plus[!,:PO4]))
 #no nutrient values for one of the samples, exclude
ii=setdiff(collect(1:size(oligo_df_plus,1)),jj)

cor(oligo_df_plus[ii,:PO4],ilr_y[ii,5])
corspearman(convert(Array{Float64},oligo_df_plus[ii,:PO4]),ilr_y[ii,5])
corspearman(convert(Array{Float64},oligo_df_plus[ii,:PO4]),ilr_y[ii,4])

plot(x=oligo_df_plus[ii,:PO4],y=ilr_y[ii,4], Geom.point)
plot(x=oligo_df_plus[ii,:PO4],y=ilr_y[ii,5], Geom.point)

## plots for Wilk's Lambda...

#Sanity check just with interpolation....
#for:
#p=5 #number of response variables
#h=1 #number of additional variables compared to reduced - 1 in our case if want to test significance of added variable

#vE = n - q- 1 #number of samples, number of vars in full model -1

#testing for light signficance:
wispE = length(wisp) - 2 - 1
suE = length(ii)- 3 -1 #had two significant vars
faE = length(fa)- 2 -2

vE=[10:30; 40; 60; 80; 100]
cval = [0.215; 0.261; 0.303; 0.341; 0.376; 0.407; 0.436; 0.462; 0.486; 0.508; 0.529; 0.548; 0.565; 0.581; 0.596; 0.610; 0.623; 0.635; 0.647; 0.658;0.668; 0.744; 0.825; 0.867; 0.893]

l1 = layer(x=vE,y=cval,Geom.point)
l2 = layer(x=vE,y=cval,Geom.line)
l3= layer(x=[wispE], y=[stepwise_df[7,:Lambda]], Geom.point, style(default_color=RGB(0,0,0),point_size=2mm),shape=[Shape.square]) #,style(default_color=RGB(0,0,0), point_shapes=[square]))
p1=plot(l1,l2,l3)


l1 = layer(x=vE,y=cval,Geom.point)
l2 = layer(x=vE,y=cval,Geom.line)
l3= layer(x=[suE], y=[stepwise_df[27,:Lambda]], Geom.point, style(default_color=RGB(0,0,0),point_size=2mm),shape=[Shape.square]) #,style(default_color=RGB(0,0,0), point_shapes=[square]))
p2=plot(l1,l2,l3)

l1 = layer(x=vE,y=cval,Geom.point)
l2 = layer(x=vE,y=cval,Geom.line)
l3= layer(x=[faE], y=[stepwise_df[40,:Lambda]], Geom.point, style(default_color=RGB(0,0,0),point_size=2mm),shape=[Shape.square]) #,style(default_color=RGB(0,0,0), point_shapes=[square]))
p3=plot(l1,l2,l3)

Gadfly.set_default_plot_size(30cm,10cm)
gridstack([p1 p2 p3])
