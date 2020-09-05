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
stepwise_df=DataFrame(season=String[], step=Int64[], base=String[], var=String[], Lambda=Float64[], F_pval=Float64[], B_hat=Array{Union{Missing, Float64},2}[])

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

        push!(stepwise_df,[season,0,string(vars[q]), "-", Λf, pvalF, B_hat])

end

#once we have a baseline, find which var explains the most variance and add from there!
r=findall(stepwise_df[!,:season] .== season)
global qbase=argmin(stepwise_df[r,:Lambda])
global qtest=setdiff(1:6,qbase)
global step=1

while step <= 6
    global qtest
    global qbase
    global step

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
        vH = length(q) #additional vars
        vE = n-(length(qbase))-1 #n - x vars in reduced model -1
        p = size(Y,2) #num y vars

        f, df1, df2 = approx_F(Λ_test,p,vE,vH)
        #chi, df = approx_Chi(Λ_test,p,vE,vH)

        #p-values:
        pvalF=1-cdf(FDist(df1,df2),f)
        #pvalchi = 1-cdf(Chisq(df),chi)
        if length(qbase) == 1
            push!(stepwise_df,[season, step, string(vars[qbase]), string(vars[q]), Λ_test, pvalF, B_hat])
        else
            push!(stepwise_df,[season, step, join(string.(vars[qbase]),"-"), string(vars[q]), Λ_test, pvalF, B_hat])
        end

    end #end qtest loop

    r=findall((stepwise_df[!,:step] .== step) .& (stepwise_df[!,:season] .== season))

    if any(stepwise_df[r,:F_pval] .< 0.05) #alright- another variable does add significant explanation!

        qbase=vcat(qbase,qtest[argmin(stepwise_df[r,:Lambda])])
        qtest=setdiff(1:6,qbase)
        step=step .+ 1

    else
        break
    end

end

end #seaons
## Sweet!  Now, backtransform for plotting if needed:


#Parameter values back to simplex:

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



# Fitted data back to simplex:

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



## if need to add to oligo_df_zero_plus dataframe:
# oligo_df_plus.ilr1 = ilr_y[:,1]
# oligo_df_plus.ilr2 = ilr_y[:,2]
# oligo_df_plus.ilr3 = ilr_y[:,3]
# oligo_df_plus.ilr4 = ilr_y[:,4]
# oligo_df_plus.ilr5 = ilr_y[:,5]

#oligo_df_plus.yearday = dayofyear.(oligo_df_plus[!,:date])
#CSV.write("/home/kristen/Documents/V6V8_analysis/scripts/rel_abnd_analysis_scripts/matlab_scripts/olgio_df_plus_with_ilr.csv", oligo_df_plus, delim = ',')




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



























## let's get some plots up in here, up in here:
Gadfly.set_default_plot_size(45cm,40cm)

Y_hat_r = B_hat_r' * Xr'
Y_hat_r=Y_hat_r'

var=:temp
varlabel="Temperature"
##
layer1 = layer(x=oligo_df_plus[ii,var], y=ilr_y[ii,1], color = dayofyear.(oligo_df_plus[ii,:date]), style(point_size=1mm))
layer1_fit=layer(x=oligo_df_plus[ii,var],y=Y_hat_r[:,1], Geom.point, style(default_color=colorant"black"))

layer2 = layer(x=oligo_df_plus[ii,var], y=ilr_y[ii,2], color = dayofyear.(oligo_df_plus[ii,:date]), style(point_size=1mm))
layer2_fit=layer(x=oligo_df_plus[ii,var],y=Y_hat_r[:,2], Geom.point, style(default_color=colorant"black"))

layer3 = layer(x=oligo_df_plus[ii,var], y=ilr_y[ii,3], color = dayofyear.(oligo_df_plus[ii,:date]), style(point_size=1mm))
layer3_fit=layer(x=oligo_df_plus[ii,var],y=Y_hat_r[:,3], Geom.point, style(default_color=colorant"black"))

layer4 = layer(x=oligo_df_plus[ii,var], y=ilr_y[ii,4], color = dayofyear.(oligo_df_plus[ii,:date]), style(point_size=1mm))
layer4_fit=layer(x=oligo_df_plus[ii,var],y=Y_hat_r[:,4], Geom.point, style(default_color=colorant"black"))

layer5 = layer(x=oligo_df_plus[ii,var], y=ilr_y[ii,5], color = dayofyear.(oligo_df_plus[ii,:date]), style(point_size=1mm))
layer5_fit=layer(x=oligo_df_plus[ii,var],y=Y_hat_r[:,5], Geom.point, style(default_color=colorant"black"))

p1=plot(layer1,layer1_fit, Guide.ylabel(lb[1]), Guide.xlabel(varlabel))
p2=plot(layer2, layer2_fit,Guide.ylabel(lb[2]), Guide.xlabel(varlabel))
p3=plot(layer3,layer3_fit,Guide.ylabel(lb[3]), Guide.xlabel(varlabel))
p4=plot(layer4,layer4_fit,Guide.ylabel(lb[4]), Guide.xlabel(varlabel))
p5=plot(layer5,layer5_fit,Guide.ylabel(lb[5]), Guide.xlabel(varlabel))

gridstack(Union{Plot,Compose.Context}[p1 p2; p3 p4; p5 Compose.context()])


## add the new and improved season aspect:
Y_hat = B_hat' * X'
Y_hat= Y_hat'

var=:PO4
layer1 = layer(x=oligo_df_plus[ii,var], y=ilr_y[ii,1], color = dayofyear.(oligo_df_plus[ii,:date]), style(point_size=1mm))
layer1_fit=layer(x=oligo_df_plus[ii,var],y=Y_hat[:,1], Geom.point, style(default_color=colorant"black"))

layer2 = layer(x=oligo_df_plus[ii,var], y=ilr_y[ii,2], color = dayofyear.(oligo_df_plus[ii,:date]), style(point_size=1mm))
layer2_fit=layer(x=oligo_df_plus[ii,var],y=Y_hat[:,1], Geom.point, style(default_color=colorant"black"))

layer3 = layer(x=oligo_df_plus[ii,var], y=ilr_y[ii,3], color = dayofyear.(oligo_df_plus[ii,:date]), style(point_size=1mm))
layer3_fit=layer(x=oligo_df_plus[ii,var],y=Y_hat[:,3], Geom.point, style(default_color=colorant"black"))

layer4 = layer(x=oligo_df_plus[ii,var], y=ilr_y[ii,4], color = dayofyear.(oligo_df_plus[ii,:date]), style(point_size=1mm))
layer4_fit=layer(x=oligo_df_plus[ii,var],y=Y_hat[:,4], Geom.point, style(default_color=colorant"black"))

layer5 = layer(x=oligo_df_plus[ii,var], y=ilr_y[ii,5], color = dayofyear.(oligo_df_plus[ii,:date]), style(point_size=1mm))
layer5_fit=layer(x=oligo_df_plus[ii,var],y=Y_hat[:,5], Geom.point, style(default_color=colorant"black"))

p1=plot(layer1,layer1_fit,Guide.ylabel(lb[1]))
p2=plot(layer2, layer2_fit,Guide.ylabel(lb[2]))
p3=plot(layer3,layer3_fit,Guide.ylabel(lb[3]))
p4=plot(layer4,layer4_fit,Guide.ylabel(lb[4]))
p5=plot(layer5,layer5_fit,Guide.ylabel(lb[5]))

gridstack(Union{Plot,Compose.Context}[p1 p2; p3 p4; p5 Compose.context()])









##
Y_hat = B_hat' * X' #B_hat is full model
seas_yhat[ii[seas],:] = Y_hat'


##
var1=:temp
var2=:lightweek
var3=:SiOH
#p1=plot(layer1,layer1_fit,Guide.ylabel(lb[1]),Scale.color_continuous(minvalue=1,maxvalue=365,colormap=p -> get(ColorSchemes.jet1, p)))
#layer2_fit=layer(x=oligo_df_plus[ii[seas],var],y=Y_hat_r[:,2], Geom.point, style(default_color=colorant"black"))

# Gadfly.set_default_plot_size(15cm,15cm)
# plot(x=dayofyear.(oligo_df_plus[!,:date]),y=oligo_df_plus[!,:lightweek],Geom.point)

cmap=Scale.color_continuous(minvalue=1,maxvalue=365,colormap=p -> get(ColorSchemes.jet1, p))

Gadfly.set_default_plot_size(50cm,30cm)
layer1v1 = layer(x=oligo_df_plus[ii[seas],var1], y=ilr_y[ii[seas],1], color = dayofyear.(oligo_df_plus[ii[seas],:date]), style(point_size=1mm))
layer2v1 = layer(x=oligo_df_plus[ii[seas],var1], y=ilr_y[ii[seas],2], color = dayofyear.(oligo_df_plus[ii[seas],:date]), style(point_size=1mm))
layer3v1 = layer(x=oligo_df_plus[ii[seas],var1], y=ilr_y[ii[seas],3], color = dayofyear.(oligo_df_plus[ii[seas],:date]), style(point_size=1mm))
layer4v1 = layer(x=oligo_df_plus[ii[seas],var1], y=ilr_y[ii[seas],4], color = dayofyear.(oligo_df_plus[ii[seas],:date]), style(point_size=1mm))
layer5v1 = layer(x=oligo_df_plus[ii[seas],var1], y=ilr_y[ii[seas],5], color = dayofyear.(oligo_df_plus[ii[seas],:date]), style(point_size=1mm))

fit1v1 = layer(x=oligo_df_plus[ii[seas],var1], y=Y_hat[:,1], style(default_color=RGB(0,0,0),point_size=1mm))
fit2v1 = layer(x=oligo_df_plus[ii[seas],var1], y=Y_hat[:,2], style(default_color=RGB(0,0,0),point_size=1mm))
fit3v1 = layer(x=oligo_df_plus[ii[seas],var1], y=Y_hat[:,3], style(default_color=RGB(0,0,0),point_size=1mm))
fit4v1 = layer(x=oligo_df_plus[ii[seas],var1], y=Y_hat[:,4], style(default_color=RGB(0,0,0),point_size=1mm))
fit5v1 = layer(x=oligo_df_plus[ii[seas],var1], y=Y_hat[:,5], style(default_color=RGB(0,0,0),point_size=1mm))

layer1v2 = layer(x=oligo_df_plus[ii[seas],var2], y=ilr_y[ii[seas],1], color = dayofyear.(oligo_df_plus[ii[seas],:date]), style(point_size=1mm))
layer2v2 = layer(x=oligo_df_plus[ii[seas],var2], y=ilr_y[ii[seas],2], color = dayofyear.(oligo_df_plus[ii[seas],:date]), style(point_size=1mm))
layer3v2 = layer(x=oligo_df_plus[ii[seas],var2], y=ilr_y[ii[seas],3], color = dayofyear.(oligo_df_plus[ii[seas],:date]), style(point_size=1mm))
layer4v2 = layer(x=oligo_df_plus[ii[seas],var2], y=ilr_y[ii[seas],4], color = dayofyear.(oligo_df_plus[ii[seas],:date]), style(point_size=1mm))
layer5v2 = layer(x=oligo_df_plus[ii[seas],var2], y=ilr_y[ii[seas],5], color = dayofyear.(oligo_df_plus[ii[seas],:date]), style(point_size=1mm))

fit1v2 = layer(x=oligo_df_plus[ii[seas],var2], y=Y_hat[:,1], style(default_color=RGB(0,0,0),point_size=1mm))
fit2v2 = layer(x=oligo_df_plus[ii[seas],var2], y=Y_hat[:,2], style(default_color=RGB(0,0,0),point_size=1mm))
fit3v2 = layer(x=oligo_df_plus[ii[seas],var2], y=Y_hat[:,3], style(default_color=RGB(0,0,0),point_size=1mm))
fit4v2 = layer(x=oligo_df_plus[ii[seas],var2], y=Y_hat[:,4], style(default_color=RGB(0,0,0),point_size=1mm))
fit5v2 = layer(x=oligo_df_plus[ii[seas],var2], y=Y_hat[:,5], style(default_color=RGB(0,0,0),point_size=1mm))

layer1v3 = layer(x=oligo_df_plus[ii[seas],var3], y=ilr_y[ii[seas],1], color = dayofyear.(oligo_df_plus[ii[seas],:date]), style(point_size=1mm))
layer2v3 = layer(x=oligo_df_plus[ii[seas],var3], y=ilr_y[ii[seas],2], color = dayofyear.(oligo_df_plus[ii[seas],:date]), style(point_size=1mm))
layer3v3 = layer(x=oligo_df_plus[ii[seas],var3], y=ilr_y[ii[seas],3], color = dayofyear.(oligo_df_plus[ii[seas],:date]), style(point_size=1mm))
layer4v3 = layer(x=oligo_df_plus[ii[seas],var3], y=ilr_y[ii[seas],4], color = dayofyear.(oligo_df_plus[ii[seas],:date]), style(point_size=1mm))
layer5v3 = layer(x=oligo_df_plus[ii[seas],var3], y=ilr_y[ii[seas],5], color = dayofyear.(oligo_df_plus[ii[seas],:date]), style(point_size=1mm))

fit1v3 = layer(x=oligo_df_plus[ii[seas],var3], y=Y_hat[:,1], style(default_color=RGB(0,0,0),point_size=1mm))
fit2v3 = layer(x=oligo_df_plus[ii[seas],var3], y=Y_hat[:,2], style(default_color=RGB(0,0,0),point_size=1mm))
fit3v3 = layer(x=oligo_df_plus[ii[seas],var3], y=Y_hat[:,3], style(default_color=RGB(0,0,0),point_size=1mm))
fit4v3 = layer(x=oligo_df_plus[ii[seas],var3], y=Y_hat[:,4], style(default_color=RGB(0,0,0),point_size=1mm))
fit5v3 = layer(x=oligo_df_plus[ii[seas],var3], y=Y_hat[:,5], style(default_color=RGB(0,0,0),point_size=1mm))

p1v1=plot(layer1v1,fit1v1,Guide.ylabel(lb[1]),cmap)
p2v1=plot(layer2v1,fit2v1,Guide.ylabel(lb[2]),cmap)
p3v1=plot(layer3v1,fit3v1,Guide.ylabel(lb[3]),cmap)
p4v1=plot(layer4v1,fit4v1,Guide.ylabel(lb[4]),cmap)
p5v1=plot(layer5v1,fit5v1,Guide.ylabel(lb[5]),cmap)

p1v2=plot(layer1v2,fit1v2,Guide.ylabel(lb[1]),cmap)
p2v2=plot(layer2v2,fit2v2,Guide.ylabel(lb[2]),cmap)
p3v2=plot(layer3v2,fit3v2,Guide.ylabel(lb[3]),cmap)
p4v2=plot(layer4v2,fit4v2,Guide.ylabel(lb[4]),cmap)
p5v2=plot(layer5v2,fit5v2,Guide.ylabel(lb[5]),cmap)

p1v3=plot(layer1v3,fit1v3,Guide.ylabel(lb[1]),cmap)
p2v3=plot(layer2v3,fit2v3,Guide.ylabel(lb[2]),cmap)
p3v3=plot(layer3v3,fit3v3,Guide.ylabel(lb[3]),cmap)
p4v3=plot(layer4v3,fit4v3,Guide.ylabel(lb[4]),cmap)
p5v3=plot(layer5v3,fit5v3,Guide.ylabel(lb[5]),cmap)
#hstack(p1,p2,p3,p4,p5)
gridstack([p1v1 p2v1 p3v1 p4v1 p5v1; p1v2 p2v2 p3v2 p4v2 p5v2;p1v3 p2v3 p3v3 p4v3 p5v3])
#gridstack(Union{Plot,Compose.Context}[p1 p2 p3 p4 p5; pl1 pl2 pl3 pl4 pl5])

## back transform - how best to plot to show the fits?

comp_yhat=invilr(seas_yhat,bal)


##### ADDITIONAL EXAMPLE CODE AND DATASET ###########################################

# # #very nice - this code is working for examples in Chapter 10, Rencher book:
# filename="/home/kristen/Desktop/Rencher_chap10_data.csv"
# test_df=CSV.read(filename,header=1) |> DataFrame
# ##
# X=convert(Array{Float64},test_df[!,[:x1,:x2,:x3]])
# X=hcat(ones(19,1),X)
# Y=convert(Array{Float64},test_df[!,[:y1,:y2,:y3]])
#
# #estimate the B's....
# B_hat= inv(X' * X) * (X' * Y) #woo hoo - this matches the example!
#
# #now, siginifance testing:
# # Overall test of significance
# # Test for variable significance
#
# #Overall:
# # Model -> b coefficients (except b0) are zero:
# y_bar = mean(Y, dims =1) #y_bar'
# n = size(Y,1)
# H_null = Y' * Y .- (n .* (y_bar * y_bar')) #E+H in book parlance
# H_full = (Y' * Y) .- (B_hat' * X' * Y) #fitted model, E
#
# Λ = det(H_full) ./ det(H_null) #compare this value to Λ(3,3,15) (#y's, #x's, n-#x's-1)
# # #reject if below a critical value from table...0.309 from look-up table -> can reject....
# #
# #to test subset - it is full model over reduced:
# H_full = (Y' * Y) .- (B_hat' * X' * Y)
# #reduced:
# Xr=hcat(X,X[:,2].^2,X[:,3].^2,X[:,4].^2,X[:,2].*X[:,3],X[:,2].*X[:,4],X[:,3].*X[:,4])
# B_hat_r = inv(Xr' * Xr) * (Xr' * Y)
# H_r = (Y' * Y) .- (B_hat_r' * Xr' * Y)
#
# Λr = det(H_r) ./ det(H_null) #overall test of signifance
# Λ_test = Λf ./ Λr
