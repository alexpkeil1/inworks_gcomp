#!/usr/bin/env julia

#= ####################################################

File usage:

julia ~/repo/inworks_gformula/simulation/inworks_appendixS2_simulation.jl _paper

- this will run the entire simulation analysis and output a CSV with simulation 
results in the output folder

Alternative file usage: 

Run interactively in julia by opening the julia and run the line (after replacing the file path):

include("<path to file>/inworks_appendixS2_simulation.jl")

This will load the functions and then you can interactively run the analyses on the 
existing simulation results that start at line 420

=# #####################################################

using LSurvival, Random, Distributions, DataFrames, CSV, GLM
# LSurvival available at https://github.com/alexpkeil1/LSurvival.jl
# can using mean vs. cum exposure give different apparent directions of hwsb

cd(expanduser("~/repo/inworks_gformula/simulation"))

####### Simulation analysis functions ###################
"""
Helper function:
inverse logit transform

```julia
expit(0.0)
```
"""
function expit(logodds::Float64)
    1 / (1 + exp(-logodds))
end


"""
Helper function:
incrementing time-varying quantities

```julia
z = float.(5:-1:1) 
z |> display
pushoff!(z, 1.1)
z
```
"""
function pushoff!(vec, x)
    pushfirst!(vec, x)
    pop!(vec)
    nothing
end


"""
Function:
Main simulation function following Young et al (cited in main text), and adapted from SAS code given at: https://github.com/CausalInference/SamplePrograms/blob/main/SNAFTM.sas

    exposure -> death
    exposure -> employment
    employment <- frailty -> death

```julia
dat = snaftOccData(Xoshiro(), N=10)
```
"""
function snaftOccData(
    rng;
    N = 10,
    HR = 2.0,
    maxtime = 20,
    shape = 1.0, #shape for weibully distributed baseline surv times (set to 1.0 for exponential)
    scale = 1 / 0.009, #scale for weibully distributed baseline surv times
    γ = [log(0.001), log(0.75), log(0.85), log(2.5), log(3.5), log(17.5)], # parameters for the L variable: intercept, time, lm1, cxm1, v, frailty
    α = [log(0.1), log(0.9), log(2.5), log(1.01), log(1.5)], #parameters for the X variable: : intercept, time, lm1, cxm1, v
    sdx = 1.85, #standard deviation of the log X variable
    expfun = "lognorm", #conditional distribution of exposure, lognormal, normal, binary(logistic)
    maxx = 10.0,
    cxsplit = 2,
)
    @assert cxsplit < maxtime
    #set random seed if none given (to prevent function errors)
    ####allocate varaible space####
    id = zeros(Int64, Int64(N * maxtime))
    v = zeros(Int64, Int64(N * maxtime))
    l = zeros(Int64, Int64(N * maxtime))
    cl = zeros(Int64, Int64(N * maxtime))
    lm1 = zeros(Int64, Int64(N * maxtime))
    clm1 = zeros(Int64, Int64(N * maxtime))
    x = zeros(Float64, Int64(N * maxtime))
    xm1 = zeros(Float64, Int64(N * maxtime))
    cx = zeros(Float64, Int64(N * maxtime))
    cxrecent = zeros(Float64, Int64(N * maxtime))
    cxolder = zeros(Float64, Int64(N * maxtime))
    mxrecent = zeros(Float64, Int64(N * maxtime))
    y = zeros(Int64(N * maxtime))
    intime = zeros(Float64, Int64(N * maxtime))
    outtime = zeros(Float64, Int64(N * maxtime))
    inage = zeros(Float64, Int64(N * maxtime))
    outage = zeros(Float64, Int64(N * maxtime))
    U = zeros(Float64, Int64(N * maxtime))
    H = zeros(Float64, Int64(N * maxtime))

    ####iterate####
    #pointer
    ptr = 0
    logtr = log(HR) / shape #log time ratio (truth)

    for _id = 1:N
        #baseline variables
        t0 = rand(rng, Weibull(shape, scale))
        frail = Int(t0 < maxtime) # iT0 in Hernan's code, t0<30, 30 = maxtime
        _agebase = rand(rng, 20:50)
        _intime = 0
        _lm1 = 0
        _xm1 = 0
        _xm2 = 0
        _xlags = zeros(maxtime) # lags
        _cxm1 = 0
        _clm1 = 0
        _chk = 0
        _v = rand(rng, Bernoulli(0.5))

        #move forward in time until censoring or failure
        while _intime < maxtime
            ptr += 1 #move up one space in dataset
            id[ptr] = _id
            #lagged, baseline variables
            lm1[ptr] = _lm1
            xm1[ptr] = _xm1
            v[ptr] = _v
            _inage = _agebase + _intime

            pushoff!(_xlags, _xm1) # add _xm1 to beginning, drop last term
            cxrecent[ptr] = sum(_xlags[1:cxsplit])
            cxolder[ptr] = sum(_xlags[(cxsplit+1):end])
            recentl = min(cxsplit, _intime) - (_intime - _clm1) # number of years worked in recent window
            # recent mean exposure during work (gives same basic set of findings as recent cumulative exposure)
            mxrecent[ptr] = recentl > 0 ? cxrecent[ptr] ./ recentl : 0
            # discrete hazard of L (leaving work)
            pl = expit(
                γ[1] +
                γ[2] * (_inage - 43) +
                γ[3] * _clm1 +
                γ[4] * (cxolder[ptr]) +  # lag between exposure and leaving work
                γ[5] * _v +
                γ[6] * frail,
            )
            l[ptr] = _lm1 == 1 ? 0 : rand(rng, Bernoulli(pl))
            cl[ptr] = _lm1 == 1 ? _clm1 : 1 + _clm1 # cumulative time at work

            # distritubion of X
            if _lm1 == 0
                mx = α[1] + α[2] * _intime + α[3] * cl[ptr] + α[4] * _xm1 + α[5] * _v
                if expfun[1:4] == "logn"
                    x[ptr] = min(maxx, rand(rng, LogNormal(mx, sdx)))
                elseif expfun[1:4] == "norm"
                    x[ptr] = rand(rng, Normal(mx, sdx)) #default is lognormal
                elseif expfun[1:4] == "bina"
                    x[ptr] = rand(rng, Bernoulli(expit(mx)))
                end
            else
                x[ptr] = 0
            end
            cx[ptr] = x[ptr] + _cxm1
            clm1[ptr] = _clm1
            #lagged variables
            _cxm1 = cx[ptr]
            _clm1 = cl[ptr]
            _lm1 = max(_lm1, l[ptr])
            _xm2 = _xm1
            _xm1 = x[ptr]

            #structual model for failure time (cumulative baseline time used up)
            _chkl1 = _chk
            _chk += exp(logtr * cxolder[ptr]) # 

            #check for failure of unit
            if t0 <= _chk # baseline time used up: event occurs
                y[ptr] = 1
                _outtime = _intime + (t0 - _chkl1) * exp(-(logtr * cx[ptr]))
                if _outtime == _intime
                    #increment a little in case remaining time is exactly 0
                    _outtime += 1e-6
                end
                intime[ptr] = _intime
                outtime[ptr] = _outtime
                inage[ptr] = _agebase + intime[ptr]
                outage[ptr] = _agebase + outtime[ptr]
                U[ptr] = t0
                H[ptr] = _chk
                break
            else  # increment by one time point
                _outtime = _intime + 1
                intime[ptr] = _intime
                outtime[ptr] = _outtime
                inage[ptr] = _agebase + intime[ptr]
                outage[ptr] = _agebase + outtime[ptr]
                U[ptr] = t0
                H[ptr] = _chk
            end
            _intime += 1 #advance time one discrete unit
        end # whileintime < maxtime & fail = 0
    end #foriid = 1:N
    dat = [
        id[1:ptr] l[1:ptr] lm1[1:ptr] cl[1:ptr] clm1[1:ptr] x[1:ptr] cx[1:ptr] mxrecent[1:ptr] cxrecent[1:ptr] cxolder[1:ptr] v[1:ptr] y[1:ptr] intime[1:ptr] outtime[1:ptr] inage[1:ptr] outage[1:ptr] U[1:ptr] H[1:ptr]
    ]
    DataFrame(
        dat,
        [
            :id,
            :l,
            :everlm1,
            :cl,
            :clm1,
            :x,
            :cx,
            :mxrecent,
            :cxrecent,
            :cxolder,
            :v,
            :y,
            :intime,
            :outtime,
            :inage,
            :outage,
            :U,
            :H,
        ],
    )

end 
snaftOccData(; kwargs...) = snaftOccData(Xoshiro(); kwargs...)

"""
Helper function:
lag a variable

```julia
id = [1,1,2,2,2]
x = rand(Xoshiro(10), length(id))
xlag = lag(id,x, default_value = 0.0)
hcat(id,x,xlag)
```
"""
function lag(id, variable; default_value = 1.0)
    idx = [findall(id .== i) for i in unique(id)]
    varlag = [vcat(default_value, variable[idlist[1:(end-1)]]) for idlist in idx]
    reduce(vcat, varlag)
end



"""
Function:
Simulate one dataset using parameters from the paper
```julia
simdata(Xoshiro(10), 10)
```
"""
function simdata(
    rng,
    N;
    γ = [log(0.1), log(1.25), log(1.25), log(1.0), log(2.5), log(5.5)],
    α = [log(0.1), log(1.8), log(1.5), log(1.1), log(2.5)],
    HR = 1.3,
)
    df = snaftOccData(
        rng;
        N = N,
        HR = HR,
        maxtime = 20,
        # mean(rand(Weibull(1, 1/0.009), 1000) .> 20) # expected number of those who survive to end of follow-up under no exposure
        shape = 1.0, #shape for weibully distributed baseline surv times (set to 1.0 for exponential)
        scale = 1 / 0.009, #scale for weibully distributed baseline surv times
        # parameters for the L variable: intercept, time, clm1, cxm1, v, frailty
        γ = γ,
        #parameters for the X variable: : intercept, time, clm1, cxm1, v
        α = α,
        sdx = 0.15, #standard deviation of the log X variable
        expfun = "lognormal", #conditional distribution of exposure, lognormal, normal, binary(logistic)
    )
    df.mx = df.cx ./ df.cl
    df
end
simdata(N; kwargs...) = simdata(Xoshiro(), N; kwargs...)


"""
Function:
Perform the full simulation for one realized dataset


```julia
b1_unadj,b1, b2, b1_unadj_signif, b1_signif, b2_signif = onerep(Xoshiro(10), 2000)
```
"""
function onerep(rng, N; kwargs...)
    df = simdata(rng, N; kwargs...)
    cr0 = glm(
        @formula(l ~ cxrecent + v + cl + inage),
        filter(x -> (x.everlm1 .== 0), df),
        Binomial(),
        contrasts = Dict(:inage => DummyCoding()),
    )
    cr1 = glm(
        @formula(l ~ cxrecent + cxolder + v + cl + inage),
        filter(x -> (x.everlm1 .== 0), df),
        Binomial(),
        contrasts = Dict(:inage => DummyCoding()),
    )
    hasint = 1
    rr0 = coeftable(cr0)
    rr1 = coeftable(cr1)
    permutedims(
        vcat(
            vcat(coef(cr0)[1+hasint], coef(cr1)[(1:2) .+ hasint]),
            (
                vcat(
                    rr0.cols[rr0.pvalcol][1+hasint],
                    rr1.cols[rr1.pvalcol][(1:2) .+ hasint],
                ) .< 0.05
            ),
        ),
    )
end
onerep(N; kwargs...) = onerep(Xoshiro(), N; kwargs...)

"""
Function:
Repeat the simulation generation/analysis over <iter> iterations (a "batch") and take the mean across the batch
"""
function repit(rng, iter; kwargs...)
    nms = [
        :agerecent, # b1_unadj
        :agerecent_adj, # b1
        :ageolder_adj, # b2
        :p_agerecent, # b1_unadj_power
        :p_agerecent_adj, # b1_power
        :p_ageolder_adj, # b2_power
    ]
    res = reduce(vcat, [onerep(rng, 2000; kwargs...) for i = 1:iter])
    DataFrame(res, nms), DataFrame(mean(res, dims = 1), nms)
end
repit(iter; kwargs...) = repit(Xoshiro(), iter; kwargs...)

"""
Function:
Repeat the batched analysis <iter> times (note this iter parameter is not equilent to the iter parameter in `repit` function)
"""
function repmore(rng, iter; fl = expanduser("temp/occsims.csv"))
    firstdone = false
    for i = 1:iter
        # parameters for the L variable: intercept, time, clm1, cx(6), v, frailty
        gvals = [
            log(0.1),
            log(1.25),
            abs.(randn(rng) * 1.0),
            abs(randn(rng) * 2.0),
            log(2.5),
            rand(rng) * 20.0,
        ]
        #parameters for the X variable: : intercept, time, clm1, cxm1, v 
        avals =
            [log(0.1), log(1.8), -abs.(randn(rng) * 0.5), abs(randn(rng) * 0.40), log(2.5)]
        hr = 1.0 + rand(rng) * 1.0
        res, ressum = try
            repit(rng, 100, γ = gvals, α = avals, HR = hr)
        catch
            nms = [
                :agerecent,
                :agerecent_adj,
                :ageolder_adj,
                :p_agerecent,
                :p_agerecent_adj,
                :p_ageolder_adj,
            ]
            "error", DataFrame(permutedims(["" for i = 1:6]), nms)
        end
        if !isnothing(res)
            parms = DataFrame(
                vcat([i, hr], gvals, avals)',
                vcat(
                    ["iter", "HR"],
                    ["γ$i" for i = 1:length(gvals)],
                    ["α$i" for i = 1:length(avals)],
                ),
            )
            CSV.write(fl, hcat(ressum, parms), append = firstdone)
            firstdone = true
        end
    end
end
repmore(iter; kwargs...) = repmore(Xoshiro(), iter; kwargs...)

###########################
# Running simulations
# programming note: this will run automatically if you just enter this program into a terminal
# if you run this interactively in a Julia session, you will need to run this portion manually
###########################
if length(ARGS) > 0
    suffix = String(ARGS[1]) # this makes it easy to run multiple times to check for consistency (need to modify the seed, however)
    # run the whole simulation - to change the seed, change "Xoshiro(1000)" to either "Xoshrio()" (no seed) or set the "1000" to another integer value for a set seed that allows replicability
    cd(expanduser("~/repo/inworks_gformula/simulation"))
    repmore(Xoshiro(1000), 10000, fl = expanduser("output/inwork_sims$suffix.csv")) # save results to an external file (this program takes a while)
end

###########################
# Analyzing simulated data
# programming note: this will run automatically if you just enter this program into a terminal
# if you run this interactively in a Julia session, you will need to run this portion manually
###########################
if length(ARGS) > 0
    using DataFrames, CSV, GLM, Loess, LaTeXStrings, Plots, StatsBase
    #suffix = String(ARGS[1]) # may need to set this explicitly if you run this 
    suffix = "_paper"
    cd(expanduser("~/repo/inworks_gformula/simulation"))
    df = CSV.read(expanduser("output/inwork_sims$suffix.csv"), DataFrame) # seed value used in paper
    df.ψ = log.(df.HR)

    # renaming parameters to fit paper's conventions (zero indexed intercept)
    rndict =
        merge(Dict("γ$i" => "γ$(i-1)" for i = 1:6), Dict("α$i" => "α$(i-1)" for i = 1:5))
    rename!(df, rndict)

    # final analytic set
    dflimit = df[
        :,
        vcat(
            [:agerecent_adj, :agerecent_adjopp, :ψ],
            [Symbol("γ$i") for i = 0:5],
            [Symbol("α$i") for i = 0:4],
        ),
    ]

    nm = names(dflimit)
    vals = reduce(
        vcat,
        [[median(cl), quantile(cl, 0.05), quantile(cl, 0.95)] for cl in eachcol(dflimit)]',
    )
    # table S5 in paper
    DataFrame(nm = nm, median = vals[:, 1], q5 = vals[:, 2], q95 = vals[:, 3])

    # now standardizing all predictors for interpretability
    stdcols = [(col .- mean(col)) / std(col) for col in eachcol(dflimit)]
    cencols = [(col .- mean(col)) for col in eachcol(dflimit)]
    stddflimit = DataFrame(reduce(hcat, stdcols), names(dflimit) .* "_std")
    cendflimit = DataFrame(reduce(hcat, cencols), names(dflimit) .* "_cen")
    dflimit2 = hcat(dflimit, stddflimit, cendflimit)

    # table S6 in paper
    res = lm(
        @formula(agerecent_adj ~ γ2_std + γ3_std + γ5_std + α2_std + α3_std + ψ_std),
        dflimit2,
    )

    #  Coefficients:
    #  ──────────────────────────────────────────────────────────────────────────────
    #                    Coef.  Std. Error        t  Pr(>|t|)   Lower 95%   Upper 95%
    #  ──────────────────────────────────────────────────────────────────────────────
    #  (Intercept)  -1.68761     0.014437   -116.89    <1e-99  -1.71591    -1.65931
    #  γ2_std        0.541404    0.0144415    37.49    <1e-99   0.513096    0.569713
    #  γ3_std        0.0162696   0.0144399     1.13    0.2599  -0.0120355   0.0445746
    #  γ5_std       -0.610195    0.0144396   -42.26    <1e-99  -0.638499   -0.58189
    #  α2_std        2.07213     0.0144416   143.48    <1e-99   2.04382     2.10044
    #  α3_std        0.0713572   0.0144394     4.94    <1e-06   0.043053    0.0996614
    #  ψ_std         0.0037981   0.0144425     0.26    0.7926  -0.0245121   0.0321083
    #  ──────────────────────────────────────────────────────────────────────────────

    #### Plots with marginal associations
    import Plots.Measures.mm
    Plots.default(
        size = [550, 430],
        fontfamily = "Times New Roman",
        markerstrokewidth = 0,
        markercolor = :black,
        markersize = 1.5,
        markeralpha = 0.5,
        bottom_margin = 1mm,
        left_margin = 3mm,
        tickfontsize = 10,
        guidefontvalign = :center,
        guidefonthalign = :center,
        tickfontvalign = :vcenter,
    )

    ## plot S3: α2 -> β1
    # loess scatterplot smoother
    lft = loess(df.α2, df.agerecent_adj)
    ord = sortperm(lft.xs[:, 1])
    x = lft.xs[ord]
    y = predict(lft)[ord]

    p = plot(
        dflimit.α2,
        dflimit.agerecent_adj,
        st = :scatter,
        ylab = "β\$_1\$: Recent exposure -> employment",
        xlab = "\nTrue α\$_2\$: # years employed -> exposure",
        label = "Simulated cohort estimate",
    ) # if exposure decreases with tenure, then you'll see negative effects of recent exposure
    plot!(p, x, y, st = :line, linewidth = 3, label = "LOESS fit")
    savefig(p, expanduser("fig/fig_S3.pdf"))


    ## plot S4: γ5 -> β1
    # loess scatterplot smoother
    lft2 = loess(df.γ5, df.agerecent_adj)
    ord2 = sortperm(lft2.xs[:, 1])
    x2 = lft2.xs[ord2]
    y2 = predict(lft2)[ord2]

    p2 = plot(
        dflimit.γ5,
        dflimit.agerecent_adj,
        st = :scatter,
        ylab = "β\$_1\$: Recent exposure -> employment",
        xlab = "\nTrue γ\$_5\$: # Frailty -> employment",
        label = "Simulated cohort estimate",
    )
    plot!(p2, x2, y2, st = :line, linewidth = 3, label = "LOESS fit")
    savefig(p2, expanduser("fig/fig_S4.pdf"))


    ## plot S5: γ2 -> β1
    # loess scatterplot smoother
    lft3 = loess(df.γ2, df.agerecent_adj)
    ord3 = sortperm(lft3.xs[:, 1])
    x3 = lft3.xs[ord3]
    y3 = predict(lft3)[ord3]

    p3 = plot(
        dflimit.γ2,
        dflimit.agerecent_adj,
        st = :scatter,
        ylab = "β\$_1\$: Recent exposure -> employment",
        xlab = "\nTrue γ\$_2\$: # Employment duration -> employment",
        label = "Simulated cohort estimate",
    )
    plot!(p3, x3, y3, st = :line, linewidth = 3, label = "LOESS fit")
    savefig(p3, expanduser("fig/fig_S5.pdf"))

    ## plot S6: γ3 -> β1
    # loess scatterplot smoother
    lft4 = loess(df.γ3, df.agerecent_adj)
    ord4 = sortperm(lft4.xs[:, 1])
    x4 = lft4.xs[ord4]
    y4 = predict(lft4)[ord4]

    p4 = plot(
        dflimit.γ3,
        dflimit.agerecent_adj,
        st = :scatter,
        ylab = "β\$_1\$: Recent exposure -> employment",
        xlab = "\nTrue γ\$_3\$: # Lagged exposure -> employment",
        label = "Simulated cohort estimate",
    )
    plot!(p4, x4, y4, st = :line, linewidth = 3, label = "LOESS fit")
    savefig(p4, expanduser("fig/fig_S6.pdf"))
end
