# To add a new cell, type '# %%'
# To add a new markdown cell, type '# %% [markdown]'
# %% 
using HDF5
using Plots
using GR
# using QuadGK
using Printf
using LaTeXStrings
using Statistics
ENV["GKSwstype"] = "100"
# %%
mode = "plane" # "plane" or "spherical".
Hmax = 1 # maximum height from galactic plane in kpc
Hstep = 1
rlim = 10 # limit of radius (active only plane mode)
kmsinkpcyr = 3.1536e7/3.085677581e16 # 1 km/sec in kpc/yr
# timestep = np.linspace(0.01,0.8,80)
time = 100 # number of snapshots
tave = 50 # minimun snapshot number to use for calculating time average
snapinterval = 0.01 # interval of snapshot in Gyr
binwidthtemp = 0.15
rootdirectory = "/home/oku/SimulationData/isogal"
tempmin = 2
tempmax = 8
hotmin = 5 # Temperature to divide gas to Cold and Hot component

H = Hstep:Hstep:Hmax
logT = tempmin:binwidthtemp:tempmax

models = [
    # ["ss_model/ver07231/fiducial", "Fiducial", "green", "solid"],
    # ["ss_model/ver07231/notstochastic", "Not stochastic", "blue", "dash"],
    # ["ss_model/ver07231/thermal", "Thermal", "red", "dot"],
    # ["ss_model/ver07231/noFB", "No FB", "purple", "dashdot"],
    # ["ss_model/ver07231/med", "High-reso", "cyan", "dash"],
    # ["original/gashalo", "Osaka", "black", "dashdotdot"],
    
    ["ss_model/ver07231/NoHalo/fiducial", "Fiducial", "green", "solid"],
    ["ss_model/ver07231/NoHalo/notstochastic", "Not stochastic", "blue", "dash"],
    ["ss_model/ver07231/NoHalo/thermal", "Thermal", "red", "dot"],
    ["ss_model/ver07231/NoHalo/noFB", "No FB", "purple", "dashdot"],
        ]
# %%
# function W3(r, h)
#   r = abs(r)/h
#   C = 8/h^3/pi
#   if r > 1
#       return 0
#   elseif r > 1/2
#       return C*2*(1-r)^3
#   else
#       return C*(1 - 6*r^2 + 6*r^3)
#   end
# end
# function W5(r, h)
#     r =abs(r)/h
#     C = 3^9/359/pi/h^3
#     if r > 1
#         return 0
#     elseif r > 2/3
#         return C*(1 - r)^5
#     elseif r > 1/3
#         return C*((1 - r)^5 - 6*(2/3 - r)^5)
#     else
#         return C*((1 - r)^5 - 6*(2/3 - r)^5 + 15*(1/3 - r)^5)
#     end
# end
function crosssection(hsml, z)
    # function integrand(x)
    #     return W5(sqrt(z^2 + x^2),hsml)*2*pi*x
    #     # return W3(sqrt(z^2 + x^2),hsml)*2*pi*x
    # end
    # return quadgk(integrand, 0, sqrt(hsml^2 - z^2))[1]
    y = abs(z)
    C = 3^9/359/pi/hsml^3
    function f(a, b)
        return a*pi*hsml*C*(b - y/hsml)^6 * (y + hsml*(b - y/hsml)/7)
    end
    r = y/hsml
    if r > 1
        return 0
    elseif r > 2/3
        return f(1/3, 1)
    elseif r > 1/3
        return f(1/3, 1) + f(-2, 2/3)
    else
        return f(1/3, 1) + f(-2, 2/3) + f(5, 1/3)
    end
end
function getcenter(sph)
    x = y = z = mtot = 0
    for i in 1:length(sph.mass)
      x += sph.mass[i] * sph.pos[1, i]
      y += sph.mass[i] * sph.pos[2, i]
      z += sph.mass[i] * sph.pos[3, i]
      mtot += sph.mass[i]
    end
    x /= mtot
    y /= mtot
    z /= mtot
    return [x, y, z]
end
function getvelocity(sph)
    x = y = z = mtot = 0
    for i in 1:length(sph.mass)
      x += sph.mass[i] * sph.vel[1, i]
      y += sph.mass[i] * sph.vel[2, i]
      z += sph.mass[i] * sph.vel[3, i]
      mtot += sph.mass[i]
    end
    x /= mtot
    y /= mtot
    z /= mtot
    return [x, y, z]
end
function padding(array)
    ret = []
    for i in 1:length(array)
        value = array[i] > 0 ? array[i] : 1e-10
        push!(ret, value)
    end
    return ret
end
function timeaverage_dMdlogT(Results, h)
    dMdlogT = Dict()
    dMdlogT["mean"] = [[0.0 for j in 1:length(logT)] for i in 1:length(models)]
    dMdlogT["std"] = [[0.0 for j in 1:length(logT)] for i in 1:length(models)]
    for i in 1:length(models)
        for j in 1:length(logT)
            tmp = Array{Float32, 1}(undef, 0)
            for k in tave:time
                val = Results["dM/dlogT"][i][h][k][j]
                if val > 0
                    push!(tmp, log10(val))
                end
            end
            if length(tmp) > 0
                dMdlogT["mean"][i][j] = Statistics.mean(tmp)
                dMdlogT["std"][i][j] = Statistics.std(tmp, corrected=false)
            else
                dMdlogT["mean"][i][j] = -Inf
                dMdlogT["std"][i][j] = 0
            end
        end
    end
    return dMdlogT
end
function timeaverage_height(Results, key)
    ret = Dict()
    ret["mean"] = [Array{Float32, 1}(undef, 0) for i in 1:length(models)]
    ret["std"] = [Array{Float32, 1}(undef, 0) for i in 1:length(models)]
    ret["H"] = [Array{Float32, 1}(undef, 0) for i in 1:length(models)]
    for i in 1:length(models)
        for j in 1:length(H)
            tmp = Array{Float32, 1}(undef, 0)
            for k in tave:time
                val = Results[key][i][j][k]
                if val > 0
                    push!(tmp, log10(val))
                end
            end
            if length(tmp) > 0
                push!(ret["mean"][i], Statistics.mean(tmp))
                push!(ret["std"][i], Statistics.std(tmp, corrected=false))
                push!(ret["H"][i], H[j])
            end
        end
    end
    return ret
end

struct Sph
    pos::Array{Float32, 2}
    vel::Array{Float32, 2}
    hsml::Array{Float32, 1}
    mass::Array{Float32, 1}
    rho::Array{Float32, 1}
    temp::Array{Float32, 1}
    u::Array{Float32, 1}
    Z::Array{Float32, 1}
end
struct Galaxy
    pos::Array{Float32, 1}
    vel::Array{Float32, 1}
    sfr::Float32
    stellarmass::Float32
end

function main(H, sph, galaxy)
    ret = Dict()
    numbintemp = length(logT)
    ret["dM/dlogT"] = zeros(Int(numbintemp))
    ret["dM"] = 0.0
    ret["dE"] = 0.0
    ret["dMz"] = 0.0
    ret["dMCold"] = 0.0
    ret["dECold"] = 0.0
    ret["dMzCold"] = 0.0
    ret["dMHot"] = 0.0
    ret["dEHot"] = 0.0
    ret["dMzHot"] = 0.0
    MS = 0.0
    MSCold = 0.0
    MSHot = 0.0
    for i in 1:length(sph.mass)
        # initialization
        dz = sph.hsml[i] + 1.0
        z = 0

        if mode == "plane"
            z = abs(sph.pos[3,i] - galaxy.pos[3,1])
            dz = abs(z - H)
            r = sqrt((sph.pos[1,i]-galaxy.pos[1,1])^2 + (sph.pos[2,i]-galaxy.pos[2,1])^2)
            if r > rlim
                continue
            end
        elseif mode == "spherical"
            z = sqrt((sph.pos[1,i]-galaxy.pos[1,1])^2 + (sph.pos[2,i]-galaxy.pos[2,1])^2 + (sph.pos[3,i]-galaxy.pos[3,1])^2)
            dz = abs(z - H)
        else
            print("set mode plane or spherical")
            break
        end

        if dz < sph.hsml[i]
            v = -1 #initialization.
            if mode == "plane"
                v = (sph.vel[3,i] - galaxy.vel[3,1])*(sph.pos[3,i] - galaxy.pos[3,1])/z
            elseif mode == "spherical"
                v = 
                 ((sph.vel[1,i]-galaxy.vel[1,1])*(sph.pos[1,i]-galaxy.pos[1,1]) 
                + (sph.vel[2,i]-galaxy.vel[2,1])*(sph.pos[2,i]-galaxy.pos[2,1]) 
                + (sph.vel[3,i]-galaxy.vel[3,1])*(sph.pos[3,i]-galaxy.pos[3,1]))/z
            end
            if v > 0
                S = crosssection(sph.hsml[i], dz)

                dM = sph.mass[i]*1e10*S*v*kmsinkpcyr
                ret["dM"] += dM

                dMS = sph.mass[i]*1e10*S*kmsinkpcyr
                MS += dMS

                energy = sph.mass[i]*(0.5*v^2 + sph.u[i])
                dE = energy*2e53*S*v*kmsinkpcyr # 2e53 is a conversion factor from 10^10 Msun km^2 s^-2 to erg
                ret["dE"] += dE

                dMz = sph.Z[i]*dM
                ret["dMz"] += dMz

                if log10(sph.temp[i]) < hotmin
                    ret["dMCold"] += dM
                    ret["dECold"] += dE
                    ret["dMzCold"] += dMz
                    MSCold += dMS
                else
                    ret["dMHot"] += dM
                    ret["dEHot"] += dE
                    ret["dMzHot"] += dMz
                    MSHot += dMS
                end

                if log10(sph.temp[i]) > tempmin && log10(sph.temp[i]) < tempmax
                    ibintemp = Int(ceil((log10(sph.temp[i]) - tempmin)/binwidthtemp))
                    ret["dM/dlogT"][ibintemp] += dM/binwidthtemp
                end
            end
        end
    end
    ret["Z"] = ret["dMz"]/ret["dM"]
    ret["ZCold"] = ret["dMzCold"]/ret["dMCold"]
    ret["ZHot"] = ret["dMzHot"]/ret["dMHot"]

    ret["V"] = ret["dM"]/MS
    ret["VCold"] = ret["dMCold"]/MSCold
    ret["VHot"] = ret["dMHot"]/MSHot

    return ret
end

## main part
Results = Dict()
Results["dM"] = [[[0.0 for k in 1:time] for j in 1:length(H)] for i in 1:length(models)]
Results["dE"] = [[[0.0 for k in 1:time] for j in 1:length(H)] for i in 1:length(models)]
Results["dMz"] = [[[0.0 for k in 1:time] for j in 1:length(H)] for i in 1:length(models)]
Results["Z"] = [[[0.0 for k in 1:time] for j in 1:length(H)] for i in 1:length(models)]
Results["MassLoadingFactor"] = [[[0.0 for k in 1:time] for j in 1:length(H)] for i in 1:length(models)]
Results["EnergyLoadingFactor"] = [[[0.0 for k in 1:time] for j in 1:length(H)] for i in 1:length(models)]
Results["MetalLoadingFactor"] = [[[0.0 for k in 1:time] for j in 1:length(H)] for i in 1:length(models)]
Results["V"] = [[[0.0 for k in 1:time] for j in 1:length(H)] for i in 1:length(models)]

Results["dMCold"] = [[[0.0 for k in 1:time] for j in 1:length(H)] for i in 1:length(models)]
Results["dECold"] = [[[0.0 for k in 1:time] for j in 1:length(H)] for i in 1:length(models)]
Results["dMzCold"] = [[[0.0 for k in 1:time] for j in 1:length(H)] for i in 1:length(models)]
Results["ZCold"] = [[[0.0 for k in 1:time] for j in 1:length(H)] for i in 1:length(models)]
Results["MassLoadingFactorCold"] = [[[0.0 for k in 1:time] for j in 1:length(H)] for i in 1:length(models)]
Results["EnergyLoadingFactorCold"] = [[[0.0 for k in 1:time] for j in 1:length(H)] for i in 1:length(models)]
Results["MetalLoadingFactorCold"] = [[[0.0 for k in 1:time] for j in 1:length(H)] for i in 1:length(models)]
Results["VCold"] = [[[0.0 for k in 1:time] for j in 1:length(H)] for i in 1:length(models)]

Results["dMHot"] = [[[0.0 for k in 1:time] for j in 1:length(H)] for i in 1:length(models)]
Results["dEHot"] = [[[0.0 for k in 1:time] for j in 1:length(H)] for i in 1:length(models)]
Results["dMzHot"] = [[[0.0 for k in 1:time] for j in 1:length(H)] for i in 1:length(models)]
Results["ZHot"] = [[[0.0 for k in 1:time] for j in 1:length(H)] for i in 1:length(models)]
Results["MassLoadingFactorHot"] = [[[0.0 for k in 1:time] for j in 1:length(H)] for i in 1:length(models)]
Results["EnergyLoadingFactorHot"] = [[[0.0 for k in 1:time] for j in 1:length(H)] for i in 1:length(models)]
Results["MetalLoadingFactorHot"] = [[[0.0 for k in 1:time] for j in 1:length(H)] for i in 1:length(models)]
Results["VHot"] = [[[0.0 for k in 1:time] for j in 1:length(H)] for i in 1:length(models)]

Results["dM/dlogT"] = [[[[] for k in 1:time] for j in 1:length(H)] for i in 1:length(models)]

Results["SFR"] = [[0.0 for k in 1:time] for i in 1:length(models)]
Results["StellarMass"] = [[0.0 for k in 1:time] for i in 1:length(models)]
Results["MetalMass"] = [[0.0 for k in 1:time] for i in 1:length(models)]
for i in 1:length(models)
    @time for j in 1:time
        ## load data
        path2snapshot = @sprintf("%s/%s/snapshot_%03d/snapshot_%03d.hdf5", rootdirectory, models[i][1], j, j) 
        sfr = 0.0
        massstars = 0.0
        h5open(path2snapshot, "r") do fp
            pos = read(fp,"PartType0/Coordinates")
            vel = read(fp,"PartType0/Velocities")
            hsml = read(fp, "PartType0/SmoothingLength")
            mass = read(fp, "PartType0/Masses")
            rho = read(fp, "PartType0/Density")
            temp = read(fp, "PartType0/Temperature")
            u = read(fp, "PartType0/InternalEnergy")
            Z = read(fp, "PartType0/MetallicitySmoothed")
            sfr = read(fp, "PartType0/StarFormationRate")
            massstars = read(fp, "PartType4/Masses")
            global sph = Sph(pos, vel, hsml, mass, rho, temp, u, Z)
        end
        galpos = getcenter(sph)
        galvel = getvelocity(sph)
        galsfr = sum(sfr)
        stellarmass = sum(massstars)
        global galaxy = Galaxy(galpos, galvel, galsfr, stellarmass)
        
        for k in 1:length(H)
            h = H[k]
            ## main calculation
            ret = main(h, sph, galaxy)

            ## save 
            for l in collect(keys(ret))
                Results[l][i][k][j] = ret[l]
            end

            Results["MassLoadingFactor"][i][k][j] = ret["dM"]/galaxy.sfr
            Results["EnergyLoadingFactor"][i][k][j] = ret["dE"]/(1e49 * galaxy.sfr)
            Results["MetalLoadingFactor"][i][k][j] = ret["dMz"]/(1e-2 * galaxy.sfr)
            
            Results["MassLoadingFactorCold"][i][k][j] = ret["dMCold"]/galaxy.sfr
            Results["EnergyLoadingFactorCold"][i][k][j] = ret["dECold"]/(1e49 * galaxy.sfr)
            Results["MetalLoadingFactorCold"][i][k][j] = ret["dMzCold"]/(1e-2 * galaxy.sfr)
            
            Results["MassLoadingFactorHot"][i][k][j] = ret["dMHot"]/galaxy.sfr
            Results["EnergyLoadingFactorHot"][i][k][j] = ret["dEHot"]/(1e49 * galaxy.sfr)
            Results["MetalLoadingFactorHot"][i][k][j] = ret["dMzHot"]/(1e-2 * galaxy.sfr)
        end
        Results["SFR"][i][j] = galaxy.sfr
        Results["StellarMass"][i][j] = galaxy.stellarmass
        for k in 1:length(sph.mass)
            Results["MetalMass"][i][j] += sph.mass[k]*sph.Z[k]
        end 
    end
end

# %%
using Plots
gr()

X = snapinterval:snapinterval:time*snapinterval

Plots.scalefontsizes(1.2)

# galaxy properties
Plots.plot(ylabel="SFR [Msun/yr]", xlabel="Time [Gyr]", ylim=(0,16))
for i in 1:length(models)
Plots.plot!(X, Results["SFR"][i], label=models[i][2], linecolor=Symbol(models[i][3]), linestyle=Symbol(models[i][4]))
end
Plots.savefig("results/SFR.png")

Plots.plot()
for i in 1:length(models)
Plots.plot!(X, Results["StellarMass"][i], label=models[i][2], linecolor=Symbol(models[i][3]), linestyle=Symbol(models[i][4]))
end
Plots.savefig("results/StellarMass.png")

Plots.plot()
for i in 1:length(models)
Plots.plot!(X, Results["MetalMass"][i], label=models[i][2], linecolor=Symbol(models[i][3]), linestyle=Symbol(models[i][4]))
end
Plots.savefig("results/MetalMass.png")

h = 1
x = H[h]
print("properties at H = $x\n")
dMdlogT = timeaverage_dMdlogT(Results, h)
Plots.plot(ylim=(-4,2), xlim=(2,8), xlabel=latexstring("\$\\log\\,(T\\,\\,[\\mathrm{K}])\$"), ylabel=latexstring("\$\\log\\,(d\\dot{M}_{\\mathrm{out}}/d\\log T\\,\\,[M_{\\odot} \\mathrm{yr}^{-1} \\mathrm{dex}^{-1}])\$"))
for i in 1:length(models)
Plots.plot!(logT, dMdlogT["mean"][i], label=models[i][2], color=Symbol(models[i][3]), linestyle=Symbol(models[i][4]), ribbon=dMdlogT["std"][i], fillalpha=.3)
end
Plots.savefig("results/dMdlogT.png")

# time evolution of outflow properties
Plots.plot(ylabel="Mass outflow rate [Msun/yr]", xlabel="Time [Gyr]",legend=:topright, yscale=:log10, ylim=(1e-3,2e3))
for i in 1:length(models)
Plots.plot!(X, Results["dM"][i][h], label=models[i][2], linecolor=Symbol(models[i][3]), linestyle=Symbol(models[i][4]))
end
Plots.savefig("results/MassOutFlowRate.png")
    
Plots.plot(xlabel="Time [Gyr]", ylabel="Mass loading factor", legend=:topleft, yscale=:log10, ylim=(5e-4,1e3))
for i in 1:length(models)
Plots.plot!(X, Results["MassLoadingFactor"][i][h], label=models[i][2], linecolor=Symbol(models[i][3]), linestyle=Symbol(models[i][4]))
end
Plots.savefig("results/MassLoadingFactor.png")


Plots.plot(ylabel="Energy outflow rate [erg/yr]", xlabel="Time [Gyr]",legend=:bottomright, yscale=:log10, ylim=(1e44,1e49))
for i in 1:length(models)
Plots.plot!(X, Results["dE"][i][h], label=models[i][2], linecolor=Symbol(models[i][3]), linestyle=Symbol(models[i][4]))
end
Plots.savefig("results/EnergyOutFlowRate.png")

Plots.plot(xlabel="Time [Gyr]", ylabel="Energy loading factor", legend=:topleft, yscale=:log10, ylim=(1e-6,2e0))
for i in 1:length(models)
Plots.plot!(X, Results["EnergyLoadingFactor"][i][h], label=models[i][2], linecolor=Symbol(models[i][3]), linestyle=Symbol(models[i][4]))
end
Plots.savefig("results/EnergyLoadingFactor.png")

Plots.plot(ylabel="Metal mass outflow rate [Msun/yr]", xlabel="Time [Gyr]",legend=:bottomright, yscale=:log10, ylim=(1e-7, 2e1))
for i in 1:length(models)
Plots.plot!(X, Results["dMz"][i][h], label=models[i][2], linecolor=Symbol(models[i][3]), linestyle=Symbol(models[i][4]))
end
Plots.savefig("results/MetalOutFlowRate.png")

Plots.plot(xlabel="Time [Gyr]", ylabel="Metal loading factor", legend=:topleft, yscale=:log10, ylim=(1e-3,1e1))
for i in 1:length(models)
Plots.plot!(X, Results["MetalLoadingFactor"][i][h], label=models[i][2], linecolor=Symbol(models[i][3]), linestyle=Symbol(models[i][4]))
end
Plots.savefig("results/MetalLoadingFactor.png")

Plots.plot(xlabel="Time [Gyr]", ylabel="Metallicity", legend=:topleft, yscale=:log10, ylim=(1e-3,1e-1))
for i in 1:length(models)
Plots.plot!(X, Results["Z"][i][h], label=models[i][2], linecolor=Symbol(models[i][3]), linestyle=Symbol(models[i][4]))
end
Plots.savefig("results/Metallicity.png")

# outflow profile against height from the galaxy

XLABEL_H = latexstring("\$H\\,\\mathrm{[kpc]}\$")
XSCALE = ""
if mode == "spherical"
    XSCALE = Symbol("log10")
elseif mode == "plane"
    XSCALE = Symbol("none")
end

AllPhase = ["MassLoadingFactor", "EnergyLoadingFactor", "MetalLoadingFactor", "V", "Z"]
ColdPhase = ["MassLoadingFactorCold", "EnergyLoadingFactorCold", "MetalLoadingFactorCold", "VCold", "ZCold"]
HotPhase = ["MassLoadingFactorHot", "EnergyLoadingFactorHot", "MetalLoadingFactorHot", "VHot", "ZHot"]
function addannotate(key, xmin, ymax)
    if key in AllPhase
        Plots.annotate!([(xmin, ymax-0.2, (latexstring("\$\\mathrm{All}\$"), 16, :black, :left))])
    elseif key in ColdPhase
        Plots.annotate!([(xmin, ymax-0.2, (latexstring("\$\\mathrm{Cold}\\,(T < 10^5\\,\\mathrm{K})\$"), 16, :black, :left))])
    elseif key in HotPhase
        Plots.annotate!([(xmin, ymax-0.2, (latexstring("\$\\mathrm{Hot}\\,(T > 10^5\\,\\mathrm{K})\$"), 16, :black, :left))])
    end
end

for key in ["MassLoadingFactor", "MassLoadingFactorCold", "MassLoadingFactorHot"]
    eta = timeaverage_height(Results, key)
    ymin = -4
    ymax = 2
    Plots.plot(xlabel=XLABEL_H, xscale=XSCALE, ylabel=latexstring("\$\\log\\,\\eta_m \$"), legend=:topright, ylim=(ymin,ymax))
    addannotate(key, Hstep, ymax)
    for i in 1:length(models)
        Plots.plot!(eta["H"][i], eta["mean"][i], label=models[i][2], color=Symbol(models[i][3]), linestyle=Symbol(models[i][4]), ribbon=eta["std"][i], fillalpha=.3)
    end
    Plots.savefig("results/"*key*"_H.png")
end

for key in ["EnergyLoadingFactor", "EnergyLoadingFactorCold", "EnergyLoadingFactorHot"]
    eta = timeaverage_height(Results, key)
    ymin = -6
    ymax = 0
    Plots.plot(xlabel=XLABEL_H, xscale=XSCALE, ylabel=latexstring("\$\\log\\,\\eta_e \$"), legend=:topright, ylim=(ymin, ymax))
    addannotate(key, Hstep, ymax)
    for i in 1:length(models)
        Plots.plot!(eta["H"][i], eta["mean"][i], label=models[i][2], color=Symbol(models[i][3]), linestyle=Symbol(models[i][4]), ribbon=eta["std"][i], fillalpha=.3)
    end
    Plots.savefig("results/"*key*"_H.png")
end

for key in ["MetalLoadingFactor", "MetalLoadingFactorCold", "MetalLoadingFactorHot"]
    eta = timeaverage_height(Results, key)
    ymin = -3
    ymax = 1
    Plots.plot(xlabel=XLABEL_H, xscale=XSCALE, ylabel=latexstring("\$\\log\\,\\eta_Z \$"), legend=:topright, ylim=(ymin, ymax))
    addannotate(key, Hstep, ymax)
    for i in 1:length(models)
        Plots.plot!(eta["H"][i], eta["mean"][i], label=models[i][2], color=Symbol(models[i][3]), linestyle=Symbol(models[i][4]), ribbon=eta["std"][i], fillalpha=.3)
    end
    Plots.savefig("results/"*key*"_H.png")
end

for key in ["V", "VCold", "VHot"]
    V = timeaverage_height(Results, key)
    ymin = 0
    ymax = 3
    Plots.plot(xlabel=XLABEL_H, xscale=XSCALE, ylabel=latexstring("\$\\log\\,(V \\mathrm{[km}\\,\\mathrm{s}^{-1}\\mathrm{]}\$"), legend=:bottomright, ylim=(ymin, ymax))
    addannotate(key, Hstep, ymax)
    for i in 1:length(models)
        Plots.plot!(V["H"][i], V["mean"][i], label=models[i][2], color=Symbol(models[i][3]), linestyle=Symbol(models[i][4]), ribbon=V["std"][i], fillalpha=.3)
    end
    Plots.savefig("results/"*key*"_H.png")
end

for key in ["Z", "ZCold", "ZHot"]
    Z = timeaverage_height(Results, key)
    ymin = -3
    ymax = -1
    Plots.plot(xlabel=XLABEL_H, xscale=XSCALE, ylabel=latexstring("\$\\log\\,Z\$"), legend=:topright, ylim=(ymin, ymax))
    addannotate(key, Hstep, ymax)
    for i in 1:length(models)
        Plots.plot!(Z["H"][i], Z["mean"][i], label=models[i][2], color=Symbol(models[i][3]), linestyle=Symbol(models[i][4]), ribbon=Z["std"][i], fillalpha=.3)
    end
    Plots.savefig("results/"*key*"_H.png")
end
