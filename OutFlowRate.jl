# To add a new cell, type '# %%'
# To add a new markdown cell, type '# %% [markdown]'
# %% 
using HDF5
using Plots
using GR
using QuadGK
using Printf
using LaTeXStrings
ENV["GKSwstype"] = "100"
# %%
Hmax = 10 # maximum height from galactic plane in kpc
Hstep = 0.5
rlim = 50 # limit of radius
kmsinkpcyr = 3.1536e7/3.085677581e16 # 1 km/sec in kpc/yr
# timestep = np.linspace(0.01,0.8,80)
time = 50 # number of snapshots in 0 -- 1 Gyr
binwidthtemp = 0.15
rootdirectory = "/home/oku/SimulationData/isogal"
tempmin = 3
tempmax = 8
hotmin = 5 # Temperature to divide gas to Cold and Hot component

H = Hstep:Hstep:Hmax
logT = tempmin:binwidthtemp:tempmax

models = [
    # ["Osaka2019_isogal", "Osaka2019 (Shimizu+19)"],
    # ["ss_model/ver03311/fiducial", "fiducial", "red", "solid"],
    # ["ss_model/ver03311/gashalo", "w. gashalo", "blue", "solid"],
    # ["ss_model/ver03311/med", "high-reso", "red", "dot"],
    # ["ss_model/ver03311/med_gashalo", "high-reso w. gashalo", "blue", "dot"],
    # ["original/fiducial", "Shimizu+19", "red", "solid"],
    # ["original/gashalo", "w. gashalo", "blue", "solid"],
    # ["original/med", "high-reso", "red", "dot"],
    # ["original/med_gashalo", "high-reso w. gashalo", "blue", "dot"],
    # ["ss_model/ver05031/fiducial", "fiducial", "red", "solid"],
    # ["ss_model/ver05031/Boost2", "Boost 2", "red", "dash"],
    # ["ss_model/ver05031/Boost4", "Boost 4", "red", "dot"],
    # ["ss_model/ver05061/fiducial", "fiducial", "blue", "dot"],
    # ["ss_model/ver05061/strongESFB", "strong ESFB", "blue", "dash"],
    # ["ss_model/ver05211/logT7.5", "logT7.5", "red", "solid"],
    # ["ss_model/ver05221/fiducial", "cummurative", "green", "solid"],
    # ["ss_model/ver05221/gashalo", "entropy w. gashalo", "purple", "solid"],
    # ["ss_model/ver05291/fiducial", "stochastic", "purple", "dash"],
    # ["ss_model/ver06061/fiducial", "Equal Weight", "red", "solid"],
    # ["ss_model/ver06181/fiducial", "Kernel Weight", "green", "solid"],
    # ["ss_model/ver07081/fiducial", "Voronoi Weight", "purple", "solid"],
    ["ss_model/ver07081/gashalo", "Voronoi Weight gashalo", "purple", "dash"],
    # ["ss_model/ver07081/med", "Voronoi Weight med", "green", "solid"],
    ["ss_model/ver07081/med_gashalo", "Voronoi Weight med gashalo", "green", "dash"],
    # ["ss_model/ver06211/med_K30T0", "Voronoi Weight med_K30T0", "purple", "dot"],
    ]
# %%
function W3(r, h)
  r = abs(r)/h
  C = 8/h^3/pi
  if r > 1
      return 0
  elseif r > 1/2
      return C*2*(1-r)^3
  else
      return C*(1 - 6*r^2 + 6*r^3)
  end
end
function W5(r, h)
    r =abs(r)/h
    C = 3^9/359/pi/h^3
    if r > 1
        return 0
    elseif r > 2/3
        return C*(1 - r)^5
    elseif r > 1/3
        return C*((1 - r)^5 - 6*(2/3 - r)^5)
    else
        return C*((1 - r)^5 - 6*(2/3 - r)^5 + 15*(1/3 - r)^5)
    end
end
function crosssection(hsml, z)
    function integrand(x)
        return W5(sqrt(z^2 + x^2),hsml)*2*pi*x
        # return W3(sqrt(z^2 + x^2),hsml)*2*pi*x
    end
    return quadgk(integrand, 0, sqrt(hsml^2 - z^2))[1]
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
    dMdlogT = [[0.0 for j in 1:length(logT)] for i in 1:length(models)]
    for i in 1:length(models)
        for j in 1:length(logT)
            for k in 1:time
                dMdlogT[i][j] += Results["dM/dlogT"][i][h][k][j]/time
            end
        end
    end
    return dMdlogT
end
function timeaverage_height(Results, key)
    ret = [[0.0 for j in 1:length(H)] for i in 1:length(models)]
    for i in 1:length(models)
        for j in 1:length(H)
            for k in 1:time
                ret[i][j] += Results[key][i][j][k]/time
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
        z = abs(sph.pos[3,i] - galaxy.pos[3,1])
        dz = abs(z - H)
        r = sqrt((sph.pos[1,i]-galaxy.pos[1,1])^2 + (sph.pos[2,i]-galaxy.pos[2,1])^2)
        if dz < sph.hsml[i] && r < rlim #&& sph.temp[i] > Trange[1] && sph.temp[i] < Trange[2]
            v = (sph.vel[3,i] - galaxy.vel[3,1])*(sph.pos[3,i] - galaxy.pos[3,1])/z
            if v > 0
                S = crosssection(sph.hsml[i], dz)

                dM = sph.mass[i]*1e10*S*v*kmsinkpcyr
                ret["dM"] += dM

                dMS = sph.mass[i]*1e10*S*kmsinkpcyr
                MS += dMS

                energy = sph.mass[i]*(0.5*v^2 + sph.u[i])
                dE = energy*2e53*S*v*kmsinkpcyr
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
Results["SFR"] = [[[0.0 for k in 1:time] for j in 1:length(H)] for i in 1:length(models)]
Results["StellarMass"] = [[[0.0 for k in 1:time] for j in 1:length(H)] for i in 1:length(models)]
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
            Results["SFR"][i][k][j] = galaxy.sfr
            Results["StellarMass"][i][k][j] = galaxy.stellarmass

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
    end
end

# %%
using Plots
gr()

X = 0.01:0.01:time*0.01

Plots.scalefontsizes(1.2)

h = 8
x = H[h]
print("properties at H = $x\n")
dMdlogT = timeaverage_dMdlogT(Results, h)
Plots.plot(ylim=(1e-4,1e2), yscale=:log10, xlim=(3,8), xlabel=latexstring("\$\\log (T\\,\\,[\\mathrm{K}])\$"), ylabel=latexstring("\$ dM_{\\mathrm{out}}/d\\log T\\,\\,[M_{\\odot} \\mathrm{yr}^{-1}]\$"))
for i in 1:length(models)
Plots.plot!(logT, padding(dMdlogT[i]), label=models[i][2], linecolor=Symbol(models[i][3]), linestyle=Symbol(models[i][4]))
end
Plots.savefig("results/dMdlogT.pdf")

Plots.plot(ylabel="Mass outflow rate [Msun/yr]", xlabel="Time [Gyr]",legend=:topright, yscale=:log10, ylim=(1e-3,2e3))
for i in 1:length(models)
Plots.plot!(X, Results["dM"][i][h], label=models[i][2], linecolor=Symbol(models[i][3]), linestyle=Symbol(models[i][4]))
end
Plots.savefig("results/MassOutFlowRate.pdf")
    
Plots.plot(xlabel="Time [Gyr]", ylabel="Mass loading factor", legend=:topleft, yscale=:log10, ylim=(5e-4,1e3))
for i in 1:length(models)
Plots.plot!(X, Results["MassLoadingFactor"][i][h], label=models[i][2], linecolor=Symbol(models[i][3]), linestyle=Symbol(models[i][4]))
end
Plots.savefig("results/MassLoadingFactor.pdf")


Plots.plot(ylabel="Energy outflow rate [erg/yr]", xlabel="Time [Gyr]",legend=:bottomright, yscale=:log10, ylim=(1e44,1e49))
for i in 1:length(models)
Plots.plot!(X, Results["dE"][i][h], label=models[i][2], linecolor=Symbol(models[i][3]), linestyle=Symbol(models[i][4]))
end
Plots.savefig("results/EnergyOutFlowRate.pdf")

Plots.plot(xlabel="Time [Gyr]", ylabel="Energy loading factor", legend=:topleft, yscale=:log10, ylim=(1e-6,2e0))
for i in 1:length(models)
Plots.plot!(X, Results["EnergyLoadingFactor"][i][h], label=models[i][2], linecolor=Symbol(models[i][3]), linestyle=Symbol(models[i][4]))
end
Plots.savefig("results/EnergyLoadingFactor.pdf")

Plots.plot(ylabel="Metal mass outflow rate [Msun/yr]", xlabel="Time [Gyr]",legend=:bottomright, yscale=:log10, ylim=(1e-7, 2e1))
for i in 1:length(models)
Plots.plot!(X, Results["dMz"][i][h], label=models[i][2], linecolor=Symbol(models[i][3]), linestyle=Symbol(models[i][4]))
end
Plots.savefig("results/MetalOutFlowRate.pdf")

Plots.plot(xlabel="Time [Gyr]", ylabel="Metal loading factor", legend=:topleft, yscale=:log10, ylim=(1e-3,1e1))
for i in 1:length(models)
Plots.plot!(X, Results["MetalLoadingFactor"][i][h], label=models[i][2], linecolor=Symbol(models[i][3]), linestyle=Symbol(models[i][4]))
end
Plots.savefig("results/MetalLoadingFactor.pdf")

Plots.plot(xlabel="Time [Gyr]", ylabel="Metallicity", legend=:topleft, yscale=:log10, ylim=(1e-3,1e-1))
for i in 1:length(models)
Plots.plot!(X, Results["Z"][i][h], label=models[i][2], linecolor=Symbol(models[i][3]), linestyle=Symbol(models[i][4]))
end
Plots.savefig("results/Metallicity.pdf")

Plots.plot(ylabel="SFR [Msun/yr]", xlabel="Time [Gyr]", ylim=(5e-1,1e1))
for i in 1:length(models)
Plots.plot!(X, Results["SFR"][i][h], label=models[i][2], linecolor=Symbol(models[i][3]), linestyle=Symbol(models[i][4]))
end
Plots.savefig("results/SFR.pdf")

Plots.plot()
for i in 1:length(models)
Plots.plot!(X, Results["StellarMass"][i][h], label=models[i][2], linecolor=Symbol(models[i][3]), linestyle=Symbol(models[i][4]))
end
Plots.savefig("results/StellarMass.pdf")

eta_M = timeaverage_height(Results, "MassLoadingFactor")
eta_E = timeaverage_height(Results, "EnergyLoadingFactor")
eta_Z = timeaverage_height(Results, "MetalLoadingFactor")

eta_MCold = timeaverage_height(Results, "MassLoadingFactorCold")
eta_ECold = timeaverage_height(Results, "EnergyLoadingFactorCold")
eta_ZCold = timeaverage_height(Results, "MetalLoadingFactorCold")

eta_MHot = timeaverage_height(Results, "MassLoadingFactorHot")
eta_EHot = timeaverage_height(Results, "EnergyLoadingFactorHot")
eta_ZHot = timeaverage_height(Results, "MetalLoadingFactorHot")

Plots.plot(xlabel="H [kpc]", ylabel="Mass loading factor", legend=:topright, yscale=:log10, ylim=(5e-4,1e3))
for i in 1:length(models)
Plots.plot!(H, eta_M[i], label=models[i][2], linecolor=Symbol(models[i][3]), linestyle=Symbol(models[i][4]))
Plots.plot!(H, eta_MCold[i], label=models[i][2]*" cold", linecolor=Symbol(models[i][3]), linestyle=Symbol(models[i][4]), linealpha=0.7)
Plots.plot!(H, eta_MHot[i], label=models[i][2]*" hot", linecolor=Symbol(models[i][3]), linestyle=Symbol(models[i][4]), linealpha=0.4)
end
Plots.savefig("results/MassLoadingFactor_H.pdf")

Plots.plot(xlabel="H [kpc]", ylabel="Energy loading factor", legend=:topright, yscale=:log10, ylim=(1e-6,2e0))
for i in 1:length(models)
Plots.plot!(H, eta_E[i], label=models[i][2], linecolor=Symbol(models[i][3]), linestyle=Symbol(models[i][4]))
Plots.plot!(H, eta_ECold[i], label=models[i][2]*" cold", linecolor=Symbol(models[i][3]), linestyle=Symbol(models[i][4]), linealpha=0.7)
Plots.plot!(H, eta_EHot[i], label=models[i][2]*" hot", linecolor=Symbol(models[i][3]), linestyle=Symbol(models[i][4]), linealpha=0.4)
end
Plots.savefig("results/EnergyLoadingFactor_H.pdf")

Plots.plot(xlabel="H [kpc]", ylabel="Metal loading factor", legend=:topright, yscale=:log10, ylim=(1e-3,1e1))
for i in 1:length(models)
Plots.plot!(H, eta_Z[i], label=models[i][2], linecolor=Symbol(models[i][3]), linestyle=Symbol(models[i][4]))
Plots.plot!(H, eta_ZCold[i], label=models[i][2]*" cold", linecolor=Symbol(models[i][3]), linestyle=Symbol(models[i][4]), linealpha=0.7)
Plots.plot!(H, eta_ZHot[i], label=models[i][2]*" hot", linecolor=Symbol(models[i][3]), linestyle=Symbol(models[i][4]), linealpha=0.4)
end
Plots.savefig("results/MetalLoadingFactor_H.pdf")

V = timeaverage_height(Results, "V")
VCold = timeaverage_height(Results, "VCold")
VHot = timeaverage_height(Results, "VHot")

Plots.plot(xlabel="H [kpc]", ylabel="V [km/s]", legend=:topright, yscale=:log10, ylim=(1e1,5e2))
for i in 1:length(models)
Plots.plot!(H, V[i], label=models[i][2], linecolor=Symbol(models[i][3]), linestyle=Symbol(models[i][4]))
Plots.plot!(H, VCold[i], label=models[i][2]*" cold", linecolor=Symbol(models[i][3]), linestyle=Symbol(models[i][4]), linealpha=0.7)
Plots.plot!(H, VHot[i], label=models[i][2]*" hot", linecolor=Symbol(models[i][3]), linestyle=Symbol(models[i][4]), linealpha=0.4)
end
Plots.savefig("results/V_H.pdf")