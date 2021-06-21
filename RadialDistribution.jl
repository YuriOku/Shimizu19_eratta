# To add a new cell, type '# %%'
# To add a new markdown cell, type '# %% [markdown]'
# %% 
using HDF5
using Plots
using Plots.PlotMeasures
using GR
using QuadGK
using Printf
using LaTeXStrings
using DelimitedFiles
ENV["GKSwstype"] = "100"
gr()
# %%
H = 0.25 # height from galactic plane in kpc
rlim = 40 # limit of radius
kmsinkpcyr = 3.1536e7/3.085677581e16 # 1 km/sec in kpc/yr
# timestep = np.linspace(0.01,0.8,80)
time = 300 # snapshot
binsize = 0.01
rootdirectory = "/home/oku/SimulationData/isogal"

models = [
    # ["Osaka2019_isogal", "Osaka2019 (Shimizu+19)"],
        #   , "geodome_model/geodome_original"\
    # ["geodome_model/ver_19.11.1", "Geodesic dome model & Cioffi+ 1988"],
    # ["geodome_model/OKU2020","Geodesic dome model & Athena fitting"],
    # ["centroid_model/ver07271_NoMomCeiling","Centroid model & Athena fitting (alpha = 0)"],
    # ["centroid_model/ver07272_nfb1","Centroid model & Athena fitting (nfb = 1)"],
    # ["centroid_model/ver07272_SFE001","Centroid model & Athena fitting (SFE = 0.01)"],
    # ["centroid_model/ver07311","Centroid model & Athena fitting (alpha = -1)"],
    # ["centroid_model/ver07311_fdens-2","Centroid model & Athena fitting (alpha = -2)"],
    # ["centroid_model/ver08041_alpha-1","Centroid model & Athena fitting (new, alpha = -1)"],
    # ["centroid_model/ver08041_alpha-2","Centroid model & Athena fitting (new, alpha = -2)"],
    # ["centroid_model/ver08041_NoThermal","Centroid model (alpha = -1)"],
    # ["centroid_model/ver08191","Centroid model (alpha = -0.24)"],
    # ["centroid_model/ver07272_CHEVALIER1974","Centroid model & Cioffi+ 1988"],
    # ["centroid_model/ver09011_n1f2a0","Centroid model"],
    # ["centroid_model/ver09011_n1f2a024","Centroid model (new, alpha = -0.24)"],
    # ["centroid_model/ver09011_n1f2a050","Centroid model (new, alpha = -0.50)"],
    # ["centroid_model/ver09011_n1f2a075","Centroid model (new, alpha = -0.75)"],
    # ["centroid_model/ver09011_n1f2a1","Centroid model (new, alpha = -1.00)"],
    # ["centroid_model/ver12151","Centroid"],
    # ["ss_model/ver01092", "Spherical superbubble"],
    # ["ss_model/ver01121", "Spherical superbubble"],
    # ["ss_model/ver01193", "Spherical superbubble"],
    # ["ss_model/ver03221", "Spherical superbubble"],
    # ["ss_model/ver12211", "Spherical superbubble"],
    # ["ss_model/ver03311", "with cool-off"],
    ["ss_model/ver03311/fiducial", "fiducial"],
    # ["ss_model/ver03311/med", "high-reso"],
    # ["ss_model/ver03311_NoCoolOff", "without cool-off"],
    # ["NoFB", "No feedback"],
    ]
# %%
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

function main(galpos)
  RDF = zeros(Int(ceil(rlim/binsize)))
  for i in 1:length(star.rshock)
    xi = star.pos[1,i] - galpos[1]
    yi = star.pos[2,i] - galpos[2]
    zi = star.pos[3,i] - galpos[3]
    for j in i+1:length(star.rshock)  
      xj = star.pos[1,j] - galpos[1]
      yj = star.pos[2,j] - galpos[2]
      zj = star.pos[3,j] - galpos[3]

      r2 = (xi - xj)^2 + (yi - yj)^2 + (zi - zj)^2
      r = sqrt(r2)
      index = Int(ceil(r/binsize))
      if index == 0
        @show r, r/binsize, Int(ceil(r/binsize)), i, j
      end
      RDF[index] += 1/(4*pi*r^2)
    end
  end
  norm = sum(RDF)
  RDF = RDF/norm/binsize
  @show RDF
  return RDF
end
struct Star
    pos::Array{Float32, 2}
    rshock::Array{Float32, 1}
    sft::Array{Float32, 1}
    sniis::Array{Float32, 1}
    sniie::Array{Float32, 1}
end
struct Sph
    pos::Array{Float32, 2}
    mass::Array{Float32, 1}
end
# %%
Volume = [[[0.0 for j in 1:time] for k in [1, 2]] for i in 1:length(models)]
for j in 1:time
  for i in 1:length(models)
        path2snapshot = @sprintf("%s/%s/snapshot_%03d/snapshot_%03d.hdf5", rootdirectory, models[i][1], j, j) 
        h5open(path2snapshot, "r") do fp
            posstar = read(fp,"PartType4/Coordinates")
            possph = read(fp,"PartType0/Coordinates")
            mass = read(fp,"PartType0/Masses")
            rshock = read(fp, "PartType4/Rshock")
            age = read(fp, "PartType4/StellarFormationTime")
            sniis = read(fp, "PartType4/SNIIStartTime")
            sniie = read(fp, "PartType4/SNIIENDTime")
            global star = Star(posstar, rshock, age, sniis, sniie)
            global sph = Sph(possph, mass)
        end
        galpos = getcenter(sph)
        
        result = main(galpos)
        Plots.plot(result, xlim=(1e-1, 1e2), ylim=(1e-3, 1e2), xscale=:log10, yscale=:log10)
        Plots.savefig("tmp.png")
  end
  @show j
end

gr()
font = Plots.font(16)
lfont = Plots.font(12)
Plots.plot(legend=:bottomleft,xlim=(0, time), ylim=(1e-8, 1e-1), yscale=:log10, xlabel="Time [Gyr]", ylabel=L"f_{\rm hot}", xtickfont=font, ytickfont=font, xlabelfont=font, ylabelfont=font, legendfont=lfont, guidefont=font, topmargin=3mm, bottommargin=3mm, rightmargin=3mm)
for i in 1:length(models)
Plots.plot!(Volume[i][1], label="hot phase", markersize=5, markercolor=:red)
Plots.plot!(Volume[i][2], label="cold phase", markersize=5, markercolor=:blue)
end
Plots.savefig("Energyloading.png")
# Plots.savefig("results/Energyloading.pdf")

