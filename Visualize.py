# %%
import yt
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid

time = 99

models = [
    # ["../ss_model/ver03311/fiducial","fiducial"],
    # ["../ss_model/ver03311/gashalo","with gashalo"],
    # ["../ss_model/ver03311/med","high-reso"],
    # ["../ss_model/ver03311/med_gashalo","high-reso with gashalo"],
    # ["../original/fiducial","fiducial"],
    # ["../original/gashalo","with gashalo"],
    # ["../original/med","high-reso"],
    # ["../original/med_gashalo","high-reso with gashalo"],
    # ["../ss_model/ver05031/fiducial","fiducial"],
    # ["../ss_model/ver05031/Boost2","Boost 2"],
    # ["../ss_model/ver05031/Boost4","Boost 4"],
    ["../ss_model/ver05061/fiducial","fiducial"],
    # ["../ss_model/ver05061/strongESFB","strong ESFB"],
    # ["../ss_model/ver05131/fiducial","thermal storage"],
    # ["../ss_model/ver05131/K0T100","K0T100"],
    # ["../ss_model/ver05211/fiducial","logT = 6.7"],
    # ["../ss_model/ver05211/logT7.5","logT = 7.5"],
    ["../ss_model/ver05221/fiducial","cummurative"],
    # ["../ss_model/ver05221/gashalo","entropy w. gashalo"],
    ["../ss_model/ver05291/fiducial","stochastic"],
    ["../ss_model/ver06061/fiducial","stochastic2"],
          ]
unit_base = {'UnitLength_in_cm'         : 3.08568e+21,
             'UnitMass_in_g'            :   1.989e+43,
             'UnitVelocity_in_cm_per_s' :        1e+5,
             'metallicity'              :   (1/0.0134, 'Zsun')
             }

bbox_lim = 1e4 #kpc

bbox = [[-bbox_lim,bbox_lim],
        [-bbox_lim,bbox_lim],
        [-bbox_lim,bbox_lim]]


ds = [0]*len(models)
center = [0]*len(models)
for i in range(len(models)):
    ds[i] = yt.load('{0}/snapshot_{1:03}/snapshot_{1:03}.hdf5'.format(models[i][0], time),unit_base=unit_base,bounding_box=bbox)
    # ad = ds[0].all_data()
    # density = ad[("PartType0","density")]
    # wdens = np.where(density == np.max(density))
    # coordinates = ad[("PartType0","Coordinates")]
    # mass_centroid = ad.quantities.weighted_average_quantity("Coordinates", "Masses")
    sp = ds[i].sphere("max", (40, 'kpc'))
    center[i] = sp.quantities.center_of_mass()
# print ('center = ',center)

# %% density
# field = "density"
# for i in range(len(models)):
#   px = yt.OffAxisProjectionPlot(ds[i], [1,0,0], ('gas', '{}'.format(field)), center[i], width=45, north_vector=[0,0,1])
#   px.set_xlabel("x")
#   px.set_ylabel("y")
#   # px = yt.ProjectionPlot(ds[i], 'z', ('gas', '{}'.format(field)), center = center[i], width=45)
#   px.set_font({'size':32})
#   # px.set_zlim(field="{}".format(field), zmin=5e-6, zmax=1e-1)
#   # px.hide_colorbar()
#   # px.hide_axes(field='x')
#   px.show()
#   # px.save("results/density-{}.pdf".format(models[i][1]))
# # %% metallicity
# field = "density"
# for i in range(len(models)):
#   px = yt.SlicePlot(ds[i], 'x', ('gas', '{}'.format(field)), center = center[i], width=400)
#   px.set_cmap(field="{}".format(field), cmap='hot')
#   # px.set_zlim(field="{}".format(field), zmin=1e-5, zmax=1e-1)
#   px.set_font({'size':32})
#   px.show()

# %%
sorted(ds[0].field_list)
# %% density
# ----------
field = "density"
fig = plt.figure()
grid = AxesGrid(fig, (0.075,0.075,0.85,0.85),
                nrows_ncols = (2, len(models)),
                axes_pad = 0.05,
                label_mode = "L",
                share_all = True,
                cbar_location="right",
                cbar_mode="single",
                cbar_size="5%",
                cbar_pad="1%")

for i in range(len(models)):
  px = yt.OffAxisProjectionPlot(ds[i], [0,0,1], ('gas', '{}'.format(field)), center = center[i], width=45, north_vector=[0,1,0])
  px.set_xlabel("x [kpc]")
  px.set_ylabel("y [kpc]")
  px.set_font({'size':20})
  px.set_zlim(field="{}".format(field), zmin=5e-6, zmax=1e-1)
  px.annotate_text([-20, 18], models[i][1], coord_system="plot")
  plot = px.plots['{}'.format(field)]
  plot.figure = fig
  plot.axes = grid[i].axes
  plot.cax = grid.cbar_axes[i]
  px.set_background_color(field)
  px._setup_plots()
for i in range(len(models)):
  px = yt.OffAxisProjectionPlot(ds[i], [0,1,0], ('gas', '{}'.format(field)), center = center[i], width=45, north_vector=[0,0,-1])
  px.set_xlabel("x [kpc]")
  px.set_ylabel("z [kpc]")
  px.set_font({'size':20})
  px.set_zlim(field="{}".format(field), zmin=5e-6, zmax=1e-1)
  plot = px.plots['{}'.format(field)]
  plot.figure = fig
  plot.axes = grid[i+len(models)].axes
  plot.cax = grid.cbar_axes[i+len(models)]
  px.set_background_color(field)
  px._setup_plots()
# plt.savefig("results/plot{}.pdf".format(field), bbox_inches="tight")
plt.show()
plt.close()
# %% temperature 20kpc
# ---------------------
field = "temperature"
fig = plt.figure()
grid = AxesGrid(fig, (0.075,0.075,0.85,0.85),
                nrows_ncols = (2, len(models)),
                axes_pad = 0.05,
                label_mode = "L",
                share_all = True,
                cbar_location="right",
                cbar_mode="single",
                cbar_size="5%",
                cbar_pad="1%")

for i in range(len(models)):
  px = yt.OffAxisProjectionPlot(ds[i], [0,0,1], ('gas', '{}'.format(field)), center = center[i], width=45, north_vector=[0,1,0], weight_field='density')
  px.set_xlabel("x [kpc]")
  px.set_ylabel("y [kpc]")
  px.set_font({'size':20})
  px.set_zlim(field="{}".format(field), zmin=1e3, zmax=1e7)
  px.set_cmap(field="{}".format(field), cmap='hot')
  px.annotate_text([-20, 18], models[i][1], coord_system="plot")
  plot = px.plots['{}'.format(field)]
  plot.figure = fig
  plot.axes = grid[i].axes
  plot.cax = grid.cbar_axes[i]
  px.set_background_color(field)
  px._setup_plots()
for i in range(len(models)):
  px = yt.OffAxisProjectionPlot(ds[i], [0,1,0], ('gas', '{}'.format(field)), center = center[i], width=45, north_vector=[0,0,-1], weight_field='density')
  px.set_xlabel("x [kpc]")
  px.set_ylabel("z [kpc]")
  px.set_font({'size':20})
  px.set_zlim(field="{}".format(field), zmin=1e3, zmax=1e7)
  px.set_cmap(field="{}".format(field), cmap='hot')
  plot = px.plots['{}'.format(field)]
  plot.figure = fig
  plot.axes = grid[i+len(models)].axes
  plot.cax = grid.cbar_axes[i+len(models)]
  px.set_background_color(field)
  px._setup_plots()
# plt.savefig("results/plot{}20kpc.pdf".format(field), bbox_inches="tight")
plt.show()


# %% temperature 200kpc
# ---------------------
field = "temperature"
fig = plt.figure()
grid = AxesGrid(fig, (0.075,0.075,0.85,0.85),
                nrows_ncols = (1, len(models)),
                axes_pad = 0.05,
                label_mode = "L",
                share_all = True,
                cbar_location="right",
                cbar_mode="single",
                cbar_size="5%",
                cbar_pad="1%")

for i in range(len(models)):
  px = yt.OffAxisProjectionPlot(ds[i], [0,1,0], ('gas', '{}'.format(field)), center = center[i], width=500, north_vector=[0,0,-1], weight_field='density')
  px.set_xlabel("x [kpc]")
  px.set_ylabel("z [kpc]")
  px.set_font({'size':20})
  px.set_zlim(field="{}".format(field), zmin=1e3, zmax=1e7)
  px.set_cmap(field="{}".format(field), cmap='hot')
  plot = px.plots['{}'.format(field)]
  plot.figure = fig
  plot.axes = grid[i].axes
  plot.cax = grid.cbar_axes[i]
  px.set_background_color(field)
  px._setup_plots()
# plt.savefig("results/plot{}200kpc.pdf".format(field), bbox_inches="tight")
plt.show()



# %% metallicity 200kpc
# ---------------------
field = "metallicity"
fig = plt.figure()
grid = AxesGrid(fig, (0.075,0.075,0.85,0.85),
                nrows_ncols = (1, len(models)),
                axes_pad = 0.05,
                label_mode = "L",
                share_all = True,
                cbar_location="right",
                cbar_mode="single",
                cbar_size="5%",
                cbar_pad="1%")

for i in range(len(models)):
  # px = yt.SlicePlot(ds[i], [0,1,0], ('gas', '{}'.format(field)), center = center[i], width=500, north_vector=[0,0,-1])
  px = yt.OffAxisProjectionPlot(ds[i], [0,1,0], ('gas', '{}'.format(field)), center = center[i], width=500, north_vector=[0,0,-1], weight_field='metallicity')
  px.set_xlabel("x [kpc]")
  px.set_ylabel("z [kpc]")
  px.set_font({'size':20})
  # px.set_zlim(field , zmin=-1e3, zmax=1e3)
  px.set_zlim(field , zmin=1e-3, zmax=1e2)
  px.set_cmap(field , cmap='Haze')
  px.set_unit(field, 'Zsun')
  px.set_background_color('{}'.format(field))
  plot = px.plots['{}'.format(field)]
  plot.figure = fig
  plot.axes = grid[i].axes
  plot.cax = grid.cbar_axes[i]
  px.set_background_color(field)
  px._setup_plots()
# plt.savefig("results/plot{}200kpc.pdf".format(field), bbox_inches="tight")
plt.show()
# %% velocity 200kpc
# ---------------------
field = "velocity_z"
fig = plt.figure()
grid = AxesGrid(fig, (0.075,0.075,0.85,0.85),
                nrows_ncols = (1, len(models)),
                axes_pad = 0.05,
                label_mode = "L",
                share_all = True,
                cbar_location="right",
                cbar_mode="single",
                cbar_size="5%",
                cbar_pad="1%")

for i in range(len(models)):
  # px = yt.SlicePlot(ds[i], [0,1,0], ('gas', '{}'.format(field)), center = center[i], width=500, north_vector=[0,0,-1])
  px = yt.SlicePlot(ds[i], "x", ('gas', '{}'.format(field)), center = center[i], width=500, north_vector=[0,0,1])
  px.set_xlabel("y [kpc]")
  px.set_ylabel("z [kpc]")
  px.set_font({'size':20})
  px.set_log("{}".format(field), False)
  # px.set_zlim(field , zmin=-1e3, zmax=1e3)
  px.set_zlim(field , zmin=-200, zmax=200)
  # px.set_cmap(field , cmap='Haze')
  px.set_unit(field, 'km/s')
  px.set_background_color('{}'.format(field))
  plot = px.plots['{}'.format(field)]
  plot.figure = fig
  plot.axes = grid[i].axes
  plot.cax = grid.cbar_axes[i]
  px.set_background_color(field)
  px._setup_plots()
# plt.savefig("results/plot{}200kpc.pdf".format(field), bbox_inches="tight")
plt.show()


# %% phase diagram
for i in range(len(models)):
  plot = yt.ParticlePlot(ds[i], ('PartType0', "Density"), ('PartType0', "temperature"),   ('PartType0', "Masses"))
  plot.set_unit(('PartType0', "Density"), 'g/cm/cm/cm')
  plot.set_unit(('PartType0', "Masses"), 'Msun')
  plot.set_log(('PartType0', "Density"), True)
  plot.set_log(('PartType0', "temperature"), True)
  plot.set_xlim(1e-32, 1e-21)
  plot.set_ylim(5e1, 1e8)
  plot.set_zlim(('PartType0', "Masses"),1e4, 1e7)
  # plot.hide_colorbar()
  plot.show()


# %% phase diagram
# ---------------------
fig = plt.figure()
grid = AxesGrid(fig, (0.075,0.075,0.85,0.85),
                nrows_ncols = (1, len(models)),
                axes_pad = 0.05,
                label_mode = "L",
                share_all = True,
                cbar_location="right",
                cbar_mode="single",
                cbar_size="5%",
                cbar_pad="1%")

for i in range(len(models)):
  print(i)
  my_sphere = ds[i].sphere("c", (1000.0, "kpc"))
  px = yt.ParticlePhasePlot(my_sphere, "density", "temperature", ["mass", ])
  px.set_unit("density", 'g/cm/cm/cm')
  px.set_unit("mass", 'Msun')
  px.set_log("density", True)
  px.set_log("temperature", True)
  px.set_font({'size':20})  
  px.set_xlim(1e-32, 1e-21)
  px.set_ylim(5e1, 1e8)
  px.set_zlim("mass" , zmin=1e4, zmax=1e7)
  plot = px.plots["mass"]
  plot.figure = fig
  plot.axes = grid[i].axes
  plot.cax = grid.cbar_axes[i]
  px._setup_plots()
  plot.axes.xaxis.set_minor_locator(plt.LogLocator(base=10.0, subs=[2.0, 5.0, 8.0]))
plt.show()
# plt.savefig("results/phasediagram.pdf", bbox_inches="tight")
# %% density 
# ----------
field = "density"
i = 0
p = yt.SlicePlot(ds[i], "x", ('gas', '{}'.format(field)), center = center[i], width=20)
# p.set_xlabel("x [kpc]")
# p.set_ylabel("y [kpc]")
p.set_font({'size':20})
p.set_zlim(field="{}".format(field), zmin=1e1, zmax=1e5)
# p.set_zlim(field="{}".format(field), zmin=1e-29, zmax=1e-22)
p.set_background_color(field)
p.annotate_contour(('gas', '{}'.format(field)), ncont=5, take_log=True, clim=(1e-28, 1e-24))
p.show()
# %%
