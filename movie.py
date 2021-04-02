# %% movie
# ---------------------------

import yt
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
import sys

tmax = 300

obj = sys.argv[1]

for time in range(tmax):

  models = [
      ["../ss_model/ver03311","with cool-off"],
      ["../ss_model/ver03311_NoCoolOff","without cool-off"]
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
  ad = [0]*len(models)
  center = [0]*len(models)
  for i in range(len(models)):
      ds[i] = yt.load('{0}/snapshot_{1:03}/snapshot_{1:03}.hdf5'.format(models[i][0], time),unit_base=unit_base,bounding_box=bbox)
    #   ad = ds[0].all_data()
    #   density = ad[("PartType0","density")]
    #   wdens = np.where(density == np.max(density))
    #   coordinates = ad[("PartType0","Coordinates")]
    #   center[i] = coordinates[wdens][0]
      sp = ds[i].sphere("max", (40, 'kpc'))
      center[i] = sp.quantities.center_of_mass()
      ad[i] = ds[i].all_data()
  # print ('center = ',center)

  fig = plt.figure()
  if obj == "movie":
    grid = AxesGrid(fig, (0.075,0.075,0.85,0.85),
                    nrows_ncols = (2,len(models)),
                    axes_pad = 0.05,
                    label_mode = "L",
                    share_all = True,
                    cbar_location="right",
                    cbar_mode="edge",
                    cbar_size="5%",
                    cbar_pad="1%")
    for i in range(len(models)):
      px = yt.OffAxisProjectionPlot(ds[i], [0,0,1], ('gas', '{}'.format("density")),  center = center[i],   width=45, north_vector=[0,1,0])
      px.set_xlabel("x [kpc]")
      px.set_ylabel("y [kpc]")
      px.set_font({'size':20})
      px.set_zlim(field="{}".format("density"), zmin=5e-6, zmax=1e-1)
      plot = px.plots['{}'.format("density")]
      plot.figure = fig
      plot.axes = grid[i].axes
      plot.cax = grid.cbar_axes[0]
      px.set_background_color("density")
      px._setup_plots()
      px = yt.OffAxisProjectionPlot(ds[i], [0,1,0], ('gas', '{}'.format ("temperature")), center = center[i],   width=45, north_vector=[0,0,1],  weight_field='density')
      px.set_xlabel("x [kpc]")
      px.set_ylabel("z [kpc]")
      px.set_font({'size':20})
      px.set_zlim(field="{}".format("temperature"), zmin=1e3, zmax=1e7)
      px.set_cmap(field="{}".format("temperature"), cmap='hot')
      plot = px.plots['{}'.format("temperature")]
      plot.figure = fig
      plot.axes = grid[i+len(models)].axes
      plot.cax = grid.cbar_axes[1]
      px.set_background_color("temperature")
      px._setup_plots()
    plt.savefig("movie/movie{:03}.png".format(time), bbox_inches="tight")
  elif obj == "phase":
    for i in range(len(models)):
      plot = yt.ParticlePlot(ds[i], ('PartType0', "Density"), ('PartType0',"temperature"), ('PartType0', "Masses"))
      plot.set_unit(('PartType0', "Density"), 'g/cm/cm/cm')
      plot.set_unit(('PartType0', "Masses"), 'Msun')
      plot.set_log(('PartType0', "Density"), True)
      plot.set_log(('PartType0', "temperature"), True)
      plot.set_xlim(1e-32, 1e-21)
      plot.set_ylim(5e1, 1e8)
      plot.set_zlim(('PartType0', "Masses"),1e4, 1e7)
      plot.show_axes()
      if i != len(models) - 1:
        plot.hide_colorbar()
      plot.save("movie/phase{}_{:03}.png".format(i, time))
  else:
    assert 0


# %%
