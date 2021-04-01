# %% movie
# ---------------------------

import yt
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid

tmax = 300

for time in range(tmax):

  models = [
      ["ss_model/ver01193","Sphericalsuperbubble"]
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
    #   ad = ds[0].all_data()
    #   density = ad[("PartType0","density")]
    #   wdens = np.where(density == np.max(density))
    #   coordinates = ad[("PartType0","Coordinates")]
    #   center[i] = coordinates[wdens][0]
      sp = ds[i].sphere("max", (40, 'kpc'))
      center[i] = sp.quantities.center_of_mass()
  # print ('center = ',center)

  fig = plt.figure()
  i = 0 # ss_model
  grid = AxesGrid(fig, (0.075,0.075,0.85,0.85),
                  nrows_ncols = (2,1),
                  axes_pad = 0.05,
                  label_mode = "L",
                  share_all = True,
                  cbar_location="right",
                  cbar_mode="each",
                  cbar_size="5%",
                  cbar_pad="1%")

  px = yt.OffAxisProjectionPlot(ds[i], [0,0,1], ('gas', '{}'.format("density")), center = center[i], width=45, north_vector=[0,1,0])
  px.set_xlabel("x [kpc]")
  px.set_ylabel("y [kpc]")
  px.set_font({'size':20})
  px.set_zlim(field="{}".format("density"), zmin=5e-6, zmax=1e-1)
  plot = px.plots['{}'.format("density")]
  plot.figure = fig
  plot.axes = grid[i].axes
  plot.cax = grid.cbar_axes[i]
  px.set_background_color("density")
  px._setup_plots()
  px = yt.OffAxisProjectionPlot(ds[i], [0,1,0], ('gas', '{}'.format("temperature")), center = center[i], width=45, north_vector=[0,0,1], weight_field='density')
  px.set_xlabel("x [kpc]")
  px.set_ylabel("z [kpc]")
  px.set_font({'size':20})
  px.set_zlim(field="{}".format("temperature"), zmin=1e3, zmax=1e7)
  px.set_cmap(field="{}".format("temperature"), cmap='hot')
  plot = px.plots['{}'.format("temperature")]
  plot.figure = fig
  plot.axes = grid[i+len(models)].axes
  plot.cax = grid.cbar_axes[i+len(models)]
  px.set_background_color("temperature")
  px._setup_plots()
  plt.savefig("movie/movie{:03}.png".format(time), bbox_inches="tight")
