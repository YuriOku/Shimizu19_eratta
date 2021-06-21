# %% movie
# ---------------------------

import yt
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
from yt.visualization.base_plot_types import get_multi_plot
from matplotlib.colors import LogNorm
import unyt
import sys

tmax = 100

obj = sys.argv[1]

for time in range(tmax):

  models = [
      ["../ss_model/ver05291/fiducial","entropy-based"]
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
      px = yt.OffAxisProjectionPlot(ds[i], [0,1,0], ('gas', '{}'.format ("temperature")), center = center[i],   width=45, north_vector=[0,0,-1],  weight_field='density')
      px.set_xlabel("x [kpc]")
      px.set_ylabel("z [kpc]")
      px.set_font({'size':20})
      px.set_zlim(field="{}".format("temperature"), zmin=1e3, zmax=1e7)
      px.set_cmap(field="{}".format("temperature"), cmap='hot')
      px.set_background_color("temperature")
      plot = px.plots['{}'.format("temperature")]
      plot.figure = fig
      plot.axes = grid[i+len(models)].axes
      plot.cax = grid.cbar_axes[1]
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
  elif obj == "oblique":
    grid = AxesGrid(fig, (0.075,0.075,0.85,0.85),
                    nrows_ncols = (len(models), 2),
                    axes_pad = 0.0,
                    label_mode = None,
                    share_all = True,
                    cbar_location="bottom",
                    cbar_mode="edge",
                    cbar_size="5%",
                    cbar_pad="1%",
                    )
    for i in range(len(models)):
      cen = center[i]
      cen[i] -= 10.*unyt.kpc
      px = yt.OffAxisProjectionPlot(ds[i], [0,1,0], ('gas', '{}'.format("density")),  center = center[i], width=20, north_vector=[0,1,-2])
      px.set_font({'size':20})
      px.set_zlim(field="{}".format("density"), zmin=5e-6, zmax=1e-1)
      px.hide_axes()
      plot = px.plots['{}'.format("density")]
      plot.figure = fig
      plot.axes = grid[i*2].axes
      plot.cax = grid.cbar_axes[0]
      px.set_background_color("density")
      px._setup_plots()
      
      cen[i] += 20.*unyt.kpc
      px = yt.OffAxisProjectionPlot(ds[i], [0,1,0], ('gas', '{}'.format ("temperature")), center = cen,   width=20,  weight_field='density', north_vector=[0,1,-2])
      px.set_font({'size':20})
      px.set_zlim(field="{}".format("temperature"), zmin=1e3, zmax=1e7)
      px.set_cmap(field="{}".format("temperature"), cmap='hot')
      px.hide_axes()
      px.set_background_color("temperature")
      plot = px.plots['{}'.format("temperature")]
      plot.figure = fig
      plot.axes = grid[i*2+1].axes
      plot.cax = grid.cbar_axes[1]
      px._setup_plots()
    plt.savefig("movie/oblique{:03}.png".format(time), bbox_inches="tight")
  
  elif obj == "dtm":
    field = "density"
    fig = plt.figure()
    grid = AxesGrid(fig, (0.075,0.075,0.85,0.85),
                    nrows_ncols = (2,len(models)),
                    axes_pad = 0.05,
                    label_mode = "L",
                    share_all = True,
                    cbar_location="right",
                    cbar_mode="single",
                    cbar_size="5%",
                    cbar_pad="1%")
    for i in range(len(models)):
      px = yt.OffAxisProjectionPlot(ds[i], [0,0,1], ('gas', '{}'.format(field)),  center = center[i],   width=45, north_vector=[0,1,0])
      px.set_xlabel("x [kpc]")
      px.set_ylabel("y [kpc]")
      px.set_font({'size':20})
      px.set_zlim(field="{}".format(field), zmin=5e-6, zmax=1e-1)
      plot = px.plots['{}'.format(field)]
      plot.figure = fig
      plot.axes = grid[0].axes
      plot.cax = grid.cbar_axes[0]
      px.set_background_color(field)
      px._setup_plots()
      px = yt.OffAxisProjectionPlot(ds[i], [0,1,0], ('gas', '{}'.format (field)), center = center[i],   width=45, north_vector=[0,0,-1])
      px.set_xlabel("x [kpc]")
      px.set_ylabel("z [kpc]")
      px.set_font({'size':20})
      px.set_zlim(field="{}".format(field), zmin=5e-6, zmax=1e-1)
      plot = px.plots['{}'.format(field)]
      plot.figure = fig
      plot.axes = grid[1].axes
      plot.cax = grid.cbar_axes[1]
      px.set_background_color(field)
      px._setup_plots()
    plt.savefig("movie/density{:03}.png".format(time), bbox_inches="tight")
    plt.close()
    fig = plt.figure()
    field = "temperature"
    grid = AxesGrid(fig, (0.075,0.075,0.85,0.85),
                    nrows_ncols = (2,len(models)),
                    axes_pad = 0.05,
                    label_mode = "L",
                    share_all = True,
                    cbar_location="right",
                    cbar_mode="single",
                    cbar_size="5%",
                    cbar_pad="1%")
    for i in range(len(models)):
      px = yt.OffAxisProjectionPlot(ds[i], [0,0,1], ('gas', '{}'.format(field)),  center = center[i],   width=45, north_vector=[0,1,0],  weight_field='density')
      px.set_xlabel("x [kpc]")
      px.set_ylabel("y [kpc]")
      px.set_font({'size':20})
      px.set_zlim(field="{}".format(field), zmin=1e3, zmax=1e7)
      px.set_cmap(field="{}".format("temperature"), cmap='hot')
      plot = px.plots['{}'.format(field)]
      plot.figure = fig
      plot.axes = grid[i].axes
      plot.cax = grid.cbar_axes[0]
      px.set_background_color(field)
      px._setup_plots()
      px = yt.OffAxisProjectionPlot(ds[i], [0,1,0], ('gas', '{}'.format (field)), center = center[i],   width=45, north_vector=[0,0,-1],  weight_field='density')
      px.set_xlabel("x [kpc]")
      px.set_ylabel("z [kpc]")
      px.set_font({'size':20})
      px.set_zlim(field="{}".format(field), zmin=1e3, zmax=1e7)
      px.set_cmap(field="{}".format("temperature"), cmap='hot')
      plot = px.plots['{}'.format(field)]
      plot.figure = fig
      plot.axes = grid[i+len(models)].axes
      plot.cax = grid.cbar_axes[1]
      px.set_background_color(field)
      px._setup_plots()
    plt.savefig("movie/temperature{:03}.png".format(time), bbox_inches="tight")
    plt.close()
    fig = plt.figure()
    field = "metallicity"
    grid = AxesGrid(fig, (0.075,0.075,0.85,0.85),
                    nrows_ncols = (2,len(models)),
                    axes_pad = 0.05,
                    label_mode = "L",
                    share_all = True,
                    cbar_location="right",
                    cbar_mode="single",
                    cbar_size="5%",
                    cbar_pad="1%")
    for i in range(len(models)):
      px = yt.OffAxisProjectionPlot(ds[i], [0,0,1], ('gas', '{}'.format(field)),  center = center[i],   width=45, north_vector=[0,1,0], weight_field="density")
      px.set_xlabel("x [kpc]")
      px.set_ylabel("y [kpc]")
      px.set_font({'size':20})
      px.set_zlim(field , zmin=1e-5, zmax=1e-1)
      px.set_unit(field, '')
      px.set_cmap(field , cmap='Haze')
      plot = px.plots['{}'.format(field)]
      plot.figure = fig
      plot.axes = grid[i].axes
      plot.cax = grid.cbar_axes[0]
      px.set_background_color(field)
      px._setup_plots()
      px = yt.OffAxisProjectionPlot(ds[i], [0,1,0], ('gas', '{}'.format (field)), center = center[i],   width=45, north_vector=[0,0,-1],  weight_field='density')
      px.set_xlabel("x [kpc]")
      px.set_ylabel("z [kpc]")
      px.set_font({'size':20})
      px.set_zlim(field , zmin=1e-5, zmax=1e-1)
      px.set_unit(field, '')
      px.set_cmap(field , cmap='Haze')
      plot = px.plots['{}'.format(field)]
      plot.figure = fig
      plot.axes = grid[i+len(models)].axes
      plot.cax = grid.cbar_axes[1]
      px.set_background_color(field)
      px._setup_plots()
    plt.savefig("movie/metallicity{:03}.png".format(time), bbox_inches="tight")
    plt.close()

  elif obj == "dtm-oblique":
    for i in range(len(models)):
      field = "density"
      px = yt.OffAxisProjectionPlot(ds[i], [0,1,0], ('gas', '{}'.format (field)), center = center[i],   width=30, north_vector=[0,1,-2])
      px.set_xlabel("x [kpc]")
      px.set_ylabel("y [kpc]")
      px.set_font({'size':20})
      px.set_zlim(field="{}".format(field), zmin=5e-6, zmax=1e-1)
      px.set_background_color(field)
      px.save("movie/{}{}_{:03}.png".format(field, i, time))

      field = "temperature"
      px = yt.OffAxisProjectionPlot(ds[i], [0,1,0], ('gas', '{}'.format (field)), center = center[i],   width=30,  weight_field='density', north_vector=[0,1,-2])
      px.set_xlabel("x [kpc]")
      px.set_ylabel("y [kpc]")
      px.set_font({'size':20})
      px.set_zlim(field="{}".format(field), zmin=1e3, zmax=1e7)
      px.set_cmap(field="{}".format("temperature"), cmap='hot')
      px.set_background_color(field)
      px.save("movie/{}{}_{:03}.png".format(field, i, time))

      field = "metallicity"
      px = yt.OffAxisProjectionPlot(ds[i], [0,1,0], ('gas', '{}'.format (field)), center = center[i],   width=30,  weight_field='density', north_vector=[0,1,-2])
      px.set_xlabel("x [kpc]")
      px.set_ylabel("y [kpc]")
      px.set_font({'size':20})
      px.set_zlim(field , zmin=1e-5, zmax=1e-1)
      px.set_cmap(field , cmap='Haze')
      px.set_unit(field, '')
      px.set_background_color(field)
      px.save("movie/{}{}_{:03}.png".format(field, i, time))
  else:
    assert 0


# %%
