#!/usr/bin/env python

import numpy as np
import re

class LoadGrace:
   """
   Simple module to load a Grace file (.agr) plot and extract
   legends, comments, and -- importantly -- data into numpy arrays.
   """
  
   def __init__(self, filename):
      self._sets = []
      self._legends = []
      self._comments = []

      f = open(filename, 'r')
      for line in f:
         if re.compile('@    s(.*) legend').search(line):
            self._legends.append( line.split('"')[1] )

         if re.compile('@    s(.*) comment').search(line):
            self._comments.append( line.split('"')[1] )

         if '@target' in line:
            tmp = []
            next(f)
            for row in f:
               if row!='&\n':
                  tmp.append( np.fromstring( row, sep=' ' ) )
               else:
                  self._sets.append( np.array(tmp) )
                  break
      f.close()

   def sets(self): return self._sets
   def legends(self): return self._legends
   def comments(self): return self._comments
   def __len__(self): return len(self._sets)

# If run as main program
if __name__ == "__main__":
   from scipy.interpolate import splev, splrep
   import matplotlib.pyplot as plt
   import argparse

   ps = argparse.ArgumentParser( description = 'Extract data from xmgrace (.agr) files')
   ps.add_argument( '-f', '--filename', type = str, help = 'xmgrace file' )
   ps.add_argument( '-o', '--outname', type = str, help = 'output file name' )
   ps.add_argument( '-m', '--multiply', type = float, help = 'output file name' )
   ## A no value argument solution.
   ps.add_argument("-g", "--grace", action='store_true', help="Load grace file by adding -g/--grace argument.")
   ps.add_argument("-s", "--pltshow", action='store_true', help="Show plot by adding -s/--pltshow argument.")
   ps.add_argument("-x", "--xmicrons2cm", action='store_true', help="Convert from microns to cm by adding -x/--xcm argument.")
   ps.add_argument("-T", "--TkeV2eV", action='store_true', help="Convert from keV to eV by adding -T/--temperature argument.")


   args = ps.parse_args()

   legends = []
   xs = []
   ys = []

   if (args.grace):
      print "Loading grace data..."
      d = LoadGrace( args.filename )
      for i in range(0, len(d)):
         legend = d.legends()[i]
         ## substitute blank spaces 
         legend = str(legend).replace(" ", "_")
         data = d.sets()[i]
         x = []
         y = []
         for j in range(0, len(data)):
            row = data[j]
            x.append(row[0])
            y.append(row[1])
         x = np.array(x)
         y = np.array(y)
         print '# legend='
         print legend 
         print
         # Append data.
         legends.append(legend)
         xs.append(x)
         ys.append(y)
   else:
      x, y = np.loadtxt(args.filename,  usecols=(0, 1), unpack=True)
      # Append data.
      legends.append(args.filename)
      xs.append(x)
      ys.append(y)

   for i in range(0, len(xs)):
      legend = legends[i]
      x = xs[i]
      y = ys[i]
      ## Obtain a given point value and gradient of temperature.
      smooth = 0 # lower less smoothing
      ## Find a spline for the temperature data.
      tck = splrep(x, y, s=smooth)
      # fine data
      N = 10000
      x_fine = np.linspace(min(x), max(x), N)
      y_fine = splev(x_fine, tck, der=0)
      #grady_fine = splev(x, tck, der=1)
      ## x in [cm]
      if (args.xmicrons2cm):
         x_fine = x_fine * 1e-4
      ## Temperature in [eV]
      if (args.TkeV2eV):
         y_fine = y_fine * 1e3 
      ## Multiply y.
      if (args.multiply):
         y_fine = y_fine * args.multiply
      ## Store fine data
      np.savetxt(args.outname+legend+'.txt', np.transpose([x_fine, y_fine]))
      if (args.pltshow):
         plt.plot(x_fine, y_fine, label=legend)
         plt.legend()
         plt.show()
