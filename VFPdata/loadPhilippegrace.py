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
   ps.add_argument( '-c', '--column', type = int, help = 'column to load' ) 
   ps.add_argument( '-r', '--row', type = int, help = 'row to load' )
   ps.add_argument( '-mx', '--xmultiply', type = float, help = 'multiply x axis' )
   ps.add_argument( '-my', '--ymultiply', type = float, help = 'multiply y axis' )
   ## A no value argument solution.
   ps.add_argument("-gf", "--GraceFile", action='store_true', help="Load grace file by adding -g/--grace argument.")
   ps.add_argument( '-mom', '--moment', type = int, help = 'moment to apply', default=0)
   ps.add_argument("-s", "--pltshow", action='store_true', help="Show plot by adding -s/--pltshow argument.")

   ## Default value of the column to be loaded.
   column = 1
   row = 1
   args = ps.parse_args()
   if (args.column):
      column = args.column
   if (args.row):
      row = args.row

   legends = []
   xs = []
   ys = []

   if (args.GraceFile):
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
         ## Multiply x.
         if (args.xmultiply):
            x = x * args.xmultiply
         ## Multiply y.
         if (args.ymultiply):
            y = y * args.ymultiply
         print '# legend='
         print legend 
         print
         # Append data.
         legends.append(legend)
         xs.append(x)
         ys.append(y)
   else:
      if (args.row):
         data = np.array(np.loadtxt(args.filename))
         x = data[0, :]
         y = data[row, :]
      else:
         x, y = np.loadtxt(args.filename,  usecols=(0, column), unpack=True)
      ## Multiply x.
      if (args.xmultiply):
         x = x * args.xmultiply
      ## Multiply y.
      if (args.ymultiply):
         y = y * args.ymultiply
      ## Apply moment on the data.
      for i in range(args.moment):
         y *= x
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
      ## Store fine data
      if (args.outname):
         np.savetxt(args.outname+legend+'.txt', np.transpose([x_fine, y_fine]))
      if (args.pltshow):
         plt.plot(x_fine, y_fine, label=legend)
         plt.legend()
         plt.show()
