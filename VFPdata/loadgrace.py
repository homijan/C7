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
   #import argparse

   #ps = argparse.ArgumentParser( description = 'Extract data from xmgrace (.agr) files')
   #ps.add_argument( 'filename', type = str, help = 'xmgrace file' )
   #args = ps.parse_args()

   #d = LoadGrace( args.filename )

   print "Loading temperature data..."
   td = LoadGrace('temp.agr')

   for i in range(0, len(td)):
      legend = td.legends()[i]
      ## substitute blank spaces 
      legend = str(legend).replace(" ", "_")
      data = td.sets()[i]
      x = []
      Te = []
      for j in range(0, len(data)):
         row = data[j]
         x.append(row[0])
         Te.append(row[1])
      x = np.array(x)
      Te = np.array(Te)
      ## Temperature in [eV]
      Te = Te * 1e3 
      ## Obtain a given point value and gradient of temperature.
      smooth = 0 # lower less smoothing
      ## Find a spline for the temperature data.
      tck = splrep(x, Te, s=smooth)
      # fine data
      N = 10000
      x_fine = np.linspace(min(x), max(x), N)
      Te_fine = splev(x_fine, tck, der=0)
      #gradTe_fine = splev(x, tck, der=1)
      #print '# legend=' + td.legends()[i] + ' comment=' + td.comments()[i]
      print '# legend='
      print legend 
      #print x_fine
      #print Te_fine 
      print
      np.savetxt('Te_'+legend+'.dat', np.transpose([x_fine, Te_fine]))
   #quit()

   print "Loading flux data..."
   fd = LoadGrace('flux.agr')

   for i in range(0, len(fd)):
      legend = fd.legends()[i]
      ## substitute blank spaces 
      legend = str(legend).replace(" ", "_")
      data = fd.sets()[i]
      x = []
      flux = []
      for j in range(0, len(data)):
         row = data[j]
         x.append(row[0])
         flux.append(row[1])
      x = np.array(x)
      flux = np.array(flux)
      ## Flux output is in microns and W/cm2
      x = x * 1e4
      flux = flux * 1e-7
      #print '# legend=' + td.legends()[i] + ' comment=' + td.comments()[i]
      print '# legend='
      print legend 
      #print flux 
      print
      np.savetxt('flux_'+legend+'.dat', np.transpose([x, flux]))
