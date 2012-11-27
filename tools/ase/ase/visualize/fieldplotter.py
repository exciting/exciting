"""plotting fields defined on atoms during a simulation."""

from ase.visualize.primiplotter import PostScriptFile, PnmFile, GifFile, JpegFile, X11Window
from ase.visualize.primiplotter import PrimiPlotter as _PrimiPlotter
import numpy
import time

class FieldPlotter(_PrimiPlotter):
    def __init__(self, atoms, datasource=None, verbose=0, timing=0,
                 interval=1, initframe=0):
        _PrimiPlotter.__init__(self, atoms, verbose=verbose, timing=timing,
                               interval=interval, initframe=initframe)
        self.datasource = datasource
        self.dims = (100,100)
        self.set_plot_plane("xy")
        self.set_data_range("plot")
        self.set_background(0.0)
        self.set_red_yellow_colors()
        
    def set_plot_plane(self, plane):
        """Set the plotting plane to xy, xz or yz (default: xy)"""
        if plane in ("xy", "xz", "yz"):
            self.plane = plane
        else:
            raise ValueError, "The argument to plotPlane must be 'xy', 'xz' or 'yz'."

    def set_data_range(self, range1, range2=None):
        """Set the range of the data used when coloring.

        This function sets the range of data values mapped unto colors
        in the final plot.
        
        Three possibilities:

        'data':        Autoscale using the data on visible atoms.
                       The range goes from the lowest to the highest
                       value present on the atoms.  If only a few atoms
                       have extreme values, the entire color range may not
                       be used on the plot, as many values may be averaged
                       on each point in the plot.

        'plot':        Autoscale using the data on the plot.  Unlike 'data'
                       this guarantees that the entire color range is used.

        min, max:      Use the range [min, max]
                       
        """
        if (range1 == "data" or range1 == "plot") and range2 == None:
            self.autorange = range1
        elif range2 != None:
            self.autorange = None
            self.range = (range1, range2)
        else:
            raise ValueError, "Illegal argument(s) to set_data_range"

    def set_background(self, value):
        """Set the data value of the background.  See also set_background_color

        Set the value of the background (parts of the plot without atoms) to
        a specific value, or to 'min' or 'max' representing the minimal or
        maximal data values on the atoms.

        Calling set_background cancels previous calls to set_background_color.
        """
        self.background = value
        self.backgroundcolor = None

    def set_background_color(self, color):
        """Set the background color.  See also set_background.

        Set the background color.  Use a single value in the range [0, 1[
        for gray values, or a tuple of three such values as an RGB color.

        Calling set_background_color cancels previous calls to set_background.
        """
        self.background = None
        self.backgroundcolor = color

    def set_red_yellow_colors(self, reverse=False):
        """Set colors to Black-Red-Yellow-White (a.k.a. STM colors)"""
        self.set_colors([(0.0, 0, 0, 0),
                        (0.33, 1, 0, 0),
                        (0.66, 1, 1, 0),
                        (1.0, 1, 1, 1)],
                       reverse)
        
    def set_black_white_colors(self, reverse=False):
        """Set the color to Black-White (greyscale)"""
        self.set_colors([(0.0, 0),  (1.0, 1)], reverse)

    def set_colors(self, colors, reverse=False):
        colors = numpy.array(colors, numpy.float)
        if len(colors.shape) != 2:
            raise ValueError, "Colors must be a 2D array."
        if reverse:
            colors[:,0] = 1 - colors[:,0]
            colors = numpy.array(colors[::-1,:])
            #print colors
        if colors[0,0] != 0.0 or colors[-1,0] != 1.0:
            raise ValueError, "First row must define the value 0 and last row must define the value 1"
        if colors.shape[1] == 2:
            self.colormode = 1
        elif colors.shape[1] == 4:
            self.colormode = 3
        else:
            raise ValueError, "Color specification must be Nx2 (grey) or Nx4 (rgb) matrix."
        self.colorfunction = InterpolatingFunction(colors[:,0], colors[:,1:])
        
    def plot(self, data=None):
        """Create a plot now.  Does not respect the interval timer.

        This method makes a plot unconditionally.  It does not look at
        the interval variable, nor is this plot taken into account in
        the counting done by the update() method if an interval
        variable was specified.

        If data is specified, it must be an array of numbers with the
        same length as the atoms.  That data will then be plotted.  If
        no data is given, the data source specified when creating the
        plotter is used.
        
        """
        if self.timing:
            self._starttimer()
        self.log("FieldPlotter: Starting plot at "
                 + time.strftime("%a, %d %b %Y %H:%M:%S"))
        if data is None:
            data = self.datasource()
        if len(data) != len(self.atoms):
            raise ValueError, ("Data has wrong length: %d instead of %d."
                               % (len(data), len(self.atoms)))
        
        invisible = self._getinvisible()
        coords = self._rotate(self._getpositions())
        radii = self._getradii()
        if self.autoscale:
            self._autoscale(coords,radii)
        scale = self.scale * self.relativescale
        coords = scale * coords
        center = self._getcenter(coords)
        offset = numpy.array(self.dims + (0.0,))/2.0 - center
        coords = coords + offset
        radii = radii * scale
        self.log("Scale is %f and size is (%d, %d)"
                 % (scale, self.dims[0], self.dims[1]))
        self.log("Physical size of plot is %f Angstrom times %f Angstrom"
                 % (self.dims[0] / scale, self.dims[1] / scale))

        # Remove invisible atoms
        selector = numpy.logical_not(invisible)
        coords = numpy.compress(selector, coords, 0)
        radii = numpy.compress(selector, radii)
        data = numpy.compress(selector, data)

        self.log("plotting data in the range [%f,%f]" %
                   (data.min(), data.max()))
        # Now create the output array
        sumarray = numpy.zeros(self.dims, numpy.float)
        weight = numpy.zeros(self.dims)

        # Loop over all atoms, and plot them
        nmiss = 0
        if self.plane == "xy":
            xy = coords[:,:2]
        elif self.plane == "xz":
            xy = coords[:,::2]
        elif self.plane == "yz":
            xy = coords[:,1:]
        else:
            raise RuntimeError, "self.plane is bogus: "+str(self.plane)
        assert xy.shape[1] == 2

        self.log("plotting %d atoms on %d * %d (= %d) grid" %
                   (len(xy), sumarray.shape[0], sumarray.shape[1],
                    len(sumarray.flat)))
                                                            
        xy = xy.astype(numpy.int)
        for i in xrange(len(xy)):
            (x, y) = xy[i]
            d = data[i]
            if (x >= 0 and x < self.dims[0] and y >= 0 and y < self.dims[1]):
                sumarray[x,y] += d
                weight[x,y] += 1
            else:
                nmiss += 1
        print "... %d atoms fell outside plot." % (nmiss,)

        datamap = self._makedatamap(sumarray, weight, data.min(), data.max())
        self.log("Range of data map: [%f, %f]" %
                   (datamap.min(), datamap.max()))
        plot = self._makeplotmap(datamap, weight)
        #self.log("Range of plot: [%f, %f]" %
        #           (min(plot.flat), max(plot.flat)))
        examinplot = plot[:]
        examinplot.shape = (plot.shape[0] * plot.shape[1],) + plot.shape[2:]
        self.log("Range of plot: %s -> %s" %
                 (str(examinplot.min(0)), str(examinplot.max(0))))
        del examinplot
        for device in self.outputdevice:
            device.inform_about_scale(scale)
            device.plotArray(self.n, numpy.swapaxes(plot,0,1))
        self.n = self.n + 1
        self.log("FieldPlotter: Finished plotting at "
                 + time.strftime("%a, %d %b %Y %H:%M:%S"))
        self.log("\n\n")

        
    def _makedatamap(self, sumarray, weight, minimum, maximum):
        background = numpy.equal(weight, 0)
        print "Number of background points:", sum(background.flat)
        datamap = sumarray / numpy.where(background, 1, weight)
        
        if self.background is not None:
            if self.background == "min":
                bg = minimum
            elif self.background == "max":
                bg = maximum
            else:
                bg = self.background
            datamap = numpy.where(background, bg, datamap)
            
        if self.autorange == "data":
            datamap = (datamap - minimum) / (maximum - minimum)
            self.log("Autorange using data.  Data range is [%f, %f]"
                     % (minimum, maximum))
        elif self.autorange == "plot":
            ma = numpy.where(background, minimum, datamap).max()
            mi = numpy.where(background, maximum, datamap).min()
            datamap = (datamap - mi) / (ma - mi)
            self.log("Autorange using plot.  Data range is [%f, %f]"
                     % (mi, ma))
        else:
            assert self.autorange == None
            datamap = (datamap - self.range[0]) / (self.range[1]
                                                   - self.range[0])
            datamap = numpy.clip(datamap, 0.0, 1.0)
            self.log("Data range specified by user: [%f, %f]" % self.range)
        datamap = numpy.where(background, bg, datamap)
        assert datamap.min() >= 0 and datamap.max() <= 1.0
        
        return datamap

    def _makeplotmap(self, datamap, weight):
        plot = numpy.zeros(self.dims + (self.colormode,), numpy.float)
        for i in range(self.dims[0]):
            for j in range(self.dims[1]):
                if self.backgroundcolor is not None and weight[i,j] == 0:
                    plot[i,j,:] = self.backgroundcolor
                else:
                    x = datamap[i,j]
                    plot[i,j,:] = self.colorfunction(x)
        return plot
    
class InterpolatingFunction:
    def __init__(self, xpoints, ypoints):
        if len(xpoints) != len(ypoints):
            raise ValueError, "Length of x and y arrays should be the same."
        idx = xpoints.argsort()
        self.xpoints = xpoints[idx]
        self.ypoints = ypoints[idx]
    def __call__(self, x):
        n = self.xpoints.searchsorted(x)
        if n == 0:
            return self.ypoints[0]
        if n == len(self.xpoints):
            return self.xpoints[-1]
        x0 = self.xpoints[n-1]
        x1 = self.xpoints[n]
        y0 = self.ypoints[n-1]
        y1 = self.ypoints[n]
        return y0 + (y1 - y0) / (x1 - x0) * (x - x0)
    
