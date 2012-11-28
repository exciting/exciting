"""An experimental package for making plots during a simulation.

A PrimiPlotter can plot a list of atoms on one or more output devices.
"""

from numpy import *
from ase.visualize.colortable import color_table
import ase.data
import sys, os, time, weakref

class PrimiPlotterBase:
    "Base class for PrimiPlotter and Povrayplotter."
    #def set_dimensions(self, dims):
    #    "Set the size of the canvas (a 2-tuple)."
    #    self.dims = dims
        
    def set_rotation(self, rotation):
        "Set the rotation angles (in degrees)."
        self.angles[:] = array(rotation) * (pi/180)
        
    def set_radii(self, radii):
        """Set the atomic radii.  Give an array or a single number."""
        self.radius = radii

    def set_colors(self, colors):
        """Explicitly set the colors of the atoms."""
        self.colors = colors

    def set_color_function(self, colors):
        """Set a color function, to be used to color the atoms."""
        if callable(colors):
            self.colorfunction = colors
        else:
            raise TypeError, "The color function is not callable."

    def set_invisible(self, inv):
        """Choose invisible atoms."""
        self.invisible = inv

    def set_invisibility_function(self, invfunc):
        """Set an invisibility function."""
        if callable(invfunc):
            self.invisibilityfunction = invfunc
        else:
            raise TypeError, "The invisibility function is not callable."

    def set_cut(self, xmin=None, xmax=None, ymin=None, ymax=None,
               zmin=None, zmax=None):
        self.cut = {"xmin":xmin, "xmax":xmax, "ymin":ymin, "ymax":ymax,
                    "zmin":zmin, "zmax":zmax}
    
    def update(self, newatoms = None):
        """Cause a plot (respecting the interval setting).

        update causes a plot to be made.  If the interval variable was
        specified when the plotter was create, it will only produce a
        plot with that interval.  update takes an optional argument,
        newatoms, which can be used to replace the list of atoms with
        a new one.
        """
        if newatoms is not None:
            self.atoms = newatoms
        if self.skipnext <= 0:
            self.plot()
            self.skipnext = self.interval
        self.skipnext -= 1
        
    def set_log(self, log):
        """Sets a file for logging.

        log may be an open file or a filename.
        """
        if hasattr(log, "write"):
            self.logfile = log
            self.ownlogfile = False
        else:
            self.logfile = open(log, "w")
            self.ownlogfile = True

    def log(self, message):
        """logs a message to the file set by set_log."""
        if self.logfile is not None:
            self.logfile.write(message+"\n")
            self.logfile.flush()
        self._verb(message)
        
    def _verb(self, txt):
        if self.verbose:
            sys.stderr.write(txt+"\n")
    
    def _starttimer(self):
        self.starttime = time.time()

    def _stoptimer(self):
        elapsedtime = time.time() - self.starttime
        self.totaltime = self.totaltime + elapsedtime
        print "plotting time %s sec (total %s sec)" % (elapsedtime,
                                                       self.totaltime)

    def _getpositions(self):
        return self.atoms.get_positions()

    def _getradii(self):
        if self.radius is not None:
            if hasattr(self.radius, "shape"):
                return self.radius   # User has specified an array
            else:
                return self.radius * ones(len(self.atoms), float)
        # No radii specified.  Try getting them from the atoms.
        try:
            return self.atoms.get_atomic_radii()
        except AttributeError:
            try:
                z = self._getatomicnumbers()
            except AttributeError:
                pass
            else:
                return ase.data.covalent_radii[z]
        # No radius available.  Defaulting to 1.0
        return ones(len(self.atoms), float)

    def _getatomicnumbers(self):
        return self.atoms.get_atomic_numbers()
    
    def _getcolors(self):
        # Try any explicitly given colors
        if self.colors is not None:
            if type(self.colors) == type({}):
                self.log("Explicit colors dictionary")
                return _colorsfromdict(self.colors,
                                       asarray(self.atoms.get_tags(),int))
            else:
                self.log("Explicit colors")
                return self.colors
        # Try the color function, if given
        if self.colorfunction is not None:
            self.log("Calling color function.")
            return self.colorfunction(self.atoms)
        # Maybe the atoms know their own colors
        try:
            c = self.atoms.get_colors()
        except AttributeError:
            c = None
        if c is not None:
            if type(c) == type({}):
                self.log("Color dictionary from atoms.get_colors()")
                return _colorsfromdict(c, asarray(self.atoms.get_tags(),int))
            else:
                self.log("Colors from atoms.get_colors()")
                return c
        # Default to white atoms
        self.log("No colors: using white")
        return ones(len(self.atoms), float)

    def _getinvisible(self):
        if self.invisible is not None:
            inv = self.invisible
        else:
            inv = zeros(len(self.atoms))
        if self.invisibilityfunction:
            inv = logical_or(inv, self.invisibilityfunction(self.atoms))
        r = self._getpositions()
        if len(r) > len(inv):
            # This will happen in parallel simulations due to ghost atoms.
            # They are invisible.  Hmm, this may cause trouble.
            i2 = ones(len(r))
            i2[:len(inv)] = inv
            inv = i2
            del i2
        if self.cut["xmin"] is not None:
            inv = logical_or(inv, less(r[:,0], self.cut["xmin"]))
        if self.cut["xmax"] is not None:
            inv = logical_or(inv, greater(r[:,0], self.cut["xmax"]))
        if self.cut["ymin"] is not None:
            inv = logical_or(inv, less(r[:,1], self.cut["ymin"]))
        if self.cut["ymax"] is not None:
            inv = logical_or(inv, greater(r[:,1], self.cut["ymax"]))
        if self.cut["zmin"] is not None:
            inv = logical_or(inv, less(r[:,2], self.cut["zmin"]))
        if self.cut["zmax"] is not None:
            inv = logical_or(inv, greater(r[:,2], self.cut["zmax"]))
        return inv        

    def __del__(self):
        if self.ownlogfile:
            self.logfile.close()
            
class PrimiPlotter(PrimiPlotterBase):
    """Primitive PostScript-based plots during a simulation.

    The PrimiPlotter plots atoms during simulations, extracting the
    relevant information from the list of atoms.  It is created using
    the list of atoms as an argument to the constructor.  Then one or
    more output devices must be attached using set_output(device).  The
    list of supported output devices is at the end.

    The atoms are plotted as circles.  The system is first rotated
    using the angles specified by set_rotation([vx, vy, vz]).  The
    rotation is vx degrees around the x axis (positive from the y
    toward the z axis), then vy degrees around the y axis (from x
    toward z), then vz degrees around the z axis (from x toward y).
    The rotation matrix is the same as the one used by RasMol.

    Per default, the system is scaled so it fits within the canvas
    (autoscale mode).  Autoscale mode is enabled and disables using
    autoscale("on") or autoscale("off").  A manual scale factor can be
    set with set_scale(scale), this implies autoscale("off").  The
    scale factor (from the last autoscale event or from set_scale) can
    be obtained with get_scale().  Finally, an explicit autoscaling can
    be triggered with autoscale("now"), this is mainly useful before
    calling get_scale or before disabling further autoscaling.
    Finally, a relative scaling factor can be set with
    SetRelativeScaling(), it is multiplied to the usual scale factor
    (from autoscale or from set_scale).  This is probably only useful in
    connection with autoscaling.

    The radii of the atoms are obtained from the first of the following
    methods which work:
    
    1.  If the radii are specified using PrimiPlotter.set_radii(r),
        they are used.  Must be an array, or a single number.

    2.  If the atoms has a get_atomic_radii() method, it is used.  This is
        unlikely.

    3.  If the atoms has a get_atomic_numbers() method, the
        corresponding covalent radii are extracted from the
        ASE.ChemicalElements module.

    4.  If all else fails, the radius is set to 1.0 Angstrom.

    The atoms are colored using the first of the following methods
    which work.

    1.  If colors are explicitly set using PrimiPlotter.set_colors(),
        they are used.

    2.  If these colors are specified as a dictionary, the tags
        (from atoms.get_tags()) are used as an index into the
        dictionary to get the actual colors of the atoms.

    3.  If a color function has been set using
        PrimiPlotter.set_color_function(), it is called with the atoms
        as an argument, and is expected to return an array of colors.

    4.  If the atoms have a get_colors() method, it is used to get the
        colors.

    5.  If these colors are specified as a dictionary, the tags
        (from atoms.get_tags()) are used as an index into the
        dictionary to get the actual colors of the atoms.

    6.  If all else fails, the atoms will be white.

    The colors are specified as an array of colors, one color per
    atom.  Each color is either a real number from 0.0 to 1.0,
    specifying a grayscale (0.0 = black, 1.0 = white), or an array of
    three numbers from 0.0 to 1.0, specifying RGB values.  The colors
    of all atoms are thus a Numerical Python N-vector or a 3xN matrix.

    In cases 1a and 3a above, the keys of the dictionary are integers,
    and the values are either numbers (grayscales) or 3-vectors (RGB
    values), or strings with X11 color names, which are then
    translated to RGB values.  Only in case 1a and 3a are strings
    recognized as colors.

    Some atoms may be invisible, and thus left out of the plot.
    Invisible atoms are determined from the following algorithm.
    Unlike the radius or the coloring, all points below are tried and
    if an atom is invisible by any criterion, it is left out of the plot.

    1.  All atoms are visible.
    
    2.  If PrimiPlotter.set_invisible() has be used to specify invisible
        atoms, any atoms for which the value is non-zero becomes invisible.

    3.  If an invisiblility function has been set with
        PrimiPlotter.set_invisibility_function(), it is called with the
        atoms as argument.  It is expected to return an integer per
        atom, any non-zero value makes that atom invisible.

    4.  If a cut has been specified using set_cut, any atom outside the
        cut is made invisible.

    Note that invisible atoms are still included in the algorithm for
    positioning and scaling the plot.

    
    The following output devices are implemented.
    
    PostScriptFile(prefix):  Create PS files names prefix0000.ps etc.

    PnmFile(prefix):  Similar, but makes PNM files.

    GifFile(prefix):  Similar, but makes GIF files.

    JpegFile(prefix):  Similar, but makes JPEG files.

    X11Window():  Show the plot in an X11 window using ghostscript.

    Output devices writing to files take an extra optional argument to
    the constructor, compress, specifying if the output file should be
    gzipped.  This is not allowed for some (already compressed) file
    formats.

    Instead of a filename prefix, a filename containing a % can be
    used.  In that case the filename is expected to expand to a real
    filename when used with the Python string formatting operator (%)
    with the frame number as argument.  Avoid generating spaces in the
    file names: use e.g. %03d instead of %3d.  
    """
    def __init__(self, atoms, verbose=0, timing=0, interval=1, initframe=0):
        """

        Parameters to the constructor:

        atoms: The atoms to be plottet.

        verbose = 0:  Write progress information to stderr.

        timing = 0:  Collect timing information.

        interval = 1: If specified, a plot is only made every
        interval'th time update() is called.  Deprecated, normally you
        should use the interval argument when attaching the plotter to
        e.g. the dynamics.

        initframe = 0: Initial frame number, i.e. the number of the
        first plot.
        
        """
        self.atoms = atoms
        self.outputdevice = []
        self.angles = zeros(3, float)
        self.dims = (512, 512)
        self.verbose = verbose
        self.timing = timing
        self.totaltime = 0.0
        self.radius = None
        self.colors = None
        self.colorfunction = None
        self.n = initframe
        self.interval = interval
        self.skipnext = 0 # Number of calls to update before anything happens.
        self.a_scale = 1
        self.relativescale = 1.0
        self.invisible = None
        self.invisibilityfunction = None
        self.set_cut()   # No cut
        self.isparallel = 0
        self.logfile = None
        self.ownlogfile = False
        
    def set_output(self, device):
        self.outputdevice.append(device)
        device.set_dimensions(self.dims)
        device.set_owner(weakref.proxy(self))

    def set_dimensions(self, dims):
        "Set the size of the canvas (a 2-tuple)."
        if self.outputdevice:
            raise RuntimeError("Cannot set dimensions after an output device has been specified.")
        self.dims = dims
        
    def autoscale(self, mode):
        if mode == "on":
            self.a_scale = 1
        elif mode == "off":
            self.a_scale = 0
        elif mode == "now":
            coords = self._rotate(self.atoms.get_positions())
            radii = self._getradii()
            self._autoscale(coords, radii)
        else:
            raise ValueError, "Unknown autoscale mode: ",+str(mode)

    def set_scale(self, scale):
        self.autoscale("off")
        self.scale = scale

    def get_scale(self):
        return self.scale

    def set_relative_scale(self, rscale = 1.0):
        self.relativescale = rscale

    def plot(self):
        """Create a plot now.  Does not respect the interval timer.

        This method makes a plot unconditionally.  It does not look at
        the interval variable, nor is this plot taken into account in
        the counting done by the update() method if an interval
        variable was specified.
        """
        if self.timing:
            self._starttimer()
        self.log("PrimiPlotter: Starting plot at "
                 + time.strftime("%a, %d %b %Y %H:%M:%S"))
        colors = self._getcolors()
        invisible = self._getinvisible()
        coords = self._rotate(self._getpositions())
        radii = self._getradii()
        if self.a_scale:
            self._autoscale(coords,radii)
        scale = self.scale * self.relativescale
        coords = scale * coords
        center = self._getcenter(coords)
        offset = array(self.dims + (0.0,))/2.0 - center
        coords = coords + offset
        self.log("Scale is %f and size is (%d, %d)"
                 % (scale, self.dims[0], self.dims[1]))
        self.log("Physical size of plot is %f Angstrom times %f Angstrom"
                 % (self.dims[0] / scale, self.dims[1] / scale))

        self._verb("Sorting.")
        order = argsort(coords[:,2])
        coords = coords[order]  ### take(coords, order)
        radii = radii[order]    ### take(radii, order)
        colors = colors[order]  ### take(colors, order)
        invisible = invisible[order]  ### take(invisible, order)
        if self.isparallel:
            id = arange(len(coords))[order] ### take(arange(len(coords)), order)
        else:
            id = None
            
        radii = radii * scale
        selector = self._computevisibility(coords, radii, invisible, id)
        coords = compress(selector, coords, 0)
        radii = compress(selector, radii)
        colors = compress(selector, colors, 0)
        self._makeoutput(scale, coords, radii, colors)
        self.log("PrimiPlotter: Finished plotting at "
                 + time.strftime("%a, %d %b %Y %H:%M:%S"))
        self.log("\n\n")
        if self.timing:
            self._stoptimer()

    def _computevisibility(self, coords, rad, invisible, id, zoom = 1):
        xy = coords[:,:2]
        typradius = sum(rad) / len(rad)
        if typradius < 4.0:
            self.log("Refining visibility check.")
            if zoom >= 16:
                raise RuntimeError, "Cannot check visibility - too deep recursion."
            return self._computevisibility(xy*2, rad*2, invisible, id, zoom*2)
        else:
            self.log("Visibility(r_typ = %.1f pixels)" % (typradius,))
        dims = array(self.dims) * zoom
        maxr = int(ceil(max(rad))) + 2
        canvas = zeros((dims[0] + 4*maxr, dims[1] + 4*maxr), int8)
        # Atoms are only invisible if they are within the canvas, or closer
        # to its edge than their radius
        visible = (greater(xy[:,0], -rad) * less(xy[:,0], dims[0]+rad)
                   * greater(xy[:,1], -rad) * less(xy[:,1], dims[1]+rad)
                   * logical_not(invisible))
        # Atoms are visible if not hidden behind other atoms
        xy = floor(xy + 2*maxr + 0.5).astype(int)
        masks = {}
        for i in xrange(len(rad)-1, -1, -1):
            if (i % 100000) == 0 and i:
                self._verb(str(i))
            if not visible[i]:
                continue
            x, y = xy[i]
            r = rad[i]
            try:
                mask, invmask, rn = masks[r]
            except KeyError:
                rn = int(ceil(r))
                nmask = 2*rn+1
                mask = (arange(nmask) - rn)**2
                mask = less(mask[:,newaxis]+mask[newaxis,:], r*r).astype(int8)
                invmask = equal(mask, 0).astype(int8)
                masks[r] = (mask, invmask, rn)
            window = logical_or(canvas[x-rn:x+rn+1, y-rn:y+rn+1], invmask)
            hidden = alltrue(window.flat)
            if hidden:
                visible[i] = 0
            else:
                canvas[x-rn:x+rn+1, y-rn:y+rn+1] = logical_or(canvas[x-rn:x+rn+1, y-rn:y+rn+1], mask)
        self.log("%d visible, %d hidden out of %d" %
                   (sum(visible), len(visible) - sum(visible), len(visible)))
        return visible
        
    def _rotate(self, positions):
        self.log("Rotation angles: %f %f %f" % tuple(self.angles))
        mat = dot(dot(_rot(self.angles[2], 2),
                      _rot(self.angles[1], 1)),
                  _rot(self.angles[0]+pi, 0))
        return dot(positions, mat)

    def _getcenter(self, coords):
        return array((max(coords[:,0]) + min(coords[:,0]),
                      max(coords[:,1]) + min(coords[:,1]), 0.0)) / 2.0

    def _autoscale(self, coords, radii):
        x = coords[:,0]
        y = coords[:,1]
        maxradius = max(radii)
        deltax = max(x) - min(x) + 2*maxradius
        deltay = max(y) - min(y) + 2*maxradius
        scalex = self.dims[0] / deltax
        scaley = self.dims[1] / deltay
        self.scale = 0.95 * min(scalex, scaley)
        self.log("Autoscale: %f" % self.scale)

    def _makeoutput(self, scale, coords, radii, colors):
        for device in self.outputdevice:
            device.inform_about_scale(scale)
            device.plot(self.n, coords, radii, colors)
        self.n = self.n + 1


class ParallelPrimiPlotter(PrimiPlotter):
    """A version of PrimiPlotter for parallel ASAP simulations.

    Used like PrimiPlotter, but only the output devices on the master
    node are used.  Most of the processing is distributed on the
    nodes, but the actual output is only done on the master.  See the
    PrimiPlotter docstring for details.
    """
    def __init__(self, *args, **kwargs):
        apply(PrimiPlotter.__init__, (self,)+args, kwargs)
        self.isparallel = 1
        import Scientific.MPI
        self.MPI = Scientific.MPI
        self.mpi = Scientific.MPI.world
        if self.mpi is None:
            raise RuntimeError, "MPI is not available."
        self.master = self.mpi.rank == 0
        self.mpitag = 42   # Reduce chance of collision with other modules.
        
    def set_output(self, device):
        if self.master:
            PrimiPlotter.set_output(self, device)

    def set_log(self, log):
        if self.master:
            PrimiPlotter.set_log(self, log)

    def _getpositions(self):
        realpos = self.atoms.get_positions()
        ghostpos = self.atoms.GetGhostCartesianPositions()
        self.numberofrealatoms = len(realpos)
        self.numberofghostatoms = len(ghostpos)
        return concatenate((realpos, ghostpos))

    def _getatomicnumbers(self):
        realz = self.atoms.get_atomic_numbers()
        ghostz = self.atoms.GetGhostAtomicNumbers()
        return concatenate((realz, ghostz))

    def _getradius(self):
        r = PrimiPlotter._getradius(self)
        if len(r) == self.numberofrealatoms + self.numberofghostatoms:
            # Must have calculated radii from atomic numbers
            return r
        else:
            assert len(r) == self.numberofrealatoms
            # Heuristic: use minimum r for the ghosts
            ghostr = min(r) * ones(self.numberofghostatoms, float)
            return concatenate((r, ghostr))

    def _getcenter(self, coords):
        # max(x) and min(x) only works for rank-1 arrays in Numeric version 17.
        maximal = maximum.reduce(coords[:,0:2])
        minimal = minimum.reduce(coords[:,0:2])
        recvmax = zeros(2, maximal.typecode())
        recvmin = zeros(2, minimal.typecode())
        self.mpi.allreduce(maximal, recvmax, self.MPI.max)
        self.mpi.allreduce(minimal, recvmin, self.MPI.min)
        maxx, maxy = recvmax
        minx, miny = recvmin
        return array([maxx + minx, maxy + miny, 0.0]) / 2.0

    def _computevisibility(self, xy, rad, invisible, id, zoom = 1):
        # Find visible atoms, allowing ghost atoms to hide real atoms.
        v = PrimiPlotter._computevisibility(self, xy, rad, invisible, id, zoom)
        # Then remove ghost atoms
        return v * less(id, self.numberofrealatoms)

    def _autoscale(self, coords, radii):
        self._verb("Autoscale")
        n = len(self.atoms)
        x = coords[:n,0]
        y = coords[:n,1]
        assert len(x) == len(self.atoms)
        maximal = array([max(x), max(y), max(radii[:n])])
        minimal = array([min(x), min(y)])
        recvmax = zeros(3, maximal.typecode())
        recvmin = zeros(2, minimal.typecode())
        self.mpi.allreduce(maximal, recvmax, self.MPI.max)
        self.mpi.allreduce(minimal, recvmin, self.MPI.min)
        maxx, maxy, maxradius = recvmax
        minx, miny = recvmin
        deltax = maxx - minx + 2*maxradius
        deltay = maxy - miny + 2*maxradius
        scalex = self.dims[0] / deltax
        scaley = self.dims[1] / deltay
        self.scale = 0.95 * min(scalex, scaley)
        self.log("Autoscale: %f" % self.scale)

    def _getcolors(self):
        col = PrimiPlotter._getcolors(self)
        nghost = len(self.atoms.GetGhostCartesianPositions())
        newcolshape = (nghost + col.shape[0],) + col.shape[1:]
        newcol = zeros(newcolshape, col.typecode())
        newcol[:len(col)] = col
        return newcol
    
    def _makeoutput(self, scale, coords, radii, colors):
        if len(colors.shape) == 1:
            # Greyscales
            ncol = 1
        else:
            ncol = colors.shape[1]  # 1 or 3.
            assert ncol == 3  # RGB values
        # If one processor says RGB, all must convert
        ncolthis = array([ncol])
        ncolmax = zeros((1,), ncolthis.typecode())
        self.mpi.allreduce(ncolthis, ncolmax, self.MPI.max)
        ncolmax = ncolmax[0]
        if ncolmax > ncol:
            assert ncol == 1
            colors = colors[:,newaxis] + zeros(ncolmax)[newaxis,:]
            ncol = ncolmax
            assert colors.shape == (len(coords), ncol)
        # Now send data from slaves to master
        data = zeros((len(coords)+1, 4+ncol), float)
        data[:-1,:3] = coords
        data[:-1,3] = radii
        data[-1,-1] = 4+ncol  # Used to communicate shape
        if ncol == 1:
            data[:-1,4] = colors
        else:
            data[:-1,4:] = colors
        if not self.master:
            self.mpi.send(data, 0, self.mpitag)
        else:
            total = [data[:-1]]  # Last row is the dimensions.
            n = len(coords)
            colsmin = colsmax = 4+ncol
            for proc in range(1, self.mpi.size):
                self._verb("Receiving from processor "+str(proc))
                fdat = self.mpi.receive(float, proc, self.mpitag)[0]
                fdat.shape = (-1, fdat[-1])
                fdat = fdat[:-1]  # Last row is the dimensions.
                total.append(fdat)
                n = n + len(fdat)
                if fdat.shape[1] < colsmin:
                    colsmin = fdat.shape[1]
                if fdat.shape[1] > colsmax:
                    colsmax = fdat.shape[1]
            self._verb("Merging data")
            # Some processors may have only greyscales whereas others
            # may have RGB.  That will cause difficulties.
            trouble = colsmax != colsmin
            data = zeros((n, colsmax), float)
            if trouble:
                assert data.shape[1] == 7
            else:
                assert data.shape[1] == 7 or data.shape[1] == 5
            i = 0
            for d in total:
                if not trouble or d.shape[1] == 7:
                    data[i:i+len(d)] = d
                else:
                    assert d.shape[1] == 5
                    data[i:i+len(d), :5] = d
                    data[i:i+len(d), 5] = d[4]
                    data[i:i+len(d), 6] = d[4]
                i = i + len(d)
            assert i == len(data)
            # Now all data is on the master
            self._verb("Sorting merged data")
            order = argsort(data[:,2])
            data = data[order]   ### take(data, order)
            coords = data[:,:3]
            radii = data[:,3]
            if data.shape[1] == 5:
                colors = data[:,4]
            else:
                colors = data[:,4:]
            PrimiPlotter._makeoutput(self, scale, coords, radii, colors)
    
class _PostScriptDevice:
    """PostScript based output device."""
    offset = (0,0)   # Will be changed by some classes
    def __init__(self):
        self.scale = 1
        self.linewidth = 1
        self.outline = 1
        
    def set_dimensions(self, dims):
        self.dims = dims

    def set_owner(self, owner):
        self.owner = owner
        
    def inform_about_scale(self, scale):
        self.linewidth = 0.1 * scale

    def set_outline(self, value):
        self.outline = value
        return self   # Can chain these calls in set_output()
        
    def plot(self, *args, **kargs):
        self.Doplot(self.PSplot, *args, **kargs)
        
    def plotArray(self, *args, **kargs):
        self.Doplot(self.PSplotArray, *args, **kargs)
        
    def PSplot(self, file, n, coords, r, colors, noshowpage=0):
        xy = coords[:,:2]
        assert(len(xy) == len(r) and len(xy) == len(colors))
        if len(colors.shape) == 1:
            gray = 1
        else:
            gray = 0
            assert(colors.shape[1] == 3)
        file.write("%!PS-Adobe-2.0\n")
        file.write("%%Creator: Primiplot\n")
        file.write("%%Pages: 1\n")        
        file.write("%%%%BoundingBox: %d %d %d %d\n" %
                   (self.offset + (self.offset[0] + self.dims[0],
                                   self.offset[1] + self.dims[1])))
        file.write("%%EndComments\n")
        file.write("\n")
        file.write("% Enforce BoundingBox\n")
        file.write("%d %d moveto %d 0 rlineto 0 %d rlineto -%d 0 rlineto\n" %
                   ((self.offset + self.dims + (self.dims[0],))))
        file.write("closepath clip newpath\n\n")
        file.write("%f %f scale\n" % (2*(1.0/self.scale,)))
        file.write("%d %d translate\n" % (self.scale * self.offset[0],
                                          self.scale * self.offset[1]))
        file.write("\n")
        if gray:
            if self.outline:
                file.write("/circ { 0 360 arc gsave setgray fill grestore stroke } def\n")
            else:
                file.write("/circ { 0 360 arc setgray fill } def\n")
        else:
            if self.outline:
                file.write("/circ { 0 360 arc gsave setrgbcolor fill grestore stroke } def\n")
            else:
                file.write("/circ { 0 360 arc setrgbcolor fill } def\n")
        file.write("%f setlinewidth 0.0 setgray\n" %
                   (self.linewidth * self.scale,))
        
        if gray:
            data = zeros((len(xy), 4), float)
            data[:,0] = colors
            data[:,1:3] = (self.scale * xy)
            data[:,3] = (self.scale * r)
            for point in data:
                file.write("%.3f %.2f %.2f %.2f circ\n" % tuple(point))
        else:
            data = zeros((len(xy), 6), float)
            data[:,0:3] = colors
            data[:,3:5] = (self.scale * xy)
            data[:,5] = (self.scale * r)
            for point in data:
                file.write("%.3f %.3f %.3f %.2f %.2f %.2f circ\n" % tuple(point))
        if not noshowpage:
            file.write("showpage\n")
            
    def PSplotArray(self, file, n, data, noshowpage=0):
        assert(len(data.shape) == 3)
        assert(data.shape[0] == self.dims[1] and data.shape[1] == self.dims[0])
        data = clip((256*data).astype(int), 0, 255)
        file.write("%!PS-Adobe-2.0\n")
        file.write("%%Creator: Fieldplotter\n")
        file.write("%%Pages: 1\n")        
        file.write("%%%%BoundingBox: %d %d %d %d\n" %
                   (self.offset + (self.offset[0] + self.dims[0],
                                   self.offset[1] + self.dims[1])))
        file.write("%%EndComments\n")
        file.write("\n")
        file.write("%d %d translate\n" % self.offset)
        file.write("%f %f scale\n" % self.dims)
        file.write("\n")
        file.write("% String holding a single line\n")
        file.write("/pictline %d string def\n" %(data.shape[1]*data.shape[2],))
        file.write("\n")
        file.write("%d %d 8\n" % self.dims)
        file.write("[%d 0 0 %d 0 0]\n" % self.dims)
        file.write("{currentfile pictline readhexstring pop}\n")
        file.write("false %d colorimage\n" % (data.shape[2],))
        file.write("\n")
        s = ""
        for d in data.flat:
            s += ("%02X" % d)
            if len(s) >= 72:
                file.write(s+"\n")
                s = ""
        file.write(s+"\n")
        file.write("\n")
        if not noshowpage:
            file.write("showpage\n")
            
class _PostScriptToFile(_PostScriptDevice):
    """Output device for PS files."""
    compr_suffix = None
    def __init__(self, prefix, compress = 0):
        self.compress = compress
        if "'" in prefix:
            raise ValueError, "Filename may not contain a quote ('): "+prefix
        if "%" in prefix:
            # Assume the user knows what (s)he is doing
            self.filenames = prefix
        else:
            self.filenames = prefix + "%04d" + self.suffix
            if compress:
                if self.compr_suffix is None:
                    raise RuntimeError, "Compression not supported."
                self.filenames = self.filenames + self.compr_suffix
        _PostScriptDevice.__init__(self)

class PostScriptFile(_PostScriptToFile):
    suffix = ".ps"
    compr_suffix = ".gz"
    offset = (50,50)
    # Inherits __init__

    def Doplot(self, plotmethod, n, *args, **kargs):
        filename = self.filenames % (n,)
        self.owner.log("Output to PostScript file "+filename)
        if self.compress:
            file = os.popen("gzip > '"+filename+"'", "w")
        else:
            file = open(filename, "w")
        apply(plotmethod, (file, n)+args, kargs)
        file.close()

class _PS_via_PnmFile(_PostScriptToFile):
    gscmd = "gs -q -sDEVICE=pnmraw -sOutputFile=- -dDEVICEWIDTH=%d -dDEVICEHEIGHT=%d - "
    # Inherits __init__

    def Doplot(self, plotmethod, n, *args, **kargs):
        filename = self.filenames % (n,)
        self.owner.log("Output to bitmapped file " + filename)
        cmd = self.gscmd + self.converter
        if self.compress:
            cmd = cmd + "| gzip "
            
        cmd = (cmd+" > '%s'") % (self.dims[0], self.dims[1], filename)
        file = os.popen(cmd, "w")
        apply(plotmethod, (file, n)+args, kargs)
        file.close()

class PnmFile(_PS_via_PnmFile):
    suffix = ".pnm"
    compr_suffix = ".gz"
    converter = ""

class GifFile(_PS_via_PnmFile):
    suffix = ".gif"
    converter = "| ppmquant -floyd 256 2>/dev/null | ppmtogif 2>/dev/null"

class JpegFile(_PS_via_PnmFile):
    suffix = ".jpeg"
    converter = "| ppmtojpeg --smooth=5"
    
class X11Window(_PostScriptDevice):
    """Shows the plot in an X11 window."""
    #Inherits __init__
    gscmd = "gs -q -sDEVICE=x11 -dDEVICEWIDTH=%d -dDEVICEHEIGHT=%d -r72x72 -"
    def Doplot(self, plotmethod, n, *args, **kargs):
        self.owner.log("Output to X11 window")
        try:
            file = self.pipe
            self.pipe.write("showpage\n")
        except AttributeError:
            filename = self.gscmd % tuple(self.dims)
            file = os.popen(filename, "w")
            self.pipe = file
        kargs["noshowpage"] = 1
        apply(plotmethod, (file, n)+args, kargs)
        file.write("flushpage\n")
        file.flush()

# Helper functions
def _rot(v, axis):
    ax1, ax2 = ((1, 2), (0, 2), (0, 1))[axis]
    c, s = cos(v), sin(v)
    m = zeros((3,3), float)
    m[axis,axis] = 1.0
    m[ax1,ax1] = c
    m[ax2,ax2] = c
    m[ax1,ax2] = s
    m[ax2,ax1] = -s
    return m

def _colorsfromdict(dict, cls):
    """Extract colors from dictionary using cls as key."""
    assert(type(dict) == type({}))
    # Allow local modifications, to replace strings with rgb values.
    dict = dict.copy()  
    isgray, isrgb = 0, 0
    for k in dict.keys():
        v = dict[k]
        if type(v) == type("string"):
            v = color_table[v]
            dict[k] = v
        try:
            if len(v) == 3:
                isrgb = 1 # Assume it is an RGB value
                if not hasattr(v, "shape"):
                    dict[k] = array(v)   # Convert to array
            else:
                raise RuntimeError, "Unrecognized color object "+repr(v)
        except TypeError:
            isgray = 1 # Assume it is a number
    if isgray and isrgb:
        # Convert all to RGB
        for k in dict.keys():
            v = dict[k]
            if not hasattr(v, "shape"):
                dict[k] = v * ones(3, float)
    # Now the dictionary is ready
    if isrgb:
        colors = zeros((len(cls),3), float)
    else:
        colors = zeros((len(cls),), float)
    for i in xrange(len(cls)):
        colors[i] = dict[cls[i]]
    return colors

