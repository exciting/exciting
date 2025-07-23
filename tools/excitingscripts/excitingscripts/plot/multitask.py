"""Please, check https://www.exciting-code.org/the-python-script-plot.multitask
to better understand how to use this script"""

import matplotlib
import matplotlib.ticker as ptk
import matplotlib.pyplot as plt
from scipy import fftpack
import numpy
import argparse
from enum import Enum

if matplotlib.__version__.split(".")[0]=="2":
    matplotlib.style.use('classic')

#-------------------------------------------------------------------------------
# Global constants
speed_light = 137.03599911
Ha_to_eV = 27.211396132

#-------------------------------------------------------------------------------

class Handle_complex(Enum):
    """
    enum to treat handle_complex
    """
    no = 1 # no need to handle complex
    full = 2 # full treatment of complex numbers
    real = 3 # only real part
    imag = 4 # only imaginary part
    abs = 5 # just consider the absolute value

class Option_preprocess(Enum):
    """
    enum to store which action to take as preprocessing
    """
    nothing = 0 # nothing
    get_efield = 1 # get electric field
    get_eps = 2 # get the dielectric function
    fourier_transform = 3
    add = 4
    sub = 5

#-------------------------------------------------------------------------------

def fft(t, f, wcut=0.0):
    """
    FFT - from time domain to (angular) frequency:
    :param t: time (array)
    :param f: function (array) f(t) to be fourier-transformed
    :param wcut: cut-off frequency (in Ha) for the low pass filter
    :return w, F: Tuple containing angular frequencies (in Ha) and the
        fourier transform
    """

    nlines = len(t)
    dt = t[1]-t[0]
    w = 2*(numpy.pi)*fftpack.fftfreq(nlines,dt)
    filter = numpy.exp(numpy.multiply( t, -wcut ))
    F = numpy.conj(numpy.multiply(fftpack.fft(numpy.multiply(filter,f)), dt))
    w,F = list(zip(*sorted(zip(w,F))))

    return numpy.array(w),numpy.array(F)

#-------------------------------------------------------------------------------
def read_file(file, columns_to_plot, nlines_skip=0, scale=[1,1],
        handle_complex=Handle_complex.no):
    """
    FFT time to (angular) frequency:
    :param str file: name of the file to be read
    :param columns_to_plot: array of 2 integers with the columns to be read from file
    :param nlines_skip: number of lines to skip
    :param scale: array of two real numbers to scale x and y
    :param handle_complex: enum with the treatment of complex numbers
    :return x, y: Tuple x and y arrays as read from file with minimal processing
    """
    # Type checking
    if not isinstance(handle_complex, Handle_complex):
        raise TypeError('Handle_complex must be an instance of Handle_complex Enum')

    ix = columns_to_plot[0]
    iy = columns_to_plot[1]

    if handle_complex == Handle_complex.no:
        x, y = numpy.loadtxt( fname=file, skiprows=nlines_skip, usecols=(ix,iy),
            dtype=float, unpack=True )
    else:
        x = numpy.loadtxt(fname=file, skiprows=nlines_skip, usecols=ix,
            dtype=float)
        y = numpy.loadtxt(fname=file, skiprows=nlines_skip, usecols=iy,
            dtype=complex)
        if handle_complex == Handle_complex.abs:
            y = abs(y)
        elif handle_complex == Handle_complex.real:
            y = y.real
        elif handle_complex == Handle_complex.imag:
            y = y.imag

    return numpy.multiply(scale[0],x), numpy.multiply(scale[1],y)


#-------------------------------------------------------------------------------
def preformat_plot():
    """
    Set up some pyplot parameters, pre-formatting the plot
    """
    params = {'ytick.minor.size': 6,
              'xtick.major.pad': 8,
              'ytick.major.pad': 4,
              'patch.linewidth': 2.,
              'axes.linewidth': 2.,
              'lines.linewidth': 1.8,
              'lines.markersize': 8.0,
              'axes.formatter.limits': (-4, 6)}

    plt.rcParams.update(params)
    plt.subplots_adjust(left=0.21, right=0.93,
                        bottom=0.18, top=0.88,
                        wspace=None, hspace=None)
    plt.grid(linestyle='--')


def format_plot(ax, xmin, xmax, ymin, ymax, xlim, ylim,
        semilog,
        xlabel, ylabel,
        hide_legend, legend_position):

    for line in ax.get_xticklines() + ax.get_yticklines():
        line.set_markersize(6)
        line.set_markeredgewidth(2)
    fonttick = 16
    plt.xticks(size=fonttick)
    plt.yticks(size=fonttick)
    plt.tick_params(direction="in", width=2, length=5)
    dxx = (xmax-xmin)/18 ; dyy = (ymax-ymin)/15
    xmin = xmin-dxx ; xmax = xmax+dxx
    ymin = ymin-dyy ; ymax = ymax+dyy
    if not(semilog):
        yfmt = ptk.ScalarFormatter()
        ax.yaxis.set_major_formatter(yfmt)
    else:
        plt.yscale('log')
    if xlim == []:
        ax.set_xlim(xmin,xmax)
    else:
        ax.set_xlim(xlim[0],xlim[1])
    if ylim == []:
        ax.set_ylim(ymin,ymax)
    else:
        ax.set_ylim(ylim[0],ylim[1])
    if not(hide_legend):
        plt.legend(loc=legend_position,borderaxespad=.8,numpoints=1,fontsize=14)
    ax.xaxis.set_major_locator(ptk.MaxNLocator(7))
    ax.set_axisbelow(True)
    fontlabel = 20
    ax.text(0.5,-0.13,xlabel,size=fontlabel,
            transform=ax.transAxes,ha='center',va='center',rotation=0)
    plt.ylabel(ylabel,size=fontlabel)

#-------------------------------------------------------------------------------
def convert_args_to_dict(args):
    """
    Parsed arguments are converted to a dictionary
    :param args: command line arguments
    """
    dict = {
        'ylabel': args.ylabel,
        'xlabel': args.xlabel,
        'files': args.f,
        'columns': args.k,
        'skip_lines': args.s,
        'captions': args.c,
        'legend_position': args.legend_position,
        'hide_legend': args.hide_legend,
        'pre_process': args.preprocess,
        'fourier_transform': args.fourier,
        'get_efield': args.get_efield,
        'get_eps': args.get_eps,
        'output': args.o,
        'scale': args.scale,
        'add': args.add,
        'sub': args.sub,
        'wcut': args.wcut,
        'handle_complex': Handle_complex[args.handle_complex],
        'xlim': args.xlim,
        'ylim': args.ylim,
        'semilog': args.semilog,
        'jind': args.jind,
        'nexc': args.nexc,
        'x': args.x,
        'y': args.y,
        'z': args.z,
        'imag_eps':args.imag_eps,
        'real_eps':args.real_eps,
        'xx': args.xx,
        'xy': args.xy,
        'xz': args.xz,
        'yx': args.yx,
        'yy': args.yy,
        'yz': args.yz,
        'zx': args.zx,
        'zy': args.zy,
        'zz': args.zz,
        'show': args.show,
        'plot_name': args.plot_name
    }
    return dict
#-------------------------------------------------------------------------------
def parse_input():
    """
    Function to parse the arguments from the command line
    """
    parser = argparse.ArgumentParser()

    # General options, common for plotting and preprocessing
    # -f: name of input files
    parser.add_argument('-f', metavar = '--files', type = str, nargs = '+',
        required=True, action='append' )
    # -s: how many lines to skip (ignore), before reading the data
    parser.add_argument('-s', metavar = '--skip', type = int, default = [],
        nargs = '*', required=False, action='append' )
    # -k: columns to read (is overwritten by options like: --x, or --xy)
    parser.add_argument('-k', metavar = '--columns', type = int, default = [],
        nargs = '*', required=False, action='append' )
    # scale: scale the data read (x and y)
    parser.add_argument('--scale', required=False, nargs = 2, default=[1.0, 1.0],
        type = float, action='store')

    # Options for preprocessing
    # preprocess: it means no plot, only preparations for a future call
    #       Options for pre-process: fourier (with wcut), get_efield, get_eps
    #       add, sub, scale, and -o for the output
    parser.add_argument('--preprocess', required=False, action='store_true')
    # fourier: obtain the fast fourier transform applying a smoothing filter,
    #       with wcut as parameter for the smoothing
    parser.add_argument('--fourier', required=False, action='store_true')
    parser.add_argument('--wcut', required=False, nargs = 1, default=[0.0],
        type = float, action='store')
    # get_efield: calculates the electric field using VECTOR_POTENTIAL.OUT, the file with
    #       the vector potential
    parser.add_argument('--get_efield', required=False, action='store_true')
    # get_eps: obtain the dielectric function, using VECTOR_POTENTIAL.OUT and CURRENT.OUT
    parser.add_argument('--get_eps', required=False, action='store_true')
    # the next options are for the tensor components of the dielectric function
    parser.add_argument('--xx', required=False, action='store_true')
    parser.add_argument('--xy', required=False, action='store_true')
    parser.add_argument('--xz', required=False, action='store_true')
    parser.add_argument('--yx', required=False, action='store_true')
    parser.add_argument('--yy', required=False, action='store_true')
    parser.add_argument('--yz', required=False, action='store_true')
    parser.add_argument('--zx', required=False, action='store_true')
    parser.add_argument('--zy', required=False, action='store_true')
    parser.add_argument('--zz', required=False, action='store_true')
    # -o: name of the file used as output
    parser.add_argument('-o', metavar = '--output', required=False,
        default = 'output.txt', type = str, action='store')
    # add: add the data read from two files (only the "y" component is added)
    parser.add_argument('--add', required=False, action='store_true')
    # sub: same as add, but subtracts, first y minus second y
    parser.add_argument('--sub', required=False, action='store_true')


    # Options for plotting
    # jind: plot the current density (this triggers the x- and y-labels)
    parser.add_argument('--jind', required=False, action='store_true')
    # --x, --y, --z: the component to be plotted (only one each time)
    parser.add_argument('--x', required=False, action='store_true')
    parser.add_argument('--y', required=False, action='store_true')
    parser.add_argument('--z', required=False, action='store_true')
    # nexc: plot the number of excitations (this triggers the x- and y-labels)
    parser.add_argument('--nexc', required=False, action='store_true')
    # imag_eps, real_eps: plot the imaginary or real part of the dielectric
    #       function (this triggers the x- and y-labels)
    parser.add_argument('--imag_eps', required=False, action='store_true')
    parser.add_argument('--real_eps', required=False, action='store_true')
    # xlim, ylim: format the plot limits for x and y axes
    parser.add_argument('--xlim', required=False, nargs = 2, default=[],
        type = float, action='store')
    parser.add_argument('--ylim', required=False, nargs = 2, default=[],
        type = float, action='store')
    # xlabel, ylabel: labels for x and y axes
    parser.add_argument('--ylabel', type = str, nargs = 1,
        default=[], required=False, action='store')
    parser.add_argument('--xlabel', type = str, nargs = 1,
        default=[], required=False, action='store')
    # -c: caption for the data read (to be shown in the legend box)
    parser.add_argument('-c', metavar = '--caption', type = str, default = [],
        nargs = '*', required=False, action='append')
    # legend_position: where to put the legend
    parser.add_argument('--legend_position', type = str,
        default = 'upper right', nargs = '?', required=False, action='store',
        choices = ['upper right','lower right','upper left','lower left'])
    # hide_legend: do not show the legend
    parser.add_argument('--hide_legend', required=False, action='store_true')
    # handle_complex: if/how to treat the complex numbers (for y-data)
    parser.add_argument('--handle_complex', type = str, required=False,
        nargs = '?', action='store', default = 'no',
        choices = ['no','full','imag','real','abs'])
    # semilog: if the plot is semilog
    parser.add_argument('--semilog', required=False, action='store_true')
    parser.add_argument('--show', required=False, action='store_true')
    parser.add_argument('--plot_name', type=str, nargs=1,
                        default=['PLOT'], required=False, action='store')

    args = parser.parse_args()

    # store the arguments parsed as a dictionary
    options = convert_args_to_dict(args)

    return options
#-------------------------------------------------------------------------------
def arrange_as_stack(list_of_lists):
    """
    Function to convert list of lists to a "stack" (only a list)
    :param list_of_lists: list with lists
    """
    stack = []
    for list_ in list_of_lists :
        for element in list_:
            stack.append( element )
    return stack

#-------------------------------------------------------------------------------
def set_implicit_options( options ):
    # Better format options['files']
    options['files'] = arrange_as_stack( options['files'] )
    # If we have only one file, don't show the legend box
    if len(options['files']) == 1:
        options['hide_legend'] = True

    # Check if we need the default options for options['columns']
    if options['columns'] == []:
        # Usual cases
        if not( options['get_eps']):
            columns = [0,1]
            if options['nexc']:
                columns = [0,2]
            elif options['x']:
                columns = [0,1] if not( options['get_efield'] ) else [0,2]
            elif options['y']:
                columns = [0,2] if not( options['get_efield'] ) else [0,4]
            elif options['z']:
                columns = [0,3] if not( options['get_efield'] ) else [0,6]
            options['columns'] = [columns]*len(options['files'])
        # Special case (get_eps)
        else:
            # xx, xy, ..., zz: first component tells which column
            # of CURRENT.OUT is to be read; the second component,
            # the one of VECTOR_POTENTIAL.OUT
            options['columns'] = [[0,2],[0,1]]
            if options['xx']:
                options['columns'] = [[0,2],[0,1]]
            if options['xy']:
                options['columns'] = [[0,4],[0,1]]
            if options['xz']:
                options['columns'] = [[0,6],[0,1]]
            if options['yx']:
                options['columns'] = [[0,2],[0,2]]
            if options['yy']:
                options['columns'] = [[0,4],[0,2]]
            if options['yz']:
                options['columns'] = [[0,6],[0,2]]
            if options['zx']:
                options['columns'] = [[0,2],[0,3]]
            if options['zy']:
                options['columns'] = [[0,4],[0,3]]
            if( options['zz'] ):
                options['columns'] = [[0,6],[0,3]]

    # Default options for options['skip_lines']
    if options['skip_lines'] == []:
        skip_default = 1  if ( options['nexc'] ) else 0
        options['skip_lines'] = [ skip_default ]*len(options['files'])
    else:
        options['skip_lines'] = arrange_as_stack( options['skip_lines'] )

    # Default options for options['captions']
    if options['captions'] == []:
        options['captions'] = options['files'][:]
    else:
        options['captions'] = arrange_as_stack( options['captions'] )

    # Default options for the axes labels
    if options['xlabel'] == []:
        options['xlabel'] = ['Time [a.u.]']
    if options['ylabel'] == []:
        if options['jind']:
            options['ylabel'] = ['Current Density [a.u.]']
        elif options['nexc']:
            options['ylabel'] = ['$N_{exc}(t)$']
        else: 
            options['ylabel'] = ['y']

    # In case we have the dielectric function
    if options['imag_eps'] or options['real_eps']:
        options['xlabel'] = ['Energy [eV]']
        options['scale'][0] = Ha_to_eV
        options['ylabel'] = ['Im($\\varepsilon$)'] if options['imag_eps'] \
            else ['Re($\\varepsilon$)']
        options['handle_complex'] = Handle_complex['imag'] if \
            options['imag_eps'] else Handle_complex['real']

#-------------------------------------------------------------------------------
def preprocess( x, y, output, option_preproc=Option_preprocess.nothing, \
        wcut=0.0, diagonal_component=True):
    """
    Preprocess x and y, depending on the desired options
    :param x: list of arrays, data to be outputed but needed for preprocessing y
    :param y: list of arrays, data to preprocess
    :param str output: name of output file
    :param option_preproc: required kind of preprocessing
    :param wcut: smoothing parameter for the fourier transform
    :param diagonal_component: to obtain the dielectric tensor, we need to
        know if the desired component belongs to the diagonal
    """
    # Obtain the electric field from the vector potential:
    # E = -1/c dA/dt
    if option_preproc == Option_preprocess.get_efield:
        y = numpy.multiply( numpy.gradient( y[0], x[0] ), -1./speed_light )
        x = x[0]
    # Obtain the dielectric function as described in
    # http://exciting-code.org/oxygen-pump-probe-spectroscopy
    elif option_preproc == Option_preprocess.get_eps:
        # electric field
        e = numpy.multiply( numpy.gradient( y[0], x[0] ), -1./speed_light )
        # fourier transform of e
        w, E = fft( x[0], e, wcut )
        # fourier transform of the current density
        _, J = fft( x[1], y[1], wcut )
        # optical conductivity
        sigma = numpy.divide( J, E, where=E!=0, \
            out=numpy.zeros(E.shape,dtype=E.dtype) )
        # dielectric tensor
        if diagonal_component:
            eps = 1 + 4*numpy.pi*1j*numpy.divide( sigma, w, \
                where=w!=0, out=numpy.zeros(w.shape,dtype=complex) )
        else:
            eps = 4*numpy.pi*1j*numpy.divide( sigma, w, \
                where=w!=0, out=numpy.zeros(w.shape,dtype=complex) )
        x = w[:]
        y = eps[:]
    # Obtain the fourier transform (with a smoothing filter)
    elif option_preproc == Option_preprocess.fourier_transform:
        x, y = fft( x[0], y[0], wcut )
    # Add two data
    elif option_preproc == Option_preprocess.add:
        x = x[0]
        y = numpy.add( y[0], y[1] )
    # Subtract to data
    elif option_preproc == Option_preprocess.sub:
        x = x[0]
        y = numpy.subtract( y[0], y[1] )
    # Now, if the only preprocessing is scaling (done when reading the input)
    else:
        x = x[0]
        y = y[0]
    # Write to output
    ff = open( output, 'w' )
    for xx, yy in zip(x,y):
        ff.write( '{0:30.20f} {1:60.20f}\n'.format(xx,yy) )
    ff.close()

#-------------------------------------------------------------------------------
def sanity_checks( options ):
    """
    Function to make some sanity checks of the command line arguments
    :param options: dictionary with the command line arguments
    """

    # positive numbers given for number of columns, scaling factor, skip_lines
    # and wcut?
    # columns
    problem = False
    for cols in options['columns']:
        if cols[0]<0 or cols[1]<0:
            problem = True
            break
    if problem:
        raise RuntimeError( 'The column indexes of data to read must >=0 ' )
    # scaling factor
    for sc in options['scale']:
        if sc <= 0:
            problem = True
            break
    if problem:
        raise RuntimeError( 'Scaling factor must be positive' )
    # skip_lines
    for sk in options['skip_lines']:
        if ( sk < 0 ):
            problem = True
            break
    if problem:
        raise RuntimeError( 'Number of lines to skip must be >= 0' )
    # wcut
    if options['wcut'][0] < 0:
        raise RuntimeError('Cutoff frequency for the smoothing must be >= 0.')


    # number of columns and files are compatible?
    error_message = 'There must be 2 columns for each file. You need to\
        specify them using -k or --x, ..., --z, --xx, ..., --zz depending \
        on the case.'
    if len( options['files'] ) != len( options['columns'] ):
        raise RuntimeError(error_message)
    for cols in options['columns']:
        if len(cols) != 2:
            raise RuntimeError( error_message )
    # number of captions and files are compatible?
    if len( options['captions'] ) != len( options['files'] ):
        raise RuntimeError('Number of files does not match number of \
            captions given')

    # number of skip_lines and files are compatible?
    if len( options['skip_lines'] ) != len( options['files'] ):
        raise RuntimeError('Number of files is not equal \
            times of lines-to-skip')

    # only one preprocess option?
    is_efield = options['get_efield']
    is_eps = options['get_eps']
    is_fft = options['fourier_transform']
    is_add = options['add']
    is_sub = options['sub']
    number_of_true_elements = sum([is_efield, is_eps, is_fft, is_add, is_sub ])
    not_ok_preproc = ( number_of_true_elements > 1 )
    error_message = 'It is not possible to perform more than one preprocessing\
        action.'
    if not_ok_preproc:
        raise RuntimeError( error_message )
    # preprocess options required with preprocess argument?
    if number_of_true_elements == 1 and not options['pre_process']:
        raise RuntimeError('Option --preprocess is required.')
    # number of files to preprocess is consistent?
    nfiles = len(options['files'])
    # only one file for get_efield and fourier_transform
    if (is_efield or is_fft) and (nfiles != 1):
        string = 'get_efield' if is_efield else 'fourier'
        raise RuntimeError('Only one file for the preprocess option: '+ string)
    # two files for get_eps, add and sub
    if (is_eps or is_add or is_sub) and (nfiles != 2):
        if is_eps:
            string = 'get_eps'
        elif is_add:
            string = 'add'
        else:
            string = 'sub'
        raise RuntimeError('Two files are needed for the preprocess option: '\
            + string )

    # no more than one of these options: x, y, z, xx, ..., zz?
    options_to_check = [ options['x'], options['y'], options['z'], \
        options['xx'], options['xy'], options['xz'], \
        options['yx'], options['yy'], options['yz'], \
        options['zx'], options['zy'], options['zz'] ]
    too_many_options = ( sum( options_to_check ) > 1 )
    error_message = 'More than one option (--x, --y, --z, --xx, --xy, ..., \
        --zz) has been required.'
    if( too_many_options ):
        raise RuntimeError( error_message )

    # no more than one option: jind, nexc, imag_eps, real_eps?
    options_to_check = [ options['jind'], options['nexc'], options['imag_eps'],\
        options['real_eps'] ]
    too_many_options = ( sum( options_to_check ) > 1 )
    if( too_many_options ):
        raise RuntimeError('More than one option among the following ones has \
            been required: jind, nexc, imag_eps, real_eps')

#-------------------------------------------------------------------------------
def main():
    args = parse_input() # get a dictionary with the parsed arguments
    set_implicit_options( args )
    sanity_checks( args )

    # Fill x and y with data read from all files
    x = [] ; y = []
    for file, columns, lines_skip in zip( args['files'], args['columns'],
            args['skip_lines'] ):
        xfile, yfile = read_file( file, columns, lines_skip,
            scale = args['scale'], handle_complex = args['handle_complex'] )
        x.append(xfile), y.append(yfile)

    # Check if we want only to preprocess or only to plot
    if args['pre_process']:
        # Convert the preprocess option to enum
        option_preproc = Option_preprocess.nothing
        if( args['get_efield'] ):
            option_preproc = Option_preprocess.get_efield
        elif( args['get_eps'] ):
            option_preproc = Option_preprocess.get_eps
        elif( args['fourier_transform'] ):
            option_preproc = Option_preprocess.fourier_transform
        elif( args['add'] ):
            option_preproc = Option_preprocess.add
        elif( args['sub'] ):
            option_preproc = Option_preprocess.sub
        # Check if we want a diagonal component of the dielectric tensor
        diagonal_component = args['xx'] or args['yy'] or args['zz']
        # Preprocess
        preprocess( x, y, args['output'], option_preproc, \
            args['wcut'][0], diagonal_component )
    else:
        fig = matplotlib.pyplot.figure( 1, figsize=(8,5.5) )
        ax = fig.add_subplot( 111 )
        preformat_plot( )

        xmin = 1.e30 ; xmax = -1.e30
        ymin = 1.e30 ; ymax = -1.e30
        for xx, yy, cc in zip( x, y, args['captions'] ):
            plt.plot( xx, yy, label=cc )
            xmin = min(min(xx),xmin) ; xmax=max(max(xx),xmax)
            ymin = min(min(yy),ymin) ; ymax=max(max(yy),ymax)
        format_plot( ax, xmin, xmax, ymin, ymax, args['xlim'], args['ylim'],
            args['semilog'], args['xlabel'][0], args['ylabel'][0],
            args['hide_legend'], args['legend_position'])

        plt.savefig(f"{args['plot_name'][0]}.png", orientation='portrait',format='png',dpi=300)

        if args['show']:
            plt.show()

#-------------------------------------------------------------------------------

if __name__ == '__main__':
    main()
