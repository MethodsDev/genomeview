import collections
import numpy

from genomeview.track import Track
from genomeview.utilities import match_chrom_format

COLORS = ["blue", "red", "green", "black"]

class Series:
    def __init__(self, x, y, color=None, label=None):
        x, y = zip(*sorted(zip(x,y)))
        self.x = x
        self.y = y
        self.color = color
        self.label = label
        

class GraphTrack(Track):
    """
    Visualizes quantitative data as a line across coordinates within the current genomic view.

    One or more datasets can be visualized (with different colors) on the same track using the
    ``add_series()`` method.
    """
    def __init__(self, name=None, x=None, y=None):
        super().__init__(name)

        self.series = collections.OrderedDict()
        
        self.min_y = 0
        self.max_y = 0
        
        if x is not None:
            self.add_series(x, y)

        self.height = 100
        self.ymargin = 5

        self.fill_coverage = False
        
    def add_series(self, x, y, color=None, label=None):
        """
        Add a dataset corresponding to a single line in the track (ie, a "series"). Note that
        while a single GraphTrack can visualize multiple datasets, they are all plotted on 
        the same y-axis and so should share the same units.

        Arguments:
            x: a list of genomic coordinates
            y: a list of data values; each y-value must correspond to a single x-value
            color: an SVG color for the line being plotted
            label: an optional text label for the graph being plotted (currently unused)
        """
        if label is None:
            label = "series_{}".format(len(self.series))
            
        assert label not in self.series

        x = numpy.asarray(x)
        y = numpy.asarray(y)

        if color is None:
            color = COLORS[len(self.series) % len(COLORS)]
            
        self.series[label] = Series(x, y, color, label)

        self.min_y = min(self.min_y, y[numpy.isfinite(y)].min())
        self.max_y = max(self.max_y, y[numpy.isfinite(y)].max())

    def ytopixels(self, yval):
        height = self.max_y - self.min_y
        return self.height - ((yval - self.min_y) / height * (self.height-2*self.ymargin) + self.ymargin)
        
    def render(self, renderer):
        for label, series in self.series.items():
            if not self.fill_coverage:
                for i in range(len(series.x)-1):
                    if any(numpy.isnan(series.x[i:i+2])) or any(numpy.isnan(series.y[i:i+2])):
                        continue
                    x1 = self.scale.topixels(series.x[i])
                    x2 = self.scale.topixels(series.x[i+1])
                    y1 = self.ytopixels(series.y[i])
                    y2 = self.ytopixels(series.y[i+1])
                    
                    yield from renderer.line(x1, y1, x2, y1, 
                        **{"stroke-width":1, "stroke":series.color, "stroke-linecap":"square", "shape-rendering":"geometricPrecision"})
                    yield from renderer.line(x2, y1, x2, y2, 
                        **{"stroke-width":1, "stroke":series.color, "stroke-linecap":"square", "shape-rendering":"geometricPrecision"})
            else:
                current_min_y = 0  # keeping track of min because higher coverage is lower y value (closer to top)
                full_path = "<path d=\"M "
                x1 = self.scale.topixels(series.x[0]) + renderer.x
                y1 = self.ytopixels(0) + renderer.y
                full_path += str(x1) + " " + str(y1)

                for i in range(len(series.x)):
                    if numpy.isnan(series.x[i]) or numpy.isnan(series.y[i]):
                        continue
                    x1 = self.scale.topixels(series.x[i])+ renderer.x
                    y1 = self.ytopixels(series.y[i]) + renderer.y
                    current_min_y = min(current_min_y, y1)
                    full_path += " L " + str(x1) + " " + str(y1)

                y1 = self.ytopixels(0) + renderer.y
                full_path += " L " + str(x1) + " " + str(y1) + "\""
                full_path += " xcenter=\"" + str((self.ytopixels(0) + current_min_y)/2) +  "\" stroke=\"" + series.color + "\" fill=\"" + series.color + "\" stroke_width=\"1\"></path>"

                yield full_path

        # since the labels are drawn at the top of the ticks, let's make sure the top tick/label is 
        # more than 12 pixels from the top of the track so it doesn't get clipped
        # TODO: this ignores the margin, as of now
        axis_max_y = self.min_y + (self.max_y - self.min_y) * (1-7/self.height)

        # ticks = get_ticks(self.min_y, axis_max_y, 4)
        ticks = numpy.linspace(self.min_y, axis_max_y, 4)

        yield from renderer.line(1, self.ytopixels(ticks[0]), 1, self.ytopixels(ticks[-1]), 
                                 **{"stroke-width":2, "stroke":"gray", "stroke-linecap":"square", "shape-rendering":"geometricPrecision"})
        for tick in ticks:
            if self.max_y > 1_000:
                label = "{:.1g}".format(tick)
            elif self.max_y < 1:
                label = "{:.1f}".format(tick)
            else:
                label = "{:,.0f}".format(tick)

            y = self.ytopixels(tick)
            yield from renderer.line(1, y, 10, y, 
                                     **{"stroke-width":2, "stroke":"gray", "stroke-linecap":"square", "shape-rendering":"geometricPrecision"})
            yield from renderer.text(14, y, label, anchor="start", fill="gray")
            
        for x in self.render_label(renderer):
            yield x


class BigWigTrack(GraphTrack):
    """
    Visualizes continuous-valued data from a bigwig file. Requires pyBigWig
    module to be installed. Supports using web URLs as well as file paths.
    """
    def __init__(self, path, nbins=1000, name=None):
        super().__init__(name)
        
        import pyBigWig
        self.bigwig = pyBigWig.open(path)
        self.nbins = 1000

    def layout(self, scale):
        super().layout(scale)
        x = []
        y = []
        binsize = max(1, int((scale.end-scale.start) / self.nbins))
        
        chrom = match_chrom_format(scale.chrom, self.bigwig.chroms().keys())

        for i in range(scale.start, scale.end, binsize):
            value = self.bigwig.stats(chrom, i, i+binsize)[0]
            value = 0.0 if value==None else value
            x.append(i+binsize/2)
            y.append(value)
        
        self.series = {"vals":Series(x, y, color="black")}
        
        self.min_y = min(y)
        self.max_y = max(y)
