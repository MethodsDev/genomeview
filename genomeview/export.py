import logging
import os
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET
import re
import math
import copy
from intervaltree import Interval, IntervalTree



def save(doc, outpath, outformat=None):
    """
    Saves document `doc` to a file at `outpath`. By default, this file 
    will be in SVG format; if it ends with .pdf or .png, or if outformat
    is specified, the document will be converted to PDF or PNG if possible.

    Conversion to PDF and PNG require rsvg-convert (provided by librsvg), 
    inkscape or webkitToPDF (PDF conversion only).

    Attributes:
        doc: the :py:class:`genomeview.Document` to be saved
        outpath: a string specifying the file to save to; file extensions of
            .pdf or .png will change the default output format
        outformat: override the file format; must be one of "pdf", "png", or
            (the default) "svg"
    """
    if isinstance(outpath, bytes):
        outpath = outpath.decode()

    if outformat is None:
        if outpath.lower().endswith(".pdf"):
            outformat = "pdf"
        elif outpath.lower().endswith(".png"):
            outformat = "png"
        else:
            outformat = "svg"

    if outformat == "svg":
        with open(outpath, "w") as outf:
            render_to_file(doc, outf)
    else:
        # render to a temporary file then convert to PDF or PNG
        with tempfile.TemporaryDirectory() as outdir:
            temp_svg_path = os.path.join(outdir, "temp.svg")
            with open(temp_svg_path, "w") as outf:
                render_to_file(doc, outf)

            tree = ET.parse(temp_svg_path)
            root_svg = tree.getroot()
            svg_splitter = svg_splitter(root_svg)
            svg_splitter.split_svg(root_svg, max_height = 5000)
            svg_slices = svg_splitter.write_splits(prefix=os.path.join(outdir, "temp_"))

            for svg_slice in svg_slices:
                convert_svg(svg_slice, outpath, outformat)



def render_to_file(doc, outf):
    """
    Renders the document as an svg to a file-like object.
    """

    for l in doc.render():
        outf.write(l + "\n")




#############################################################################
######################### low-level functionality ###########################
#############################################################################

class SvgSplitter:
    # section_regex = re.compile(".*clipPath$")
    # read_regex = re.compile(".*path$")
    # text_regex = re.compile(".*text$")
    # path_dot = re.compile("L [0-9\.]+ ([0-9\.]+)")
    path_y_regex = re.compile(r'(?:M|L) (?:[-0-9\.]+) ([-0-9\.]+)')
    rect_y_regex = re.compile(r'y="([-0-9\.]+)" .* height="([-0-9\.]+)')
    line_y_regex = re.compile(r'y[12]="([-0-9\.]+)')
    text_y_regex = re.compile(r'y="([-0-9\.]+)')
    
    def __init__(self, root):
        ET.register_namespace("", "http://www.w3.org/2000/svg")  # Register the SVG namespace the avoid "ns0:" or ":ns0" being inserted in tags
        self.root = root
        self.split_template = self.make_template_svg_split(root)
        self.split_svgs = []

    def make_template_svg_split(self, root):
        return ET.Element('svg', {
                'version': root.attrib['version'],
                'baseProfile': root.attrib['baseProfile'],
                'width': root.attrib['width'],
                'height': "0"
                # 'height': str(min(max_height, original_height - i * max_height)) #,
                # 'xmlns': root.attrib['xmlns']
            })
    
    def get_new_svg_split(self):
        return copy.deepcopy(self.split_template)

    # def update_splits_height(self, element, current_height = 0):
    #     for subelement in element:
    #         current_height += update_splits_height(subelement, current_height)
    #     current_height += element.attrib['height']

    def update_splits_height(self, element, current_height=0):
        for subelement in element:
            current_height += self.update_splits_height(subelement, current_height)
        current_height += float(element.attrib.get('height', 0))
        return current_height
    
    def write_splits(self, prefix):
        output_files = []
        for i, split_svg in enumerate(self.split_svgs):
            split_svg.write(f'{prefix}_slice_{i + 1}.svg')
            output_files.append(f'{prefix}_slice_{i + 1}.svg')
        return(output_files)

    def is_read_block(self, element):
        if self.has_children(element):
            return False
        path_found = False
        for subelement in element:
            if self.is_path(subelement) and len(subelement.attrib['d'].split(" ")) == 16:
                return True
        return False

    def is_clipPath(self, element):
        return True if element.tag[-8:] == "clipPath" else False

    def is_path(self, element):
        return True if element.tag[-4:] == "path" else False

    def is_rect(self, element):
        return True if element.tag[-4:] == "rect" else False

    def is_text(self, element):
        return True if element.tag[-4:] == "text" else False

    def is_line(self, element):
        return True if element.tag[-4:] == "line" else False
    
    def is_defs(self, element):
        return True if element.tag[-4:] == "defs" else False

    def is_g(self, element):
        return True if element.tag[-1] == "g" else False
 
    def has_children(self, element):
        for el in element.iter():
            if el == element:
                continue
            if self.is_clipPath(el):  # or self.is_g(el):
                return True
        return False
   
    def get_path_max_y(self, element):
        max_y = 0

        matches = self.path_y_regex.finditer(element.attrib['d'])
        for match in matches:
            if float(match.group(1)) > max_y:
                max_y = float(match.group(1))

        return max_y

    def get_path_y_range(self, element):
        min_y = math.inf
        max_y = -math.inf

        matches = self.path_y_regex.finditer(element.attrib['d'])
        for match in matches:
            if float(match.group(1)) > max_y:
                max_y = float(match.group(1))
            if float(match.group(1)) < min_y:
                min_y = float(match.group(1))

        return (min_y, max_y)  # could return an Interval(min_y, max_y)
        
    # def get_line_y_range(self, element):
    #     min_y = -math.inf
    #     max_y = math.inf

    #     matches = line_y_regex.finditer(element.tag)
    #     for match in matches:
    #         if float(match.group(1)) > max_y:
    #             max_y = float(match.group(1))
    #         if float(match.group(1)) < min_y:
    #             min_y = float(match.group(1))
        
    #     return (min_y, max_y)  # could return an Interval(min_y, max_y)

    def get_line_y_range(self, element):
        min_y = min(float(element.attrib['y1']), float(element.attrib['y2']))
        max_y = max(float(element.attrib['y1']), float(element.attrib['y2']))
        return (min_y, max_y)  # could return an Interval(min_y, max_y)

    def get_line_y(self, element):  # assumes only horizontal lines
        return float(element.attrib['y1'])    
    
    # def get_rect_y_range(self, element):
    #     match = rect_regex.search(element.tag)
    #     return (float(matfch.group(1)), float(match.group(1)) + float(match.group(2))) # could return an Interval(min_y, max_y)

    def get_rect_y_range(self, element):
        return (float(element.attrib['y']), float(element.attrib['y']) + float(element.attrib['height']))
    
    def get_text_y_range(self, element):
        return (float(element.attrib['y']) - float(element.attrib['font-size']), float(element.attrib['y']))

    # can probably just keep the min/max y based on is_path elements alone
    def get_element_min_y(self, element):
        if self.is_path(element):
            return self.get_path_y_range(element)[0]
        elif self.is_rect(element):
            return self.get_rect_y_range(element)[0]
        elif self.is_text(element):
            return self.get_text_y_range(element)[0]
        # elif self.is_line(element):
        #     return self.get_line_y(element)
        # elif self.is_defs(element):
        #    return None

    def get_element_max_y(self, element):
        if self.is_path(element):
            return self.get_path_y_range(element)[1]
        elif self.is_rect(element):
            return self.get_rect_y_range(element)[1]
        elif self.is_text(element):
            return self.get_text_y_range(element)[1]
        # elif self.is_line(element):
        #     return self.get_line_y(element)

    def find_or_create_binterval(self, y_range, bins_tree, max_height):
        # Search for overlapping bins
        overlapping_bins = bins_tree.at(y_range[0])
        if len(overlapping_bins) == 1:
            return list(overlapping_bins)[0]

        if len(bins_tree) >= 1:
            current_binterval = sorted(bins_tree)[0]
            if y_range[0] < current_binterval.begin:
                while True:
                    previous_binterval = bins_tree.at(current_binterval.begin - 1)
                    if not previous_binterval:
                        new_bin = {'elements': [], 'min_y': current_binterval.begin - max_height, 'max_y': current_binterval.begin - 0.01}
                        new_binterval = Interval(current_binterval.begin - max_height, current_binterval.begin - 0.01, new_bin)
                        bins_tree.add(new_binterval)
                        current_binterval = new_binterval
                    else:
                        current_binterval = list(previous_binterval)[0]
                    if current_binterval.begin < y_range[0] < current_binterval.end:
                        return current_binterval
            else:
                while True:
                    next_binterval = bins_tree.at(current_binterval.end + 1)
                    if not next_binterval:
                        new_bin = {'elements': [], 'min_y': current_binterval.end + 0.01, 'max_y': current_binterval.end + max_height}
                        new_binterval = Interval(current_binterval.end + 0.01, current_binterval.end + max_height, new_bin)
                        bins_tree.add(new_binterval)
                        current_binterval = new_binterval
                    else:
                        current_binterval = list(next_binterval)[0]
                    if current_binterval.begin <= y_range[0] < current_binterval.end:
                        return current_binterval
    
        else:
            new_bin = {'elements': [], 'min_y': y_range[0], 'max_y': y_range[0] + max_height}
            # new_interval = Interval(y_range[0], y_range[1] + max_height, new_bin)
            new_binterval = Interval(y_range[0], y_range[0] + max_height - 0.01, new_bin)
            bins_tree.add(new_binterval)
            return new_binterval

    def create_bins(self, elements, max_height):
        bins_tree = IntervalTree()
        defs_elements = []
        current_bin = None
        initial_intron_line = None

        for element in elements:
            if self.is_defs(element):
                defs_elements.append(element)
                current_bin = None  # Reset current bin for elements following 'defs'
                continue

            # For 'path' elements
            if self.is_path(element):
                y_range = self.get_path_y_range(element)

                current_binterval = self.find_or_create_binterval(y_range, bins_tree, max_height)
                current_bin = current_binterval.data

                if initial_intron_line:
                    current_bin['elements'].append(initial_intron_line)
                    initial_intron_line = None

            # Add element to the current bin
            if current_bin is not None:
                current_bin['elements'].append(element)
            else:
                # Handle 'text' elements following 'defs' elements
                if self.is_line(element):
                    initial_intron_line = element
                elif self.is_text(element):
                    y_range = self.get_text_y_range(element)
                    bin_for_text = self.find_or_create_binterval(y_range, bins_tree, max_height).data
                    bin_for_text['elements'].append(element)

        # Convert IntervalTree to list of bins
        bins = [interval.data for interval in sorted(bins_tree)]

        # Append 'defs' elements to each bin after creation
        for bin in bins:
            bin['elements'] = defs_elements + bin['elements']

        return bins
    
    def split_svg(self, current_element, current_split = None, max_height = 32000, current_height = 0, root = True):
        last_clip_path_height = 0
        
        if not current_split:
            current_split = self.get_new_svg_split()
        
        for child_element in current_element:
            if child_element.tag.endswith('clipPath'):
                last_clip_path_height = float(child_element.find('.//{http://www.w3.org/2000/svg}rect').attrib['height'])
                continue
            elif self.is_g(child_element):
                child_element_height = last_clip_path_height
            else:
                child_element_height = float(child_element.attrib['height'])

            if float(current_split.attrib['height']) + child_element_height < max_height:
                # keep as whole with wtv is before, and can keep going
                current_split.append(child_element)
                current_split.attrib['height'] = str(float(current_split.attrib['height']) + child_element_height + 10) # +10 because there is a 10 between each group at the same level
            elif child_element_height < max_height:
                # finish previous part and start new one with this
                self.split_svgs.append(ET.ElementTree(current_split))
    
                current_split = self.get_new_svg_split()
                current_split.append(child_element)
                current_split.attrib['height'] = str(float(current_split.attrib['height']) + child_element_height + 10) # +10 because there is a 10 between each group at the same level
            else:  # needs to be split somehow, adding however much can be added first, then make new splits with the rest
                # add this level branch
                if self.has_children(child_element):
                    # go deeper
                    self.split_svg(child_element, current_split, max_height, root = False)
                    # update the height of this, maybe do at the end of the whole recursion instead
                # below here, can probable close the split after this element is processed before going further to make subsections easier to handle
                elif self.is_read_block(child_element):
                    self.split_svgs.append(ET.ElementTree(current_split))
                    
                    elements = list(child_element)
                    bins = self.create_bins(elements, max_height)
                    # print("generated " + str(len(bins)) + " bins")
                    for bin in bins:
                        new_svg = self.get_new_svg_split()
                        new_svg.attrib['height'] = str(bin['max_y'] - bin['min_y'])
                        if bin == bins[0]:
                            y_offset = math.inf
                            y_trim = -math.inf
                            for elem in bin['elements']:
                                current_min_y = self.get_element_min_y(elem)
                                if current_min_y and current_min_y < y_offset:
                                    y_offset = current_min_y
                                current_max_y = self.get_element_max_y(elem)
                                if current_max_y and current_max_y > y_trim:
                                    y_trim = current_max_y
                            # print("first bin offset = ", str(y_offset))
                            # print("first bin trim at = ", str(y_trim))
                            new_svg.attrib['height'] = str(y_trim - y_offset)
                        else:
                            y_offset = bin['min_y']
                            # print("next bin offset = ", str(y_offset))
                        if bin == bins[-1]:
                            y_trim = -math.inf
                            for elem in bin['elements']:
                                current_y = self.get_element_max_y(elem)
                                if current_y and current_y > y_trim:
                                    y_trim = current_y
                            # print("last bin trim at = ", str(y_trim))
                            new_svg.attrib['height'] = str(y_trim - y_offset + 10)
                        for elem in bin['elements']:
                            # need to subtract min_y from all height/y coordinates
                            self.offset_y(elem, y_offset)
                            new_svg.append(elem)
                        self.split_svgs.append(ET.ElementTree(new_svg))
                    current_split = None
        if root and current_split:
            self.split_svgs.append(ET.ElementTree(current_split))

    def offset_y(self, element, offset):
        if self.is_path(element):
            self.offset_path_y(element, offset)
        elif self.is_line(element):
            self.offset_line_y(element, offset)
        elif self.is_rect(element):
            self.offset_rect_y(element, offset)
        elif self.is_text(element):
            self.offset_text_y(element, offset)
        # self.is_defs(element) doesn't need updating

    def offset_path_y(self, element, offset):
        groups = element.attrib['d'].split(" ")
        for i in range(2, len(groups), 3):
            if groups[i].replace('.', '', 1).isdigit():  # Check if it's a number
                groups[i] = str(float(groups[i]) - offset)
        element.attrib['d'] = ' '.join(groups)
        
    def offset_line_y(self, element, offset):  # assumes only horizontal lines
        element.attrib['y1'] = str(float(element.attrib['y1']) - offset)
        element.attrib['y2'] = str(float(element.attrib['y2']) - offset)

    def offset_rect_y(self, element, offset):
        element.attrib['y'] = str(float(element.attrib['y']) - offset)
    
    def offset_text_y(self, element, offset):
        element.attrib['y'] = str(float(element.attrib['y']) - offset)



def convert_svg(inpath, outpath, outformat):
    converter = _getExportConverter(outformat)

    if converter == "webkittopdf":
        exportData = _convertSVG_webkitToPDF(inpath, outpath, outformat)
    elif converter == "librsvg":
        exportData = _convertSVG_rsvg_convert(inpath, outpath, outformat)
    elif converter == "inkscape":
        exportData = _convertSVG_inkscape(inpath, outpath, outformat)

    return exportData


def _getExportConverter(exportFormat, requested_converter=None):
    if requested_converter == "webkittopdf" and exportFormat=="png":
        logging.error("webkitToPDF does not support export to PNG; use librsvg or inkscape instead, or "
            "export to PDF")
        sys.exit(1)

    if exportFormat == "png" and requested_converter is None:
        return "librsvg"

    if requested_converter == "rsvg-convert":
        return "librsvg"

    if requested_converter in [None, "webkittopdf"]:
        if _checkWebkitToPDF():
            return "webkittopdf"

    if requested_converter in [None, "librsvg"]:
        if _checkRSVGConvert():
            return "librsvg"

    if requested_converter in [None, "inkscape"]:
        if _checkInkscape():
            return "inkscape"

    raise Exception("No converter found for conversion to {}".format(exportFormat))
    return None



def _checkWebkitToPDF():
    try:
        subprocess.check_call("webkitToPDF", stderr=subprocess.PIPE, shell=True)
        return True
    except subprocess.CalledProcessError:
        return False

def _checkRSVGConvert():
    try:
        subprocess.check_call("rsvg-convert -v", stdout=subprocess.PIPE, shell=True)
        return True
    except subprocess.CalledProcessError:
        return False

def _checkInkscape():
    try:
        subprocess.check_call("inkscape --version", stdout=subprocess.PIPE, shell=True)
        return True
    except subprocess.CalledProcessError:
        return False



def _convertSVG_webkitToPDF(inpath, outpath, outformat):
    if outformat.lower() != "pdf":
        return None

    try:
        cmd = "webkitToPDF {} {}".format(inpath, outpath)
        subprocess.check_call(cmd, shell=True)#, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError:
        return None

    return open(outpath, "rb").read()

def _convertSVG_inkscape(inpath, outpath, outformat):
    options = ""
    outformat = outformat.lower()
    if outformat == "png":
        options = "--export-dpi 150 --export-background white"

    try:
        subprocess.check_call("inkscape {} {} --export-{}={}".format(options, inpath, outformat, outpath), 
            shell=True)
    except subprocess.CalledProcessError as e:
        print("EXPORT ERROR:", str(e))

    return open(outpath, "rb").read()


def _convertSVG_rsvg_convert(inpath, outpath, outformat):
    options = ""
    outformat = outformat.lower()
    if outformat == "png":
        options = "-a --background-color white"

    try:
        subprocess.check_call("rsvg-convert -f {} {} -o {} {}".format(outformat, options, outpath, inpath), shell=True)
    except subprocess.CalledProcessError as e:
        print("EXPORT ERROR:", str(e))

    return open(outpath, "rb").read()