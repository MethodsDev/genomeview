import collections
import gzip
import math
import os
import pysam
import inspect
import ipywidgets as widgets
import re

from intervaltree import Interval, IntervalTree

import genomeview
from genomeview import utilities


def visualize_data(file_paths, chrom, start, end, reference_path=None, 
                   width=900, axis_on_top=False):
    """
    Creates a GenomeView document to display the data in the specified
    files (eg bam, bed, etc).

    Args:
        file_paths: this specifies the file paths to be rendered. It must be 
            either a list/tuple of the paths, or a dictionary mapping 
            {track_name:path}. (If you are using a python version prior to 3.6, 
            use collections.ordereddict to ensure the order remains the same.)
            Currently supports files ending in .bam, .cram, .bed, .bed.gz, 
            .bigbed, or .bigwig (or .bw). Most of these file types require a
            separate index file to be present (eg a .bam.bai or a .bed.gz.tbi 
            file must exist).
        chrom: chromosome (or contig) to be rendered
        start: start coordinate of region to be rendered
        end: end coordinate of region to be rendered
        reference_path: path to fasta file specifying reference genomic 
            sequence. This is required in order to display mismatches
            in bam tracks.
        width: the pixel width of the document
        axis_on_top: specifies whether the axis should be added at the bottom
            (default) or at the top
    """
    if reference_path is not None:
        source = genomeview.FastaGenomeSource(reference_path)
    else:
        source = None

    doc = genomeview.Document(width)
    
    view = genomeview.GenomeView(chrom, start, end, "+", source)
    doc.add_view(view)

    def add_axis():
        axis_track = genomeview.Axis("axis")
        view.add_track(axis_track)

    if axis_on_top:
        add_axis()

    if isinstance(file_paths, collections.abc.Mapping):
        names = file_paths.keys()
        file_paths = [file_paths[name] for name in names]
    else:
        names = [None] * len(file_paths)
        file_paths = file_paths
        
    for name, path in zip(names, file_paths):
        if path.lower().endswith(".bam") or path.lower().endswith(".cram"):
            if utilities.is_paired_end(path):
                cur_track = genomeview.PairedEndBAMTrack(path, name=name)
            else:
                cur_track = genomeview.SingleEndBAMTrack(path, name=name)
                if utilities.is_long_frag_dataset(path):
                    cur_track.min_indel_size = 5

        elif path.lower().endswith(".bed") or path.lower().endswith(".bed.gz") or path.lower().endswith(".bigbed") or path.lower().endswith(".bb"):
            cur_track = genomeview.BEDTrack(path, name=name)

        elif path.lower().endswith(".bigwig") or path.lower().endswith(".bw"):
            cur_track = genomeview.BigWigTrack(path, name=name)

        else:
            suffix =  os.path.basename(path)
            raise ValueError("Unknown file suffix: {}".format(suffix))

        view.add_track(cur_track)

    if not axis_on_top:
        add_axis()

    return doc



### newly added wrappers

def get_regions_by_read_id(bam_file, read_id):
    regions = []

    with pysam.AlignmentFile(bam_file, "rb") as bam_in:
        for read in bam_in.fetch():
            if read.query_name != read_id:
                continue

            regions.append(Interval(read.reference_start, read.reference_end, bam_in.get_reference_name(read.reference_id) + "+" if read.is_forward else "-"))

    return(regions)


def find_read_in_bam(read_id, bam, opener_fn=pysam.AlignmentFile):
    regions = []
    virtual_bams = []

    with opener_fn(bam) as bam_in:
        bam_refs = bam_in.references
        for read in bam_in.fetch():
            if read.query_name != read_id:
                continue
            regions.append(Interval(read.reference_start, read.reference_end, read.reference_name + ("+" if read.is_forward else "-")))
            virtual_bams.append(genomeview.VirtualBAM([read], bam_refs))
    
    if len(regions) == 0:
        print("Error: read either not found or not aligned")
        return -1

    return (regions, virtual_bams)


def get_virtualbam_max_coverage(coverage_bam):
    all_starts = []
    all_ends = []
    for inter in sorted(coverage_bam.reads_interval_tree):
        all_starts.append(inter.begin)
        all_ends.append(inter.end)
    all_ends.sort()
    current_max = current_coverage = i = j = 0
    while i < len(all_starts) or j < len(all_ends):
        if i < len(all_starts) and all_starts[i] < all_ends[j]:
            current_coverage += 1
            i += 1
            current_max = max(current_max, current_coverage)
        elif i < len(all_starts):
            current_coverage -= 1
            j += 1
        elif j < len(all_ends):
            current_coverage -= 1
            j += 1
    return current_max


class TighterSingleEndBAMTrack(genomeview.SingleEndBAMTrack):
    def __init__(self, *args, **kwdargs):
        super().__init__(*args, **kwdargs)
        self.row_height = 3
        self.margin_y = 2

def color_from_bed(interval):
    if interval.tx.color and len(interval.tx.color.split(",")) == 3:
        # hex_colors = interval.tx.color.split(",")
        # return "#{0:02x}{1:02x}{2:02x}".format(int(hex_colors[0]), int(hex_colors[1]), int(hex_colors[2]))
        return "rgb(" + interval.tx.color + ")"
    else:
        return genomeview.color_by_strand(interval)


# adding chrom and strand accessors for Interval from data slot based on the usage made in the code below
@property
def interval_chrom(self):
    return self.data[:-1]

Interval.chrom = interval_chrom

@property
def interval_strand(self):
    return self.data[-1:]

Interval.strand = interval_strand


def interval_data_reduce(current_data, new_data):
    if current_data == new_data:
        return current_data
    else:
        return None



class Configuration():
    def __init__(self, genome_fasta, bed_annotation, gtf_annotation = None, bed_color_fn=color_from_bed):
        self.source = genomeview.genomesource.FastaGenomeSource(genome_fasta)
        self.bed_annotation = bed_annotation
        if gtf_annotation:
            self.index_gtf(gtf_annotation)
        self.bed_color_fn = bed_color_fn

    def index_gtf(self, gtf_annotation):
        gene_id_regex = re.compile('gene_id "([a-zA-Z0-9\._]+)";')
        gene_name_regex = re.compile('gene_name "([a-zA-Z0-9\.]+)";')
        transcript_id_regex = re.compile('transcript_id "([a-zA-Z0-9\._\^\-]+=?)";')
        exon_id_regex = re.compile('exon_id "([a-zA-Z0-9\._]+)";')
        # transcript_name_regex = re.compile('transcript_name "([a-zA-Z0-9\.]+)";')
        # protein_id_regex = re.compile('protein_id "([a-zA-Z0-9\.]+)";')

        self.gene_name_to_gene_id = {}
        self.gene_to_transcripts = {}
        self.transcript_to_gene = {}
        self.gene_to_exons = {}
        self.transcript_to_exons = {}
        self.id_to_coordinates = {}

        with gzip.open(gtf_annotation, "r") as gtf_file:
            # current_gene_interval = None
            for entry in pysam.tabix_iterator(gtf_file, pysam.asGTF()):
                if entry.feature == "gene":
                    # gene_id = (entry.attributes.split(";")[0]).split(" ")[1].strip('"')
                    res = gene_id_regex.search(entry.attributes)
                    if res:
                        gene_id = res.group(1)
                    else:
                        print("missing gene_id in a gene entry, skipping entry")


                    #gene_name = (entry.attributes.split(";")[2]).split(" ")[2].strip('"')
                    res = gene_name_regex.search(entry.attributes)
                    if res:
                        gene_name = res.group(1)

                    self.gene_name_to_gene_id[gene_name] = gene_id
                    self.id_to_coordinates[gene_id] = Interval(entry.start, entry.end, entry.contig + entry.strand)
                    self.gene_to_transcripts[gene_id] = []
                    self.gene_to_exons[gene_id] = IntervalTree()

                elif entry.feature == "transcript":
                    # gene_id = (entry.attributes.split(";")[0]).split(" ")[1].strip('"')
                    gene_id = None
                    res = gene_id_regex.search(entry.attributes)
                    if res:
                        gene_id = res.group(1)

                        if gene_id not in self.id_to_coordinates:
                            self.id_to_coordinates[gene_id] = Interval(entry.start, entry.end, entry.contig + entry.strand)
                            self.gene_to_transcripts[gene_id] = []
                            self.gene_to_exons[gene_id] = IntervalTree()

                    # transcript_id = (entry.attributes.split(";")[1]).split(" ")[2].strip('"')
                    res = transcript_id_regex.search(entry.attributes)
                    if res:
                        transcript_id = res.group(1)
                    else:
                        print("missing transcript_id in a transcript entry, skipping entry:")
                        print(entry)

                    self.id_to_coordinates[transcript_id] = Interval(entry.start, entry.end, entry.contig + entry.strand)
                    self.transcript_to_gene[transcript_id] = gene_id
                    self.transcript_to_exons[transcript_id] = IntervalTree()
                    if gene_id:
                        self.gene_to_transcripts[gene_id].append(transcript_id)

                elif entry.feature == "exon":
                    # gene_id = None
                    res = gene_id_regex.search(entry.attributes)
                    if res:
                        gene_id = res.group(1)
                        self.gene_to_exons[gene_id].add(Interval(entry.start, entry.end, entry.contig + entry.strand))

                    transcript_id = None
                    res = transcript_id_regex.search(entry.attributes)
                    if res:
                        transcript_id = res.group(1)
                        self.transcript_to_exons[transcript_id].add(Interval(entry.start, entry.end, entry.contig + entry.strand))

                    exon_id = None
                    res = exon_id_regex.search(entry.attributes)
                    if res:
                        exon_id = res.group(1)
                        if exon_id not in self.id_to_coordinates:
                            self.id_to_coordinates[exon_id] = Interval(entry.start, entry.end, entry.contig + entry.strand)
                    #else:
                    #    print("missing exon_id in a exon entry, skipping entry")

                    # gene_id = (entry.attributes.split(";")[0]).split(" ")[1].strip('"')
                    # transcript_id = (entry.attributes.split(";")[1]).split(" ")[2].strip('"')
                    # exon_id = (entry.attributes.split(";")[7]).split(" ")[2].strip('"')

    def update_bed(self, bed_annotation):
        self.bed_annotation = bed_annotation

    def add_bed_tracks_to_view(self, view):
        if self.bed_annotation:
            if type(self.bed_annotation) is list:
                for bed_path in bed_annotation:
                    bed_track = genomeview.BEDTrack(bed_path)
                    bed_track.color_fn = self.bed_color_fn
                    view.add_track(bed_track)
            elif type(self.bed_annotation) is dict:
                for bed_name, bed_path in self.bed_annotation.items():
                    bed_track = genomeview.BEDTrack(bed_path, name=bed_name)
                    bed_track.color_fn = self.bed_color_fn
                    view.add_track(bed_track)
            else:
                bed_track = genomeview.BEDTrack(self.bed_annotation)
                bed_track.color_fn = self.bed_color_fn
                view.add_track(bed_track)


    def add_virtualbed_tracks_to_view(self, view, chrom, start, end, vertical_layout=True, use_names=True):
        if self.bed_annotation:
            if type(self.bed_annotation) is list:
                for bed_path in self.bed_annotation:
                    virtual_bed = genomeview.bedtrack.VirtualBEDTrack(bed_path)
                    virtual_bed.index(chrom, start, end, field_defs=None)
                    virtual_bed.color_fn = self.bed_color_fn
                    virtual_bed.vertical_layout = vertical_layout
                    view.add_track(virtual_bed)
            elif type(self.bed_annotation) is dict:
                for bed_name, bed_path in self.bed_annotation.items():
                    if use_names:
                        view.add_track(genomeview.track.TrackLabel(bed_name))
                    else:
                        view.add_track(genomeview.track.TrackLabel(""))
                    virtual_bed = genomeview.bedtrack.VirtualBEDTrack(bed_path, name="")
                    virtual_bed.index(chrom, start, end, field_defs=None)
                    virtual_bed.color_fn = self.bed_color_fn
                    virtual_bed.vertical_layout = vertical_layout
                    view.add_track(virtual_bed)
            else:
                virtual_bed = genomeview.bedtrack.VirtualBEDTrack(self.bed_annotation)
                virtual_bed.index(chrom, start, end, field_defs=None)
                virtual_bed.color_fn = self.bed_color_fn
                virtual_bed.vertical_layout = vertical_layout
                view.add_track(virtual_bed)


    def make_genomeview_row(self, start, end, chrom, strand, bams_list, 
                            padding_perc = 0.1, 
                            with_coverage = True,
                            include_secondary = False,
                            row = None, 
                            tighter_track = False):
        padding = math.ceil((end - start) * padding_perc)

        if row is None:
            row = genomeview.ViewRow("row")
        gene_view = genomeview.GenomeView(chrom, max(0, start - padding), end + padding, "+", self.source)
        # gene_view = genomeview.GenomeView(chrom, start - padding, end + padding, "+", self.source)
        gene_view.add_track(genomeview.track.TrackLabel(chrom + " : " + str(start - padding) + " - " + str(end + padding)))

        self.add_bed_tracks_to_view(gene_view)

        gene_view.add_track(genomeview.Axis())
        for key, value in bams_list.items():
            if isinstance(value, genomeview.VirtualBAM):
                opener_kwargs = {'opener_fn': lambda x: x}
            else:
                opener_kwargs = {}


            if with_coverage:
                coverage_track = genomeview.BAMCoverageTrack(value, name=key, **opener_kwargs)
                gene_view.add_track(coverage_track)
            if tighter_track:
                bam_track = TighterSingleEndBAMTrack(value, name=key, **opener_kwargs)
            else:
                bam_track = genomeview.SingleEndBAMTrack(value, name=key, **opener_kwargs)
            if include_secondary:
                coverage_track.include_secondary = True
                bam_track.include_secondary = True
            gene_view.add_track(bam_track)
        row.add_view(gene_view)
        return row


    def add_single_view_row_to_plot(self, doc, bams_list, 
                                    interval = None, 
                                    data = None, 
                                    padding_perc = 0.1,
                                    with_coverage = True,
                                    include_secondary = False,
                                    tighter_track = False):
        if data is not None:
            chrom = data[2]
            strand = data[3] if len(data) > 2 else True 
            start = data[1].begin
            end = data[1].end
        elif interval is not None:
            chrom = interval.chrom
            strand = interval.strand
            start = interval.begin
            end  = interval.end
        else:
            raise("Neither an Interval or data structure has been provided.")

        row = self.make_genomeview_row(start, end, chrom, strand, bams_list, 
                                       padding_perc = padding_perc, 
                                       with_coverage = with_coverage,
                                       include_secondary = include_secondary,
                                       row = None, 
                                       tighter_track = tighter_track)
        doc.elements.append(row)
        return doc


    def add_multi_view_row_to_plot(self, doc, bams_list, 
                                   interval_list = None, 
                                   data_list = None,
                                   padding_perc = 0.1,
                                   with_coverage = True,
                                   include_secondary = False,
                                   tighter_track = False):
        row = genomeview.ViewRow("row")

        if interval_list is not None:
            for interval in interval_list:
                chrom = interval.chrom
                strand = interval.strand
                start = interval.begin
                end  = interval.end
                
                self.make_genomeview_row(start, end, chrom, strand, bams_list, 
                                         padding_perc = padding_perc, 
                                         with_coverage = with_coverage,
                                         include_secondary = include_secondary,
                                         row = row, 
                                         tighter_track = tighter_track)
        elif data_list is not None:
            for data in data_list:
                chrom = data[2]
                strand = data[3] if len(data) > 2 else True 
                start = data[1].begin
                end = data[1].end
                
                self.make_genomeview_row(start, end, chrom, strand, bams_list, 
                                         padding_perc = padding_perc,
                                         with_coverage = with_coverage,
                                         include_secondary = include_secondary,
                                         row = row, 
                                         tighter_track = tighter_track)
        doc.elements.append(row)
        return doc


    def plot_interval(self, bams_list, 
                      interval = None, 
                      data = None, 
                      padding_perc = 0.1,
                      with_coverage = True,
                      include_secondary = False,
                      view_width = 1600, 
                      tighter_track = False):
        doc = genomeview.Document(view_width)
        return self.add_single_view_row_to_plot(doc, bams_list, 
                                                interval = interval, 
                                                data = data, 
                                                padding_perc = padding_perc,
                                                with_coverage = with_coverage,
                                                include_secondary = include_secondary,
                                                tighter_track = tighter_track)


    def plot_intervals(self, bams_list, 
                       interval_list = None, 
                       data_list = None, 
                       N_per_row = 1, 
                       padding_perc = 0.1,
                       with_coverage = True,
                       include_secondary = False,
                       view_width = 1600,
                       tighter_track = False):
        doc = genomeview.Document(view_width)

        if interval_list is not None:
            for i in range(0, len(interval_list), N_per_row):
                doc = self.add_multi_view_row_to_plot(doc, bams_list, 
                                                      interval_list = interval_list[i:i+N_per_row], 
                                                      data_list = None, 
                                                      padding_perc = padding_perc,
                                                      with_coverage = with_coverage,
                                                      include_secondary = include_secondary,
                                                      tighter_track = tighter_track)
        elif data_list is not None:
            for i in range(0, len(data_list), N_per_row):  
                doc = self.add_multi_view_row_to_plot(doc, bams_list, 
                                                      interval_list = None, 
                                                      data_list = data_list[i:i+N_per_row], 
                                                      padding_perc = padding_perc, 
                                                      with_coverage = with_coverage,
                                                      include_secondary = include_secondary,
                                                      tighter_track = tighter_track)
        return doc

    def get_feature_info(self, feature):
        feature_id = None
        feature_type = None
        if feature in self.gene_name_to_gene_id:
            feature_id = self.gene_name_to_gene_id[feature]
            feature_type = "gene"
        elif feature in self.gene_to_exons:
            feature_id = feature
            feature_type = "gene"
        elif feature in self.transcript_to_exons:
            feature_id = feature
            feature_type = "transcript"
        elif feature in self.id_to_coordinates:
            feature_id = feature
            feature_type = "exon"
        
        return(feature_id, feature_type)



    def plot_feature(self, feature, bams_list, 
                     padding_perc = 0.1,
                     with_coverage = True,
                     include_secondary = False,
                     view_width = 1600, 
                     tighter_track = False):

        (feature_id, feature_type) = self.get_feature_info(feature)

        return self.plot_interval(bams_list, 
                                  interval = self.id_to_coordinates[feature_id], 
                                  padding_perc = padding_perc,
                                  with_coverage = with_coverage,
                                  include_secondary = include_secondary,
                                  view_width = view_width, 
                                  tighter_track = tighter_track)


    # plot_feature for a list in tabs
    def plot_features(self, features, bams_list, output_format="svg",
                      **kwargs):

        params = inspect.signature(self.plot_feature).parameters
        filtered_kwargs = {k: v for k, v in kwargs.items() if k in params}

        features_tab = widgets.Tab()
        tab_contents = []

        for feature in features:
            tab_contents.append(self.plot_feature(feature, bams_list, **filtered_kwargs).get_widget(output_format))

        features_tab.children = tab_contents
        features_tab.titles = features

        return features_tab



    # double check no risk of conflict with passing an opener_fn through kwargs to find_read_in_bam that would also affect some other things down the line with VirtualBAM
    def plot_read(self, read_id, bam, output_format="svg", **kwargs):

        (regions, virtual_bams) = find_read_in_bam(read_id, bam, **kwargs)

        # Get the parameter names of plot_interval
        params = inspect.signature(self.plot_interval).parameters

        # Filter kwargs to only include keys that are in plot_interval's parameters
        filtered_kwargs = {k: v for k, v in kwargs.items() if k in params}

        all_widgets = []
        for region, virtual_bam in zip(regions, virtual_bams):
            all_widgets.append(self.plot_interval(bams_list={"read": virtual_bam}, interval=region).get_widget(output_format), **filtered_kwargs)

        return widgets.VBox(all_widgets)



    def plot_exons(self, feature, bams_list, merge_exons = True,
                   padding_perc = 0.05, 
                   with_coverage = True,
                   include_secondary = False, 
                   view_width = 1600,
                   normalize_interval_width = False,
                   N_per_row = 99999,
                   as_widget = False,
                   tighter_track = False):

        (feature_id, feature_type) = self.get_feature_info(feature)
        
        if feature_type == "exon":
            return self.plot_feature(feature_id, bams_list=bams_list, padding_perc=padding_perc, with_coverage=with_coverage, include_secondary=include_secondary, view_width=view_width, tighter_track=tighter_track)

        elif feature_type == "transcript":
            exons_list = sorted(self.transcript_to_exons[feature_id])
        elif feature_type == "gene":
            if merge_exons:
                tmp_exons = self.gene_to_exons[feature_id].copy()
                tmp_exons.merge_overlaps(data_reducer = interval_data_reduce)
                exons_list = sorted(tmp_exons)
            else:
                exons_list = sorted(self.gene_to_exons[feature_id])

        if as_widget:
            all_views = []
            all_titles = []
            for exon in exons_list:
                doc = self.make_intervals_row_from_virtual(genomeview.Document(view_width), [exon], bams_list=bams_list, padding_perc=padding_perc, with_coverage=with_coverage, include_secondary=include_secondary, normalize_interval_width=normalize_interval_width, tighter_track=tighter_track)
                all_views.append(widgets.HTML(doc._repr_svg_()))
                # all_views.append(doc.get_widget())
                all_titles.append("Exon:: " + exon.data + " : " + str(exon.begin) + " - " + str(exon.end))

            stack = widgets.Stack(all_views, selected_index=0)
            dropdown = widgets.Dropdown(options=all_titles)
            widgets.jslink((dropdown, 'index'), (stack, 'selected_index'))
            return(widgets.VBox([dropdown, stack]))

        else:
            doc = genomeview.Document(view_width)
            for i in range(0, len(exons_list), N_per_row):
               doc = self.make_intervals_row_from_virtual(doc, exons_list[i:i+N_per_row], bams_list=bams_list, padding_perc=padding_perc, with_coverage=with_coverage, include_secondary=include_secondary, normalize_interval_width=normalize_interval_width, tighter_track=tighter_track)
            return doc



    def plot_splice_junctions(self, feature, bams_list,
                              padding_perc = 0.05, 
                              with_coverage = True,
                              include_secondary = False, 
                              view_width = 1600,
                              normalize_interval_width = False,
                              as_widget = False,
                              tighter_track = False):

        (feature_id, feature_type) = self.get_feature_info(feature)
        
        if feature_type == "exon":
            print("Error, feature type provided is an exon, hence there is no splice junction within it.")
            return

        exons_pairs = []
        if feature_type == "transcript":
            all_exons = sorted(self.transcript_to_exons[feature_id])
            for i in range(1, len(all_exons)):
                exons_pairs.append((all_exons[i-1], all_exons[i]))

        elif feature_type == "gene":
            for transcript_id in self.gene_to_transcripts[feature_id]:
                all_exons = sorted(self.transcript_to_exons[transcript_id])
                for i in range(1, len(all_exons)):
                    exons_pairs.append((all_exons[i-1], all_exons[i]))

        if len(exons_pairs) == 0:
            print("No splice junctions for the requested feature, this probably means it only has 1 exon")
            return

        seen = set()
        seen_add = seen.add
        exons_pairs = [x for x in exons_pairs if not (x in seen or seen_add(x))]

        if as_widget:
            all_views = []
            all_titles = []
            for pair in exons_pairs:
                doc = self.make_intervals_row_from_virtual(genomeview.Document(view_width), pair, bams_list=bams_list, padding_perc=padding_perc, with_coverage=with_coverage, include_secondary=include_secondary, normalize_interval_width=normalize_interval_width, tighter_track=tighter_track)
                all_views.append(widgets.HTML(doc._repr_svg_()))
                # all_views.append(doc.get_widget())
                all_titles.append("Splice junction btw: exon:: " + pair[0].data + ":" + str(pair[0].begin) + "-" + str(pair[0].end) + " and exon::" + pair[1].data + " : " + str(pair[1].begin) + " - " + str(pair[1].end))

            stack = widgets.Stack(all_views, selected_index=0)
            dropdown = widgets.Dropdown(options=all_titles)
            widgets.jslink((dropdown, 'index'), (stack, 'selected_index'))
            return(widgets.VBox([dropdown, stack]))

        else:
            doc = genomeview.Document(view_width)
            for pair in exons_pairs:
               doc = self.make_intervals_row_from_virtual(doc, pair, bams_list=bams_list, padding_perc=padding_perc, with_coverage=with_coverage, include_secondary=include_secondary, normalize_interval_width=normalize_interval_width, tighter_track=tighter_track)
            return doc
     

    def make_intervals_row_from_virtual(self, doc, intervals_list, bams_list,
                                        padding_perc = 0.1, 
                                        with_coverage = True,
                                        include_secondary = False,
                                        row = None, 
                                        normalize_interval_width = False,
                                        tighter_track = False):


        if row is None:
            row = genomeview.ViewRow("row")

        total_interval_size = 0
        left_bound = math.inf
        right_bound = -math.inf

        smallest_interval_size = math.inf
        for interval in intervals_list:
            total_interval_size += interval.end - interval.begin
            smallest_interval_size = min(smallest_interval_size, interval.end - interval.begin)
            left_bound = min(left_bound, interval.begin)
            right_bound = max(right_bound, interval.end)

        reserved_width = (len(intervals_list) - 1) * row.space_between + doc.margin_x * 2
        if not normalize_interval_width:
            padding = math.ceil(smallest_interval_size * padding_perc)
            total_interval_size = total_interval_size + (padding * len(intervals_list))
            per_base_size = (doc.width - reserved_width)/total_interval_size

        max_coverage_dict = {}
        virtual_bams_dict = {}
        for key, value in bams_list.items():

            virtual_bams_dict[key] = []
            all_reads_for_coverage = set()

            bam_refs = None
            with pysam.AlignmentFile(value, "rb") as bam:
                bam_refs = bam.references
                for interval in intervals_list:
                    for read in bam.fetch(interval.chrom, interval.begin - padding, interval.end + padding):
                        if not include_secondary and read.is_secondary:
                            continue
                        all_reads_for_coverage.add(read)

                coverage_bam = genomeview.bamtrack.VirtualBAM(all_reads_for_coverage, bam_refs)
                coverage_bam.index()

                for read in coverage_bam.fetch():
                    bam_track = genomeview.bamtrack.VirtualBAM([read], bam_refs)
                    bam_track.quick_consensus = False
                    virtual_bams_dict[key].append(bam_track)
                
                max_coverage_dict[key] = get_virtualbam_max_coverage(coverage_bam)
        
        first_interval = True

        for interval in intervals_list:
            start = interval.begin
            end = interval.end
            chrom = interval.chrom
            strand = interval.strand
            
            if normalize_interval_width:
                padding = math.ceil((end - start) * padding_perc)
                interval_width = (doc.width - reserved_width)/len(intervals_list)
            else:
                interval_width = math.floor((interval.end - interval.begin + padding) * per_base_size)
            interval_view = genomeview.GenomeView(chrom, start - padding, end + padding, "+", source=self.source)
            interval_view.add_track(genomeview.track.TrackLabel(chrom + (" +" if strand else " -") + " : " + str(start - padding) + " - " + str(end + padding)))
            
            # BED type features
            self.add_virtualbed_tracks_to_view(interval_view, chrom, left_bound, right_bound, use_names=first_interval)

            interval_view.add_track(genomeview.Axis())
            for key, value in bams_list.items():
                if first_interval:
                    interval_view.add_track(genomeview.track.TrackLabel(key))
                else:
                    interval_view.add_track(genomeview.track.TrackLabel(""))  # for spacing
                if with_coverage:
                    coverage_track = genomeview.BAMCoverageTrack(value, name="")
                    coverage_track.max_y = max_coverage_dict[key]
                    interval_view.add_track(coverage_track)
                for virtual_bam in virtual_bams_dict[key]:
                    if tighter_track:
                        bam_track = TighterSingleEndBAMTrack(virtual_bam, name=None, opener_fn=lambda x: x)
                    else:
                        bam_track = genomeview.bamtrack.SingleEndBAMTrack(virtual_bam, name=None, opener_fn=lambda x: x)
                    if include_secondary:
                        coverage_track.include_secondary = True
                        bam_track.include_secondary = True
                    interval_view.add_track(bam_track)
            interval_view.pixel_width = interval_width
            interval_view.margin_y = 0
            row.add_view(interval_view)

            first_interval = False

        doc.elements.append(row)
        return doc



