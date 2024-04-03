import collections
import gzip
import math
import os
import pysam
import inspect
import ipywidgets as widgets

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

            regions.append(Interval(read.reference_start, read.reference_end, bam_in.get_reference_name(read.reference_id)))

    return(regions)


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


class Configuration():
    def __init__(self, genome_path, annotation_path, gtf_path = None, tss_path = None):
        self.source = genomeview.genomesource.FastaGenomeSource(genome_path)
        self.annotation_path = annotation_path
        self.tss_path = tss_path
        if gtf_path:
            self.index_gtf(gtf_path)

    def index_gtf(self, gtf_path):
        self.gene_name_to_gene_id = {}
        self.gene_to_transcripts = {}
        self.transcript_to_gene = {}
        self.gene_to_exons = {}
        self.transcript_to_exons = {}
        self.id_to_coordinates = {}

        with gzip.open(gtf_path, "r") as gtf_file:
            current_gene_interval = None
            for entry in pysam.tabix_iterator(gtf_file, pysam.asGTF()):
                if entry.feature == "gene":
                    gene_id = (entry.attributes.split(";")[0]).split(" ")[1].strip('"')
                    gene_name = (entry.attributes.split(";")[2]).split(" ")[2].strip('"')

                    self.gene_name_to_gene_id[gene_name] = gene_id
                    self.id_to_coordinates[gene_id] = Interval(entry.start, entry.end, entry.contig)
                    self.gene_to_transcripts[gene_id] = []
                    self.gene_to_exons[gene_id] = IntervalTree()

                elif entry.feature == "transcript":
                    gene_id = (entry.attributes.split(";")[0]).split(" ")[1].strip('"')
                    transcript_id = (entry.attributes.split(";")[1]).split(" ")[2].strip('"')
                    self.id_to_coordinates[transcript_id] = Interval(entry.start, entry.end, entry.contig)
                    self.gene_to_transcripts[gene_id].append(transcript_id)
                    self.transcript_to_gene[transcript_id] = gene_id
                    self.transcript_to_exons[transcript_id] = IntervalTree()

                elif entry.feature == "exon":
                    gene_id = (entry.attributes.split(";")[0]).split(" ")[1].strip('"')
                    transcript_id = (entry.attributes.split(";")[1]).split(" ")[2].strip('"')
                    exon_id = (entry.attributes.split(";")[7]).split(" ")[2].strip('"')
                    if exon_id not in self.id_to_coordinates:
                        self.id_to_coordinates[exon_id] = Interval(entry.start, entry.end, entry.contig)
                    self.gene_to_exons[gene_id].add(Interval(entry.start, entry.end, entry.contig))
                    self.transcript_to_exons[transcript_id].add(Interval(entry.start, entry.end, entry.contig))


    def make_genomeview_row(self, start, end, chrom, strand, bams_list, 
                            padding_perc = 0.1, 
                            with_coverage = True,
                            with_TSS = True, 
                            include_secondary = False,
                            row = None, 
                            tighter_track = False):
        padding = math.ceil((end - start) * padding_perc)

        if row is None:
            row = genomeview.ViewRow("row")
        gene_view = genomeview.GenomeView(chrom, max(0, start - padding), end + padding, strand, self.source)
        gene_view.add_track(genomeview.track.TrackLabel(chrom + (" +" if strand else " -") + " : " + str(start - padding) + " - " + str(end + padding)))
        if self.annotation_path:
            gene_view.add_track(genomeview.BEDTrack(self.annotation_path, name="annot"))
        if with_TSS and self.tss_path:
            gene_view.add_track(genomeview.BEDTrack(self.tss_path, name="ref_TSS"))
        
        gene_view.add_track(genomeview.Axis())
        for key, value in bams_list.items():
            if with_coverage:
                coverage_track = genomeview.BAMCoverageTrack(value, name=key)
                gene_view.add_track(coverage_track)
            if tighter_track:
                bam_track = TighterSingleEndBAMTrack(value, name=key)
            else:
                bam_track = genomeview.SingleEndBAMTrack(value, name=key)
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
                                    with_TSS=True, 
                                    include_secondary = False,
                                    tighter_track = False):
        if data is not None:
            chrom = data[2]
            strand = data[3] if len(data) > 2 else True 
            start = data[1].begin
            end = data[1].end
        elif interval is not None:
            chrom = interval.data
            strand = True
            start = interval.begin
            end  = interval.end
        else:
            raise("Neither an Interval or data structure has been provided.")

        row = self.make_genomeview_row(start, end, chrom, strand, bams_list, 
                                       padding_perc = padding_perc, 
                                       with_coverage = with_coverage,
                                       with_TSS = with_TSS, 
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
                                   with_TSS = True, 
                                   include_secondary = False,
                                   tighter_track = False):
        row = genomeview.ViewRow("row")

        if interval_list is not None:
            for interval in interval_list:
                chrom = interval.data
                strand = True
                start = interval.begin
                end  = interval.end
                
                self.make_genomeview_row(start, end, chrom, strand, bams_list, 
                                         padding_perc = padding_perc, 
                                         with_coverage = with_coverage,
                                         with_TSS = with_TSS, 
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
                                         with_TSS = with_TSS, 
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
                      with_TSS = True, 
                      include_secondary = False,
                      view_width = 1600, 
                      tighter_track = False):
        doc = genomeview.Document(view_width)
        return self.add_single_view_row_to_plot(doc, bams_list, 
                                           interval = interval, 
                                           data = data, 
                                           padding_perc = padding_perc,
                                           with_coverage = with_coverage,
                                           with_TSS = with_TSS, 
                                           include_secondary = include_secondary,
                                           tighter_track = tighter_track)


    def plot_intervals(self, bams_list, 
                       interval_list = None, 
                       data_list = None, 
                       N_per_row = 1, 
                       padding_perc = 0.1,
                       with_coverage = True,
                       with_TSS = True, 
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
                                                      with_TSS = with_TSS, 
                                                      include_secondary = include_secondary,
                                                      tighter_track = tighter_track)
        elif data_list is not None:
            for i in range(0, len(data_list), N_per_row):  
                doc = self.add_multi_view_row_to_plot(doc, bams_list, 
                                                      interval_list = None, 
                                                      data_list = data_list[i:i+N_per_row], 
                                                      padding_perc = padding_perc, 
                                                      with_coverage = with_coverage,
                                                      with_TSS = with_TSS, 
                                                      include_secondary = include_secondary,
                                                      tighter_track = tighter_track)
        return doc


    def plot_feature(self, feature, bams_list, 
                     padding_perc = 0.1,
                     with_coverage = True,
                     with_TSS = True,
                     include_secondary = False,
                     view_width = 1600, 
                     tighter_track = False):
        return self.plot_interval(bams_list, 
                                  interval = self.id_to_coordinates[feature], 
                                  padding_perc = padding_perc,
                                  with_coverage = with_coverage,
                                  with_TSS = with_TSS, 
                                  include_secondary = include_secondary,
                                  view_width = view_width, 
                                  tighter_track = tighter_track)


    # plot_feature for a list in tabs
    def plot_features(self, features, bams_list, outformat="png",
                      **kwargs):

        params = inspect.signature(self.plot_feature).parameters
        filtered_kwargs = {k: v for k, v in kwargs.items() if k in params}

        features_tab = widgets.Tab()
        tab_contents = []

        for feature in features:
            tab_contents.append(self.plot_feature(feature, bams_list, **filtered_kwargs).get_widget(outformat))

        features_tab.children = tab_contents
        features_tab.titles = features

        return features_tab


    def plot_read(self, read_id, bam, **kwargs):

        regions = []

        with pysam.AlignmentFile(bam, "rb") as bam_in:
            for read in bam_in.fetch():
                if read.query_name != read_id:
                    continue
                regions.append(Interval(read.reference_start, read.reference_end, bam_in.get_reference_name(read.reference_id)))
            
        if len(regions) == 0:
            print("Error: read either not found or not aligned")
            return -1

        # Get the parameter names of plot_interval
        params = inspect.signature(self.plot_interval).parameters

        # Filter kwargs to only include keys that are in plot_interval's parameters
        filtered_kwargs = {k: v for k, v in kwargs.items() if k in params}

        all_widgets = []
        for region in regions:
            all_widgets.append(self.plot_interval(bams_list={"bam": bam}, interval=region).get_widget(), **filtered_kwargs)

        return widgets.VBox(all_widgets)



    def plot_exons(self, feature_to_plot, bams_list, merge_exons = True,
                   padding_perc = 0.05, 
                   with_coverage = True, 
                   with_TSS = False, 
                   include_secondary = False, 
                   view_width = 1600,
                   equal_size_exons = False,
                   N_per_row = 99999,
                   tighter_track = False):

        feature_id = None
        feature_type = None
        if feature_to_plot in self.gene_name_to_gene_id:
            feature_id = self.gene_name_to_gene_id[feature_to_plot]
            feature_type = "gene"
        elif feature_to_plot in self.gene_to_exons:
            feature_id = feature_to_plot
            feature_type = "gene"
        elif feature_to_plot in self.transcript_to_exons:
            feature_id = feature_to_plot
            feature_type = "transcript"
        elif feature_to_plot in self.id_to_coordinates:
            feature_id = feature_to_plot
            feature_type = "exon"
        
        if feature_type == "exon":
            return self.plot_feature(feature_id, bams_list=bams_list, padding_perc=padding_perc, with_coverage=with_coverage, with_TSS=with_TSS, include_secondary=include_secondary, view_width=view_width, tighter_track=tighter_track)
        elif feature_type == "transcript":

            #chrom = self.id_to_coordinates[self.transcript_to_exons[feature_id][0].data].data
            #chrom = None
            #for exon in self.transcript_to_exons[feature_id]:
            #    chrom = self.id_to_coordinates[exon.data].data
            #    break

            if equal_size_exons:
                # self.plot_intervals(bams_list, interval_list = [self.id_to_coordinates[exon_interval.data] for exon_interval in sorted(self.transcript_to_exons[feature_id])], N_per_row = N_per_row, padding_perc=padding_perc, with_coverage=with_coverage, with_TSS=with_TSS, include_secondary=include_secondary, tighter_track=tighter_track)
                self.plot_intervals(bams_list, interval_list = sorted(self.transcript_to_exons[feature_id]), N_per_row = N_per_row, padding_perc=padding_perc, with_coverage=with_coverage, with_TSS=with_TSS, include_secondary=include_secondary, tighter_track=tighter_track)
            else:
                doc = genomeview.Document(view_width)
                exons_list = sorted(self.transcript_to_exons[feature_id])
                for i in range(0, len(exons_list), N_per_row):
                    # exons_row = self.make_exons_row_from_virtual(exons_list[i:i+N_per_row], bams_list=bams_list, view_width=view_width, padding_perc=padding_perc, with_coverage=with_coverage, with_TSS=with_TSS, include_secondary=include_secondary, tighter_track=tighter_track)
                    # doc.elements.append(exons_row)
                    doc = self.make_exons_row_from_virtual(doc, exons_list[i:i+N_per_row], bams_list=bams_list, view_width=view_width, padding_perc=padding_perc, with_coverage=with_coverage, with_TSS=with_TSS, include_secondary=include_secondary, tighter_track=tighter_track)
                return doc

        elif feature_type == "gene":
            doc = genomeview.Document(view_width)
            if merge_exons:
                chrom = self.id_to_coordinates[feature_id].data
                tmp_exons = self.gene_to_exons[feature_id].copy()
                
                tmp_exons.merge_overlaps()
                exons_list_without_chroms = sorted(tmp_exons)
                exons_list = []
                for exon in exons_list_without_chroms:
                    if exon.data is None:
                        exons_list.append(Interval(exon.begin, exon.end, chrom))
                    else:
                        exons_list.append(exon)
            else:
                exons_list = sorted(self.gene_to_exons[feature_id])
                # exons_list = [self.id_to_coordinates[exon_interval.data] for exon_interval in sorted(self.gene_to_exons[feature_id])]
            
            #chrom = self.id_to_coordinates[self.gene_to_exons[feature_id][0].data].data
            #chrom = None
            #for exon in self.gene_to_exons[feature_id]:
            #    chrom = self.id_to_coordinates[exon.data].data
            #    break

            for i in range(0, len(exons_list), N_per_row):
                # exons_row = self.make_exons_row_from_virtual(exons_list[i:i+N_per_row], bams_list=bams_list, view_width=view_width, padding_perc=padding_perc, with_coverage=with_coverage, with_TSS=with_TSS, include_secondary=include_secondary, tighter_track=tighter_track)
                # doc.elements.append(exons_row)
                doc = self.make_exons_row_from_virtual(doc, exons_list[i:i+N_per_row], bams_list=bams_list, view_width=view_width, padding_perc=padding_perc, with_coverage=with_coverage, with_TSS=with_TSS, include_secondary=include_secondary, tighter_track=tighter_track)
            return doc
  


     

    def make_exons_row_from_virtual(self, doc, exons_list, bams_list, view_width,
                                    strand = True,
                                    padding_perc = 0.1, 
                                    with_coverage = True,
                                    with_TSS = True, 
                                    include_secondary = False,
                                    row = None, 
                                    tighter_track = False):


        if row is None:
            row = genomeview.ViewRow("row")

        total_exon_size = 0
        left_bound = math.inf
        right_bound = -math.inf
        # chrom = None

        smallest_exon_size = math.inf
        for exon_interval in exons_list:
            # chrom = self.id_to_coordinates[exon_interval.data].data
            total_exon_size += exon_interval.end - exon_interval.begin
            smallest_exon_size = min(smallest_exon_size, exon_interval.end - exon_interval.begin)
            left_bound = min(left_bound, exon_interval.begin)
            right_bound = max(right_bound, exon_interval.end)

        padding = math.ceil(smallest_exon_size * padding_perc)
        total_exon_size = total_exon_size + (padding * len(exons_list))

        reserved_width = (len(exons_list) - 1) * row.space_between + doc.margin_x * 2
        per_base_size = (view_width - reserved_width)/total_exon_size

        max_coverage_dict = {}
        virtual_bams_dict = {}
        for key, value in bams_list.items():

            virtual_bams_dict[key] = []
            all_reads_for_coverage = set()

            bam_refs = None
            with pysam.AlignmentFile(value, "rb") as bam:
                bam_refs = bam.references
                for exon_interval in exons_list:
                    for read in bam.fetch(exon_interval.data, exon_interval.begin, exon_interval.end):
                        all_reads_for_coverage.add(read)

                coverage_bam = genomeview.bamtrack.VirtualBAM(all_reads_for_coverage, bam_refs)
                coverage_bam.index()

                for read in coverage_bam.fetch():
                    virtual_bams_dict[key].append(genomeview.bamtrack.VirtualBAM([read], bam_refs))
                
                max_coverage_dict[key] = get_virtualbam_max_coverage(coverage_bam)
        
        for exon_interval in exons_list:
            start = exon_interval.begin
            end = exon_interval.end
            chrom = exon_interval.data

            exon_width = math.floor((exon_interval.end - exon_interval.begin + padding) * per_base_size)
            exon_view = genomeview.GenomeView(chrom, start - padding, end + padding, strand, source=self.source)
            exon_view.add_track(genomeview.track.TrackLabel(chrom + (" +" if strand else " -") + " : " + str(start - padding) + " - " + str(end + padding)))
            # BED type features
            if self.annotation_path:
                virtual_bed = genomeview.bedtrack.VirtualBEDTrack(self.annotation_path)
                virtual_bed.index(chrom, left_bound, right_bound, field_defs=None)
                exon_view.add_track(virtual_bed)
            if with_TSS and self.tss_path:
                virtual_bed = genomeview.bedtrack.VirtualBEDTrack(self.tss_path)
                virtual_bed.index(chrom, left_bound, right_bound, fields_defs=None)
                exon_view.add_track(virtual_bed)
            
            exon_view.add_track(genomeview.Axis())
            for key, value in bams_list.items():
                
                if with_coverage:
                    coverage_track = genomeview.BAMCoverageTrack(value, name=key)
                    coverage_track.max_y = max_coverage_dict[key]
                    exon_view.add_track(coverage_track)
                for virtual_bam in virtual_bams_dict[key]:
                    if tighter_track:
                        bam_track = TighterSingleEndBAMTrack(virtual_bam, name=None, opener_fn=lambda x: x)
                    else:
                        bam_track = genomeview.bamtrack.SingleEndBAMTrack(virtual_bam, name=None, opener_fn=lambda x: x)
                    if include_secondary:
                        coverage_track.include_secondary = True
                        bam_track.include_secondary = True
                    exon_view.add_track(bam_track)
            exon_view.pixel_width = exon_width
            exon_view.margin_y = 0
            row.add_view(exon_view)

        doc.elements.append(row)
        return doc



