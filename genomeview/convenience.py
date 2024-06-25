import collections
from abc import ABC, abstractmethod
from collections.abc import MutableMapping
import gzip
import math
import os
import pysam
import inspect
import ipywidgets as widgets
import re
from functools import partial

from intervaltree import Interval, IntervalTree
from Bio.Seq import Seq

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



# new general helpers
def my_hook_compressed(filename, mode):
    if 'b' not in mode:
        mode += 't'
    ext = os.path.splitext(filename)[1]
    if ext == '.gz':
        import gzip
        return gzip.open(filename, mode)
    elif ext == '.bz2':
        import bz2
        return bz2.open(filename, mode)
    else:
        return open(filename, mode)



# def flatten(dictionary, parent_key='', separator='_'):
#     items = []
#     for key, value in dictionary.items():
#         new_key = parent_key + separator + key if parent_key else key
#         if isinstance(value, MutableMapping):
#             items.extend(flatten(value, new_key, separator=separator).items())
#         else:
#             items.append((new_key, value))
#    return dict(items)

def flatten(dictionary, separator='_'):
    items = []
    for key, value in dictionary.items():
        if isinstance(value, MutableMapping):
            for el in flatten(value, separator=separator):
                items.append(key + separator + el)
        elif isinstance(value, list):
            for el in value:
                items.append(key + separator + el)
        else:
            items.append(key + separator + value)
    return items


### newly added wrappers

def get_regions_by_read_id(bam_file, read_id, opener_fn=pysam.AlignmentFile):
    regions = []

    with opener_fn(bam) as bam_in:
        for read in bam_in.fetch():
            if read.query_name != read_id:
                continue

            regions.append(Interval(read.reference_start, read.reference_end, bam_in.get_reference_name(read.reference_id) + "+" if read.is_forward else "-"))

    return(regions)


def find_read_in_bam(read_id, bams_dict, silence_error=False):
    regions = []
    virtual_bams = []

    for key, value in bams_dict.items():
        if isinstance(value, genomeview.VirtualBAM):
            opener_fn = lambda x: x
        else:
            opener_fn = pysam.AlignmentFile

        with opener_fn(value) as bam_in:
            bam_refs = bam_in.references
            for read in bam_in.fetch():
                if read.query_name != read_id:
                    continue
                regions.append(Interval(read.reference_start, read.reference_end, read.reference_name + ("+" if read.is_forward else "-")))
                virtual_bams.append(genomeview.VirtualBAM([read], bam_refs))
    
    if len(regions) == 0 and not silence_error:
        print("Error: read either not found or not aligned")
        return (-1, -1)

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


def get_bam_opener(bam):
    if isinstance(bam, genomeview.VirtualBAM):
        opener_fn = lambda x: x
    else:
        opener_fn = pysam.AlignmentFile
    return opener_fn


def make_bed_track(bed, name=None): # , chrom=None, start=None, end=None):
    if isinstance(bed, genomeview.VirtualBEDTrack):
        # bed_track = genomeview.VirtualBEDTrack(bed, name)
        bed_track = bed
        # bed_track.index(chrom, start, end)
    else:
        bed_track = genomeview.BEDTrack(bed, name)
    return bed_track


# similar to regular pysam.AlignmentSegment.get_tag() but instead of raising a KeyError when tag is missing, just return None
def get_read_tag(read, tag):
    if not read.has_tag(tag):
        return None

    return read.get_tag(tag)


class CellBarcode(ABC):
    @abstractmethod
    def get_barcode(read):
    # def get_barcode(self, read):
        pass


class HaasStyleCellBarcode(CellBarcode):
    def get_barcode(read):
    # def get_barcode(self, read):
        return Seq(read.query_name.split("^")[0]).reverse_complement()


class ONTCellBarcode(CellBarcode):
    def get_barcode(read):
    # def get_barcode(self, read):
        return get_read_tag("BC")


# 10X, PipSeq
class StandardCellBarcode(CellBarcode):
    def get_barcode(read):
    # def get_barcode(self, read):
        return get_read_tag("CB")



def is_in_whitelist(ident, whitelist):
    if isinstance(whitelist, dict):
        for key, whitelist in whitelist.items():
            if ident in whitelist:
                return key
    elif isinstance(whitelist, list) or isinstance(whitelist, set):
        if ident in whitelist:
            return "all"
    return None


# whitelist_dict in the form of {condition: whitelisted_barcodes} or a plain list/set that gives the name "all"
def split_bam_by_cellbarcode_whitelist(bam_name,
                                       bam_file,
                                       interval,
                                       cellbarcode_whitelist=None,
                                       cellbarcode_from=StandardCellBarcode):
    if cellbarcode_whitelist is None:
        return {bam_name + "_all": bam_file}

    opener_fn = get_bam_opener(bam_file)

    tmp_reads = {}
    with opener_fn(bam_file) as bam:
        refs = bam.references
        for read in bam.fetch(interval.chrom, interval.begin, interval.end):
            if read.is_unmapped or read.is_secondary:
                continue

            cell_barcode = cellbarcode_from.get_barcode(read)
            if not cell_barcode:
                continue

            key = is_in_whitelist(cell_barcode, cellbarcode_whitelist)
            if key:
                if key not in tmp_reads:
                    tmp_reads[key] = []
                tmp_reads[key].append(read)
    
    virtual_bams_dict = {}
    for key, reads in tmp_reads.items():
        virtual_bams_dict[bam_name + "_" + key] = genomeview.VirtualBAM(reads, refs)

    return virtual_bams_dict


class Classification(ABC):
    @abstractmethod
    def get_classification(self, read):
        pass
        # should return a list of classifcations (because of possible ambiguous)


# read_assignments.tsv(.gz)
class IsoQuantClassification(Classification):
    ISOQUANT_READ_ASSIGNMENTS_DEFS = {
        "read_id":0,
        "chr":1, 
        "strand": 2,
        "isoform_id": 3,
        "gene_id":4, 
        "assignment_type":5,
        "assignment_events":6,
        "exons":7,
        "additional_info":8
    }

    def __init__(self, file_path, ambiguous_classification = True):
        self.read_to_gene_id_to_isoform_id = {}
        self.read_to_assignment_type = {}
        self.ambiguous_classification = ambiguous_classification
        self.index(file_path)

    # call index() on file after creating the object so that the different Isoquant files can be indexed on the same object
    def index(self, file_path):
        with my_hook_compressed(file_path, "rt") as f:
            for line in f:
                if line[0] == "#":
                    continue

                fields = line.rstrip().split("\t")
                
                values =  dict.fromkeys(self.ISOQUANT_READ_ASSIGNMENTS_DEFS)
                for field_name, field in self.ISOQUANT_READ_ASSIGNMENTS_DEFS.items():
                    cur_value = None
                    if len(fields) > field:
                        cur_value = fields[field]
                    values[field_name] = cur_value

                self.read_to_assignment_type[values['read_id']] = values['assignment_type']
        
                if values['isoform_id'] is not None and values['isoform_id'] != ".":
                    if values['read_id'] not in self.read_to_gene_id_to_isoform_id:
                        self.read_to_gene_id_to_isoform_id[values['read_id']] = {}
                    if values['gene_id'] not in self.read_to_gene_id_to_isoform_id[values['read_id']]:
                        self.read_to_gene_id_to_isoform_id[values['read_id']][values['gene_id']] = []
                    self.read_to_gene_id_to_isoform_id[values['read_id']][values['gene_id']].append(values['isoform_id'])

    def get_classification(self, read):
        if read.query_name not in self.read_to_gene_id_to_isoform_id:
            return None

        if self.read_to_assignment_type[read.query_name] == "ambiguous":
            return ["ambiguous_isoform"]

        if self.read_to_assignment_type[read.query_name] == "inconsistent":
            return ["inconsistent_gene"]

        return flatten(self.read_to_gene_id_to_isoform_id[read.query_name])



class BAMtagClassification(Classification):
    def __init__(tag):
        self.tag = tag

    def get_classification(self, read):
        return get_read_tag(read, self.tag)


def split_bam_by_classification(bam_file,
                                bam_name,
                                interval,
                                classification_from,
                                cellbarcode_whitelist = None,
                                cellbarcode_from = StandardCellBarcode):
    if classification_from is None:
        print("No way of getting classification provided")
        return -1

    opener_fn = get_bam_opener(bam_file)

    tmp_reads = {}
    with opener_fn(bam_file) as bam:
        refs = bam.references
        for read in bam.fetch(interval.chrom, interval.begin, interval.end):
            if read.is_unmapped or read.is_secondary:
                continue

            if cellbarcode_whitelist is not None:
                cell_barcode = cellbarcode_getter.get_barcode(read)
                if not cell_barcode:
                    continue
                if is_in_whitelist(cell_barcode, cellbarcode_whitelist) is None:
                    continue

            classifications = classification_from.get_classification(read)
            if classifications is None:
                if bam_name + "_unclassified" not in tmp_reads:
                    tmp_reads[bam_name + "_unclassified"] = []
                tmp_reads[bam_name + "_unclassified"].append(read)

            else:
                for classification in classifications:
                    if bam_name + "_" + classification not in tmp_reads:
                        tmp_reads[bam_name + "_" + classification] = []
                    tmp_reads[bam_name + "_" + classification].append(read)

    
    virtual_bams_dict = {}
    for key, reads in tmp_reads.items():
        virtual_bams_dict[key] = genomeview.VirtualBAM(reads, refs)

    return virtual_bams_dict



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



class Configuration:
    def __init__(self, genome_fasta, bed_annotation, gtf_annotation = None, bed_color_fn=color_from_bed):
        self.genome_fasta = genomeview.genomesource.FastaGenomeSource(genome_fasta)
        self.bed_annotation = bed_annotation
        self.bed_color_fn = bed_color_fn

        self.gene_name_to_gene_id = {}
        self.gene_id_to_gene_name = {}
        self.gene_to_transcripts = {}
        self.transcript_to_gene = {}
        self.gene_to_exons = {}
        self.transcript_to_exons = {}
        self.id_to_coordinates = {}

        if gtf_annotation:
            self.index_gtf(gtf_annotation)


    def shallow_copy(self):
        return Configuration(self.genome_fasta, None, bed_color_fn = self.bed_color_fn)


    def index_gtf(self, gtf_annotation):
        gene_id_regex = re.compile('gene_id "([a-zA-Z0-9\.\-_]+)";')
        gene_name_regex = re.compile('gene_name "([a-zA-Z0-9\.\-_]+)";')
        transcript_id_regex = re.compile('transcript_id "([a-zA-Z0-9\._\^\-]+=?)";')
        exon_id_regex = re.compile('exon_id "([a-zA-Z0-9\.\-_]+)";')
        # transcript_name_regex = re.compile('transcript_name "([a-zA-Z0-9\.]+)";')
        # protein_id_regex = re.compile('protein_id "([a-zA-Z0-9\.]+)";')

        with gzip.open(gtf_annotation, "r") as gtf_file:
            for entry in pysam.tabix_iterator(gtf_file, pysam.asGTF()):
                if entry.feature == "gene":
                    res = gene_id_regex.search(entry.attributes)
                    if res:
                        gene_id = res.group(1)
                    else:
                        print("missing gene_id in a gene entry, skipping entry")


                    res = gene_name_regex.search(entry.attributes)
                    if res:
                        gene_name = res.group(1)
                        self.gene_name_to_gene_id[gene_name] = gene_id
                        self.gene_id_to_gene_name[gene_id] = gene_name

                    self.id_to_coordinates[gene_id] = Interval(entry.start, entry.end, entry.contig + entry.strand)
                    self.gene_to_transcripts[gene_id] = []
                    self.gene_to_exons[gene_id] = IntervalTree()

                elif entry.feature == "transcript":
                    gene_id = None
                    res = gene_id_regex.search(entry.attributes)
                    if res:
                        gene_id = res.group(1)

                        if gene_id not in self.id_to_coordinates:
                            self.id_to_coordinates[gene_id] = Interval(entry.start, entry.end, entry.contig + entry.strand)
                            self.gene_to_transcripts[gene_id] = []
                            self.gene_to_exons[gene_id] = IntervalTree()

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


    def update_bed(self, bed_annotation):
        self.bed_annotation = bed_annotation


    def get_bed_entries(self, interval):
        all_known_annotations = {}
        if self.bed_annotation:
            if type(self.bed_annotation) is list:
                for bed_path in self.bed_annotation:
                    for bed_entry in genomeview.bedtrack.bed_fetch(bed_path, interval.chrom, interval.begin, interval.end):
                        all_known_annotations[bed_entry.name] = genomeview.VirtualBEDTrack(transcripts = [bed_entry])
            elif type(self.bed_annotation) is dict:
                for bed_name, bed_path in self.bed_annotation.items():
                    for bed_entry in genomeview.bedtrack.bed_fetch(bed_path, interval.chrom, interval.begin, interval.end):
                        all_known_annotations[bed_entry.name] = genomeview.VirtualBEDTrack(transcripts = [bed_entry])
            else:
                for bed_entry in genomeview.bedtrack.bed_fetch(self.bed_annotation, interval.chrom, interval.begin, interval.end):
                    all_known_annotations[bed_entry.name] = genomeview.VirtualBEDTrack(transcripts = [bed_entry])
        return all_known_annotations


    def add_bed_tracks_to_view(self, view, vertical_layout=True, use_names=True):
        if self.bed_annotation:
            if type(self.bed_annotation) is list:
                for bed_path in self.bed_annotation:
                    bed_track = make_bed_track(bed_path)
                    bed_track.color_fn = self.bed_color_fn
                    bed_track.vertical_layout = vertical_layout
                    view.add_track(bed_track)
            elif type(self.bed_annotation) is dict:
                for bed_name, bed_path in self.bed_annotation.items():
                    if use_names:
                        view.add_track(genomeview.track.TrackLabel(bed_name))
                    else:
                        view.add_track(genomeview.track.TrackLabel(""))
                    virtual_bed = make_bed_track(bed_path, name="")
                    virtual_bed.color_fn = self.bed_color_fn
                    virtual_bed.vertical_layout = vertical_layout
                    view.add_track(virtual_bed)
            else:
                virtual_bed = make_bed_track(self.bed_annotation)
                virtual_bed.color_fn = self.bed_color_fn
                virtual_bed.vertical_layout = vertical_layout
                view.add_track(virtual_bed)


    def make_genomeview_row(self, start, end, chrom, strand, bams_dict,
                            padding_perc = 0.1,
                            add_track_label = "auto",
                            add_reads_label = "auto",
                            add_coverage_label = "auto",
                            with_reads = True,
                            with_axis = True,
                            with_coverage = True,
                            with_bed = True,
                            with_bed_label = False,
                            vertical_layout_reads = False,
                            include_secondary = False,
                            quick_consensus = True,
                            row = None, 
                            view_width = None,
                            view_margin_y = None,
                            fill_coverage = False,
                            coverage_track_max_y = None,
                            tighter_track = False,
                            **kwargs):

        padding = math.ceil((end - start) * padding_perc)

        if row is None:
            row = genomeview.ViewRow("row")
        gene_view = genomeview.GenomeView(chrom, max(0, start - padding), end + padding, "+", self.genome_fasta)
        # gene_view = genomeview.GenomeView(chrom, start - padding, end + padding, "+", self.source)
        if add_track_label:
            if add_track_label == "auto":
                gene_view.add_track(genomeview.track.TrackLabel(chrom + " : " + str(start - padding) + " - " + str(end + padding)))
            else:
                gene_view.add_track(genomeview.track.TrackLabel(add_track_label))

        if with_bed:
            self.add_bed_tracks_to_view(gene_view, use_names = with_bed_label)

        if with_axis:
            gene_view.add_track(genomeview.Axis())
        for key, value in bams_dict.items():
            if isinstance(value, genomeview.VirtualBAM):
                opener_kwargs = {'opener_fn': lambda x: x}
            else:
                opener_kwargs = {}


            if with_coverage:
                if add_coverage_label:
                    if add_coverage_label == "auto":
                        add_coverage_label = key
                else:
                    add_coverage_label = ""
                coverage_track = genomeview.BAMCoverageTrack(value, name=add_coverage_label, **opener_kwargs)

                if fill_coverage:
                    coverage_track.fill_coverage = True

                if coverage_track_max_y:
                    if type(coverage_track_max_y) is dict:
                        coverage_track.max_y = coverage_track_max_y[key]
                    else:
                        coverage_track.max_y = coverage_track_max_y
                gene_view.add_track(coverage_track)
            if with_reads:
                if add_reads_label:
                    if add_reads_label == "auto":
                        add_reads_label = key
                else:
                    add_reads_label = ""


                if tighter_track:
                    bam_track = TighterSingleEndBAMTrack(value, name=add_reads_label, **opener_kwargs)
                else:
                    bam_track = genomeview.SingleEndBAMTrack(value, name=add_reads_label, **opener_kwargs)
                if include_secondary:
                    coverage_track.include_secondary = True
                    bam_track.include_secondary = True
                bam_track.quick_consensus = quick_consensus
                bam_track.vertical_layout = vertical_layout_reads
                gene_view.add_track(bam_track)

        if view_width:
            gene_view.pixel_width = view_width
        if view_margin_y:
            gene_view.margin_y = view_margin_y

        row.add_view(gene_view)
        return row


    def add_single_view_row_to_plot(self, doc,
                                    interval = None, 
                                    data = None, **kwargs):

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

        row = self.make_genomeview_row(start, end, chrom, strand, **kwargs)
        doc.elements.append(row)
        return doc


    def add_multi_view_row_to_plot(self, doc,
                                   interval_list = None, 
                                   data_list = None, **kwargs):

        row = genomeview.ViewRow("row")

        if interval_list is not None:
            for interval in interval_list:
                chrom = interval.chrom
                strand = interval.strand
                start = interval.begin
                end  = interval.end
                
                self.make_genomeview_row(start, end, chrom, strand,
                                         row = row,  **kwargs)

        elif data_list is not None:
            for data in data_list:
                chrom = data[2]
                strand = data[3] if len(data) > 2 else True 
                start = data[1].begin
                end = data[1].end
                
                self.make_genomeview_row(start, end, chrom, strand,
                                         row = row,  **kwargs)

        doc.elements.append(row)
        return doc


    def plot_interval(self, view_width = 1600, **kwargs):

        doc = genomeview.Document(view_width)
        return self.add_single_view_row_to_plot(doc, **kwargs)



    def plot_intervals(self,
                       interval_list = None, 
                       data_list = None, 
                       N_per_row = 1,  **kwargs):

        doc = genomeview.Document(view_width)

        if interval_list is not None:
            for i in range(0, len(interval_list), N_per_row):
                doc = self.add_multi_view_row_to_plot(doc,
                                                      interval_list = interval_list[i:i+N_per_row], 
                                                      data_list = None,
                                                       **kwargs)

        elif data_list is not None:
            for i in range(0, len(data_list), N_per_row):  
                doc = self.add_multi_view_row_to_plot(doc,
                                                      interval_list = None, 
                                                      data_list = data_list[i:i+N_per_row],  **kwargs)
        return doc

    def get_feature_info(self, feature):
        feature_id = None
        feature_type = None
        if feature in self.gene_name_to_gene_id:
            feature_id = self.gene_name_to_gene_id[feature]
            feature_type = "gene"
        elif feature in self.gene_id_to_gene_name:
            feature_id = feature
            feature_type = "gene"
        elif feature in self.transcript_to_exons:
            feature_id = feature
            feature_type = "transcript"
        elif feature in self.id_to_coordinates:
            feature_id = feature
            feature_type = "exon"
        
        return(feature_id, feature_type)


    def get_interval_from_feature(self, feature):
        if isinstance(feature, tuple):
            (feature_id, feature_type) = feature
        else:
            (feature_id, feature_type) = self.get_feature_info(feature)
        return self.id_to_coordinates[feature_id]


    def get_gene_name(self, feature_info):
        if isinstance(feature_info, tuple):
            (feature_id, feature_type) = feature_info
        else:
            (feature_id, feature_type) = self.get_feature_info(feature_info)

        if feature_type == "transcript":
            if feature_id not in self.transcript_to_gene:
                print("transcript not associated with any gene")
                return feature_id
            feature_id = self.transcript_to_gene[feature_id]
            feature_type = "gene"

        if feature_type == "gene":
            if feature_id in self.gene_name_to_gene_id:
                return feature_id
            elif feature_id in self.gene_id_to_gene_name:
                return self.gene_id_to_gene_name[feature_id]
            else:
                print("no gene name associated with this gene id")
                return feature_id

        print("unknown feature type")
        return feature_id


    def plot_feature(self, feature, **kwargs):

        (feature_id, feature_type) = self.get_feature_info(feature)

        return self.plot_interval(interval = self.id_to_coordinates[feature_id], **kwargs)

    # plot_feature for a list in tabs
    def plot_features(self, 
                      features,
                      output_format="svg",
                      **kwargs):

        # params = inspect.signature(self.plot_feature).parameters
        # filtered_kwargs = {k: v for k, v in kwargs.items() if k in params}

        features_tab = widgets.Tab()
        tab_contents = []

        for feature in features:
            tab_contents.append(self.plot_feature(feature, **kwargs).get_widget(output_format))

        features_tab.children = tab_contents
        features_tab.titles = features

        return features_tab


    def plot_read(self, read_id, bams_dict, interval="auto", output_format="svg", silence_error=False, **kwargs):

        # Get the parameter names of plot_interval
        # params = inspect.signature(find_read_in_bam).parameters
        # Filter kwargs to only include keys that are in plot_interval's parameters
        # filtered_kwargs = {k: v for k, v in kwargs.items() if k in params}

        (regions, virtual_bams) = find_read_in_bam(read_id, bams_dict, silence_error)

        if regions == -1:
            return None


        all_widgets = []
        for region, virtual_bam in zip(regions, virtual_bams):
            if isinstance(interval, str) and interval == "auto":
                tmp = self.plot_interval(bams_dict={"read": virtual_bam}, interval=region, **kwargs).get_widget(output_format)
            else:
                tmp = self.plot_interval(bams_dict={"read": virtual_bam}, interval=interval, **kwargs).get_widget(output_format)
            if tmp is not None:
                all_widgets.append(tmp)

        return widgets.VBox(all_widgets)


    def plot_reads(self, read_ids, bams_dict, interval, output_format="svg", **kwargs):
        if not isinstance(interval, Interval):
            print("Error, Interval() required but not provided")
            return -1

        first = False
        all_widgets = []

        all_widgets.append(self.plot_interval(interval=interval, bams_dict={}, **kwargs).get_widget(output_format))

        for bam_name, bam_file in bams_dict.items():
            for read_id in read_ids:
                all_widgets.extend(self.plot_read(read_id, {bam_name: bam_file}, interval, output_format,
                                                                                            silence_error=True,
                                                                                            with_coverage = False,
                                                                                            with_axis = False,
                                                                                            with_bed = False,
                                                                                            add_track_label = False,
                                                                                            add_reads_label = False,
                                                                                            **kwargs).children)
        return widgets.VBox(all_widgets)


    def plot_exons(self, 
                   feature,
                   merge_exons = True,
                   N_per_row = 99999,
                   view_width = 1600,
                   as_widget = False,
                   **kwargs):

        (feature_id, feature_type) = self.get_feature_info(feature)
        
        if feature_type == "exon":
            return self.plot_feature(feature_id, **kwargs)

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
                doc = self.make_intervals_row_through_virtual(genomeview.Document(view_width), [exon], **kwargs)
                all_views.append(widgets.HTML(doc._repr_svg_()))
                all_titles.append("Exon:: " + exon.data + " : " + str(exon.begin) + " - " + str(exon.end))

            stack = widgets.Stack(all_views, selected_index=0)
            dropdown = widgets.Dropdown(options=all_titles)
            widgets.jslink((dropdown, 'index'), (stack, 'selected_index'))
            return(widgets.VBox([dropdown, stack]))

        else:
            doc = genomeview.Document(view_width)
            for i in range(0, len(exons_list), N_per_row):
               self.make_intervals_row_through_virtual(doc, exons_list[i:i+N_per_row], **kwargs)
            return doc



    def plot_splice_junctions(self, 
                              feature,
                              view_width = 1600,
                              as_widget = False, 
                              **kwargs):


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
                doc = self.make_intervals_row_through_virtual(genomeview.Document(view_width), pair, **kwargs)
                all_views.append(widgets.HTML(doc._repr_svg_()))
                all_titles.append("Splice junction btw: exon:: " + pair[0].data + ":" + str(pair[0].begin) + "-" + str(pair[0].end) + " and exon::" + pair[1].data + " : " + str(pair[1].begin) + " - " + str(pair[1].end))

            stack = widgets.Stack(all_views, selected_index=0)
            dropdown = widgets.Dropdown(options=all_titles)
            widgets.jslink((dropdown, 'index'), (stack, 'selected_index'))
            return(widgets.VBox([dropdown, stack]))

        else:
            doc = genomeview.Document(view_width)
            for pair in exons_pairs:
               self.make_intervals_row_through_virtual(doc, pair, **kwargs)
            return doc
     

    def make_intervals_row_through_virtual(self,
                                           doc,
                                           intervals_list,
                                           bams_dict,
                                           padding_perc = 0.1, 
                                           add_track_label = "auto",
                                           add_reads_label = "auto",
                                           add_coverage_label = "auto",
                                           with_coverage = True,
                                           include_secondary = False,
                                           row = None, 
                                           normalize_interval_width = False,
                                           shared_max_coverage = False,
                                           **kwargs):


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
        for key, value in bams_dict.items():

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

        if shared_max_coverage:
            if type(shared_max_coverage) is bool:  # is True but not a value
                max_coverage_dict = max(max_coverage_dict.values())
            else:
                max_coverage_dict = shared_max_coverage

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

            row = self.make_genomeview_row(start = start, 
                                           end = end,
                                           chrom = chrom,
                                           strand = strand,
                                           bams_dict = bams_dict, 
                                           padding_perc = 0, 
                                           add_track_label = add_track_label,
                                           add_reads_label = add_reads_label,
                                           add_coverage_label = add_coverage_label,
                                           row = row, 
                                           view_width = interval_width, 
                                           view_margin_y = 0,
                                           coverage_track_max_y = max_coverage_dict,
                                           **kwargs)

            if add_track_label:
                add_track_label = "\n"
            add_reads_label = None
            add_coverage_label = None

        doc.elements.append(row)
        return doc


    def get_gene_tab_title(self, feature_name):

        (feature_id, feature_type) = self.get_feature_info(feature_name)
        gene_name = self.get_gene_name((feature_id, feature_type))
        if gene_name != feature_name:
            feature_name = gene_name  # + "_" + feature_name
        return feature_name
                

    # custom bed dict accepts a dict with the same keys are bams_dict only
    def organize_tab_section(self,
                             bams_dict,
                             interval,
                             tab_name,
                             custom_bed_dict = None,
                             with_coverage = True,
                             **kwargs):

        bed_config = self.shallow_copy()

        if custom_bed_dict is not None:
            all_bed_entries = genomeview.VirtualBEDTrack()
            for virtual_bed in custom_bed_dict.values():
                all_bed_entries.transcripts.extend(virtual_bed.transcripts)
            bed_config.update_bed(all_bed_entries)
        else:
            bed_config.update_bed(self.bed_annotation)

        shared_static_svg = bed_config.plot_interval(bams_dict={},
                                                     interval = interval,
                                                     with_bed = True,
                                                     with_reads = False,
                                                     with_coverage = False,
                                                     add_track_label = False,
                                                     **kwargs
                                                    )._repr_svg_()
        
        tab_sections = []
        for key, bam in bams_dict.items():
            if bam is None:
                continue
            unique_id = f"{tab_name}_{key}"

            static_svg = ""
            resizable_svg = ""
            if custom_bed_dict is not None and key in custom_bed_dict:
                bed_config.update_bed(custom_bed_dict[key])
                static_svg += bed_config.plot_interval(bams_dict={},
                                                       interval = interval,
                                                       with_bed = True,
                                                       with_reads = False,
                                                       with_coverage = False,
                                                       add_track_label = False,
                                                       **kwargs
                                                      )._repr_svg_()

            if with_coverage:
                static_svg += self.plot_interval(bams_dict = {key: bam},
                                                 interval = interval,
                                                 with_reads = False,
                                                 with_coverage = True,
                                                 with_bed = False,
                                                 add_track_label = False,
                                                 fill_coverage = True,
                                                 **kwargs
                                                )._repr_svg_() + "</br>"

            resizable_svg += self.plot_interval(bams_dict = {key: bam},
                                                interval = interval,
                                                with_reads = True,
                                                with_coverage = False,
                                                with_axis = False,
                                                with_bed = False,
                                                add_track_label = False,
                                                add_reads_label = False,
                                                vertical_layout_reads = True,
                                                **kwargs
                                               )._repr_svg_() + "</br>"

            tab_sections.append({
                'unique_id': unique_id,
                # 'name': f"{tab_name}_{key}",
                'static_svg': static_svg,
                'resizable_svg': resizable_svg,
                'expended': True
            })

        return {'tab_name': tab_name,
                'shared_static_svg': shared_static_svg,
                'tab_sections': tab_sections}



    def organize_tabs(self,
                      bams_dict_dict,
                      tab_title_fn = None,  # function that takes a feature_name/id or tuple(feature_id, feature_type) as input
                      custom_bed_dict_dict = None,
                      **kwargs):

        if tab_title_fn is None:  # workaround because can't set a self.method as default parameter
            tab_title_fn = self.get_gene_tab_title

        tabs = []
        for feature_name, bams_dict in bams_dict_dict.items():
            print(feature_name)
            (feature_id, feature_type) = self.get_feature_info(feature_name)
            interval = self.id_to_coordinates[feature_id]

            custom_bed_dict = None
            if custom_bed_dict_dict is not None and feature_name in custom_bed_dict_dict:
                custom_bed_dict = custom_bed_dict_dict[feature_name]
            # if custom_bed_dict_dict is None or feature_name not in custom_bed_dict_dict:
            #     custom_bed_dict = None
            # else:
            #     custom_bed_dict = custom_bed_dict_dict[feature_name]

            tabs.append(self.organize_tab_section(bams_dict = bams_dict,
                                                  interval = interval, 
                                                  tab_name = tab_title_fn(feature_id),
                                                  custom_bed_dict = custom_bed_dict,
                                                  **kwargs))

        return tabs


    def plot_by_features_as_tab(self,
                                bam_dict,
                                features_list,
                                page_title = "Plot by feature", # and split by barcode whitelist category if provided
                                cellbarcode_whitelist = None,
                                cellbarcode_from = StandardCellBarcode,  # can provide a (partial) callable or "tag:__" to use a BAM tag
                                **kwargs):

        virtual_bams_dict = {}
        if cellbarcode_from is not None and cellbarcode_whitelist is not None:

            for feature in features_list:

                interval = self.get_interval_from_feature(feature)
                virtual_bams_dict[feature] = {}

                for bam_name, bam_file in bam_dict.items():

                    #for bam_name, bam_file in bams_dict.items():
                    virtual_bams_dict[feature].update(split_bam_by_cellbarcode_whitelist(bam_name,
                                                                                         bam_file,
                                                                                         interval,
                                                                                         cellbarcode_whitelist = cellbarcode_whitelist,
                                                                                         cellbarcode_from = cellbarcode_from))
        else:
            for feature in features_list:
                virtual_bams_dict[feature] = {}

                for bam_name, bam_file in bam_dict.items():
                    virtual_bams_dict[feature][bam_name] = bam_file


        tabs = self.organize_tabs(virtual_bams_dict, **kwargs)

        return genomeview.templates.render_tab_titles(tabs, page_title)



    #  def plot_single_isoform()


    # # single gene, single bam
    # def split_bam_by_classification(self,
    #                                 bam_file,
    #                                 interval,
    #                                 classification_from = IsoQuantClassification,
    #                                 whitelist = None,
    #                                 cellbarcode_from = None,
    #                                 **kwargs):
        
    #     virtual_bams_dict = sort_bam_reads_by_classification(bam_file, interval, classification_from, whitelist, cellbarcode_from)

    #     return virtual_bams_dict



    # returns two dict with same keys, one with dicts of split BAMs, and one with associated BED entries
    # here whitelist is just that, a whitelist, there is no split based on which dict key a barcode is found in
    def split_bams_dict_by_classification(self,
                                          bams_dict,
                                          gene,
                                          classification_from,
                                          # cellbarcode_whitelist = None,
                                          # cellbarcode_from = None,
                                          **kwargs):
        
        (feature_id, feature_type) = self.get_feature_info(gene)
        if feature_type != "gene":
            print("feature provided is not a gene")
            return -1

        interval = self.get_interval_from_feature((feature_id, feature_type))

        virtual_bams_dict = {}
        for bam_name, bam_file in bams_dict.items():
            virtual_bams_dict.update(split_bam_by_classification(bam_file = bam_file,
                                                                 bam_name = bam_name,
                                                                 interval = interval,
                                                                 classification_from = classification_from,
                                                                 **kwargs))

        # parse known annotations
        all_known_annotations = self.get_bed_entries(interval)


        virtual_bed_dict = {}
        for known_annotation, virtual_bed in all_known_annotations.items():
            for classification, virtual_bam in virtual_bams_dict.items():
                if "ambiguous" in classification or "unclassified" in classification:
                    continue
                else:
                    if known_annotation in classification or known_annotation.split("|")[0] in classification:
                        if classification not in virtual_bed_dict:
                            virtual_bed_dict[classification] = virtual_bed
                        else:
                            virtual_bed_dict[classification].transcripts.extend(virtual_bed.transcripts)

#        # assign known annotations to classifications by substring matching, could make this into a method parameter provided as input
#        for classification, virtual_bam in virtual_bams_dict.items():
#            if "ambiguous" in classification or "unclassified" in classification:  # special classifications that are preserved
#                continue  # making a virtualBED with None in it just means need to handle that case as well elsewhere
#                # virtual_bed_dict[classification] = None
#            else:
#                for known_annotation, virtual_bed in all_known_annotations.items():
#                    if classification in known_annotation:
#                        if classification not in annotation_matched_virtual_bams_dict:
#                            virtual_bed_dict[classification] = virtual_bed
#                        else:
#                            virtual_bed_dict[classification].transcripts.append(virtual_bed.transcripts)

        return (virtual_bams_dict, virtual_bed_dict)
        


    def plot_by_classification_over_features(self,
                                             bams_dict,
                                             features_list,
                                             classification_from,
                                             page_title = "Split by Classification, Plot by feature",
                                             **kwargs):

        virtual_bams_dict_dict = {}
        virtual_beds_dict_dict = {}
        for feature in features_list:
            (virtual_bams_dict_dict[feature], virtual_beds_dict_dict[feature]) = \
                            self.split_bams_dict_by_classification(bams_dict = bams_dict,
                                                                   gene = feature,
                                                                   classification_from = classification_from,
                                                                   **kwargs
                                                                   )

        tabs = self.organize_tabs(bams_dict_dict = virtual_bams_dict_dict,
                                  custom_bed_dict_dict = virtual_beds_dict_dict,
                                  **kwargs)

        return genomeview.templates.render_tab_titles(tabs, page_title)











