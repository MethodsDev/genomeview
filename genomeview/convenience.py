import collections
import gzip
import math
import os
import pysam

from intervaltree import Interval

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


class TighterSingleEndBAMTrack(genomeview.SingleEndBAMTrack):
    def __init__(self, *args, **kwdargs):
        super().__init__(*args, **kwdargs)
        self.row_height = 3
        self.margin_y = 1.5


class Configuration():
    def __init__(self, genome_path, annotation_path, gtf_path = None, tss_path = None):
        self.source = genomeview.genomesource.FastaGenomeSource(genome_path)
        self.annotation_path = annotation_path
        self.tss_path = tss_path
        if gtf_path:
            self.index_gtf(gtf_path)

    def index_gtf(self, gtf_path):
        self.transcript_to_gene = {}
        self.gene_to_transcripts = {}
        self.id_to_coordinates = {}

        with gzip.open(gtf_path, "r") as gtf_file:
            current_gene_interval = None
            for entry in pysam.tabix_iterator(gtf_file, pysam.asGTF()):
                if entry.feature == "gene":
                    gene_id = (entry.attributes.split(";")[0]).split(" ")[1].strip('"')
                    self.id_to_coordinates[gene_id] = Interval(entry.start-1, entry.end-1, entry.contig)
                    self.gene_to_transcripts[gene_id] = []
                elif entry.feature == "transcript":
                    gene_id = (entry.attributes.split(";")[0]).split(" ")[1].strip('"')
                    transcript_id = (entry.attributes.split(";")[1]).split(" ")[2].strip('"')
                    self.id_to_coordinates[transcript_id] = Interval(entry.start-1, entry.end-1, entry.contig)
                    self.transcript_to_gene[transcript_id] = gene_id
                    self.gene_to_transcripts[gene_id].append(transcript_id)


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
        gene_view.add_track(genomeview.track.TrackLabel(chrom + (" +" if strand else " -") + " : " + str(start) + " - " + str(end)))
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
            for interval in intervals_list:
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
                       view_width = 1600):
        doc = genomeview.Document(view_width)

        if interval_list is not None:
            for i in range(0, len(interval_list), N_per_row):  
                doc = self.add_multi_view_row_to_plot(doc, bams_list, 
                                                      interval = interval_list[i:i+N_per_row], 
                                                      data = None, 
                                                      padding_perc = padding_perc,
                                                      with_coverage = with_coverage,
                                                      with_TSS = with_TSS, 
                                                      include_secondary = include_secondary,
                                                      tighter_track = tighter_track)
        elif data_list is not None:
            for i in range(0, len(data_list), N_per_row):  
                doc = self.add_multi_view_row_to_plot(doc, bams_list, 
                                                      interval = None, 
                                                      data = data_list[i:i+N_per_row], 
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
                     squish_reads = False):
        return self.plot_interval(bams_list, 
                                  interval = self.id_to_coordinates[feature], 
                                  padding_perc = padding_perc,
                                  with_coverage = with_coverage,
                                  with_TSS = with_TSS, 
                                  include_secondary = include_secondary,
                                  view_width = view_width, 
                                  tighter_track = squish_reads)



