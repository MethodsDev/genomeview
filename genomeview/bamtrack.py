import collections
import logging
import pysam
import numpy as np
import pandas as pd

import os
from dataclasses import dataclass

from intervaltree import IntervalTree

from genomeview.track import Track
from genomeview.intervaltrack import Interval, IntervalTrack
from genomeview import MismatchCounts
from genomeview.utilities import match_chrom_format
from genomeview.graphtrack import GraphTrack, BINNED_COLORS


def allreads(read):
    return True

def color_by_strand(interval):
    # brightness = 0.2 + (cur_reads[0].mapq/40.0*0.8)

    if interval.strand == "-":
        color = "#8C8FCE"
        if interval.read.is_secondary:
            color = "#8CB0CE"
    else:
        color = "#E89E9D"
        if interval.read.is_secondary:
            color = "#E8C49D"
    return color

    
class SingleEndBAMTrack(IntervalTrack):
    """
    Displays bam as single-ended reads
    
    Attributes:
        nuc_colors (dict): defines the SVG colors used to display mismatched nucleotides
        insertion_color, deletion_color, clipping_color (str): SVG colors for insertions, 
             deletions and soft/hard clipping

        quick_consensus (bool): specify whether the quick consensus mode should be used. When activated, 
            mismatches wrt the reference genome are only shown when at least several reads support
            a variant at that position (useful when displaying high-error rate data types eg 
            pacbio). Only relevant if draw_mismatches is also True. (default: True)
        draw_mismatches (bool): whether to show mismatches with respect to the reference genome.
            (default: True).

        include_secondary (bool): whether to draw alignments specified as "secondary" in the BAM flags 
            (default: True).
        
        include_read_fn: callback function used to specify which reads should be included in 
            the display. The function takes as its only argument a read (pysam.AlignedSegment) 
            and returns True (yes, display the read) or False (no, don't display). If this 
            function is not specified, by default all reads are shown.
    """
    def __init__(self, bam_path, name=None, opener_fn=pysam.AlignmentFile):
        """
        Args:
            bam_path (str): path of the bam file to display
            name (str): name of the track (optional - use None if you don't want to specify a name)

        """
        super().__init__([], name=name)

        self.bam_path = bam_path
        # with pysam.AlignmentFile(bam_path) as bam:
        self.opener_fn = opener_fn
        with self.opener_fn(bam_path) as bam:
            self.bam_references = bam.references
        # self.bam = pysam.AlignmentFile(bam_path)
        self.intervals = self
        self.mismatch_counts = None
        
        self.nuc_colors = {"A":"blue", "C":"red", "G":"green", "T":"black", "N":"gray"}
        self.insertion_color = "purple"
        self.clipping_color = "cyan"
        self.deletion_color = "cyan"

        self.quick_consensus = True
        self.draw_mismatches = True
        self.include_secondary = False

        self.min_indel_size = 0
        self.min_insertion_label_size = 5
        self.min_cigar_line_width = 2
        
        self.draw_read_labels = False

        self.include_read_fn = allreads
        self.color_fn = color_by_strand
        
    def fetch(self):
        """
        Iterator over reads from the bam file

        Overload this method in subclasses to feed this track reads from a different source
        (for example, reads that are already in memory, rather than being read from a file).
        """
        chrom = self.match_chrom_format(self.scale.chrom)
        start, end = self.scale.start, self.scale.end
        
        # with pysam.AlignmentFile(self.bam_path) as bam:
        with self.opener_fn(self.bam_path) as bam:
            for read in bam.fetch(chrom, start, end):
                if not self.include_read_fn or self.include_read_fn(read):
                    yield read
        
    def __iter__(self):
        c = 0
        for i, read in enumerate(self.fetch()):
            c += 1
            if read.is_unmapped: continue
            if read.is_secondary and not self.include_secondary: continue
            id_ = read.query_name + str(i)
            interval = Interval(id_, self.scale.chrom, read.reference_start, read.reference_end, 
                                not read.is_reverse)
            interval.read = read
            if self.draw_read_labels:
                interval.label = read.query_name
            yield interval

    def match_chrom_format(self, chrom):
        """
        Ensures that the input argument `chrom` matches the chromosome name formatting in
        the bam file being visualized (ie "chr14" vs "14").
        """
        # return match_chrom_format(chrom, self.bam.references)
        return match_chrom_format(chrom, self.bam_references)
        
    def layout(self, scale):
        super().layout(scale)
        self.reset_mismatch_counts()

    def reset_mismatch_counts(self):
        if self.quick_consensus and self.draw_mismatches:
            self.mismatch_counts = MismatchCounts(
                self.scale.chrom, self.scale.start, self.scale.end)

            # workaround for some quirk of pysam with crams and large cigars
            # (or something like that, opening a fresh file handle seems to fix the issue)
            #bam = pysam.AlignmentFile(self.bam_path)
            # with pysam.AlignmentFile(bam_path) as bam:
            with self.opener_fn(self.bam_path) as bam:
                self.mismatch_counts.tally_reads(bam)

    def layout_interval(self, interval):
        super().layout_interval(interval)

    def draw_interval(self, renderer, interval):
        """
        Draw a read and then, if ``self.draw_mismatches`` is True, draw mismatches/indels 
        on top.
        """
        extra_args = {"title": "read id: " + interval.read.query_name + "\nCIGAR: " + interval.read.cigarstring}
        yield from super().draw_interval(renderer, interval, extra_args=extra_args)

        if self.draw_mismatches:
            yield from self._draw_cigar(renderer, interval)

    def _draw_mismatch(self, renderer, length, genome_position, sequence_position, yoffset, alnseq):
        extras = {"stroke":"none"}

        try:
            refseq = self.scale.get_seq(genome_position, genome_position + length)
        except AssertionError:
            logging.debug("Unable to get reference sequence; will not draw mismatches")
            return

        for i in range(length):
            if genome_position+i < self.scale.start: continue
            if genome_position+i >= self.scale.end: break

            alt = alnseq[sequence_position+i]
            ref = refseq[i]

            if alt != ref:
                curstart = self.scale.topixels(genome_position+i)
                curend = self.scale.topixels(genome_position+i+1)

                color = self.nuc_colors[alt]  # self.nuc_colors[alnseq[sequence_position+i]]

                if not self.mismatch_counts or alt=="N" or self.mismatch_counts.query(alt, genome_position+i):
                    width = max(curend-curstart, self.min_cigar_line_width)
                    midpoint = (curstart+curend)/2
                    extras["title"] = ref + " -> " + alt + "\nreference position: " + str(genome_position + i)
                    yield from renderer.rect(midpoint-width/2, yoffset, width, self.row_height, fill=color, 
                                             **extras)

    def _draw_deletion(self, renderer, length, genome_position, yoffset):
        extras = {"stroke":"none"}

        if length > self.min_indel_size:
            curstart = self.scale.topixels(genome_position)
            curend = self.scale.topixels(genome_position+length)

            if genome_position > self.scale.end: return
            if genome_position+length < self.scale.start: return
            
            width = max(curend-curstart, self.min_cigar_line_width*2)
            midpoint = (curstart+curend)/2
            ymid = yoffset+self.row_height/2

            yield from renderer.rect(midpoint-width/2, yoffset, width, self.row_height, fill="white", 
                                     **extras)
            yield from renderer.line(midpoint-width/2, ymid, midpoint+width/2, ymid, 
                                     stroke="black", **{"stroke-width":1})

    def _draw_insertion(self, renderer, length, genome_position, yoffset):
        if length > self.min_indel_size:
            curstart = self.scale.topixels(genome_position-0.5)
            curend = self.scale.topixels(genome_position+0.5)

            if genome_position > self.scale.end: return
            if genome_position < self.scale.start: return

            midpoint = (curstart+curend)/2

            stroke_width = 1
            width = stroke_width
            ibeam_extension = 1.0
            
            font_size = self.row_height * 0.8
            draw_label = False
            length_string = str(length)
            label_width = len(length_string) * font_size * 0.9

            if length >= self.min_insertion_label_size:
                draw_label = True
            else:
                if label_width < self.scale.relpixels(length*1.5):
                    draw_label = True
            if draw_label:
                width = label_width

            yield from renderer.line(
               midpoint-width/2-ibeam_extension, yoffset+stroke_width/2, 
               midpoint+width/2+ibeam_extension, yoffset+stroke_width/2, stroke=self.insertion_color, 
               **{"stroke-width":stroke_width})

            yield from renderer.rect(
               midpoint-width/2, yoffset, width, self.row_height, 
               fill=self.insertion_color, **{"stroke":"none"})

            yield from renderer.line(
               midpoint-width/2-ibeam_extension, yoffset+self.row_height-stroke_width/2, 
               midpoint+width/2+ibeam_extension, yoffset+self.row_height-stroke_width/2, 
               stroke=self.insertion_color, **{"stroke-width":stroke_width})

            if draw_label:
                yield from renderer.text(midpoint, yoffset+self.row_height*0.75, length_string,
                    size=font_size, fill="white", **{"font-weight":"bold"})

    def _draw_clipping(self, renderer, length, genome_position, yoffset, side, strand):
        extras = {"stroke":"none"}

        if length >= 5:
            # always draw clipping, irrespective of consensus sequence or mode
            if side == "left":
                curstart = self.scale.topixels(genome_position - length - 0.5)
                curend = self.scale.topixels(genome_position + 0.5)
                if strand == "+": # don't make an arrow, just a rectangle
                    strand = None
            elif side == "right":
                curstart = self.scale.topixels(genome_position - 0.5)
                curend = self.scale.topixels(genome_position+length + 0.5)
                if strand == "-": # don't make an arrow, just a rectangle
                    strand = None
            else:
                curstart = self.scale.topixels(genome_position - 0.5)
                curend = self.scale.topixels(genome_position + 0.5)

            width = max(curend-curstart, self.min_cigar_line_width*2)
            midpoint = (curstart+curend)/2

            left = midpoint-width/2
            top = yoffset
            if strand is None:
                yield from renderer.rect(midpoint-width/2, yoffset, width, self.row_height, fill=self.clipping_color,
                                         **extras)
            else:
                arrow_width = min(self.row_height / 2, self.margin_x * 0.7, self.scale.relpixels(30))
                direction = "right" if strand == "+" else "left"

                yield from renderer.block_arrow(left, top, width, self.row_height,
                    arrow_width=arrow_width, direction=direction,
                    fill=self.clipping_color, **{"stroke":"none"})


    def _draw_cigar(self, renderer, interval):
        """
        draw mismatches/insertions/deletions and clipping
        """
        read = interval.read
        # if read.is_secondary: return ## here
        
        # min_width = 2

        row = self.intervals_to_rows[interval.id]
        yoffset = row*(self.row_height+self.margin_y)

        genome_position = read.reference_start
        sequence_position = 0
        alnseq = read.query_sequence

        for code, length in read.cigartuples:
            length = int(length)
            if code == 0 and alnseq is not None: #"M":
                yield from self._draw_mismatch(renderer, length, genome_position, sequence_position, yoffset, alnseq)

                sequence_position += length
                genome_position += length
            elif code in [2,3]: #in "D" "N":
                # if not self.mismatch_counts or self.mismatch_counts.query("DEL", genome_position, genome_position+length+1):
                yield from self._draw_deletion(renderer, length, genome_position, yoffset)

                genome_position += length
            elif code == 1: # I
                # if not self.mismatch_counts or self.mismatch_counts.query("INS", genome_position-2, genome_position+2):
                yield from self._draw_insertion(renderer, length, genome_position, yoffset)

                sequence_position += length
            elif code in [4, 5]: #"HS":
                if sequence_position == 0:
                    yield from self._draw_clipping(renderer, length, genome_position, yoffset, "left", interval.strand)
                else:
                    yield from self._draw_clipping(renderer, length, genome_position, yoffset, "right", interval.strand)

                if code == 4:
                    sequence_position += length



@dataclass
class PileupRead:
    alignment: pysam.AlignedSegment
    query_position: int
    is_del: bool
    is_refskip: bool
    indel: int

@dataclass
class PileupColumn:
    pos: int
    pileups: list[PileupRead]
    
    def __post_init__(self):
        self.n = len(self.pileups)


class VirtualBAM():
    def __init__(self, reads, references):
        self.reads = reads
        self.references = references
        self.is_indexed = False
        self.reads_interval_tree = None
        self.full_reads_interval_tree = None

        # toggle this to not return reads that only align around the region but not actually inside the region when fetching/piluping
        self.aligned_chunks_only = False

        # toggle to return all reads stored in order, regardless of start/end
        self.dumb_fetch = False

    def __enter__(self):
        return self
 
    def __exit__(self, *args):
        return

    def index(self):
        if self.is_indexed:
            return
        self.reads_interval_tree = IntervalTree()
        self.full_reads_interval_tree = IntervalTree()
        for read in self.reads:
            interval_start = current_position = read.reference_start
            for cigar_code, length in read.cigartuples:
                if cigar_code in [0, 2, 7, 8]:  # M, D, =, X
                    current_position += length
                elif cigar_code == 3:  # N: skipped region from the reference
                    if interval_start != current_position:
                        self.reads_interval_tree.addi(interval_start, current_position, read)
                    current_position += length
                    interval_start = current_position
                # else irrelevant
            if interval_start != current_position:
                self.reads_interval_tree.addi(interval_start, current_position, read)
            self.full_reads_interval_tree.addi(read.reference_start, read.reference_end, read)
        self.is_indexed = True

    def fetch(self, chrom=None, start=None, end=None):
        seen = set()
        if self.dumb_fetch:
            for read in self.reads:
                yield read
            return

        if chrom is None:
            if self.is_indexed:
                if self.aligned_chunks_only:
                    for interval in sorted(self.reads_interval_tree[start:end]):
                        if interval.data not in seen:
                            seen.add(interval.data)
                            yield interval.data
                else:
                    for interval in sorted(self.full_reads_interval_tree[start:end]):
                        yield interval.data
            else:
                for read in self.reads:
                    yield read
        else:
            chrom = match_chrom_format(chrom, self.references)
            if self.is_indexed:
                if self.aligned_chunks_only:
                    for interval in sorted(self.reads_interval_tree[start:end]):
                        if interval.data.reference_name == chrom:
                            if interval.data not in seen:
                                seen.add(interval.data)
                                yield interval.data
                else:
                    for interval in sorted(self.full_reads_interval_tree[start:end]):
                        if interval.data.reference_name == chrom:
                            yield interval.data
            else:
                for read in self.reads:
                    if read.reference_name == chrom and read.reference_start < end and read.reference_end > start:
                        yield read


    def point_fetch(self, chrom=None, start=None, end=None):
        if chrom is None:
            if self.is_indexed:
                for interval in sorted(self.reads_interval_tree[start:end]):
                    yield interval.data
            else:
                for read in self.reads:
                    yield read
        else:
            chrom = match_chrom_format(chrom, self.references)
            if self.is_indexed:
                for interval in sorted(self.reads_interval_tree[start:end]):
                #for interval in self.reads_interval_tree[start]:
                    if interval.data.reference_name == chrom:
                        yield interval.data
            else:
                for read in self.reads:
                    if read.reference_name == chrom and read.reference_start < end and read.reference_end > start:
                        yield read


    # truncate is always True for the implementation, but have the argument existing for compatibility
    def pileup(self, chrom, start, end, truncate=True, min_base_quality=13, step_size=100):
        chrom = match_chrom_format(chrom, self.references)
    
        # for ref_pos in range(start, end, step_size):
        if step_size > end - start:
            step_size = end - start

        for window_start in range(start, end, step_size):

            window_end = window_start + step_size
            window_end = end if window_end > end else window_end

            pileups = np.empty(window_end - window_start, dtype=object)
            pileups[...] = [[] for _ in range(pileups.shape[0])]
                   
            for read in self.fetch(chrom, window_start, window_end):
                ref_position = read.reference_start
                query_position = 0  # Initialize query_position to the start of the read
                current_window_index = 0
    
                for cigar_code, length in read.cigartuples:
                    overlap_length = 0

                    if ref_position >= window_end: # likely only when window end happens at the exact same position as a cigar change
                        break
                    
                    elif cigar_code == 4 or (cigar_code == 1 and ref_position < window_start): # S: Soft clipping, I: insertion to the reference
                        query_position += length
                        continue

                    elif ref_position + length <= window_start:
                        if cigar_code in [0, 2, 3, 7, 8]:  # M, D, N, =, X consume reference sequence
                            ref_position += length
                        if cigar_code in [0, 7, 8]:  # M, =, X consume query sequence
                            query_position += length
                        continue

                    elif ref_position < window_start: # start of current cigar does not overlap window, but part of it does
                        overlap_length = length - (window_start - ref_position)
                        overlap_length = window_end - window_start if overlap_length > window_end - window_start else overlap_length # in case the cigar operation extends beyond the end of the window
                        current_window_index = 0

                    else: # should always be within window at this point
                        overlap_length = window_end - ref_position if ref_position + length > window_end else length
                        current_window_index = ref_position - window_start

                    if cigar_code in [0, 7, 8]: # M, =, X
                        for current_window_offset in range(current_window_index, current_window_index + overlap_length):
                            # ref_pos = window_start + current_window_offset
                            if ref_position < window_start:
                                final_query_position = query_position + (window_start - ref_position) + current_window_offset - current_window_index
                            else:
                                final_query_position = query_position + current_window_offset - current_window_index
                            base_qual = read.query_qualities[final_query_position] if read.query_qualities else 0
                            if base_qual >= min_base_quality:
                                pileups[current_window_offset].append(PileupRead(read, final_query_position, False, False, 0))

                        ref_position += length
                        query_position += length
                    elif cigar_code == 1: # I: insertion to the reference
                        #pileups[ref_position - window_start].append(PileupRead(read, None, False, False, length))
                        query_position += length

                    elif cigar_code == 2: # D: deletion from the reference
                        for current_window_offset in range(current_window_index, current_window_index + overlap_length):
                            pileups[current_window_offset].append(PileupRead(read, None, True, False, length))
                        ref_position += length
                    elif cigar_code == 3:  # N: skipped region from the reference
                        for current_window_offset in range(current_window_index, current_window_index + overlap_length):
                            pileups[current_window_offset].append(PileupRead(read, None, False, True, 0))
                        ref_position += length                    
            
            i = 0
            for ref_pos in range(window_start, window_end):
                if pileups[i]:
                    yield PileupColumn(ref_pos, pileups[i])
                i += 1


class PairedEndBAMTrack(SingleEndBAMTrack):
    """
    Displays paired-end reads together (otherwise, same as :py:class:`genomeview.SingleEndBAMTrack`).

    Attributes:
        overlap_color: color used to highlight portions of read pairs that are overlapping one another
    """
    def __init__(self, bam_path, name=None):
        super().__init__(bam_path, name)

        self.overlap_color = "lime"

    def layout(self, scale):
        if scale == self.scale: return
        
        self.scale = scale
        self.reset_mismatch_counts()
        
        # for chrom, start, end in self.scale.regions():
        chrom, start, end = self.scale.chrom, self.scale.start, self.scale.end
        cur_read_coords = collections.defaultdict(list)

        for read in self.fetch():
            if read.is_unmapped: continue
            cur_read_coords[read.query_name].append(
                (read.reference_start, read.reference_end, read.next_reference_start, read.is_proper_pair))
        
        # a bit of hocus-pocus to deal with reads whose mates map outside of our region of interest
        for pair in cur_read_coords.values():
            if len(pair) == 1:
                read_end = pair[0]
                if read_end[3]:
                    pair.append((read_end[2], read_end[2]))
                pair.sort()

        for read_name, coords in sorted(cur_read_coords.items(), key=lambda x: x[1]):
            pair_start = coords[0][0]
            pair_end = coords[-1][1]
            interval = Interval(read_name, chrom, pair_start, pair_end)
            if self.draw_read_labels:
                interval.label = read_name
            self.layout_interval(interval)
            #self.intervals.append(interval)
                
        self.height = (len(self.rows)+1) * (self.row_height+self.margin_y)
    
    
    def draw_read_pair(self, renderer, reads):
        reads = [read for read in reads if not read.is_unmapped]
        if len(reads) == 0: return

        chrom = reads[0].reference_name
        row = self.intervals_to_rows[reads[0].query_name]
        
        pair_start = None
        if len(reads) > 1:
            pair_start = reads[0].reference_end
            pair_end = reads[-1].reference_start
        elif reads[0].is_proper_pair:
            # some more hocus-pocus to deal with reads whose mates map outside of our region of interest
            if reads[0].next_reference_start < reads[0].reference_start:
                pair_start = reads[0].next_reference_start
                pair_end = reads[0].reference_start
            else:
                pair_start = reads[0].reference_start
                pair_end = reads[0].next_reference_start

        if pair_start is not None:
            x1 = self.scale.topixels(pair_start)
            x2 = self.scale.topixels(pair_end)
            y = row*(self.row_height+self.margin_y) + self.row_height/2 # refactor

            yield from renderer.line(x1, y, x2, y, **{"stroke-width":1, "stroke":"gray"})
        
        for i, read_end in enumerate(reads):
            interval = Interval(read_end.query_name, chrom, read_end.reference_start,
                                read_end.reference_end, not read_end.is_reverse)
            
            if self.draw_read_labels:
                interval.label = "{}_{}".format(read_end.query_name, 1 if read_end.is_read1 else 2)

            interval.read = read_end

            yield from self.draw_interval(renderer, interval)

            if i == 1 and self.draw_read_labels:
                end = self.scale.topixels(read_end.reference_end)
                top = row*(self.row_height+self.margin_y)

                yield from renderer.text(end+self.label_distance, top+self.row_height,
                                         read_end.query_name, anchor="start")

            
    def render(self, renderer):
        read_buffer = {}
        for read in self.fetch():
            if read.is_unmapped: continue
            if read.query_name in read_buffer:
                other_read = read_buffer.pop(read.query_name)
                cur_reads = [other_read, read]
            else:
                read_buffer[read.query_name] = read
                continue
            
            for read_end in cur_reads:
                yield from self.draw_read_pair(renderer, cur_reads)
                
        for read_name, read in read_buffer.items():
            yield from self.draw_read_pair(renderer, [read])
        
        for x in self.render_label(renderer):
            yield x



def _get_filter_fn(keyfn, value):
    def filter_fn(read):
        return keyfn(read) == value
    return filter_fn

def get_group_by_tag_fn(tag):
    """
    creates a grouping function based on the values of "tag", to be used by GroupedBAMTrack
    for example, use tag="HP" with 10x genomics data to split the view into reads from 
    haplotype 1, haplotype 2, and those missing haplotype information
    """
    def group_by_tag(read):
        if not read.has_tag(tag):
            return "missing"
        return str(read.get_tag(tag))
    return group_by_tag

class GroupedBAMTrack(Track):
    """
    Displays reads from a BAM, separated out into groups based on a feature of the reads. 
    For example, group reads based on the value of tag.

    Attributes:
        keyfn: the function used to specify the groupings of reads. Takes as input a read 
            (:py:class:`pysam.AlignedSegment`).
        bam_track_class: the class used to display each group of reads, should probably be
            either :class:`genomeview.bamtrack.PairedEndBAMTrack` or 
            :class:`genomeview.bamtrack.SingleEndBAMTrack`

        space_between (float): the amount of space (pixels) between groups. (Default: 10)
        category_label_fn: a function that nicely formats the category labels. Takes as argument
            the result of the keyfn and should return a string. (Default: render as string)
    """
    def __init__(self, bam_path, keyfn, bam_track_class, name=None):
        """
        """
        super().__init__(name)
        self.keyfn = keyfn
        self.bam_track_class = bam_track_class
        self.bam_path = bam_path
        self.bam = pysam.AlignmentFile(bam_path)
        self.subtracks = []
        
        self.space_between = 10
        self.category_label_fn = str
        
    def layout(self, scale):
        self.scale = scale
        
        categories = set()
        chrom = match_chrom_format(self.scale.chrom, self.bam.references)
        
        for read in self.bam.fetch(chrom, self.scale.start, self.scale.end):
            category = self.keyfn(read)
            categories.add(category)
        
        categories = sorted(categories)
        self.height = 0
        for category in categories:
            cur_track = self.bam_track_class(self.bam_path, name=self.category_label_fn(category))
            cur_track.include_read_fn = _get_filter_fn(self.keyfn, category)
            cur_track.layout(scale)
            self.height += cur_track.height + self.space_between
            
            self.subtracks.append(cur_track)
            
    def render(self, renderer):
        cury = 0
        for subtrack in self.subtracks:
            subrenderer = renderer.subrenderer(y=cury, height=subtrack.height)
            yield from subrenderer.render(subtrack)
            cury += subtrack.height + self.space_between


class BAMCoverageTrack(GraphTrack):
    def __init__(self, bam_path, name=None, opener_fn=pysam.AlignmentFile):
        if name is None and isinstance(name, str):
            name = os.path.basename(bam_path.split(".")[0])
        super().__init__(name=name)
        
        self.bam_path = bam_path
        # self.bam = pysam.AlignmentFile(bam_path)
        self.opener_fn = opener_fn
        with self.opener_fn(bam_path) as bam:
            self.bam_references = bam.references
        self.include_secondary = False
        self.bin_size = 0
        self.priming_orientation = "3p"
        self.cached_series = False
        
    def layout(self, scale):
        super().layout(scale)

        if len(self.series) > 0 and self.cached_series:
            return

        chrom = match_chrom_format(scale.chrom, self.bam_references)

        if self.bin_size == 0:
            self._add_single_coverage(scale, chrom)
        else:
            self._add_binned_coverage(scale, chrom)

    def _add_single_coverage(self, scale, chrom):

        counts = collections.defaultdict(int)
        
        with self.opener_fn(self.bam_path) as bam:
            for read in bam.fetch(chrom, scale.start, scale.end):
                if read.is_secondary and not self.include_secondary:
                    continue
                for i in read.get_reference_positions():
                    counts[i] += 1
        
        x = np.arange(scale.start, scale.end+1)
        y = np.empty(scale.end - scale.start + 1, dtype=int)
        for i, curx in enumerate(x):
            y[i] = counts[curx]
        
        s = pd.Series(y, index=x).sort_index()
        s = s[(s!=s.shift(-1))|(s!=s.shift(1))]
        
        x = s.index
        y = s.values
        
        if len(x):
            self.add_series(x, y)

    def _add_binned_coverage(self, scale, chrom):
        # Initialize list of counters for each bin's start and end positions
        # +2 because +1 for included and +1 for reads that are aligned since upstream
        num_bins = ((scale.end - scale.start) // self.bin_size) + 2

        # array of coverage tracks
        # coverage = np.zeros((num_bins, scale.end - scale.start), dtype=int)
        coverage = collections.defaultdict(lambda: np.zeros(num_bins, dtype=int))

        # flag to indicate if reads start from left or right side of the figure
        is_fwd = (
            (scale.strand == "+" and self.priming_orientation == "5p")
            or (scale.strand == "-" and self.priming_orientation == "3p")
        )

        # Parse the BAM file and bin by alignment start position, then keep track of all covered blocks' start and ends
        with self.opener_fn(self.bam_path) as bam:
            for read in bam.fetch(chrom, scale.start, scale.end):
                if read.is_secondary and not self.include_secondary:
                    continue

                # get all the reference coordinates that are aligned to the read
                aligned_pairs = read.get_aligned_pairs(matches_only=True)

                # this shouldn't happen, should it?
                # could also check if any coordinates are in the region of interest
                if not aligned_pairs:
                    continue

                if is_fwd:
                    start_pos = read.reference_start
                    if start_pos < scale.start:
                        bin_index = 0
                    else:
                        bin_index = ((start_pos - scale.start) // self.bin_size) + 1
                else:
                    end_pos = read.reference_end
                    if end_pos > scale.end:
                        bin_index = 0
                    else:
                        bin_index = ((scale.end - end_pos) // self.bin_size) + 1

                for _, j in aligned_pairs:
                    if scale.start <= j < scale.end:
                        coverage[j - scale.start][bin_index] += 1

        x = np.array(sorted(coverage))
        # flag when x coordinate is changing, always keep this
        xdiff = np.diff(x, prepend=-1) > 1
        # use cumsum to calculate total coverage at all locations
        cumulative_coverage = np.vstack([np.cumsum(coverage[i]) for i in x])

        # reverse this because the tracks overlap, need shortest in front
        for bin_index in reversed(range(num_bins)):
            color = BINNED_COLORS[bin_index % len(BINNED_COLORS)]
            y = cumulative_coverage[bin_index, :]
            # calculate diff to figure out when coverage is changing
            ix = xdiff | (np.diff(y, prepend=-1) != 0)

            # only plot changing coordinates to save space
            self.add_series(x[ix], y[ix], color=color)
