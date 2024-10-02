from genomeview.utilities import my_hook_compressed, flatten
from genomeview.bam_read_operations import get_read_tag

from abc import ABC, abstractmethod

import pysam



class Classification(ABC):
    """
    Abstract Class to return an isoform classification for a given read.
    Should be inherited from and have the .get_classification(self, read, gene_id) method implemented where read is a pysam.AlignmentSegment object
    """

    @abstractmethod
    def get_classification(self, read, gene_id):
        pass
        # should return a list of classifcations (because of possible ambiguous)


class IsoQuantClassification(Classification):
    """
    Implementation of Classification class to use with a read_assignments.tsv(.gz) output from IsoQuant.
    Works by internally indexing the file provided. (may be able to use a read name bgzip index for random access)
    """

#    ISOQUANT_READ_ASSIGNMENTS_DEFS = {
#        "read_id":0,
#        "chr":1, 
#        "strand": 2,
#        "isoform_id": 3,
#        "gene_id":4, 
#        "assignment_type":5,
#        "assignment_events":6,
#        "exons":7,
#        "additional_info":8
#    }

    def __init__(self, file_path, ambiguous_classification = True):
        self.read_to_gene_id_to_isoform_id = {}
        self.read_to_assignment_type = {}
        self.ambiguous_classification = ambiguous_classification

        if type(file_path) is list:
            for fp in file_path:
                self.index(fp)
        elif type(file_path) is dict:
            for fp in file_path.values():
                self.index(fp)
        else:
            self.index(file_path)

    # call index() on file after creating the object so that the different Isoquant files can be indexed on the same object
    def index(self, file_path):
        with my_hook_compressed(file_path, "rt") as f:
            for line in f:
                if line[0] == "#":
                    continue

                fields = line.rstrip().split("\t")
                
                # values =  dict.fromkeys(self.ISOQUANT_READ_ASSIGNMENTS_DEFS)
                # for field_name, field in self.ISOQUANT_READ_ASSIGNMENTS_DEFS.items():
                #     cur_value = None
                #     if len(fields) > field:
                #         cur_value = fields[field]
                #     values[field_name] = cur_value

                # self.read_to_assignment_type[values['read_id']] = values['assignment_type']
                self.read_to_assignment_type[fields[0]] = fields[5]
        
                # if values['isoform_id'] is not None and values['isoform_id'] != ".":
                #     if values['read_id'] not in self.read_to_gene_id_to_isoform_id:
                #         self.read_to_gene_id_to_isoform_id[values['read_id']] = {}
                #     if values['gene_id'] not in self.read_to_gene_id_to_isoform_id[values['read_id']]:
                #         self.read_to_gene_id_to_isoform_id[values['read_id']][values['gene_id']] = []
                #     self.read_to_gene_id_to_isoform_id[values['read_id']][values['gene_id']].append(values['isoform_id'])
                if fields[3] is not None and fields[3] != ".":
                    if fields[0] not in self.read_to_gene_id_to_isoform_id:
                        self.read_to_gene_id_to_isoform_id[fields[0]] = {}
                    if fields[4] not in self.read_to_gene_id_to_isoform_id[fields[0]]:
                        self.read_to_gene_id_to_isoform_id[fields[0]][fields[4]] = []
                    self.read_to_gene_id_to_isoform_id[fields[0]][fields[4]].append(fields[3])

    def get_classification(self, read, gene_id):
        if read.query_name not in self.read_to_gene_id_to_isoform_id:
            return None

        if gene_id not in self.read_to_gene_id_to_isoform_id[read.query_name]:
            return ["other_gene"]

        if "unique" in self.read_to_assignment_type[read.query_name]:
            return flatten(self.read_to_gene_id_to_isoform_id[read.query_name])

        return [self.read_to_assignment_type[read.query_name]]


class BAMtagClassification(Classification):
    """
    Implementation of Classification class to use the information stored in a given BAM tag.
    """

    def __init__(self, tag):
        self.tag = tag

    def get_classification(self, read, gene_id):
        return [get_read_tag(read, self.tag)]



