"""
This module provides managing classes and methods for the tool.
"""
import matplotlib.colors
import numpy as np
import collections
import distinctipy
import subprocess
import tempfile
import typing
import seaborn
import random
import shutil
import sys
import os

import Bio.Align.AlignInfo
import Bio.Data.CodonTable
import Bio.SeqUtils
import Bio.SeqRecord
import Bio.GenBank
import Bio.AlignIO
import Bio.Align
import Bio.SeqIO
import Bio.Seq

import scipy.cluster.hierarchy
import scipy.spatial.distance

import pandas as pd

import BCBio.GFF

import lovis4u.Manager
import lovis4u.Methods
import lovis4u.Drawing


class Feature:
    """A Feature object represents a locus' feature and its properties.

    Attributes:
        feature_id (str): Feature identifier.
        feature_type (str): Type of element (e.g. CDS or tRNA).
        start (int): 1-based start genomic coordinate.
        end (int): 1-based end genomic coordinates
        strand (int): Genomic strand (1: plus strand, -1: minus strand).
        name (str): Name of the feature which will be used as a label.
        sequence (Bio.Seq.Seq): Feature's sequence.
        record (Bio.SeqRecord.SeqRecord): SeqRecord object of corresponding feature sequence.
        group (str): Feature group that defines feature's colour and meant to represent a set of homologous features.
            Can be set with Loci.mmseqs_cluster() method that uses mmseqs clustering to define CDS feature groups.
        group_type (str): Type of feature group that allow to visualise different set of feature groups differently
            (e.g. plot labels only for features which group_type is "variable" or "labeled").
        category (str): Feature category which initially is built to handle phrogs category annotation.
            In visualisation it defines the "category" colour annotation under features.
            Supposed to represent clusters on locus or any second layer of feature properties.
        vis_prms (dict): Visualisation parameters that holds colours, label and other info for Drawing methods.
        overlapping (bool): Whether feature overlaps with visualised region or not.
        prms (lovis4u.Manager.Parameters): Parameters' class object that holds config and cmd arguments.

    """

    def __init__(self, feature_id: str, feature_type: str, start: int, end: int, strand: int, name: str,
                 sequence: Bio.Seq.Seq, group: str, group_type: str, category: str, vis_prms: dict, overlapping: bool,
                 parameters: lovis4u.Manager.Parameters):
        """Create a Feature object.

        Arguments:
            feature_id (str): Feature identifier.
            feature_type (str): Type of element (e.g. CDS or tRNA). Currently only CDS are supported.
            start (int): 1-based start genomic coordinate.
            end (int): 1-based end genomic coordinates
            strand (int): Genomic strand (1: plus strand, -1: minus strand).
            name (str): Name of the feature which will be used as a label.
            sequence (Bio.Seq.Seq): Feature's sequence.
            group (str): Feature group that defines feature's colour and meant to represent a set of homologous features.
                Can be set with Loci.mmseqs_cluster() method that uses mmseqs clustering to define CDS feature groups.
            group_type (str): Type of feature group that allow to visualise different set of feature groups differently
                (e.g. plot labels only for features which group_type is "variable" or "labeled").
            category (str): Feature category which initially is built to handle phrogs category annotation.
                In visualisation it defines the "category" colour annotation under features.
                Supposed to represent clusters on locus or any second layer of feature properties.
            vis_prms (dict): Visualisation parameters that holds colours, label and other info for Drawing methods.
            prms (lovis4u.Manager.Parameters): Parameters' class object that holds config and cmd arguments.

        """
        self.feature_type = feature_type
        self.feature_id = feature_id
        self.start = start
        self.end = end
        self.strand = strand
        self.name = name
        self.sequence = sequence
        self.record = Bio.SeqRecord.SeqRecord(seq=self.sequence, id=self.feature_id)
        self.group = group  # maybe rename later
        self.group_type = group_type
        self.category = category
        self.vis_prms = vis_prms
        self.vis_prms["label"] = str(name)
        self.overlapping = overlapping
        self.prms = parameters


class Locus:
    """A Locus object represents a particular locus that will be one of the sequence tracks on final figure.

    Attributes:
        seq_id (str): Sequence identifier. Can be used to label locus.
        coordinates (list): List of regions to be shown. Each region format: dict with keys: start, end, strand and
            corresponding 1-based start, end coordinates and strand (1: plus strand, -1: minus strand).
        size (int): total length of regions to be plotted.
        description (str): Sequence description that can be used to label locus.
        length (int): full length of the locus independent on region that should be plotted.
        circular (bool): Bool value whether locus is circular or not. It defines whether you have gap or not passing
            1 value on the final figure.
        features (list): list of Feature objects that overlapped with coordinates.
        order (int): index on ordered list of loci visualisation. Can be pre-defined or found based on loci
            hierarchical clustering with cluster_sequences category.
        group (int | str): locus group that defines set of closely-related loci.
        category_colours (dict): colours for locus' features categories.
        prms (lovis4u.Manager.Parameters): Parameters' class object that holds config and cmd arguments.

    """

    def __init__(self, seq_id: str, coordinates: list, description: str, sequence: Bio.Seq.Seq, length: int,
                 circular: bool, features: list, order: int, parameters: lovis4u.Manager.Parameters,
                 group: typing.Union[int, str] = 1):
        """Create a Locus object.

        Arguments:
            seq_id (str): Sequence identifier. Can be used to label locus.
            coordinates (list): List of regions to be shown. Each region format: dict with keys: start, end, strand and
                corresponding 1-based start, end coordinates and strand (1: plus strand, -1: minus strand).
            description (str): Sequence description that can be used to label locus.
            sequence (Bio.Seq.Seq): Locus Sequence.
            length (int): full length of the locus independent on region that should be plotted.
            circular (bool): Bool value whether locus is circular or not. It defines whether you have gap or not passing
                1 value on the final figure.
            features (list): list of Feature objects that overlapped with coordinates.
            order (int): index on ordered list of loci visualisation. Can be pre-defined or found based on loci
            hierarchical clustering with cluster_sequences category.
            parameters (lovis4u.Manager.Parameters): Parameters' class object that holds config and cmd arguments.
            group (int | str): locus group that defines set of closely-related loci [1].

        """
        self.seq_id = seq_id
        self.coordinates = coordinates
        self.size = 0
        taken_coordinates = []
        for coordinate in coordinates:
            self.size += abs(coordinate["end"] - coordinate["start"] + 1)
            taken_coordinates += list(range(coordinate["start"], coordinate["end"] + 1))
        if len(taken_coordinates) != len(set(taken_coordinates)):
            raise lovis4u.Manager.lovis4uError(f"Specified coordinates seems to be overlapped"
                                               f" or not in 0-based format.")
        for coordinate in coordinates:
            if coordinate["start"] < 1:
                raise lovis4u.Manager.lovis4uError(f"Coordinates for {seq_id}: {coordinate} is in 0-based format"
                                                   f" while input should be in 1-based.")
            if coordinate["end"] > length:
                raise lovis4u.Manager.lovis4uError(f"Coordinates for {seq_id}: {coordinate} is out of sequence length"
                                                   f" ({length} nt).")

        self.description = description
        self.sequence = sequence
        self.length = length
        self.circular = circular
        self.features = features
        self.order = order
        self.group = group
        self.category_colours = dict()
        self.prms = parameters

    def calculate_sequence_property(self, window_size: int = 29, property_name: str = "gc") -> list:
        """Method to calculate gc content or gc skew

        Returns:
            list: list with corresponding values for each nucleotide

        """
        try:
            result_list = []
            locus_length = len(self.sequence)
            for i in range(locus_length):
                start = (i - window_size // 2) % locus_length
                end = (i + window_size // 2 + 1) % locus_length
                if start < end:
                    window_seq = self.sequence[start:end]
                else:
                    window_seq = self.sequence[start:] + self.sequence[:end]
                if property_name == "gc":
                    result_list.append(Bio.SeqUtils.gc_fraction(window_seq) * 100)
                elif property_name == "gc_skew":
                    result_list.append(Bio.SeqUtils.GC_skew(window_seq, len(window_seq))[0])
            return result_list
        except Exception as error:
            raise lovis4u.Manager.lovis4uError("Unable to calculate gc content for a sequence.") from error


class Loci:
    """A Loci object holds information about all loci to be plotted and methods for data preparation.

    Attributes:
        loci (list): List of Locus objects.
        locus_annotation (pd.DataFrame): Table with information about each locus that defines visualisation
            (e.g. coordinates for visualisation, description, etc).
        feature_annotation (pd.DataFrame): Table with information about each feature that defines visualisation
            (e.g. group, name, category, etc).
        prms (lovis4u.Manager.Parameters): Parameters' class object that holds config and cmd arguments.

    """

    def __init__(self, parameters: lovis4u.Manager.Parameters):
        """Create a Loci object.

        Arguments:
            parameters (lovis4u.Manager.Parameters): Parameters' class object that holds config and cmd arguments.
        """
        self.loci = []
        self.locus_annotation = pd.DataFrame(columns=["sequence_id", "length", "coordinates", "circular", "description",
                                                      "order", "group"]).set_index("sequence_id")
        self.feature_annotation = pd.DataFrame(columns=["feature_id", "locus_id", "coordinates", "feature_type", "name",
                                                        "group", "group_type", "category", "fill_colour",
                                                        "stroke_colour",
                                                        "show_label"]).set_index("feature_id")
        self.prms = parameters

    def __load_annotation_file(self, file_path: str, annotation_columns: list, index_column: str) -> pd.DataFrame:
        """Private method to load an annotation file.

        Arguments:
            file_path (str): File path for an annotation file to be loaded.
            annotation_columns (list): List of columns that should be considered.
            index_column (str): Column name to be considered as index.

        Returns:
              pd.DataFrame: Preprocessed annotation file.
        """
        annotation_table = pd.read_table(file_path)
        found_allowed_columns = [i for i in annotation_columns if i in annotation_table.columns]
        not_found_allowed_columns = [i for i in annotation_columns if i not in annotation_table.columns]
        annotation_table = annotation_table[found_allowed_columns].set_index(index_column)
        annotation_table[not_found_allowed_columns] = None
        return annotation_table

    def load_locus_annotation_file(self, file_path: str) -> None:
        """Load loci annotation file.

        Arguments:
            file_path (str): File path for a loci annotation file to be loaded.

        Returns:
            None

        """
        annotation_columns = ["sequence_id", "length", "coordinates", "circular", "description", "order", "group"]
        self.locus_annotation = self.__load_annotation_file(file_path, annotation_columns, "sequence_id")
        return None

    def load_feature_annotation_file(self, file_path: str) -> None:
        """Load features annotation file.

        Arguments:
            file_path (str): File path for a features annotation file to be loaded.

        Returns:
            None

        """
        annotation_columns = ["feature_id", "locus_id", "coordinates", "feature_type", "name", "group", "group_type",
                              "category", "fill_colour", "stroke_colour", "show_label"]
        self.feature_annotation = self.__load_annotation_file(file_path, annotation_columns, "feature_id")
        return None

    def __update_locus_annotation(self, record_id: str, record_description: str, record_length: int) -> None:
        """Private method for updating loci annotation.

        Arguments:
            record_id (str): Sequence identifier.
            record_description (str): Sequence description.
            record_length (int): Sequence length.

        Returns:
            None

        """
        if record_id not in self.locus_annotation.index:
            self.locus_annotation.loc[record_id] = {col: None for col in self.locus_annotation.columns}
        if self.prms.args["windows"]:
            coordinates = []
            for window in self.prms.args["windows"]:
                subwindows = window.strip().split(",")
                for subwindow in subwindows:
                    subwindow_split = subwindow.strip().split(":")
                    if len(subwindow_split) == 4 and subwindow_split[0] == record_id:
                        coordinate = ":".join(subwindow_split[1:])
                        coordinates.append(coordinate)
            coordinates = ",".join(coordinates)
        else:
            coordinates = f"1:{record_length}:1"
        default_values = dict(length=record_length, coordinates=coordinates,
                              description=record_description, circular=1, order=len(self.loci), group=1)
        self.locus_annotation.loc[record_id] = self.locus_annotation.loc[record_id].fillna(default_values)
        return None

    def __update_feature_annotation(self, feature_id: str, locus_id: str, coordinates: str, feature_type: str,
                                    category: str, name: str) -> None:
        """Private method for updating feature annotation.

        Arguments:
            feature_id (str): Feature identifier.
            locus_id (str): Sequence description.
            coordinates (str): Feature coordinates.
            category (str): Feature type.
            name (str): Feature name.


        Returns:
            None

        """
        if feature_id not in self.feature_annotation.index:
            self.feature_annotation.loc[feature_id] = {col: None for col in self.feature_annotation.columns}

        if self.feature_annotation.loc[feature_id]["group_type"] in self.prms.args["feature_group_types_to_show_label"] \
                or "all" in self.prms.args["feature_group_types_to_show_label"]:
            show_label = 1
        else:
            show_label = 0
        stroke_colour = "default"
        if self.prms.args["set_feature_stroke_colour_based_on_fill_colour"] and \
                self.feature_annotation.loc[feature_id]["fill_colour"] and \
                self.feature_annotation.loc[feature_id]["fill_colour"] != "default":
            stroke_colour = lovis4u.Methods.scale_lightness(self.feature_annotation.loc[feature_id]["fill_colour"],
                                                            self.prms.args["feature_stroke_colour_relative_lightness"])
        default_values = dict(locus_id=locus_id, coordinates=coordinates, feature_type=feature_type,
                              name=name, group="", group_type="undefined", category=category, fill_colour="default",
                              stroke_colour=stroke_colour,
                              show_label=show_label)
        self.feature_annotation.loc[feature_id] = self.feature_annotation.loc[feature_id].fillna(default_values)
        return None

    def load_loci_from_extended_gff(self, input_f: str, ilund4u_mode: bool = False) -> None:
        """Load loci from the folder with gff files. Each GFF file also should contain corresponding nucleotide
            sequence. Such files are produced for example by pharokka annotation tool.

        All files with extension other than .gff (not case-sensitive) will be ignored.

        Arguments:
            input_folder: folder name with gff files.

        Returns:
            None

        """

        try:
            if isinstance(input_f, str):
                input_folder = input_f
                if not os.path.exists(input_folder):
                    raise lovis4u.Manager.lovis4uError(f"Folder/file {input_folder} does not exist.")
                if os.path.isfile(input_f):
                    gff_files = [input_f]
                else:
                    gff_files = [os.path.join(input_folder, f) for f in os.listdir(input_folder)]
            elif isinstance(input_f, list):
                gff_files = input_f
            else:
                raise lovis4u.Manager.lovis4uError(f"The input for the GFF parsing function must be either a folder or "
                                                   f"a list of files.")
            if not gff_files:
                raise lovis4u.Manager.lovis4uError(f"Folder {input_f} does not contain files.")
            if self.prms.args["verbose"]:
                print(f"○ Reading gff file{'s' if len(gff_files) > 1 else ''}...", file=sys.stdout)
            for gff_file in gff_files:
                try:
                    gff_records = list(BCBio.GFF.parse(gff_file, limit_info=dict(gff_type=["CDS", "tRNA", "tmRNA",
                                                                                           "RNA",
                                                                                           "pseudogene", "rRNA",
                                                                                           "misc_RNA",
                                                                                           "ncRNA", "lncRNA"])))
                    if len(gff_records) != 1:
                        print(f"○ Warning: gff file {gff_file} contains information for more than 1 "
                              f"sequence. File will be skipped.")
                        continue
                    gff_record = gff_records[0]
                    try:
                        record_locus_sequence = gff_record.seq
                    except Bio.Seq.UndefinedSequenceError:
                        print(f"○ Warning: gff file {gff_file} doesn't contain corresponding sequences.")
                        continue
                    if self.prms.args["gff_description_source"] in gff_record.annotations:
                        record_description = gff_record.annotations[self.prms.args["gff_description_source"]][0]
                        if isinstance(record_description, tuple):
                            record_description = " ".join(record_description)
                    else:
                        record_description = ""
                    if self.prms.args["use_filename_as_contig_id"]:
                        gff_record.id = os.path.splitext(os.path.basename(gff_file))[0]
                    self.__update_locus_annotation(gff_record.id, record_description, len(record_locus_sequence))
                    locus_annotation_row = self.locus_annotation.loc[gff_record.id]
                    coordinates = [dict(zip(["start", "end", "strand"], map(int, c.split(":")))) for c in
                                   locus_annotation_row["coordinates"].split(",")]
                    record_locus = Locus(seq_id=gff_record.id, coordinates=coordinates,
                                         description=locus_annotation_row["description"],
                                         sequence=Bio.SeqRecord.SeqRecord(seq=record_locus_sequence, id=gff_record.id),
                                         circular=locus_annotation_row["circular"],
                                         length=locus_annotation_row["length"], parameters=self.prms, features=[],
                                         order=locus_annotation_row["order"])
                    features_ids = [i.id for i in gff_record.features]
                    if len(features_ids) != len(set(features_ids)):
                        raise lovis4u.Manager.lovis4uError(f"Gff file {gff_file} contains duplicated feature ids while"
                                                           f" only unique are allowed.")
                    for gff_feature in gff_record.features:
                        feature_id = gff_feature.id
                        if ilund4u_mode or self.prms.args["add_locus_id_prefix"]:
                            if gff_record.id not in feature_id:
                                feature_id = f"{gff_record.id}-{feature_id}"
                        if gff_feature.type == "CDS":
                            transl_table = self.prms.args["default_transl_table"]
                            if "transl_table" in gff_feature.qualifiers.keys():
                                transl_table = int(gff_feature.qualifiers["transl_table"][0])
                            name = ""
                            if self.prms.args["gff_CDS_name_source"] in gff_feature.qualifiers:
                                name = gff_feature.qualifiers[self.prms.args["gff_CDS_name_source"]][0]
                            category = ""
                            if self.prms.args["gff_CDS_category_source"] in gff_feature.qualifiers:
                                category = ",".join(gff_feature.qualifiers[self.prms.args["gff_CDS_category_source"]])
                            sequence = gff_feature.translate(record_locus_sequence, table=transl_table,
                                                             cds=False)[:-1]
                        else:
                            name, category = "", ""
                            if self.prms.args["gff_noncoding_name_source"] in gff_feature.qualifiers:
                                name = self.prms.args["gff_noncoding_name_source"][0]
                            else:
                                for ans in self.prms.args["gff_noncoding_name_alternative_source"]:
                                    if ans in gff_feature.qualifiers:
                                        name = gff_feature.qualifiers[ans][0]
                                        break
                            if not name:
                                name = feature_id
                            category = gff_feature.type
                            sequence = gff_feature.extract(record_locus_sequence)
                        for coordinate in record_locus.coordinates:
                            overlapping = False
                            start, end = coordinate["start"], coordinate["end"]
                            if start <= gff_feature.location.start + 1 <= end or start <= gff_feature.location.end <= end:
                                overlapping = True
                                break
                        if not overlapping and not self.prms.args["cluster_all_proteins"]:
                            continue
                        self.__update_feature_annotation(feature_id, record_locus.seq_id,
                                                         f"{int(gff_feature.location.start) + 1}:"
                                                         f"{int(gff_feature.location.end)}:{gff_feature.location.strand}",
                                                         gff_feature.type, category, name)
                        feature_annotation_row = self.feature_annotation.loc[feature_id]
                        feature = Feature(feature_type=feature_annotation_row["feature_type"],
                                          feature_id=feature_id, start=int(gff_feature.location.start) + 1,
                                          end=int(gff_feature.location.end), strand=gff_feature.location.strand,
                                          name=feature_annotation_row["name"],
                                          sequence=sequence,
                                          group=feature_annotation_row["group"],
                                          group_type=feature_annotation_row["group_type"],
                                          category=feature_annotation_row["category"],
                                          vis_prms=dict(fill_colour=feature_annotation_row["fill_colour"],
                                                        stroke_colour=feature_annotation_row["stroke_colour"],
                                                        show_label=feature_annotation_row["show_label"],
                                                        hmmscan_hit=0),
                                          overlapping=overlapping, parameters=self.prms)
                        record_locus.features.append(feature)
                    self.loci.append(record_locus)
                except:
                    print(f"○ Warning: gff file {gff_file} was not read properly and skipped")
                    if self.prms.args["parsing_debug"]:
                        self.prms.args["debug"] = True
                        raise lovis4u.Manager.lovis4uError()
            seq_id_to_order = self.locus_annotation["order"].to_dict()
            loci_ids = [l.seq_id for l in self.loci]
            if len(loci_ids) != len(set(loci_ids)):
                raise lovis4u.Manager.lovis4uError(f"The input gff files have duplicated contig ids. "
                                                   f"You can use `--use-filename-as-id` parameter to use file name "
                                                   f"as contig id which can help to fix the problem.")
            feature_ids = [feature.feature_id for locus in self.loci for feature in locus.features]
            if len(feature_ids) != len(set(feature_ids)):
                raise lovis4u.Manager.lovis4uError(f"The input gff files have duplicated features ids. "
                                                   f"You can use `--add-locus-id-prefix` (`-alip`) parameter to add "
                                                   f"contig id prefix to each feature id which can help to "
                                                   f"fix the problem.")
            self.loci.sort(key=lambda locus: seq_id_to_order[locus.seq_id])
            if self.prms.args["verbose"]:
                print(f"⦿ {len(self.loci)} {'locus was' if len(self.loci) == 1 else 'loci were'} loaded from the gff "
                      f"files folder", file=sys.stdout)
            return None
        except Exception as error:
            raise lovis4u.Manager.lovis4uError("Unable to load loci from gff folder.") from error

    def load_loci_from_gb(self, input_folder: str) -> None:
        """Load loci from the folder with genbank files.

        All files with extension other than .gb (not case-sensitive) will be ignored.

        Arguments:
            input_folder: folder name with gb files.

        Returns:
            None

        """
        if not os.path.exists(input_folder):
            raise lovis4u.Manager.lovis4uError(f"Folder/file {input_folder} does not exist.")
        try:
            if os.path.isfile(input_folder):
                gb_files = [input_folder]
            else:
                gb_files = [os.path.join(input_folder, f) for f in os.listdir(input_folder)]
            if not gb_files:
                raise lovis4u.Manager.lovis4uError(f"Folder {input_folder} does not contain files.")
            if self.prms.args["verbose"]:
                print(f"○ Reading gb file{'s' if len(gb_files) > 1 else ''}...", file=sys.stdout)
            for gb_file in gb_files:
                try:
                    gb_records = list(Bio.SeqIO.parse(gb_file, "genbank"))
                    if len(gb_records) != 1:
                        print(f"○ Warning: gb file {gb_file} contains information for more than 1 "
                              f"sequence. File will be skipped.")
                        continue
                    gb_record = gb_records[0]
                    record_locus_sequence = gb_record.seq
                    if self.prms.args["genbank_description_source"] == "description":
                        record_description = gb_record.description
                    elif "annotations:" in self.prms.args["genbank_description_source"]:
                        feature_description_key = self.prms.args["genbank_description_source"].split(":")[1]
                        record_description = gb_record.annotations[feature_description_key]
                    else:
                        record_description = ""
                    if self.prms.args["use_filename_as_contig_id"]:
                        gb_record.id = os.path.splitext(os.path.basename(gb_file))[0]
                    self.__update_locus_annotation(gb_record.id, record_description, len(record_locus_sequence))
                    locus_annotation_row = self.locus_annotation.loc[gb_record.id]
                    coordinates = [dict(zip(["start", "end", "strand"], map(int, c.split(":")))) for c in
                                   locus_annotation_row["coordinates"].split(",")]
                    record_locus = Locus(seq_id=gb_record.id, coordinates=coordinates,
                                         description=locus_annotation_row["description"],
                                         sequence=Bio.SeqRecord.SeqRecord(seq=record_locus_sequence, id=gb_record.id),
                                         circular=locus_annotation_row["circular"],
                                         length=locus_annotation_row["length"], parameters=self.prms, features=[],
                                         order=locus_annotation_row["order"])

                    gb_CDSs = [i for i in gb_record.features if i.type == "CDS"]
                    first_CDS_record = gb_CDSs[0]
                    id_source = self.prms.args["genbank_id_source"]
                    if self.prms.args["genbank_id_source"] not in first_CDS_record.qualifiers:
                        for alternative_id_source in self.prms.args["genbank_id_alternative_source"]:
                            if alternative_id_source in first_CDS_record.qualifiers:
                                id_source = alternative_id_source
                                if self.prms.args["verbose"]:
                                    print(f"○ Warning: there is no <{self.prms.args['genbank_id_source']}> attribute "
                                          f"for CDS records in {gb_file}. Alternative <{id_source}> was used instead.",
                                          file=sys.stdout)
                                break
                        if id_source == self.prms.args["genbank_id_source"]:
                            print(f"There is no <{self.prms.args['genbank_id_source']}> "
                                  f"attribute for CDS record found in {gb_file}. We tried to"
                                  f" find any from the alternative list: "
                                  f"{','.join(self.prms.args['genbank_id_alternative_source'])}"
                                  f", but they also weren't found.")  # add about cmd parameter
                    features_ids = [i.qualifiers[id_source][0] for i in gb_CDSs if id_source in i.qualifiers]
                    if len(features_ids) != len(set(features_ids)):
                        print(f"GB file {gb_record} contains duplicated feature ids while"
                              f" only unique are allowed.")
                    for gb_feature in gb_record.features:
                        if id_source not in gb_feature.qualifiers:
                            print(f"    ○ Warning: genbank feature for {gb_file} located at {gb_feature.location} "
                                  f"was skipped\n    since it does not have id qualifier {id_source} found for other "
                                  f"features.\n    it could be a case of zero length ORF. ", file=sys.stdout)
                            continue
                        feature_id = gb_feature.qualifiers[id_source][0].replace("|", "_")
                        if self.prms.args["add_locus_id_prefix"]:
                            if gb_record.id not in feature_id:
                                feature_id = f"{gb_record.id}-{feature_id}"
                        transl_table = self.prms.args["default_transl_table"]
                        if "transl_table" in gb_feature.qualifiers.keys():
                            transl_table = int(gb_feature.qualifiers["transl_table"][0])
                        name = ""
                        if gb_feature.type == "CDS":
                            if self.prms.args["genbank_CDS_name_source"] in gb_feature.qualifiers:
                                name = gb_feature.qualifiers[self.prms.args["genbank_CDS_name_source"]][0]
                            category = ""
                            if self.prms.args["genbank_CDS_category_source"] in gb_feature.qualifiers:
                                category = ",".join(
                                    gb_feature.qualifiers[self.prms.args["genbank_CDS_category_source"]])
                            sequence = gb_feature.translate(record_locus_sequence, table=transl_table, cds=False)[:-1]
                        else:
                            name, category = "", ""
                            if self.prms.args["gff_noncoding_name_source"] in gb_feature.qualifiers:
                                name = self.prms.args["gff_noncoding_name_source"][0]
                            else:
                                for ans in self.prms.args["gff_noncoding_name_alternative_source"]:
                                    if ans in gb_feature.qualifiers:
                                        name = gb_feature.qualifiers[ans][0]
                                        break
                            if not name:
                                name = feature_id
                            category = gb_feature.type
                            sequence = gb_feature.extract(record_locus_sequence)
                        for coordinate in record_locus.coordinates:
                            overlapping = False
                            start, end = coordinate["start"], coordinate["end"]
                            if start <= gb_feature.location.start + 1 <= end or start <= gb_feature.location.end <= end:
                                overlapping = True
                                break
                        if not overlapping and not self.prms.args["cluster_all_proteins"]:
                            continue
                        self.__update_feature_annotation(feature_id, record_locus.seq_id,
                                                         f"{int(gb_feature.location.start) + 1}:"
                                                         f"{int(gb_feature.location.end)}:"
                                                         f"{gb_feature.location.strand}", gb_feature.type, category,
                                                         name)
                        feature_annotation_row = self.feature_annotation.loc[feature_id]
                        feature = Feature(feature_type=feature_annotation_row["feature_type"],
                                          feature_id=feature_id, start=int(gb_feature.location.start) + 1,
                                          end=int(gb_feature.location.end),
                                          strand=gb_feature.location.strand,
                                          name=feature_annotation_row["name"],
                                          sequence=sequence,
                                          group=feature_annotation_row["group"],
                                          group_type=feature_annotation_row["group_type"],
                                          category=feature_annotation_row["category"],
                                          vis_prms=dict(fill_colour=feature_annotation_row["fill_colour"],
                                                        stroke_colour=feature_annotation_row["stroke_colour"],
                                                        show_label=feature_annotation_row["show_label"],
                                                        hmmscan_hit=0),
                                          overlapping=overlapping,
                                          parameters=self.prms)

                        record_locus.features.append(feature)
                    self.loci.append(record_locus)
                except:
                    print(f"○ Warning: gb file {gb_file} was not read properly and skipped")
                    if self.prms.args["parsing_debug"]:
                        self.prms.args["debug"] = True
                        raise lovis4u.Manager.lovis4uError()
            seq_id_to_order = self.locus_annotation["order"].to_dict()
            loci_ids = [l.seq_id for l in self.loci]
            if len(loci_ids) != len(set(loci_ids)):
                raise lovis4u.Manager.lovis4uError(f"The input gb files have duplicated contig ids. "
                                                   f"You can use `--use-filename-as-id` parameter to use file name "
                                                   f"as contig id which can help to fix the problem.")
            feature_ids = [feature.feature_id for locus in self.loci for feature in locus.features]
            if len(feature_ids) != len(set(feature_ids)):
                raise lovis4u.Manager.lovis4uError(f"The input gb files have duplicated features ids. "
                                                   f"You can use `--add-locus-id-prefix` (`-alip`) parameter to add "
                                                   f"contig id prefix to each feature id which can help to "
                                                   f"fix the problem.")
            self.loci.sort(key=lambda locus: seq_id_to_order[locus.seq_id])
            if self.prms.args["verbose"]:
                print(f"⦿ {len(self.loci)} {'locus was' if len(self.loci) == 1 else 'loci were'} loaded from the "
                      f"genbank files folder", file=sys.stdout)
            return None
        except Exception as error:
            raise lovis4u.Manager.lovis4uError("Unable to load loci from gb folder.") from error

    def save_locus_annotation_table(self) -> None:
        """Save loci annotation table to the output folder.

        Output file name is locus_annotation_table.tsv

        Returns:
            None

        """
        try:
            if not os.path.exists(self.prms.args["output_dir"]):
                os.mkdir(self.prms.args["output_dir"])
            locus_ids = [locus.seq_id for locus in self.loci]
            self.locus_annotation = self.locus_annotation.loc[self.locus_annotation.index.isin(locus_ids)]
            file_path = os.path.join(self.prms.args["output_dir"], "locus_annotation_table.tsv")
            self.locus_annotation.to_csv(file_path, sep="\t", index_label="sequence_id")
            if self.prms.args["verbose"]:
                print(f"⦿ Loci annotation table was saved to {file_path}", file=sys.stdout)
            return None
        except Exception as error:
            raise lovis4u.Manager.lovis4uError("Unable to save loci annotation table.") from error

    def save_feature_annotation_table(self) -> None:
        """Save feature annotation table to the output folder.

        Output file name is feature_annotation_table.tsv

        Returns:
            None

        """
        try:
            if not os.path.exists(self.prms.args["output_dir"]):
                os.mkdir(self.prms.args["output_dir"])
            feature_ids = [feature.feature_id for locus in self.loci for feature in locus.features]
            self.feature_annotation = self.feature_annotation.loc[self.feature_annotation.index.isin(feature_ids)]
            file_path = os.path.join(self.prms.args["output_dir"], "feature_annotation_table.tsv")
            self.feature_annotation.to_csv(file_path, sep="\t", index_label="feature_id")
            if self.prms.args["verbose"]:
                print(f"⦿ Feature annotation table was saved to {file_path}", file=sys.stdout)
            return None
        except Exception as error:
            raise lovis4u.Manager.lovis4uError("Unable to save feature annotation table.") from error

    def mmseqs_cluster(self) -> pd.DataFrame:
        """Cluster all proteins using mmseqs in order to define groups of homologues.

        Returns:
            pd.DataFrame: parsed mmseqs table (pandas dataframe) with columns: cluster, protein_id; where cluster is
                defined by representative sequence id within a corresponding cluster.

        """
        if self.prms.args["verbose"]:
            print(f"○ Running mmseqs for protein clustering...", file=sys.stdout)
        try:
            feature_records = [feature.record for locus in self.loci for feature in locus.features if
                               feature.feature_type == "CDS"]
            temp_input = tempfile.NamedTemporaryFile()
            Bio.SeqIO.write(feature_records, temp_input.name, "fasta")
            if not os.path.exists(self.prms.args["output_dir"]):
                os.mkdir(self.prms.args["output_dir"])
            mmseqs_output_folder = os.path.join(self.prms.args["output_dir"], "mmseqs")
            if os.path.exists(mmseqs_output_folder):
                shutil.rmtree(mmseqs_output_folder)
            os.mkdir(mmseqs_output_folder)
            Bio.SeqIO.write(feature_records, os.path.join(mmseqs_output_folder, "input_proteins.fa"), "fasta")
            mmseqs_output_folder_db = os.path.join(mmseqs_output_folder, "DB")
            os.mkdir(mmseqs_output_folder_db)
            mmseqs_stdout = open(os.path.join(mmseqs_output_folder, "mmseqs_stdout.txt"), "w")
            mmseqs_stderr = open(os.path.join(mmseqs_output_folder, "mmseqs_stderr.txt"), "w")
            subprocess.run([self.prms.args["mmseqs_binary"], "createdb", temp_input.name,
                            os.path.join(mmseqs_output_folder_db, "sequencesDB")], stdout=mmseqs_stdout,
                           stderr=mmseqs_stderr)
            subprocess.run([self.prms.args["mmseqs_binary"], "cluster",
                            os.path.join(mmseqs_output_folder_db, "sequencesDB"),
                            os.path.join(mmseqs_output_folder_db, "clusterDB"),
                            os.path.join(mmseqs_output_folder_db, "tmp"),
                            "--cluster-mode", str(self.prms.args["mmseqs_cluster_mode"]),
                            "--cov-mode", str(self.prms.args["mmseqs_cov_mode"]),
                            "--min-seq-id", str(self.prms.args["mmseqs_min_seq_id"]),
                            "-c", str(self.prms.args["mmseqs_c"]),
                            "-s", str(self.prms.args["mmseqs_s"])], stdout=mmseqs_stdout, stderr=mmseqs_stderr)
            subprocess.run([self.prms.args["mmseqs_binary"], "createtsv",
                            os.path.join(mmseqs_output_folder_db, "sequencesDB"),
                            os.path.join(mmseqs_output_folder_db, "sequencesDB"),
                            os.path.join(mmseqs_output_folder_db, "clusterDB"),
                            os.path.join(mmseqs_output_folder, "mmseqs_clustering.tsv")],
                           stdout=mmseqs_stdout, stderr=mmseqs_stderr)
            mmseqs_clustering_results = pd.read_table(os.path.join(mmseqs_output_folder, "mmseqs_clustering.tsv"),
                                                      sep="\t", header=None, names=["cluster", "protein_id"])
            mmseqs_clustering_results = mmseqs_clustering_results.set_index("protein_id")

            num_of_unique_clusters = len(set(mmseqs_clustering_results["cluster"].to_list()))
            num_of_proteins = len(mmseqs_clustering_results.index)
            if self.prms.args["verbose"]:
                print(f"⦿ {num_of_unique_clusters} clusters for {num_of_proteins} proteins were found with mmseqs\n"
                      f"\tmmseqs clustering results were saved to "
                      f"{os.path.join(mmseqs_output_folder, 'mmseqs_clustering.tsv')}", file=sys.stdout)
            return mmseqs_clustering_results
        except Exception as error:
            raise lovis4u.Manager.lovis4uError("Unable to run mmseqs clustering. In case you use a linux machine"
                                               ",  have you run a post-install 'lovis4u --linux` command to switch to"
                                               " the linux mmseqs binary?") from error

    def define_feature_groups(self, dataframe: pd.DataFrame, group_column_name: str = "cluster") -> None:
        """Set features attribute "group" based on input dataframe.

        By default it is designed to use mmseqs_cluster() function results as input. If you already have precomputed
            feature groups you can set them with feature table.

        Arguments:
            dataframe (pd.DataFrame): dataframe with feature id - group pairs. Its index column should
                represent all loci CDS features.
            group_column_name (str): column name of the dataframe that represent corresponding group to each feature.

        Returns:
            None

        """
        try:

            for locus in self.loci:
                for feature in locus.features:
                    if feature.group and self.prms.args["keep_predefined_groups"]:
                        continue
                    if feature.feature_type == "CDS":
                        feature.group = dataframe.loc[feature.feature_id, group_column_name]
                    else:
                        feature.group = feature.name
                    self.feature_annotation.loc[feature.feature_id, "group"] = feature.group
            return None
        except Exception as error:
            raise lovis4u.Manager.lovis4uError("Unable to define protein features groups.") from error

    def remove_non_overlapping_features(self) -> None:
        """Removes features that are not overlapping with visualisation window.

        Returns:
            None
        """
        try:
            ids_of_non_overlapping_objects = []
            for locus in self.loci:
                ids_of_non_overlapping_objects += [obj.feature_id for obj in locus.features if not obj.overlapping]
                filtered_objects = [obj for obj in locus.features if obj.overlapping]
                locus.features = filtered_objects
            if ids_of_non_overlapping_objects:
                print("○ Warning message: LoVis4u clusters all proteins by default to define their classes"
                      "\n\t('variable' or 'conserved'), including those outside the visualisation window."
                      "\n\tTo cluster only proteins within the visualised area, use the -cl-owp parameter.")
            self.feature_annotation = self.feature_annotation.drop(ids_of_non_overlapping_objects)
        except Exception as error:
            raise lovis4u.Manager.lovis4uError("Unable to clean non overlapping features.") from error

    def define_labels_to_be_shown(self):
        """Set feature visaulisation attribute "show_label" based on feature groups.

        controlled by feature_labels_to_ignore, feature_group_types_to_show_label, and
            feature_group_types_to_show_label_on_first_occurrence parameters.

        Returns:
            None

        """
        try:
            added_first_occurrence_labels = []
            for locus in self.loci:
                for feature in locus.features:
                    if self.prms.args["show_all_feature_labels"]:
                        feature.vis_prms["show_label"] = 1
                        continue
                    if feature.vis_prms["label"] not in self.prms.args["feature_labels_to_ignore"]:
                        if "any" in self.prms.args["feature_group_types_to_show_label"]:
                            feature.vis_prms["show_label"] = 1
                        elif feature.group_type in self.prms.args["feature_group_types_to_show_label"] or \
                                (feature.vis_prms["hmmscan_hit"] and self.prms.args["show_all_label_for_query_"
                                                                                    "proteins"]):
                            feature.vis_prms["show_label"] = 1
                        elif (feature.group_type in \
                              self.prms.args["feature_group_types_to_show_label_on_first_occurrence"]) or \
                                (feature.vis_prms["hmmscan_hit"] and self.prms.args["show_label_on_first_occurrence"
                                                                                    "_for_query_proteins"]):
                            if feature.group not in added_first_occurrence_labels or len(self.loci) == 1:
                                feature.vis_prms["show_label"] = 1
                                added_first_occurrence_labels.append(feature.group)
                    else:
                        feature.vis_prms["show_label"] = 0
        except Exception as error:
            raise lovis4u.Manager.lovis4uError("Unable define feature labels to be shown.") from error

    def cluster_sequences(self, dataframe: pd.DataFrame, one_cluster: bool) -> None:
        """Define loci order and clusters with proteome similarity based hierarchical clustering.
            This function changes the order of loci that are plotted and also updates corresponding to each loci group
            attribute which defines homologues groups of proteomes.

        It's designed to use as input mmseqs_cluster() function results. However, if you have obtained homologues
            groups by other method you can also build pandas dataframe based on that with index corresponding to
            feature id and column "cluster" corresponding to the group.

        Arguments:
              dataframe (pd.DataFrame): dataframe with feature id - group pairs. Its index column should
                represent all loci CDS features.
                one_cluster (bool): consider all sequences to be members of one cluster, but still define the
                optimal order.

        Returns:
            None

        """
        try:
            proteins_loci_dict = collections.defaultdict(collections.deque)
            loci_clusters_dict = dict()
            number_of_loci = len(self.loci)
            if number_of_loci < 2:
                return None
            proteome_sizes = pd.Series(np.zeros(number_of_loci, dtype=int))
            for locus_index in range(number_of_loci):
                locus = self.loci[locus_index]
                loci_clusters = [dataframe.loc[feature.feature_id, "cluster"] for feature in locus.features if
                                 feature.feature_type == "CDS"]
                loci_clusters_dict[locus_index] = list(set(loci_clusters))
                proteome_sizes.iloc[locus_index] = len(set(loci_clusters))
                for l_cl in loci_clusters:
                    proteins_loci_dict[l_cl].append(locus_index)

            loci_ids = [locus.seq_id for locus in self.loci]
            similarity_matrix = pd.DataFrame(0.0, index=loci_ids, columns=loci_ids)
            for locus_index in range(number_of_loci):
                counts = pd.Series(np.zeros(number_of_loci, dtype=int))
                for cluster in loci_clusters_dict[locus_index]:
                    js = proteins_loci_dict[cluster]
                    counts.iloc[js] += 1
                locus_size = proteome_sizes[locus_index]
                norm_factors = pd.Series(0.5 * (locus_size + proteome_sizes) / (locus_size * proteome_sizes),
                                         index=counts.index)
                weights = counts.mul(norm_factors)
                similarity_matrix.iloc[locus_index] = weights
            symmetric_distance_matrix = 1 - similarity_matrix
            np.fill_diagonal(symmetric_distance_matrix.values, 0)
            linkage_matrix = scipy.cluster.hierarchy.linkage(
                scipy.spatial.distance.squareform(symmetric_distance_matrix),
                method="average")
            dendrogram = scipy.cluster.hierarchy.dendrogram(linkage_matrix, no_plot=True)
            if not one_cluster:
                clusters = pd.Series(scipy.cluster.hierarchy.fcluster(linkage_matrix,
                                                                      self.prms.args["clustering_h_value"],
                                                                      criterion="distance"),
                                     index=loci_ids)
                for locus in self.loci:
                    locus.group = clusters[locus.seq_id]
                    self.locus_annotation.loc[locus.seq_id, "group"] = locus.group
            order = dendrogram["leaves"][::-1]
            self.locus_annotation["initial_order"] = self.locus_annotation["order"]
            for locus_index in range(number_of_loci):
                locus = self.loci[locus_index]
                self.locus_annotation.loc[locus.seq_id, "order"] = order.index(locus_index)
            self.locus_annotation.sort_values(by="order", inplace=True)
            seq_id_to_order = self.locus_annotation["order"].to_dict()
            self.loci.sort(key=lambda locus: seq_id_to_order[locus.seq_id])

            reordered_similarity_matrix = similarity_matrix.reindex(index=self.locus_annotation.index,
                                                                    columns=self.locus_annotation.index)
            if not os.path.exists(self.prms.args["output_dir"]):
                os.mkdir(self.prms.args["output_dir"])
            file_path = os.path.join(self.prms.args["output_dir"], "proteome_similarity_matrix.tsv")
            reordered_similarity_matrix.to_csv(file_path, sep="\t")
            num_of_loci_groups = len(set(self.locus_annotation["group"].to_list()))
            if self.prms.args["verbose"]:
                if num_of_loci_groups == 1:
                    print(f"⦿ Loci order and {num_of_loci_groups} cluster was defined with proteome similarity based "
                          f"hierarchical clustering", file=sys.stdout)
                elif num_of_loci_groups > 1:
                    print(f"⦿ Loci order and {num_of_loci_groups} clusters were defined with proteome similarity based "
                          f"hierarchical clustering", file=sys.stdout)
                print(f"⦿ Proteome similarity matrix of loci was saved to {file_path}", file=sys.stdout)
            return None
        except Exception as error:
            raise lovis4u.Manager.lovis4uError("Unable to cluster loci sequences.") from error

    def find_variable_feature_groups(self, mmseqs_results: pd.DataFrame) -> None:
        """Define feature group type attributes (variable or conserved) based on their conservation in corresponding
            loci group feature.

        It's designed to use as input mmseqs_cluster() function results. However, if you have obtained homologues
            groups by other method you can also build pandas dataframe based on that with index corresponding to
            feature id and column "cluster" corresponding to the group.

        Arguments:
              mmseqs_results (pd.DataFrame): dataframe with feature id - group pairs. Its index column should
                represent all loci CDS features.

        Returns:
            None

        """
        try:
            loci_clusters_sizes = self.locus_annotation["group"].value_counts()
            loci_clusters_cutoff_v = np.round(self.prms.args["CDS_is_variable_cutoff"] * loci_clusters_sizes).astype(
                int)
            loci_clusters_cutoff_c = np.round(self.prms.args["CDF_is_conserved_cutoff"] * loci_clusters_sizes).astype(
                int)
            loci_clusters_cutoff_v[loci_clusters_cutoff_v == 0] = 1
            cluster_types = collections.defaultdict(dict)
            for cluster in set(mmseqs_results["cluster"].to_list()):
                cluster_proteins = mmseqs_results[mmseqs_results["cluster"] == cluster].index
                cluster_loci = [locus for locus in self.loci if
                                any(feature.feature_id in cluster_proteins for feature in locus.features)]
                cluster_loci_groups = [locus.group for locus in cluster_loci]
                for cluster_locus_group in cluster_loci_groups:
                    current_group_cluster_loci = [locus.seq_id for locus in cluster_loci if
                                                  locus.group == cluster_locus_group]
                    current_group_cluster_size = len(set(current_group_cluster_loci))
                    if loci_clusters_sizes[cluster_locus_group] > 1 and \
                            current_group_cluster_size <= loci_clusters_cutoff_v[cluster_locus_group]:
                        cluster_types[cluster_locus_group][cluster] = "variable"
                    elif loci_clusters_sizes[cluster_locus_group] > 1 and \
                            (loci_clusters_cutoff_v[cluster_locus_group] < current_group_cluster_size <
                             loci_clusters_cutoff_c[cluster_locus_group]):
                        cluster_types[cluster_locus_group][cluster] = "intermediate"
                    else:
                        cluster_types[cluster_locus_group][cluster] = "conserved"
            for locus in self.loci:
                locus_group = locus.group
                for feature in locus.features:
                    if feature.feature_type == "CDS":
                        if feature.group_type and self.prms.args["keep_predefined_groups"]:
                            continue
                        feature.group_type = cluster_types[locus_group][feature.group]
                    else:
                        feature.group_type = "noncoding"
                    self.feature_annotation.loc[feature.feature_id, "group_type"] = feature.group_type
            return None
        except Exception as error:
            raise lovis4u.Manager.lovis4uError("Unable to define variable feature groups.") from error

    def set_feature_colours_based_on_groups(self) -> None:
        """Define features fill colour based on corresponding feature group and group types.

        Returns:
            None

        """
        try:
            feature_groups = set([feature.group for locus in self.loci for feature in locus.features
                                  if (feature.group and feature.feature_type == "CDS")])
            if self.prms.args["feature_group_types_to_set_colour"] and \
                    "all" not in self.prms.args["feature_group_types_to_set_colour"]:
                feature_groups = set([feature.group for locus in self.loci for feature in locus.features
                                      if feature.group and feature.feature_type == "CDS" and feature.group_type in
                                      self.prms.args["feature_group_types_to_set_colour"]])
            number_of_unique_feature_groups = len(feature_groups)
            if self.prms.args["groups_fill_colour_palette_lib"] == "seaborn":
                colours_rgb = seaborn.color_palette(self.prms.args["groups_fill_colour_seaborn_palette"],
                                                    number_of_unique_feature_groups,
                                                    desat=self.prms.args["groups_fill_colour_seaborn_desat"])
                random.shuffle(colours_rgb)
            elif self.prms.args["groups_fill_colour_palette_lib"] == "distinctipy":
                colours_rgb = distinctipy.get_colors(number_of_unique_feature_groups,
                                                     exclude_colours=[(1, 1, 1), (0, 0, 0)],
                                                     pastel_factor=self.prms.args["groups_fill_colours_pastel_factor"])
            colours = list(map(lambda x: matplotlib.colors.rgb2hex(x), colours_rgb))
            colours_dict = {g: c for g, c in zip(list(feature_groups), colours)}
            for locus in self.loci:
                for feature in locus.features:
                    if feature.group in feature_groups and feature.feature_type == "CDS":
                        if self.prms.args["keep_predefined_colours"] and feature.vis_prms["fill_colour"] != "default":
                            continue
                        feature.vis_prms["fill_colour"] = colours_dict[feature.group]
                        self.feature_annotation.loc[feature.feature_id, "fill_colour"] = feature.vis_prms["fill_colour"]
            return None
        except Exception as error:
            raise lovis4u.Manager.lovis4uError("Unable to set feature colours based on groups.") from error

    def set_category_colours(self, use_table: bool = True) -> None:
        """Define colours for each category.

        Arguments:
            use_table (bool): Bool value whether table with predefined colours should be used or not.

        Returns:
            None

        """
        try:
            colours_dict = dict()
            if use_table:
                colours_dict.update(
                    pd.read_table(self.prms.args["category_colours"]).set_index("category")["colour"].to_dict())

            feature_categories = list(set([feature.category for locus in self.loci for feature in locus.features
                                           if feature.category and feature.category]))
            if not feature_categories:
                if self.prms.args["verbose"]:
                    print("○ Warning: there are no feature categories to set colours", file=sys.stdout)
            colours_dict = {cat: col for cat, col in colours_dict.items() if cat in feature_categories}

            feature_categories = [ff for ff in feature_categories if ff not in colours_dict.keys()]
            number_of_unique_feature_functions = len(feature_categories)
            colours_rgb = seaborn.color_palette(self.prms.args["category_colour_seaborn_palette"],
                                                number_of_unique_feature_functions,
                                                desat=self.prms.args["category_colour_seaborn_desat"])
            random.shuffle(colours_rgb)
            colours = list(map(lambda x: matplotlib.colors.rgb2hex(x), colours_rgb))
            colours_dict.update({g: c for g, c in zip(list(feature_categories), colours)})
            for locus in self.loci:
                locus.category_colours = colours_dict
            return None
        except Exception as error:
            raise lovis4u.Manager.lovis4uError("Unable to set category colours.") from error

    def pyhmmer_annotation(self) -> None:
        """Run pyhhmmer hmmscan against a set of databases for additional annotation of hotspot proteins.

        Arguments:
            proteomes (Proteomes): Proteomes object.

        Returns:
            None

        """
        try:
            if self.prms.args["verbose"]:
                print(f"○ Preparing data for additional protein annotation with pyhmmer hmmscan...",
                      file=sys.stdout)

            cds_records = [feature.record for locus in self.loci for feature in locus.features if
                           feature.feature_type == "CDS"]
            hmmscan_output_folder = os.path.join(self.prms.args["output_dir"], "hmmscan")
            if os.path.exists(hmmscan_output_folder):
                shutil.rmtree(hmmscan_output_folder)
            os.mkdir(hmmscan_output_folder)
            hmmscan_input_fasta_file_path = os.path.join(hmmscan_output_folder, "input_proteins.fa")
            Bio.SeqIO.write(cds_records, hmmscan_input_fasta_file_path, "fasta")

            alignment_table = lovis4u.Methods.run_pyhmmer(hmmscan_input_fasta_file_path, len(cds_records), self.prms)
            if not alignment_table.empty:
                found_hits_for = alignment_table.index.to_list()
                for locus in self.loci:
                    for feature in locus.features:
                        if feature.feature_id in found_hits_for:
                            alignment_table_row = alignment_table.loc[feature.feature_id]
                            if self.prms.args["update_category_with_database_name"]:
                                feature.category = alignment_table_row["target_db"]
                            if self.prms.args["update_protein_name_with_target_name"]:
                                feature.name = alignment_table_row["target"]
                                feature.vis_prms["label"] = feature.name
                            feature.vis_prms["hmmscan_hit"] = 1
            return None
        except Exception as error:
            raise lovis4u.Manager.lovis4uError("Unable to run pyhmmer hmmscan annotation.") from error

    def reorient_loci(self, ilund4u_mode: bool = False) -> None:
        """Auto re-orient loci (reset strands) of loci if they are not matched.

        Function tries to maximise co-orientation of homologous features.

        Returns:
            None

        """
        try:
            count_of_changed_strands = 0
            loci = [locus for locus in self.loci]
            for locus_index in range(1, len(loci)):
                p_locus = loci[locus_index - 1]
                c_locus = loci[locus_index]
                p_locus_strands = list(set([c["strand"] for c in p_locus.coordinates]))
                c_locus_strands = list(set([c["strand"] for c in c_locus.coordinates]))
                if len(p_locus_strands) == 1 and len(c_locus_strands) == 1:
                    if not ilund4u_mode:
                        pr_locus_features_groups = set([f.group for f in p_locus.features])
                        c_locus_features_groups = set([f.group for f in c_locus.features])
                    else:
                        pr_locus_features_groups = set(
                            [f.group for f in p_locus.features if f.group_type == "conserved"])
                        c_locus_features_groups = set(
                            [f.group for f in c_locus.features if f.group_type == "conserved"])
                    overlapped_f_groups = pr_locus_features_groups & c_locus_features_groups
                    prl_strand, cl_strand = p_locus_strands[0], c_locus_strands[0]
                    pr_locus_features_strands = {f.group: f.strand * prl_strand for f in p_locus.features if
                                                 f.group in overlapped_f_groups}
                    c_locus_features_strands = {f.group: f.strand * cl_strand for f in c_locus.features if
                                                f.group in overlapped_f_groups}
                    codirection_score = 0

                    for ovg in overlapped_f_groups:
                        codirection_score += pr_locus_features_strands[ovg] * c_locus_features_strands[ovg]
                    if codirection_score < 0:
                        count_of_changed_strands += 1
                        annot_coordinates = []
                        for cc in loci[locus_index].coordinates:
                            cc["strand"] *= -1
                            annot_coordinates.append(f"{cc['start']}:{cc['end']}:{cc['strand']}")
                        loci[locus_index].coordinates = loci[locus_index].coordinates[::-1]
                        annot_coordinates = annot_coordinates[::-1]
                        self.locus_annotation.loc[c_locus.seq_id, "coordinates"] = ",".join(annot_coordinates)
                else:
                    if self.prms.args["verbose"]:
                        print("○ Warning: loci reorientation cannot be applied for loci that have both strands in"
                              " pre-defined coordinates for visualisation")
            if self.prms.args["verbose"]:
                if count_of_changed_strands == 0:
                    print(f"⦿ Orientation was not changed for any locus", file=sys.stdout)
                elif count_of_changed_strands == 1:
                    print(f"⦿ Orientation was changed for 1 locus", file=sys.stdout)
                elif count_of_changed_strands > 1:
                    print(f"⦿ Orientation was changed for {count_of_changed_strands} loci", file=sys.stdout)

            return None
        except Exception as error:
            raise lovis4u.Manager.lovis4uError("Unable to define variable feature groups.") from error

    def get_loci_lengths_and_n_of_regions(self) -> list[list[int]]:
        """Get loci lengths and number of regions.

        Returns:
            list: list each element of each contains locus size and number of breaks for visualisation track.

        """
        try:
            loci_sizes = []
            for locus in self.loci:
                number_of_gaps = len(locus.coordinates) - 1
                if locus.circular:
                    for i in range(number_of_gaps):
                        if locus.coordinates[i]["end"] == locus.length and locus.coordinates[i + 1]["start"] == 1:
                            number_of_gaps -= 1
                loci_sizes.append([locus.size, number_of_gaps])
            return loci_sizes
        except Exception as error:
            raise lovis4u.Manager.lovis4uError("Unable to get loci lengths.") from error

    def get_contig_sizes(self):
        """Helper method to get a dictionary with contig sizes for bedgraph file parsing

        Returns:
            dict: dictionary with contig sizes
        """
        try:
            contig_sizes = dict()
            for index, row in self.locus_annotation.iterrows():
                contig_sizes[index] = row["length"]
            return contig_sizes
        except Exception as error:
            raise lovis4u.Manager.lovis4uError("Unable to get contig sizes dict.") from error


class BedGraph:
    """A BedGraph object holds information and process bedgraph files with coverage.

    Attributes:
        filepath (str): path to a bedgraph file
        contig_sizes (dict): Dictionary with contig sizes
        coverage (dict):  Dictionary with format: key - contig id, value - pd.Series with coverage counts.
        parameters (lovis4u.Manager.Parameters): Parameters' class object that holds config and cmd arguments.

    """

    def __init__(self, filepath: str, contig_sizes: dict, parameters: lovis4u.Manager.Parameters, isbigWig=False):
        """Create a BedGraph object

        Arguments:
            filepath (str): path to a bedgraph file.
            contig_sizes:  Dictionary with format: key - contig id, value - full length of the corresponding contig.
            parameters (lovis4u.Manager.Parameters): Parameters' class object that holds config and cmd arguments.

        """
        self.filepath = filepath
        self.prms = parameters
        self.contig_sizes = contig_sizes
        self.coverage = self.__read_file(isbigWig)

    def __read_file(self, isbigWig):
        """Read correspnding bedgraph file.

        Returns:
            coverage (dict):  Dictionary with format: key - contig id, value - pd.Series with coverage counts.
        """
        try:
            if isbigWig:
                temp_file_bedgraph = tempfile.NamedTemporaryFile()
                os.system(f"{self.prms.args['bigWigToBedGraph_binary']} {self.filepath} {temp_file_bedgraph.name}")
                old_name = self.filepath
                self.filepath = temp_file_bedgraph.name
            coverage = dict()
            for contig_id, length in self.contig_sizes.items():
                coverage[contig_id] = pd.Series(np.zeros(length, dtype=int))
            with open(self.filepath, "r") as bedgraph:
                for line in bedgraph:
                    contig, start, end, counts = line.strip().split("\t")
                    if contig in coverage.keys():
                        coverage[contig][int(start):int(end)] = int(counts)
            self.filepath = old_name
            ids_to_remove = [contig_id for contig_id, cov in coverage.items() if cov.sum() == 0]
            for contig_id in ids_to_remove:
                del coverage[contig_id]
            return coverage
        except Exception as error:
            raise lovis4u.Manager.lovis4uError(f"Unable to read {self.filepath} bedgraph file.") from error

    def get_window_coverage(self, contig, start, end):
        """Get coverage for a particular window.

        Arguments:
            contig (str): contig id
            start (int): start coordinate (0-based)
            end (int): end coordinate (0-based)
        """
        try:
            window_coverage = self.coverage[contig][start:end]
            return window_coverage
        except Exception as error:
            raise lovis4u.Manager.lovis4uError(f"Unable to extract window coverage {contig, start, end} for "
                                               f"{self.filepath} bedgraph file.") from error


class CoverageProfiles:
    """A CoverageProfiles object holds information about a set of bedgraph profiles

    Attributes:
            bedgraph_files (list): list with paths to bedgraph files.
            contig_sizes:  Dictionary with format: key - contig id, value - full length of the corresponding contig.
            bedgraphs (list): list with BedGraph objects
            parameters (lovis4u.Manager.Parameters): Parameters' class object that holds config and cmd arguments.

    """

    def __init__(self, contig_sizes: dict, parameters: lovis4u.Manager.Parameters, bedgraph_files: list = [],
                 bigwig_files: list = []):
        """Create a CoverageProfiles object

        Arguments:
            bedgraph_files (list): list with paths to bedgraph files.
            bigwig_files (list): list with paths to bigwig files.
            parameters (lovis4u.Manager.Parameters): Parameters' class object that holds config and cmd arguments.
            contig_sizes:  Dictionary with format: key - contig id, value - full length of the corresponding contig.

        """
        self.bedgraph_files = bedgraph_files
        self.bigwig_files = bigwig_files
        self.prms = parameters
        self.contig_sizes = contig_sizes
        self.bedgraphs = []
        self.bedgraphs += [BedGraph(filepath, contig_sizes, parameters, isbigWig = False) for filepath in bedgraph_files]
        self.bedgraphs += [BedGraph(filepath, contig_sizes, parameters, isbigWig = True) for filepath in bigwig_files]
