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
    """A Feature object represents a locus' feature (currently only CDS) and its properties.

    Attributes:
        feature_id (str): Feature identifier.
        feature_type (str): Type of element (e.g. CDS or tRNA). Currently only CDS are supported.
        start (int): 1-based start genomic coordinate.
        end (int): 1-based end genomic coordinates
        strand (int): Genomic strand (1: plus strand, -1: minus strand).
        name (str): Name of the feature which will be used as a label.
        sequence (Bio.Seq.Seq): Feature's sequence.
        record (Bio.SeqRecord.SeqRecord): SeqRecord object of corresponding feature sequence.
        group (str): Feature group that defines feature's color and meant to represent a set of homologous features.
            Can be set with Loci.mmseqs_cluster() method that uses mmseqs clustering to define CDS feature groups.
        group_type (str): Type of feature group that allow to visualise different set of feature groups differently
            (e.g. plot labels only for features which group_type is "variable" or "labeled").
        category (str): Feature category which initially is built to handle phrogs category annotation.
            In visualisation it defines the "category" color annotation under features.
            Supposed to represent clusters on locus or any second layer of feature properties.
        vis_prms (dict): Visualisation parameters that holds colors, label and other info for Drawing methods.
        prms (lovis4u.Manager.Parameters): Parameters' class object that holds config and cmd arguments.

    """

    def __init__(self, feature_id: str, feature_type: str, start: int, end: int, strand: int, name: str,
                 sequence: Bio.Seq.Seq, group: str, group_type: str, category: str, vis_prms: dict,
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
            group (str): Feature group that defines feature's color and meant to represent a set of homologous features.
                Can be set with Loci.mmseqs_cluster() method that uses mmseqs clustering to define CDS feature groups.
            group_type (str): Type of feature group that allow to visualise different set of feature groups differently
                (e.g. plot labels only for features which group_type is "variable" or "labeled").
            category (str): Feature category which initially is built to handle phrogs category annotation.
                In visualisation it defines the "category" color annotation under features.
                Supposed to represent clusters on locus or any second layer of feature properties.
            vis_prms (dict): Visualisation parameters that holds colors, label and other info for Drawing methods.
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
        category_colors (dict): colors for locus' features categories.
        prms (lovis4u.Manager.Parameters): Parameters' class object that holds config and cmd arguments.

    """

    def __init__(self, seq_id: str, coordinates: list, description: str, length: int, circular: bool,
                 features: list, order: int, parameters: lovis4u.Manager.Parameters, group: typing.Union[int, str] = 1):
        """Create a Locus object.

        Arguments:
            seq_id (str): Sequence identifier. Can be used to label locus.
            coordinates (list): List of regions to be shown. Each region format: dict with keys: start, end, strand and
                corresponding 1-based start, end coordinates and strand (1: plus strand, -1: minus strand).
            description (str): Sequence description that can be used to label locus.
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
        self.length = length
        self.circular = circular
        self.features = features
        self.order = order
        self.group = group
        self.category_colors = dict()
        self.prms = parameters


class Loci:
    """A Loci object holds information about all loci to be plotted and methods for data preparation.

    Attributes:
        loci (list): List of Locus objects.
        loci_annotation (pd.DataFrame): Table with information about each locus that defines visualisation
            (e.g. coordinates for visualisation, description, etc).
        feature_annotation (pd.DataFrame): Table with information about each feature that defines visualisation
            (e.g. group, name, category, etc).
        prms (lovis4u.Manager.Parameters): Parameters' class object that holds config and cmd arguments.

    """

    def __init__(self, parameters=lovis4u.Manager.Parameters):
        """Create a Loci object.

        Arguments:
            parameters (lovis4u.Manager.Parameters): Parameters' class object that holds config and cmd arguments.
        """
        self.loci = []
        self.loci_annotation = pd.DataFrame(columns=["sequence_id", "length", "coordinates", "circular", "description",
                                                     "order", "group"]).set_index("sequence_id")
        self.feature_annotation = pd.DataFrame(columns=["feature_id", "locus_id", "coordinates", "feature_type", "name",
                                                        "group", "group_type", "category", "fill_color", "stroke_color",
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

    def load_loci_annotation_file(self, file_path: str) -> None:
        """Load loci annotation file.

        Arguments:
            file_path (str): File path for a loci annotation file to be loaded.

        Returns:
            None

        """
        annotation_columns = ["sequence_id", "length", "coordinates", "circular", "description", "order", "group"]
        self.loci_annotation = self.__load_annotation_file(file_path, annotation_columns, "sequence_id")
        return None

    def load_features_annotation_file(self, file_path: str) -> None:
        """Load features annotation file.

        Arguments:
            file_path (str): File path for a features annotation file to be loaded.

        Returns:
            None

        """
        annotation_columns = ["feature_id", "locus_id", "coordinates", "feature_type", "name", "group", "group_type",
                              "category", "fill_color", "stroke_color", "show_label"]
        self.feature_annotation = self.__load_annotation_file(file_path, annotation_columns, "feature_id")
        return None

    def __update_loci_annotation(self, record_id: str, record_description: str, record_length: int) -> None:
        """Private method for updating loci annotation.

        Arguments:
            record_id (str): Sequence identifier.
            record_description (str): Sequence description.
            record_length (int): Sequence length.

        Returns:
            None

        """
        if record_id not in self.loci_annotation.index:
            self.loci_annotation.loc[record_id] = {col: None for col in self.loci_annotation.columns}

        default_values = dict(length=record_length, coordinates=f"1:{record_length}:1",
                              description=record_description, circular=1, order=len(self.loci), group=1)
        self.loci_annotation.loc[record_id] = self.loci_annotation.loc[record_id].fillna(default_values)
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
        stroke_color = "default"
        if self.prms.args["set_feature_stroke_color_based_on_fill_color"] and \
                self.feature_annotation.loc[feature_id]["fill_color"] and \
                self.feature_annotation.loc[feature_id]["fill_color"] != "default":
            stroke_color = lovis4u.Methods.scale_lightness(self.feature_annotation.loc[feature_id]["fill_color"],
                                                           self.prms.args["feature_stroke_color_relative_lightness"])
        default_values = dict(locus_id=locus_id, coordinates=coordinates, feature_type=feature_type,
                              name=name, group="", group_type="", category=category, fill_color="default",
                              stroke_color=stroke_color,
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
                    raise lovis4u.Manager.lovis4uError(f"Folder {input_folder} does not exist.")
                gff_files = [os.path.join(input_folder, f) for f in os.listdir(input_folder) if
                             os.path.splitext(f)[-1].lower() == ".gff"]
            elif isinstance(input_f, list):
                gff_files = input_f
            else:
                raise lovis4u.manager.lovis4uError(f"The input for the GFF parsing function must be either a folder or "
                                                   f"a list of files.")
            if not gff_files:
                raise lovis4u.manager.lovis4uError(f"Folder {input_f} does not contain .gff(.GFF) files.")
            for gff_file_path in gff_files:
                gff_file = gff_file_path
                gff_records = list(BCBio.GFF.parse(gff_file_path, limit_info=dict(gff_type=["CDS"])))
                if len(gff_records) != 1:
                    raise lovis4u.Manager.lovis4uError(f"Gff file {gff_file} does not contain information for only"
                                                       f" 1 sequence.")
                gff_record = gff_records[0]
                try:
                    record_locus_sequence = gff_record.seq
                except Bio.Seq.UndefinedSequenceError as error:
                    raise lovis4u.Manager.lovis4uError(f"gff file doesn't contain corresponding sequences.") from error
                if self.prms.args["gff_description_source"] in gff_record.annotations:
                    record_description = gff_record.annotations[self.prms.args["gff_description_source"]][0]
                    if isinstance(record_description, tuple):
                        record_description = " ".join(record_description)
                else:
                    record_description = ""

                self.__update_loci_annotation(gff_record.id, record_description, len(record_locus_sequence))
                locus_annotation_row = self.loci_annotation.loc[gff_record.id]
                coordinates = [dict(zip(["start", "end", "strand"], map(int, c.split(":")))) for c in
                               locus_annotation_row["coordinates"].split(",")]
                record_locus = Locus(seq_id=gff_record.id, coordinates=coordinates,
                                     description=locus_annotation_row["description"],
                                     circular=locus_annotation_row["circular"],
                                     length=locus_annotation_row["length"], parameters=self.prms, features=[],
                                     order=locus_annotation_row["order"])
                features_ids = [i.id for i in gff_record.features]
                if len(features_ids) != len(set(features_ids)):
                    raise lovis4u.Manager.lovis4uError(f"Gff file {gff_file} contains duplicated feature ids while"
                                                       f" only unique are allowed.")
                for gff_feature in gff_record.features:
                    feature_id = gff_feature.id
                    if ilund4u_mode:
                        if gff_record.id not in feature_id:
                            feature_id = f"{gff_record.id}-{feature_id}"
                    transl_table = self.prms.args["default_transl_table"]
                    if "transl_table" in gff_feature.qualifiers.keys():
                        transl_table = int(gff_feature.qualifiers["transl_table"][0])
                    name = ""
                    if self.prms.args["gff_CDS_name_source"] in gff_feature.qualifiers:
                        name = gff_feature.qualifiers[self.prms.args["gff_CDS_name_source"]][0]
                    category = ""
                    if self.prms.args["gff_CDS_category_source"] in gff_feature.qualifiers:
                        category = ",".join(gff_feature.qualifiers[self.prms.args["gff_CDS_category_source"]])
                    for coordinate in record_locus.coordinates:
                        start, end = coordinate["start"], coordinate["end"]
                        if start <= gff_feature.location.start + 1 <= end or start <= gff_feature.location.end <= end:
                            self.__update_feature_annotation(feature_id, record_locus.seq_id,
                                                             f"{int(gff_feature.location.start) + 1}:"
                                                             f"{int(gff_feature.location.end)}:{gff_feature.location.strand}",
                                                             "CDS", category, name)
                            feature_annotation_row = self.feature_annotation.loc[feature_id]
                            feature = Feature(feature_type=feature_annotation_row["feature_type"],
                                              feature_id=feature_id, start=int(gff_feature.location.start) + 1,
                                              end=int(gff_feature.location.end), strand=gff_feature.location.strand,
                                              name=feature_annotation_row["name"],
                                              sequence=gff_feature.translate(record_locus_sequence, table=transl_table,
                                                                             cds=False)[:-1],
                                              group=feature_annotation_row["group"],
                                              group_type=feature_annotation_row["group_type"],
                                              category=feature_annotation_row["category"],
                                              vis_prms=dict(fill_color=feature_annotation_row["fill_color"],
                                                            stroke_color=feature_annotation_row["stroke_color"],
                                                            show_label=feature_annotation_row["show_label"]),
                                              parameters=self.prms)
                            record_locus.features.append(feature)
                self.loci.append(record_locus)
            seq_id_to_order = self.loci_annotation["order"].to_dict()
            self.loci.sort(key=lambda locus: seq_id_to_order[locus.seq_id])
            if self.prms.args["verbose"]:
                if len(self.loci) == 1:
                    print(f"📥 {len(self.loci)} locus was loaded from extended gff files folder", file=sys.stdout)
                else:
                    print(f"📥 {len(self.loci)} loci were loaded from extended gff files folder", file=sys.stdout)
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
            raise lovis4u.Manager.lovis4uError(f"Folder {input_folder} does not exist.")
        try:
            gb_files = [f for f in os.listdir(input_folder) if os.path.splitext(f)[-1].lower() == ".gb"]
            if not gb_files:
                raise lovis4u.Manager.lovis4uError(f"Folder {input_folder} does not contain .gb(.GB) files.")
            for gb_file in gb_files:
                gb_file_path = os.path.join(input_folder, gb_file)
                gb_records = list(Bio.SeqIO.parse(gb_file_path, "genbank"))
                if len(gb_records) != 1:
                    raise lovis4u.Manager.lovis4uError(f"gb file {gb_file} does not contain information for only"
                                                       f" 1 sequence.")
                gb_record = gb_records[0]
                record_locus_sequence = gb_record.seq
                if self.prms.args["genbank_description_source"] == "description":
                    record_description = gb_record.description
                elif "annotations:" in self.prms.args["genbank_description_source"]:
                    feature_description_key = self.prms.args["genbank_description_source"].split(":")[1]
                    record_description = gb_record.annotations[feature_description_key]
                else:
                    record_description = ""
                self.__update_loci_annotation(gb_record.id, record_description, len(record_locus_sequence))
                locus_annotation_row = self.loci_annotation.loc[gb_record.id]
                coordinates = [dict(zip(["start", "end", "strand"], map(int, c.split(":")))) for c in
                               locus_annotation_row["coordinates"].split(",")]
                record_locus = Locus(seq_id=gb_record.id, coordinates=coordinates,
                                     description=locus_annotation_row["description"],
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
                                print(f"❗ Warning: there is no <{self.prms.args['genbank_id_source']}> attribute "
                                      f"for CDS records in {gb_file}. Alternative <{id_source}> was used instead.",
                                      file=sys.stdout)
                            break
                    if id_source == self.prms.args["genbank_id_source"]:
                        raise lovis4u.Manager.lovis4uError(f"There is no <{self.prms.args['genbank_id_source']}> "
                                                           f"attribute for CDS record found in {gb_file}. We tried to"
                                                           f" find any from the alternative list: "
                                                           f"{','.join(self.prms.args['genbank_id_alternative_source'])}"
                                                           f", but they also weren't found.")  # add about cmd parameter
                features_ids = [i.qualifiers[id_source][0] for i in gb_CDSs]
                if len(features_ids) != len(set(features_ids)):
                    raise lovis4u.Manager.lovis4uError(f"GB file {gb_record} contains duplicated feature ids while"
                                                       f" only unique are allowed.")
                for gb_feature in gb_CDSs:
                    feature_id = gb_feature.qualifiers[id_source][0].replace("|", "_")
                    transl_table = self.prms.args["default_transl_table"]
                    if "transl_table" in gb_feature.qualifiers.keys():
                        transl_table = int(gb_feature.qualifiers["transl_table"][0])
                    name = ""
                    if self.prms.args["genbank_CDS_name_source"] in gb_feature.qualifiers:
                        name = gb_feature.qualifiers[self.prms.args["genbank_CDS_name_source"]][0]
                    category = ""
                    if self.prms.args["genbank_CDS_category_source"] in gb_feature.qualifiers:
                        category = ",".join(gb_feature.qualifiers[self.prms.args["genbank_CDS_category_source"]])
                    for coordinate in record_locus.coordinates:
                        start, end = coordinate["start"], coordinate["end"]
                        if start <= gb_feature.location.start + 1 <= end or start <= gb_feature.location.end <= end:
                            self.__update_feature_annotation(feature_id, record_locus.seq_id,
                                                             f"{int(gb_feature.location.start) + 1}:"
                                                             f"{int(gb_feature.location.end)}:"
                                                             f"{gb_feature.location.strand}", "CDS", category, name)
                            feature_annotation_row = self.feature_annotation.loc[feature_id]
                            feature = Feature(feature_type=feature_annotation_row["feature_type"],
                                              feature_id=feature_id, start=int(gb_feature.location.start) + 1,
                                              end=int(gb_feature.location.end),
                                              strand=gb_feature.location.strand, name=feature_annotation_row["name"],
                                              sequence=gb_feature.translate(record_locus_sequence, table=transl_table,
                                                                            cds=False)[:-1],
                                              group=feature_annotation_row["group"],
                                              group_type=feature_annotation_row["group_type"],
                                              category=feature_annotation_row["category"],
                                              vis_prms=dict(fill_color=feature_annotation_row["fill_color"],
                                                            stroke_color=feature_annotation_row["stroke_color"],
                                                            show_label=feature_annotation_row["show_label"]),
                                              parameters=self.prms)

                            record_locus.features.append(feature)
                self.loci.append(record_locus)
            seq_id_to_order = self.loci_annotation["order"].to_dict()
            self.loci.sort(key=lambda locus: seq_id_to_order[locus.seq_id])
            if self.prms.args["verbose"]:
                if len(self.loci) == 1:
                    print(f"📥 {len(self.loci)} locus was loaded from genbank files folder", file=sys.stdout)
                else:
                    print(f"📥 {len(self.loci)} loci were loaded from genbank files folder", file=sys.stdout)

            return None
        except Exception as error:
            raise lovis4u.Manager.lovis4uError("Unable to load loci from gb folder.") from error

    def save_loci_annotation_table(self) -> None:
        """Save loci annotation table to the output folder.

        Output file name is loci_annotation_table.tsv

        Returns:
            None

        """
        try:
            if not os.path.exists(self.prms.args["output_dir"]):
                os.mkdir(self.prms.args["output_dir"])
            file_path = os.path.join(self.prms.args["output_dir"], "loci_annotation_table.tsv")
            self.loci_annotation.to_csv(file_path, sep="\t", index_label="sequence_id")
            if self.prms.args["verbose"]:
                print(f"💌 Loci annotation table was saved to {file_path}", file=sys.stdout)
            return None
        except Exception as error:
            raise lovis4u.Manager.lovis4uError("Unable to save loci annotation table.") from error

    def save_features_annotation_table(self) -> None:
        """Save feature annotation table to the output folder.

        Output file name is features_annotation_table.tsv

        Returns:
            None

        """
        try:
            if not os.path.exists(self.prms.args["output_dir"]):
                os.mkdir(self.prms.args["output_dir"])
            file_path = os.path.join(self.prms.args["output_dir"], "features_annotation_table.tsv")
            self.feature_annotation.to_csv(file_path, sep="\t", index_label="feature_id")
            if self.prms.args["verbose"]:
                print(f"💌 Feature annotation table was saved to {file_path}", file=sys.stdout)
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
            print(f"🚀 Running mmseqs for protein clustering...", file=sys.stdout)
        try:
            feature_records = [feature.record for locus in self.loci for feature in locus.features]
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
                print(f"✅ {num_of_unique_clusters} clusters for {num_of_proteins} proteins were found with mmseqs\n"
                      f"💌 mmseqs clustering results were saved to "
                      f"{os.path.join(mmseqs_output_folder, 'mmseqs_clustering.tsv')}", file=sys.stdout)
            return mmseqs_clustering_results
        except Exception as error:
            raise lovis4u.Manager.lovis4uError("Unable to run mmseqs clustering.") from error

    def define_features_groups(self, dataframe: pd.DataFrame, group_column_name: str = "cluster") -> None:
        """Set features attribute "group" based on input dataframe.

        By default is designed to use mmseqs_cluster() function results as input. If you already have precomputed
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
                    feature.group = dataframe.loc[feature.feature_id, group_column_name]
                    self.feature_annotation.loc[feature.feature_id, "group"] = feature.group
            return None
        except Exception as error:
            raise lovis4u.Manager.lovis4uError("Unable to define protein features groups.") from error

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
                        elif feature.group_type in self.prms.args["feature_group_types_to_show_label"]:
                            feature.vis_prms["show_label"] = 1
                        elif feature.group_type in \
                                self.prms.args["feature_group_types_to_show_label_on_first_occurrence"]:
                            if feature.group not in added_first_occurrence_labels:
                                feature.vis_prms["show_label"] = 1
                                added_first_occurrence_labels.append(feature.group)
                    else:
                        feature.vis_prms["show_label"] = 0
        except Exception as error:
            raise lovis4u.Manager.lovis4uError("Unable define feature labels to be shown.") from error

    def cluster_sequences(self, dataframe: pd.DataFrame, same_cluster=False) -> None:
        """Define loci order and clusters with proteome similarity based hierarchical clustering.
            This function changes the order of loci that are plotted and also updates corresponding to each loci group
            attribute which defines homologues groups of proteomes.

        It's designed to use as input mmseqs_cluster() function results. However, if you have obtained homologues
            groups by other method you can also build pandas dataframe based on that with index corresponding to
            feature id and column "cluster" corresponding to the group.

        Arguments:
              dataframe (pd.DataFrame): dataframe with feature id - group pairs. Its index column should
                represent all loci CDS features.

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
                loci_clusters = [dataframe.loc[feature.feature_id, "cluster"] for feature in locus.features]
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
            if not same_cluster:
                clusters = pd.Series(scipy.cluster.hierarchy.fcluster(linkage_matrix, 0.35, criterion="distance"),
                                     index=loci_ids)
                for locus in self.loci:
                    locus.group = clusters[locus.seq_id]
                    self.loci_annotation.loc[locus.seq_id, "group"] = locus.group
            order = dendrogram["leaves"][::-1]
            self.loci_annotation["initial_order"] = self.loci_annotation["order"]
            for locus_index in range(number_of_loci):
                locus = self.loci[locus_index]
                self.loci_annotation.loc[locus.seq_id, "order"] = order.index(locus_index)
            self.loci_annotation.sort_values(by="order", inplace=True)
            seq_id_to_order = self.loci_annotation["order"].to_dict()
            self.loci.sort(key=lambda locus: seq_id_to_order[locus.seq_id])

            reordered_similarity_matrix = similarity_matrix.reindex(index=self.loci_annotation.index,
                                                                    columns=self.loci_annotation.index)
            if not os.path.exists(self.prms.args["output_dir"]):
                os.mkdir(self.prms.args["output_dir"])
            file_path = os.path.join(self.prms.args["output_dir"], "proteome_similarity_matrix.tsv")
            reordered_similarity_matrix.to_csv(file_path, sep="\t")
            num_of_loci_groups = len(set(self.loci_annotation["group"].to_list()))
            if self.prms.args["verbose"]:
                if num_of_loci_groups == 1:
                    print(f"🫧 Loci order and {num_of_loci_groups} cluster was defined with proteome similarity based "
                          f"hierarchical clustering", file=sys.stdout)
                elif num_of_loci_groups > 1:
                    print(f"🫧 Loci order and {num_of_loci_groups} clusters were defined with proteome similarity based "
                          f"hierarchical clustering", file=sys.stdout)
                print(f"💌 Proteome similarity matrix of loci was saved to {file_path}", file=sys.stdout)
            return None
        except Exception as error:
            raise lovis4u.Manager.lovis4uError("Unable to cluster loci sequences.") from error

    def find_variable_feature_groups(self, mmseqs_results: pd.DataFrame) -> None:
        """Define feature group type attributes (variable or shell/core) based on their conservation in corresponding
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
            loci_clusters_sizes = self.loci_annotation["group"].value_counts()
            loci_clusters_cutoffs = np.round(self.prms.args["CDS_is_variable_cutoff"] * loci_clusters_sizes).astype(int)
            loci_clusters_cutoffs[loci_clusters_cutoffs == 0] = 1
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
                            current_group_cluster_size <= loci_clusters_cutoffs[cluster_locus_group]:
                        cluster_types[cluster_locus_group][cluster] = "variable"
                    else:
                        cluster_types[cluster_locus_group][cluster] = "shell/core"
            for locus in self.loci:
                locus_group = locus.group
                for feature in locus.features:
                    if feature.group_type and self.prms.args["keep_predefined_groups"]:
                        continue
                    feature.group_type = cluster_types[locus_group][feature.group]
                    self.feature_annotation.loc[feature.feature_id, "group_type"] = feature.group_type
            return None
        except Exception as error:
            raise lovis4u.Manager.lovis4uError("Unable to define variable feature groups.") from error

    def set_feature_colors_based_on_groups(self) -> None:
        """Define features fill color based on corresponding feature group and group types.

        Returns:
            None

        """
        try:
            feature_groups = set([feature.group for locus in self.loci for feature in locus.features if feature.group])
            if self.prms.args["feature_group_types_to_set_color"] and \
                    "all" not in self.prms.args["feature_group_types_to_set_color"]:
                feature_groups = set([feature.group for locus in self.loci for feature in locus.features
                                      if feature.group and feature.group_type in
                                      self.prms.args["feature_group_types_to_set_color"]])
            number_of_unique_feature_groups = len(feature_groups)
            if self.prms.args["groups_fill_color_palette_lib"] == "seaborn":
                colors_rgb = seaborn.color_palette(self.prms.args["groups_fill_color_seaborn_palette"],
                                                   number_of_unique_feature_groups,
                                                   desat=self.prms.args["groups_fill_color_seaborn_desat"])
                random.shuffle(colors_rgb)
            elif self.prms.args["groups_fill_color_palette_lib"] == "distinctipy":
                colors_rgb = distinctipy.get_colors(number_of_unique_feature_groups,
                                                    exclude_colors=[(1, 1, 1), (0, 0, 0)],
                                                    pastel_factor=self.prms.args["groups_fill_colors_pastel_factor"])
            colors = list(map(lambda x: matplotlib.colors.rgb2hex(x), colors_rgb))
            colors_dict = {g: c for g, c in zip(list(feature_groups), colors)}
            for locus in self.loci:
                for feature in locus.features:
                    if feature.group in feature_groups:
                        if self.prms.args["keep_predefined_colors"] and feature.vis_prms["fill_color"] != "default":
                            continue
                        feature.vis_prms["fill_color"] = colors_dict[feature.group]
                        self.feature_annotation.loc[feature.feature_id, "fill_color"] = feature.vis_prms["fill_color"]
            return None
        except Exception as error:
            raise lovis4u.Manager.lovis4uError("Unable to set feature colors based on groups.") from error

    def set_category_colors(self, use_table: bool = True) -> None:
        """Define colors for each category.

        Arguments:
            use_table (bool): Bool value whether table with predefined colors should be used or not.

        Returns:
            None

        """
        try:
            colors_dict = dict()
            if use_table:
                colors_dict.update(
                    pd.read_table(self.prms.args["category_colors"]).set_index("category")["color"].to_dict())

            feature_categories = list(set([feature.category for locus in self.loci for feature in locus.features
                                           if feature.category and feature.category]))
            if not feature_categories:
                if self.prms.args["verbose"]:
                    print("❗ Warning: there are no feature categories to set colors", file=sys.stdout)
            colors_dict = {cat: col for cat, col in colors_dict.items() if cat in feature_categories}

            feature_categories = [ff for ff in feature_categories if ff not in colors_dict.keys()]
            number_of_unique_feature_functions = len(feature_categories)
            colors_rgb = seaborn.color_palette(self.prms.args["category_color_seaborn_palette"],
                                               number_of_unique_feature_functions,
                                               desat=self.prms.args["category_color_seaborn_desat"])
            random.shuffle(colors_rgb)
            colors = list(map(lambda x: matplotlib.colors.rgb2hex(x), colors_rgb))
            colors_dict.update({g: c for g, c in zip(list(feature_categories), colors)})
            for locus in self.loci:
                locus.category_colors = colors_dict
            return None
        except Exception as error:
            raise lovis4u.Manager.lovis4uError("Unable to set category colors.") from error

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
                            [f.group for f in p_locus.features if f.group_type == "shell/core"])
                        c_locus_features_groups = set(
                            [f.group for f in c_locus.features if f.group_type == "shell/core"])
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
                        self.loci_annotation.loc[c_locus.seq_id, "coordinates"] = ",".join(annot_coordinates)
                else:
                    if self.prms.args["verbose"]:
                        print("❗ Warning: loci reorientation cannot be applied for loci that have both strands in"
                              " pre-defined coordinates for visualisation")
            if self.prms.args["verbose"]:
                if count_of_changed_strands == 0:
                    print(f"🔁 Orientation was not changed for any locus", file=sys.stdout)
                elif count_of_changed_strands == 1:
                    print(f"🔁 Orientation was changed for 1 locus", file=sys.stdout)
                elif count_of_changed_strands > 1:
                    print(f"🔁 Orientation was changed for {count_of_changed_strands} loci", file=sys.stdout)

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
