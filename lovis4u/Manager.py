"""
This module provides managing classes and methods for the tool.
"""
import collections
import statistics
import argparse
import configs
import time
import math
import sys
import os
import re

import reportlab.pdfbase.pdfmetrics as pdfmetrics
from reportlab.lib.units import cm, mm
import reportlab.pdfbase.ttfonts
import reportlab.pdfgen.canvas
import reportlab.rl_config
import reportlab.pdfgen

import lovis4u


class lovis4uError(Exception):
    """A class for exceptions parsing inherited from the Exception class.

    """
    pass


class Parameters:
    """A Parameters object holds and parse command line's and config's arguments.

    A Parameters object have to be created in each script since it's used almost by each
        class of the tool as a mandatory argument.

    Attributes:
        args (dict): dictionary that holds all arguments.
    """

    def __init__(self):
        """Create a Parameters object.

        """
        self.args = dict(debug=False, verbose=False)
        self.cmd_arguments = dict()

    def parse_cmd_arguments(self) -> None:
        """Parse command-line args

        Returns:
            None

        """
        parser = argparse.ArgumentParser(prog="lovis4u", add_help=False,
                                         usage="lovis4u [-gff gff_folder | -gb gb_folder] [optional args]")
        parser.add_argument("-data", "--data", dest="lovis4u_data", action="store_true")
        parser.add_argument("-linux", "--linux", dest="linux", action="store_true", default=None)
        parser.add_argument("-mac", "--mac", dest="mac", action="store_true", default=None)
        parser.add_argument("-get-hmms", "--get-hmms", dest="get_hmms", action="store_true")
        mutually_exclusive_group = parser.add_mutually_exclusive_group()
        mutually_exclusive_group.add_argument("-gff", "--gff", dest="gff", type=str, default=None)
        mutually_exclusive_group.add_argument("-gb", "--gb", dest="gb", type=str, default=None)
        parser.add_argument("-ufid", "--use-filename-as-id", dest="use_filename_as_contig_id", action="store_true",
                            default=None)
        parser.add_argument("-alip", "--add-locus-id-prefix", dest="add_locus_id_prefix", action="store_true",
                            default=None)
        parser.add_argument("-laf", "--locus-annotation-file", dest="locus-annotation", type=str, default=None)
        parser.add_argument("-faf", "--feature-annotation-file", dest="feature-annotation", type=str, default=None)
        parser.add_argument("-mmseqs-off", "--mmseqs-off", dest="mmseqs", action="store_false")
        parser.add_argument("-hmmscan", "--run-hmmscan", dest="run_hmmscan_search", action="store_true", default=None)
        parser.add_argument("-dm", "--defence-models", dest="defence_models",
                            choices=["both", "DefenseFinder", "PADLOC"], default=None)
        parser.add_argument("-hmm", "--add-hmm-models", dest="hmm_models", action='append', nargs='*', default=[])
        parser.add_argument("-omh", "--only-mine-hmms", dest="only_mine_hmms", action="store_true", default=None)
        parser.add_argument("-kdn", "--keep-default-name", dest="update_protein_name_with_target_name",
                            action="store_false", default=None)
        parser.add_argument("-kdc", "--keep-default-category", dest="update_category_with_database_name",
                            action="store_false", default=None)
        parser.add_argument("-salq", "--show-all-labels-for-query", dest="show_all_label_for_query_proteins",
                            action="store_true", default=None)
        parser.add_argument("-cl-owp", "--cluster-only-window-proteins", dest="cluster_all_proteins",
                            action="store_false", default=None)
        parser.add_argument("-fv-off", "--find-variable-off", dest="find-variable", action="store_false")
        parser.add_argument("-cl-off", "--clust_loci-off", dest="clust_loci", action="store_false")
        parser.add_argument("-oc", "--one-cluster", dest="one_cluster", action="store_true", default=None)
        parser.add_argument("-reorient_loci", "--reorient_loci", dest="reorient_loci", action="store_true")
        parser.add_argument("-lls", "--locus-label-style", dest="locus_label_style",
                            choices=["id", "description", "full"], default=None)
        parser.add_argument("-llp", "--locus-label-position", dest="locus_label_position",
                            choices=["left", "bottom"], default=None)
        parser.add_argument("-sgc-off", "--set-group-colour-off", dest="set-group-colour", action="store_false")
        parser.add_argument("-sgcf", "--set-group-colour-for", dest="feature_group_types_to_set_colour", nargs="*",
                            type=str, default=None)
        parser.add_argument("-scc", "--set-category-colour", dest="set-category-colour", action="store_true")
        parser.add_argument("-cct", "--category-colour-table", dest="category_colours", type=str, default=None)
        parser.add_argument("-safl", "--show-all-feature-labels", dest="show_all_feature_labels",
                            action="store_true")
        parser.add_argument("-sflf", "--show-feature-label-for", dest="feature_group_types_to_show_label", nargs="*",
                            type=str, default=None)
        parser.add_argument("-sfflf", "--show-first-feature-label-for",
                            dest="feature_group_types_to_show_label_on_first_occurrence", nargs="*",
                            type=str, default=None)
        parser.add_argument("-ifl", "--ignored-feature-labels", dest="feature_labels_to_ignore", nargs="*",
                            type=str, default=None)
        parser.add_argument("-snl", "--show-noncoding-labels", dest="show_noncoding_labels", action="store_true",
                            default=None)
        parser.add_argument("-sfnl", "--show-first-noncoding-label", dest="show_first_noncoding_label",
                            action="store_true", default=None)
        parser.add_argument("-hl", "--homology-links", dest="homology-track", action="store_true")
        parser.add_argument("-slt", "--scale-line-track", dest="draw_scale_line_track", action="store_true",
                            default=None)
        parser.add_argument("-sxa", "--show-x-axis", dest="draw_individual_x_axis", action="store_true", default=None)
        parser.add_argument("-hix", "--hide-x-axis", dest="draw_individual_x_axis", action="store_false", default=None)
        parser.add_argument("-dml", "--draw-middle-line", dest="draw_middle_line", action="store_true", default=None)
        parser.add_argument("-mm-per-nt", "--mm-per-nt", dest="mm_per_nt", type=float, default=None)
        parser.add_argument("-fw", "--figure-width", dest="figure_width", type=float, default=None)
        parser.add_argument("-o", dest="output_dir", type=str, default=None)
        parser.add_argument("--pdf-name", dest="pdf-name", type=str, default="lovis4u.pdf")
        parser.add_argument("-c", dest="config_file", type=str, default="standard")
        parser.add_argument("-v", "--version", action="version", version="%(prog)s 0.0.14")
        parser.add_argument("-q", "--quiet", dest="verbose", default=True, action="store_false")
        parser.add_argument("--parsing-debug", "-parsing-debug", dest="parsing_debug", action="store_true")
        parser.add_argument("--debug", "-debug", dest="debug", action="store_true")
        parser.add_argument("-h", "--help", dest="help", action="store_true")
        args = parser.parse_args()
        args = vars(args)
        if len(sys.argv[1:]) == 0:
            args["help"] = True
        args_to_keep = ["locus-annotation", "feature-annotation", "gb", "gff"]
        filtered_args = {k: v for k, v in args.items() if v is not None or k in args_to_keep}
        self.cmd_arguments = filtered_args
        if args["lovis4u_data"]:
            lovis4u.Methods.copy_package_data()
            sys.exit()
        if args["linux"]:
            lovis4u.Methods.adjust_paths("linux")
            sys.exit()
        if args["mac"]:
            lovis4u.Methods.adjust_paths("mac")
            sys.exit()
        if args["get_hmms"]:
            self.load_config()
            lovis4u.Methods.get_HMM_models(self.args)
            sys.exit()
        if args["help"]:
            help_message_path = os.path.join(os.path.dirname(__file__), "lovis4u_data", "help.txt")
            with open(help_message_path, "r") as help_message:
                print(help_message.read(), file=sys.stdout)
                sys.exit()
        if not args["gff"] and not args["gb"]:
            raise lovis4u.Manager.lovis4uError("-gff or -gb parameter with folder path should be provided")
        return None

    def load_config(self, path: str = "standard") -> None:
        """Load configuration file.

        Arguments
            path (str): path to a config file or name (only standard available at this moment).

        Returns:
            None

        """
        try:
            initial_path = ""
            if path in ["standard", "A4p1", "A4p2", "A4L"]:
                initial_path = path
                path = os.path.join(os.path.dirname(__file__), "lovis4u_data", f"{path}.cfg")
            config = configs.load(path).get_config()
            internal_dir = os.path.dirname(__file__)
            for key in config["root"].keys():
                if isinstance(config["root"][key], str):
                    if config["root"][key] == "None":
                        config["root"][key] = None
                if isinstance(config["root"][key], str) and "{internal}" in config["root"][key]:
                    config["root"][key] = config["root"][key].replace("{internal}",
                                                                      os.path.join(internal_dir, "lovis4u_data"))
            config["root"]["output_dir"] = config["root"]["output_dir"].replace("{current_date}",
                                                                                time.strftime("%Y_%m_%d-%H_%M"))
            keys_to_transform_to_list = ["feature_group_types_to_set_colour", "feature_group_types_to_show_label",
                                         "genbank_id_alternative_source", "feature_labels_to_ignore",
                                         "feature_group_types_to_show_label_on_first_occurrence",
                                         "gff_noncoding_name_alternative_source", "hmm_config_names", "database_names"]
            for ktl in keys_to_transform_to_list:
                if isinstance(config["root"][ktl], str):
                    if config["root"][ktl] != "None":
                        config["root"][ktl] = [config["root"][ktl]]
                    else:
                        config["root"][ktl] = []
                config["root"][ktl] = list(config["root"][ktl])
            self.args.update(config["root"])
            self.load_palette()
            self.load_fonts()
            if self.cmd_arguments:
                self.args.update(self.cmd_arguments)

            if os.path.exists(self.args["output_dir"]):
                if self.args["verbose"]:
                    if initial_path:
                        print(f"○ Loaded configuration file: '{initial_path}'. List of available: "
                              f"standard (auto-size),\n\tA4p1 (A4 page one-column [90mm]), A4p2 (A4 page two-column"
                              f" [190mm]), and A4L (A4 landscape [240mm]).\n\tUse -c/--config <name> parameter "
                              f"for choosing.", file=sys.stdout)
                    print("○ Warning: the output folder already exists. Results will be rewritten (without removal "
                          "other files in this folder)", file=sys.stdout)
            self.args["locus_label_description_font_face"] = self.args[f"locus_label_description_font_face_" \
                                                                       f"{self.args['locus_label_position']}"]
            self.args["locus_label_id_font_face"] = self.args[f"locus_label_id_font_face_" \
                                                              f"{self.args['locus_label_position']}"]
            if self.args["show_noncoding_labels"]:
                self.args["feature_group_types_to_show_label"].append("noncoding")
            if self.args["show_first_noncoding_label"]:
                self.args["feature_group_types_to_show_label_on_first_occurrence"].append("noncoding")

            if self.args["hmm_models"]:
                if self.args["only_mine_hmms"]:
                    self.args["hmm_config_names"] = []
                    self.args["database_names"] = []
                for hmm_model in self.args["hmm_models"]:
                    hmm_folder_path = hmm_model[0]
                    if not os.path.exists(hmm_folder_path):
                        raise lovis4uError(f"HMM profiles folder {hmm_folder_path} does not exists.")
                    if len(hmm_model) == 1:
                        database_name = os.path.basename(hmm_folder_path)
                    else:
                        database_name = hmm_model[1]
                    self.args[f"users_hmm_{database_name}"] = hmm_folder_path
                    self.args["hmm_config_names"].append(f"users_hmm_{database_name}")
                    self.args["database_names"].append(database_name)
            # Check conflicts
            if self.args["draw_individual_x_axis"] and self.args["locus_label_position"] == "bottom":
                raise lovis4uError("Individual x-axis cannot be plotted when locus label position"
                                   " set as 'bottom'.")
            return None
        except Exception as error:
            raise lovis4uError("Unable to parse the specified config file. Please check your config file "
                               "or written name.") from error

    def load_palette(self) -> None:
        """Load palette file.

        Returns:
            None

        """
        try:
            palette_path = self.args[f"palette"]
            self.args[f"palette"] = configs.load(palette_path).get_config()["root"]
            return None
        except Exception as error:
            raise lovis4uError("Unable to load palette.") from error

    def load_fonts(self) -> None:
        """Load fonts.

        Returns:
            None

        """
        try:
            font_pattern = re.compile(r"^font_(.*)$")
            font_subdict = {font_pattern.match(key).group(1): value for key, value in self.args.items() if
                            font_pattern.match(key)}
            for font_type, font_path in font_subdict.items():
                pdfmetrics.registerFont(reportlab.pdfbase.ttfonts.TTFont(font_type, font_path))
            return None
        except Exception as error:
            raise lovis4uError("Unable to load fonts.") from error


class CanvasManager:
    """Canvas manager object responsible for preprocessing data for visualisation and interaction between visualisation
        and raw data.

    Attributes:
         layout (dict): Defines size and coordinate system of a canvas.
         tracks (list): List containing Track objects each of them represents visualisation unit (e.g. particular locus).
         cross_tracks (list): List containing CrossTrack objects each of them represents visualisation unit that
            interacts with multiple regular Track objects.
        prms (Parameters): Parameters' class object that holds config and cmd arguments.

    """

    def __init__(self, parameters):
        """Create a CanvasManager object.

        Arguments:
            parameters (Parameters): Parameters' class object that holds config and cmd arguments.

        """
        self.layout = dict()
        self.tracks = []
        self.cross_tracks = []
        self.prms = parameters

    def define_layout(self, loci) -> None:
        """Define canvas' layout based on input loci.

        Arguments:
            loci (lovis4u.DataProcessing.Loci): Loci object with information about sequences and features.

        Returns:
            None

        """
        try:
            annotated_descriptions = [i.description for i in loci.loci if i.description]
            if not annotated_descriptions and self.prms.args["locus_label_style"] != "id":
                if self.prms.args["verbose"]:
                    print("○ Warning message: the annotation lacks description. Locus label style is "
                          "changed to 'id'")
                self.prms.args["locus_label_style"] = "id"
            if self.prms.args["locus_label_position"] == "left":
                if self.prms.args["locus_label_style"] == "full":
                    label_height = self.prms.args["feature_height"] * mm * 0.4
                else:
                    label_height = min(1, self.prms.args["locus_label_size"]) * self.prms.args["feature_height"] * mm
                label_font_size = lovis4u.Methods.str_height_to_size(label_height,
                                                                     self.prms.args["locus_label_id_font_face"])
            elif self.prms.args["locus_label_position"] == "bottom":
                label_font_size = self.prms.args["bottom_locus_label_font_size"]
                label_height = lovis4u.Methods.str_font_size_to_height(label_font_size,
                                                                       self.prms.args["locus_label_id_font_face"])
            else:
                raise lovis4u.Manager.lovis4uError("Locus label position parameter should be either 'bottom' "
                                                   "or 'left'.")
            self.prms.args["locus_label_height"] = label_height
            self.prms.args["locus_label_font_size"] = label_font_size
            max_id_string_width = max([pdfmetrics.stringWidth(i.seq_id, self.prms.args["locus_label_id_font_face"],
                                                              self.prms.args["locus_label_font_size"]) for i in
                                       loci.loci])
            if self.prms.args["locus_label_style"] != "id":
                max_descr_string_width = \
                    max([pdfmetrics.stringWidth(i, self.prms.args["locus_label_description_font_face"],
                                                self.prms.args["locus_label_font_size"])
                         for i in annotated_descriptions])

            if self.prms.args["locus_label_style"] == "full":
                max_label_string_width = max(max_id_string_width, max_descr_string_width)
            elif self.prms.args["locus_label_style"] == "id":
                max_label_string_width = max_id_string_width
            elif self.prms.args["locus_label_style"] == "description":
                max_label_string_width = max_descr_string_width
            loci_lengths = loci.get_loci_lengths_and_n_of_regions()
            self.layout["total_nt_width"] = max(i[0] for i in loci_lengths)
            if self.prms.args["figure_width"]:
                if self.prms.args["locus_label_position"] == "left":
                    figure_width_for_loci = self.prms.args["figure_width"] * mm - max_label_string_width - \
                                            2 * self.prms.args["margin"] * mm \
                                            - self.prms.args["gap_after_locus_label"] * mm
                elif self.prms.args["locus_label_position"] == "bottom":
                    figure_width_for_loci = self.prms.args["figure_width"] * mm - 2 * self.prms.args["margin"] * mm
                self.prms.args["mm_per_nt"] = mm * figure_width_for_loci / self.layout["total_nt_width"]
            else:
                width_per_nt = self.prms.args["mm_per_nt"] * mm
                figure_width_for_loci = width_per_nt * self.layout["total_nt_width"]
                if figure_width_for_loci < 8 * cm:
                    width_per_nt = 8 * cm / self.layout["total_nt_width"]
                self.prms.args["mm_per_nt"] = width_per_nt / mm
            self.layout["width_per_nt"] = self.prms.args["mm_per_nt"] / mm
            self.layout["x_gap_between_regions"] = self.prms.args["gap_between_regions"] * mm
            each_loci_region_width = [
                (i[0] * self.layout["width_per_nt"]) + (i[1] * self.layout["x_gap_between_regions"])
                for i in loci_lengths]
            max_loci_region_length = max(each_loci_region_width)
            self.layout["locus_label_left_border"] = self.prms.args["margin"] * mm
            if self.prms.args["locus_label_position"] == "left":
                self.layout["locus_label_right_border"] = self.layout[
                                                              "locus_label_left_border"] + max_label_string_width
                self.layout["loci_tracks_left_border"] = self.layout["locus_label_right_border"] + \
                                                         self.prms.args["gap_after_locus_label"] * mm
                self.layout["loci_tracks_right_border"] = self.layout["loci_tracks_left_border"] + \
                                                          max_loci_region_length
            elif self.prms.args["locus_label_position"] == "bottom":
                self.layout["loci_tracks_left_border"] = self.prms.args["margin"] * mm
                self.layout["loci_tracks_right_border"] = self.layout["loci_tracks_left_border"] + \
                                                          max_loci_region_length
            self.layout["figure_width"] = self.layout["loci_tracks_right_border"] + self.prms.args["margin"] * mm
            self.layout["figure_height"] = self.prms.args["margin"] * mm
            return None
        except Exception as error:
            raise lovis4u.Manager.lovis4uError("Unable to define layout from loci data.") from error

    def add_loci_tracks(self, loci) -> None:
        """Add loci tracks to your canvas.

        Arguments:
            loci (lovis4u.DataProcessing.Loci): Loci object with information about sequences and features.

        Returns:
            None

        """
        try:
            for locus in loci.loci:
                locus_loader = LocusLoader(self.prms)
                locus_loader.prepare_track_specific_data(locus, self.layout.copy())
                locus_track_height = locus_loader.calculate_track_height()
                self.layout["figure_height"] += locus_track_height + self.prms.args["gap"] * mm
                locus_track = locus_loader.create_track()
                self.tracks.append(locus_track)
            if self.prms.args["verbose"]:
                print(f"⦿ {len(loci.loci)} loci tracks were added to the canvas", file=sys.stdout)
            return None
        except Exception as error:
            raise lovis4u.Manager.lovis4uError("Unable to add loci tracks to the canvas.") from error

    def add_categories_colour_legend_track(self, loci) -> None:
        """Add categories colour legend tracks to your canvas.

        Arguments:
            loci (lovis4u.DataProcessing.Loci): Loci object with information about sequences and features.

        Returns:
            None

        """
        try:
            colour_legend_loader = CategoriesColorLegendLoader(self.prms)
            colour_legend_loader.prepare_track_specific_data(self.layout.copy(), loci)
            colour_legend_track_height = colour_legend_loader.calculate_track_height()
            if colour_legend_track_height != 0:
                self.layout["figure_height"] += colour_legend_track_height + self.prms.args["gap"] * mm
            colour_legend_track = colour_legend_loader.create_track()
            if isinstance(colour_legend_track, lovis4u.Drawing.ColorLegendVis):
                self.tracks.append(colour_legend_track)
                if self.prms.args["verbose"]:
                    print(f"⦿ Categories colour legend track was added to the canvas", file=sys.stdout)
                return None
        except Exception as error:
            raise lovis4u.Manager.lovis4uError("Unable to add categories colour legend track to the canvas.") from error

    def add_scale_line_track(self) -> None:
        """Add scale line tracks to your canvas.

        Returns:
            None

        """
        try:
            scale_loader = ScaleLoader(self.prms)
            scale_loader.prepare_track_specific_data(self.layout.copy())
            scale_track_height = scale_loader.calculate_track_height()
            self.layout["figure_height"] += scale_track_height + self.prms.args["gap"] * mm
            scale_track = scale_loader.create_track()
            self.tracks.append(scale_track)
            if self.prms.args["verbose"]:
                print(f"⦿ Scale line track was added to the canvas", file=sys.stdout)
        except Exception as error:
            raise lovis4u.Manager.lovis4uError("Unable to add scale line track to the canvas.") from error

    def add_homology_track(self) -> None:
        """Add homology track to your canvas.

        You should add this track after you added loci tracks.

        Returns:
            None

        """
        try:
            loci_tracks = [i for i in self.tracks.copy() if isinstance(i, lovis4u.Drawing.LocusVis)]
            if not loci_tracks:
                raise lovis4u.Manager.lovis4uError("Unable to create homology track if no loci track was added.")
            for lt in loci_tracks:
                lt.track_data["clean_features_coordinates"] = True
            self.cross_tracks.append(lovis4u.Drawing.HomologyTrack(self.layout, loci_tracks, self.prms))
            if self.prms.args["verbose"]:
                print(f"⦿ Homology track was added to the canvas", file=sys.stdout)
            return None
        except Exception as error:
            raise lovis4u.Manager.lovis4uError("Unable to add scale line track to the canvas.") from error

    def plot(self, filename: str) -> None:
        """Plot all added tracks and save the plot as pdf.

        Arguments:
            filename (str): filename for the output pdf.

        Returns:
            None

        """
        try:
            file_path = os.path.join(self.prms.args["output_dir"], filename)
            self.layout["figure_height"] += (self.prms.args["margin"] - self.prms.args["gap"]) * mm
            plot = Canvas(file_path, self.layout["figure_width"], self.layout["figure_height"], self.prms)

            for cross_track in self.cross_tracks:
                cross_track.draw(plot.canvas)

            current_y_coordinate = self.layout["figure_height"] - self.prms.args["margin"] * mm
            for track in self.tracks:
                track.layout["current_y_coordinate"] = current_y_coordinate
                track.layout["figure_height"] = self.layout["figure_height"]
                track.draw(plot.canvas)
                current_y_coordinate -= track.track_data["track_height"] + self.prms.args["gap"] * mm
            plot.save()
            if self.prms.args["verbose"]:
                print(f"⦿ lovis4u plot was saved as: {file_path}", file=sys.stdout)
            return None
        except Exception as error:
            raise lovis4u.Manager.lovis4uError("Unable to plot the canvas and save the figure.") from error


class Loader:
    """Parent class for tracks loaders.

    Attributes:
        prms (Parameters): Parameters' class object.
        layout (dict): Layout built by CanvasManager's define_layout() method.
        track_data (dict): Track specific data that will be sent to the Drawing module.

    """

    def __init__(self, parameters: Parameters):
        """Parent's constructor for creating a Loader class object.

        Arguments:
            parameters (Parameters): Parameters' class object that holds config and cmd arguments.

        """
        self.prms = parameters
        self.layout = None
        self.track_data = None

    def prepare_track_specific_data(self) -> None:
        """Empty parent's method for data preparation.

        Returns:
            None

        """
        pass

    def calculate_track_height(self) -> None:
        """Empty parent's method for track height calculation.

        Returns:
            None

        """
        pass

    def create_track(self) -> None:
        """Empty parent's method for track initialisation.

        Returns:
            None

        """
        pass


class LocusLoader(Loader):
    """A LocusLoader object prepares data for a Locus track Drawing object.

    Attributes:
        prms (Parameters): Parameters' class object.
        layout (dict): Layout built by CanvasManager's define_layout() method.
        track_data (dict): Track specific data that will be sent to the Drawing module.


    """

    def __init__(self, parameters):
        """Create a LocusLoader object.

        Arguments:
            parameters (Parameters): Parameters' class object that holds config and cmd arguments.

        """
        super().__init__(parameters)

    def prepare_track_specific_data(self, locus, layout: dict) -> None:
        """Prepare LocusLoader specific data.

        Attributes:
            locus (lovis4u.DataProcessing.Locus): corresponding locus object.
            layout (dict): Layout built by CanvasManager's define_layout() method.

        Returns:
            None

        """
        try:
            self.layout = layout
            layout["inverse_y_coordinate"] = layout["figure_height"]
            track_data = dict()
            self.track_data = track_data
            track_data["locus_id_width"] = pdfmetrics.stringWidth(locus.seq_id,
                                                                  self.prms.args["locus_label_id_font_face"],
                                                                  self.prms.args["locus_label_font_size"])
            two_space_width = pdfmetrics.stringWidth("  ", self.prms.args["locus_label_id_font_face"],
                                                     self.prms.args["locus_label_font_size"])
            track_data["two_space_width"] = two_space_width
            track_data["locus_id"] = locus.seq_id
            track_data["locus_description"] = locus.description
            if locus.description:
                track_data["locus_description_width"] = pdfmetrics.stringWidth(locus.description,
                                                                               self.prms.args["locus_label_description"
                                                                                              "_font_face"],
                                                                               self.prms.args["locus_label_font_size"])
            if self.prms.args["locus_label_position"] == "bottom":
                txt_coordinates = []
                for i in range(len(locus.coordinates)):
                    cc = locus.coordinates[i].copy()
                    cc_txt = f"{cc['start']}:{cc['end']}{'(+)' if cc['strand'] == 1 else '(-)'}"
                    if i > 0 and locus.circular:
                        if cc["start"] == 1 and locus.coordinates[i - 1]["end"] == locus.length:
                            txt_coordinates[-1] = f"{locus.coordinates[i - 1]['start']}:{cc['end']}" \
                                                  f"{'(+)' if cc['strand'] == 1 else '(-)'}"
                            continue
                    txt_coordinates.append(cc_txt)
                track_data["text_coordinates"] = ", ".join(txt_coordinates)
                track_data["text_coordinates_width"] = pdfmetrics.stringWidth(locus.description,
                                                                              self.prms.args[
                                                                                  "locus_label_id_font_face"],
                                                                              self.prms.args["locus_label_font_size"])

            track_data["f_label_font_size"] = self.prms.args["feature_label_font_size"]
            track_data["f_label_height"] = lovis4u.Methods.str_font_size_to_height(
                track_data["f_label_font_size"], self.prms.args["feature_label_font_face"])
            track_data["feature_label_gap"] = self.prms.args["feature_label_gap"] * track_data["f_label_height"]

            # Managing features positions and parameters
            track_data["clean_features_coordinates"] = False
            track_data["features"] = []
            features_taken_nt_coordinates = []
            for feature in locus.features:
                features_taken_nt_coordinates.append([feature.start, feature.end])
                feature.vis_prms["type"] = feature.feature_type
                if feature.feature_type == "CDS":
                    if feature.vis_prms["fill_colour"] == "default":
                        feature.vis_prms["fill_colour"] = lovis4u.Methods.get_colour("feature_default_fill_colour",
                                                                                     self.prms)
                    if self.prms.args["set_feature_stroke_colour_based_on_fill_colour"] and \
                            feature.vis_prms["stroke_colour"] == "default":
                        scale_l = self.prms.args["feature_stroke_colour_relative_lightness"]
                        feature.vis_prms["stroke_colour"] = lovis4u.Methods.scale_lightness(
                            feature.vis_prms["fill_colour"],
                            scale_l)
                    elif not self.prms.args["set_feature_stroke_colour_based_on_fill_colour"] and \
                            feature.vis_prms["stroke_colour"] == "default":
                        feature.vis_prms["stroke_colour"] = lovis4u.Methods.get_colour("feature_default_stroke_colour",
                                                                                       self.prms)
                else:
                    feature.vis_prms["fill_colour"] = None
                    feature.vis_prms["stroke_colour"] = lovis4u.Methods.get_colour("noncoding_feature_default_"
                                                                                   "stroke_colour", self.prms)
                f_label_width = 0
                if feature.vis_prms["show_label"] and feature.vis_prms["label"]:
                    f_label_width = pdfmetrics.stringWidth(feature.vis_prms["label"],
                                                           self.prms.args["feature_label_font_face"],
                                                           track_data["f_label_font_size"])
                feature_vis_data = feature.vis_prms
                feature_vis_data["coordinates"] = lovis4u.Methods.feature_nt_to_x_transform(feature.start, feature.end,
                                                                                            feature.strand, locus,
                                                                                            layout)
                feature_vis_data["label_width"] = f_label_width
                feature_vis_data["feature_width"] = feature_vis_data["coordinates"]["end"] - \
                                                    feature_vis_data["coordinates"]["start"]
                feature_vis_data["group"] = feature.group
                track_data["features"].append(feature_vis_data)
            # Managing category visualisation
            feature_cateregories = set([feature.category for feature in locus.features if feature.category])
            track_data["functions_coordinates"] = None
            if feature_cateregories and locus.category_colours:
                track_data["category_colours"] = locus.category_colours
                track_data["functions_coordinates"] = dict()
                for ff in feature_cateregories:
                    ff_features = [feature for feature in locus.features if feature.category == ff]
                    ff_coordinates = [[f.vis_prms["coordinates"]["start"], f.vis_prms["coordinates"]["end"]] for f in
                                      ff_features]
                    track_data["functions_coordinates"][ff] = ff_coordinates
            # Managing feature labels' positions:
            taken_label_coordinates = collections.defaultdict(list)
            if sum([fvd["show_label"] for fvd in track_data["features"]]) > 0:
                space_width = pdfmetrics.stringWidth(" ", self.prms.args["feature_label_font_face"],
                                                     track_data["f_label_font_size"])
                sorted_features = sorted(track_data["features"], key=lambda x: x["feature_width"] - x["label_width"],
                                         reverse=True)
                for fvd in [fvd for fvd in sorted_features if fvd["label_width"]]:
                    feature_center = fvd["coordinates"]["center"]
                    width_diff = fvd["label_width"] - fvd["feature_width"]
                    left_position = [fvd["coordinates"]["start"] - width_diff, fvd["coordinates"]["end"]]
                    centered_position = [feature_center - fvd["label_width"] / 2,
                                         feature_center + fvd["label_width"] / 2]
                    right_position = [fvd["coordinates"]["start"], fvd["coordinates"]["end"] + width_diff]
                    for pos in [left_position, centered_position, right_position]:
                        overlap = 0
                        if pos[0] < layout["loci_tracks_left_border"]:
                            overlap = layout["loci_tracks_left_border"] - pos[0]
                        if pos[1] > layout["loci_tracks_right_border"]:
                            overlap = layout["loci_tracks_right_border"] - pos[1]
                        if overlap:
                            pos[0] += overlap
                            pos[1] += overlap
                    label_position = centered_position
                    if width_diff > 0:
                        left_pos_overlap = [of for of in track_data["features"] if
                                            (of != fvd and left_position[0] < of["coordinates"]["end"] <
                                             fvd["coordinates"]["end"] and of["show_label"])]
                        right_pos_overlap = [of for of in track_data["features"] if
                                             (of != fvd and right_position[1] > of["coordinates"]["start"] >
                                              fvd["coordinates"]["start"] and of["show_label"])]
                        if left_pos_overlap and not right_pos_overlap:
                            min_n_distance = min([fvd["coordinates"]["start"] - nf["coordinates"]["end"]
                                                  for nf in left_pos_overlap])
                            if min_n_distance < space_width:
                                right_position[0] += space_width + min(0, min_n_distance)
                                right_position[1] += space_width + min(0, min_n_distance)
                            label_position = right_position
                        elif not left_pos_overlap and right_pos_overlap:
                            min_n_distance = min([nf["coordinates"]["start"] - fvd["coordinates"]["end"]
                                                  for nf in right_pos_overlap])
                            label_position = left_position
                            if min_n_distance < space_width and left_position[0] != layout["loci_tracks_left_border"]:
                                left_position[0] -= space_width
                                left_position[1] -= space_width
                            else:
                                label_position = centered_position
                        else:
                            label_position = centered_position
                    fvd["label_position"] = label_position
                    for label_row in range(0, len(track_data["features"])):
                        overlapped = False
                        for taken_coordinate in taken_label_coordinates[label_row]:
                            if taken_coordinate[0] <= label_position[0] <= taken_coordinate[1] or \
                                    taken_coordinate[0] <= label_position[1] <= taken_coordinate[1] or \
                                    label_position[0] <= taken_coordinate[0] <= label_position[1] or \
                                    label_position[0] <= taken_coordinate[1] <= label_position[1]:
                                overlapped = True
                        if not overlapped:
                            fvd["label_row"] = label_row
                            taken_label_coordinates[label_row].append(label_position)
                            break
                    fvd["label_y_bottom"] = track_data["feature_label_gap"] + \
                                            (fvd["label_row"] * track_data["f_label_height"]) + \
                                            (fvd["label_row"] * track_data["feature_label_gap"])
                for fvd in track_data["features"]:
                    if fvd["label_width"]:
                        if fvd["label_row"] > 0:
                            taken_middle_rows = [i for i in range(fvd["label_row"] - 1, -1, -1) if
                                                 any(taken_coordinate[0] <= fvd["coordinates"]["center"] <=
                                                     taken_coordinate[1]
                                                     for taken_coordinate in taken_label_coordinates[i])]
                            label_line_upper = fvd["label_y_bottom"] - track_data["feature_label_gap"] / 2
                            label_line_bottom = track_data["feature_label_gap"] / 2
                            label_line_coordinates = []
                            ll_start = label_line_bottom
                            for tmr in sorted(taken_middle_rows):
                                tmr_start = 0.5 * track_data["feature_label_gap"] + \
                                            (tmr * track_data["f_label_height"]) + \
                                            (tmr * track_data["feature_label_gap"])
                                tmr_end = tmr_start + track_data["f_label_height"] + \
                                          0.5 * track_data["feature_label_gap"]
                                label_line_coordinates.append([ll_start, tmr_start])
                                ll_start = tmr_end
                            label_line_coordinates.append([ll_start, label_line_upper])
                            fvd["label_line_coordinates"] = label_line_coordinates
            track_data["n_label_rows"] = sum([1 for k, v in taken_label_coordinates.items() if v])
            # Managing middle line indicating locus borders
            if self.prms.args["draw_middle_line"]:
                regions_for_middle_line = [[c["start"], c["end"]] for c in locus.coordinates]
                for ftc in features_taken_nt_coordinates:
                    ftcs, ftse = ftc
                    for added_region in regions_for_middle_line:
                        new_regions = []
                        to_remove = False
                        if added_region[0] <= ftcs <= added_region[1]:
                            to_remove = True
                            if ftcs > added_region[0]:
                                new_regions.append([added_region[0], ftcs - 1])
                        if added_region[0] <= ftse <= added_region[1]:
                            to_remove = True
                            if ftse < added_region[1]:
                                new_regions.append([ftse + 1, added_region[1]])
                        if to_remove:
                            regions_for_middle_line.remove(added_region)
                        regions_for_middle_line += new_regions
                middle_line_coordinates = [lovis4u.Methods.region_nt_to_x_transform(rml[0], rml[1], locus, layout)
                                           for rml in regions_for_middle_line]
                track_data["middle_line_coordinates"] = middle_line_coordinates
            track_data["proteome_size"] = len(locus.features)  # to change if we get other features
            # Managing individual x axis
            if self.prms.args["draw_individual_x_axis"]:
                track_data["x_axis_annotation"] = dict()
                self.prms.args["x_axis_ticks_labels_height"] = lovis4u.Methods.str_font_size_to_height(
                    self.prms.args["x_axis_ticks_labels_font_size"],
                    self.prms.args["x_axis_ticks_labels_font_face"]) / mm
                track_data["x_axis_annotation"]["label_size"] = self.prms.args["x_axis_ticks_labels_font_size"]
                axis_regions = []
                axis_tics_coordinates = []
                for coordinate in locus.coordinates:
                    axis_regions.append(lovis4u.Methods.region_nt_to_x_transform(coordinate["start"], coordinate["end"],
                                                                                 locus, layout))
                    current_tics_coordinates = [coordinate["start"], coordinate["end"]]
                    axis_tics_coordinates += current_tics_coordinates
                if len(axis_tics_coordinates) == 2:
                    nt_range = coordinate["end"] - coordinate["start"]
                    axis_tics_coordinates.append(coordinate["start"] + int(0.25 * nt_range))
                    axis_tics_coordinates.append(coordinate["start"] + int(0.5 * nt_range))
                    axis_tics_coordinates.append(coordinate["end"] - int(0.25 * nt_range))

                axis_tics_positions = list(map(lambda x: lovis4u.Methods.nt_to_x_transform(x, locus, layout, "center"),
                                               axis_tics_coordinates))
                sorted_tics = sorted(zip(axis_tics_positions, axis_tics_coordinates))
                sorted_positions, sorted_coordinates = zip(*sorted_tics)
                axis_tics_positions = list(sorted_positions)
                axis_tics_coordinates = list(sorted_coordinates)
                axis_tics_labels = list(map(str, axis_tics_coordinates))
                axis_tics_label_size = self.prms.args["x_axis_ticks_labels_font_size"]
                axis_tics_label_width = list(map(lambda x: pdfmetrics.stringWidth(x,
                                                                                  self.prms.args[
                                                                                      "x_axis_ticks_labels_font_face"],
                                                                                  axis_tics_label_size),
                                                 axis_tics_labels))

                space_width = pdfmetrics.stringWidth(" ", self.prms.args["x_axis_ticks_labels_font_face"],
                                                     axis_tics_label_size)
                tics_labels_coordinates = axis_tics_positions.copy()
                # Shitty, refactor later
                for t_i in range(len(axis_tics_coordinates)):
                    label_width = axis_tics_label_width[t_i]
                    tick_position = axis_tics_positions[t_i]
                    if t_i == 0:
                        tics_labels_coordinates[t_i] += label_width * 0.5
                    elif t_i != len(axis_tics_coordinates) - 1:
                        center_label_coordinates = [tick_position - 0.5 * label_width,
                                                    tick_position + 0.5 * label_width]
                        if center_label_coordinates[0] - space_width <= axis_tics_positions[t_i - 1]:
                            tics_labels_coordinates[t_i] += label_width * 0.5
                            if axis_tics_positions[t_i] - axis_tics_positions[t_i - 1] < space_width:
                                tics_labels_coordinates[t_i] += space_width
                        if center_label_coordinates[1] + space_width >= axis_tics_positions[t_i + 1]:
                            tics_labels_coordinates[t_i] -= label_width * 0.5
                            if axis_tics_positions[t_i + 1] - axis_tics_positions[t_i] < space_width:
                                tics_labels_coordinates[t_i] -= space_width
                    if tics_labels_coordinates[t_i] + label_width * 0.5 > layout["loci_tracks_right_border"]:
                        tics_labels_coordinates[t_i] -= ((tics_labels_coordinates[t_i] + label_width * 0.5) -
                                                         layout["loci_tracks_right_border"])
                track_data["x_axis_annotation"]["axis_regions"] = axis_regions
                track_data["x_axis_annotation"]["axis_tics_position"] = axis_tics_positions
                track_data["x_axis_annotation"]["axis_tics_labels"] = axis_tics_labels
                track_data["x_axis_annotation"]["tics_labels_coordinates"] = tics_labels_coordinates
                return None
        except Exception as error:
            raise lovis4u.Manager.lovis4uError("Unable to prepare Locus specific data.") from error

    def calculate_track_height(self) -> float:
        """Calculate track height to define layout.

        Returns:
            float: track height.

        """
        try:
            track_height = self.prms.args["feature_height"] * mm + \
                           (self.track_data["n_label_rows"] * self.track_data["f_label_height"] * \
                            (1 + self.prms.args["feature_label_gap"]))
            if self.prms.args["draw_individual_x_axis"]:
                track_height += (self.prms.args["feature_bottom_gap"] + self.prms.args["x_axis_ticks_height"] * 1.3 + \
                                 self.prms.args["x_axis_ticks_labels_height"]) * mm
            elif not self.prms.args["draw_individual_x_axis"] and self.track_data["functions_coordinates"]:
                track_height += (self.prms.args["feature_bottom_gap"] + \
                                 self.prms.args["category_annotation_line_width"]) * mm
            if self.prms.args["locus_label_position"] == "bottom":
                track_height += self.prms.args["locus_label_height"] + self.prms.args["feature_bottom_gap"] * mm

            self.track_data["track_height"] = track_height
            return track_height
        except Exception as error:
            raise lovis4u.Manager.lovis4uError("Unable to calculate a locus track height.") from error

    def create_track(self):
        """Initialise a LocusVis track object.

        Returns:
            lovis4u.Drawing.LocusVis: visualisation track.

        """
        try:
            return lovis4u.Drawing.LocusVis(self.layout, self.track_data, self.prms)
        except Exception as error:
            raise lovis4u.Manager.lovis4uError("Unable to create a locus track object.") from error


class ScaleLoader(Loader):
    """A LocusLoader object prepares data for a Scale track Drawing object.

    Attributes:
        prms (Parameters): Parameters' class object.
        layout (dict): Layout built by CanvasManager's define_layout() method.
        track_data (dict): Track specific data that will be sent to the Drawing module.

    """

    def __init__(self, parameters):
        """Create a ScaleLoader object.

        Arguments:
            parameters (Parameters): Parameters' class object that holds config and cmd arguments.

        """
        super().__init__(parameters)

    def prepare_track_specific_data(self, layout: dict) -> None:
        """Prepare ScaleLoader specific data.

        Attributes:
            layout (dict): Layout built by CanvasManager's define_layout() method.

        Returns:
            None

        """
        try:
            self.layout = layout
            total_nt_width = self.layout["total_nt_width"]
            raw_scale_line_nt_width = round(total_nt_width * self.prms.args["scale_line_relative_size"])
            raw_scale_line_nt_width_pow = int(math.log(raw_scale_line_nt_width, 10))
            scale_line_nt_width = round(raw_scale_line_nt_width // math.pow(10, raw_scale_line_nt_width_pow) *
                                        math.pow(10, raw_scale_line_nt_width_pow))
            track_data = dict()
            self.track_data = track_data
            track_data["scale_line_nt_width"] = scale_line_nt_width
            track_data["coordinates"] = [self.layout["loci_tracks_left_border"],
                                         self.layout["loci_tracks_left_border"] +
                                         layout["width_per_nt"] * scale_line_nt_width]
            track_data["scale_line_width"] = track_data["coordinates"][1] - track_data["coordinates"][0]
            track_data["scale_label"] = f"{scale_line_nt_width} nt"
            track_data["scale_line_label_font_size"] = self.prms.args["scale_line_label_font_size"]
            track_data["scale_line_label_height"] = lovis4u.Methods.str_font_size_to_height(
                track_data["scale_line_label_font_size"], self.prms.args["scale_line_label_font_face"])
            track_data["scale_line_label_width"] = pdfmetrics.stringWidth(track_data["scale_label"],
                                                                          self.prms.args["scale_line_label_font_face"],
                                                                          track_data["scale_line_label_font_size"])
            if (track_data["scale_line_width"] - track_data["scale_line_label_width"]) / \
                    track_data["scale_line_width"] > 0.1:
                track_data["style"] = "fancy"
                track_data["space_width"] = pdfmetrics.stringWidth(" ", self.prms.args["scale_line_label_font_face"],
                                                                   track_data["scale_line_label_font_size"])
            else:
                track_data["style"] = "default"
            return None
        except Exception as error:
            raise lovis4u.Manager.lovis4uError("Unable to prepare scale track specific data.") from error

    def calculate_track_height(self):
        """Calculate track height to define layout.

        Returns:
            float: track height.

        """
        try:
            if self.track_data["style"] == "fancy":
                track_height = self.track_data["scale_line_label_height"]
            else:
                track_height = (1.2 * self.prms.args["scale_line_tics_height"] +
                                self.track_data["scale_line_label_height"])
            self.track_data["track_height"] = track_height
            return track_height
        except Exception as error:
            raise lovis4u.Manager.lovis4uError("Unable to calculate a scale track height.") from error

    def create_track(self):
        """Initialise a ScaleVis track object.

        Returns:
            lovis4u.Drawing.LocusVis: visualisation track.

        """
        try:
            return lovis4u.Drawing.ScaleVis(self.layout, self.track_data, self.prms)
        except Exception as error:
            raise lovis4u.Manager.lovis4uError("Unable to create a scale track object.") from error


class CategoriesColorLegendLoader(Loader):
    """A CategoriesColorLegendLoader object prepares data for a categories colour legend track Drawing object.

    Attributes:
        prms (Parameters): Parameters' class object.
        layout (dict): Layout built by CanvasManager's define_layout() method.
        track_data (dict): Track specific data that will be sent to the Drawing module.

    """

    def __init__(self, parameters):
        """Create a CategoriesColorLegendLoader object.

        Arguments:
            parameters (Parameters): Parameters' class object that holds config and cmd arguments.

        """
        super().__init__(parameters)

    def prepare_track_specific_data(self, layout: dict, loci) -> None:
        """Prepare ScaleLoader specific data.

        Attributes:
            layout (dict): Layout built by CanvasManager's define_layout() method.
            loci (lovis4u.DataProcessing.Loci): Loci object with information about sequences and features.

        Returns:
            None

        """
        try:
            self.layout = layout
            self.track_data = dict()
            label_height = lovis4u.Methods.str_font_size_to_height(
                self.prms.args["colour_legend_label_font_size"], self.prms.args["colour_legend_font_face"])
            line_width = self.prms.args["colour_legend_line_width"] * mm
            line_gap = 0.4 * label_height
            self.track_data["line_width"] = line_width
            self.track_data["colour_legend_label_size"] = self.prms.args["colour_legend_label_font_size"]
            left_border = self.layout["loci_tracks_left_border"]
            right_border = self.layout["loci_tracks_right_border"]
            x_gap = pdfmetrics.stringWidth(" " * 5, self.prms.args["colour_legend_font_face"],
                                           self.track_data["colour_legend_label_size"])
            y_gap = 0.5 * label_height
            colour_dict = dict()
            for locus in loci.loci:
                colour_dict.update(locus.category_colours)
            self.track_data["labels"] = []
            current_x = left_border
            n_of_rows = 0
            current_y = - (label_height + line_gap + line_width)
            for label, colour in colour_dict.items():
                label_dict = dict()
                label_width = pdfmetrics.stringWidth(label, self.prms.args["colour_legend_font_face"],
                                                     self.track_data["colour_legend_label_size"])
                label_end = current_x + label_width
                if label_end > right_border:
                    current_x = left_border
                    current_y -= (label_height + line_gap + line_width + y_gap)
                    n_of_rows += 1
                label_dict["label_x"] = current_x
                label_dict["label_width"] = label_width
                label_dict["relative_y"] = current_y
                label_dict["relative_y_text"] = current_y + line_gap + line_width
                label_dict["colour"] = colour
                label_dict["label"] = label
                self.track_data["labels"].append(label_dict)
                current_x = label_dict["label_x"] + label_width + x_gap
            return None
        except Exception as error:
            raise lovis4u.Manager.lovis4uError("Unable to prepare categories colour legend track specific data.") \
                from error

    def calculate_track_height(self):
        """Calculate track height to define layout.

        Returns:
            float: track height.

        """
        try:
            if self.track_data["labels"]:
                min_relative_y = min([i["relative_y"] for i in self.track_data["labels"]])
                track_height = abs(min_relative_y)
            else:
                track_height = 0
            self.track_data["track_height"] = track_height
            return track_height
        except Exception as error:
            raise lovis4u.Manager.lovis4uError("Unable to calculate a categories colour legend track height.") \
                from error

    def create_track(self):
        """Initialise a ColorLegendVis track object.

        Returns:
            lovis4u.Drawing.LocusVis | None: visualisation track.

        """
        try:
            if self.track_data["labels"]:
                return lovis4u.Drawing.ColorLegendVis(self.layout, self.track_data, self.prms)
            else:
                if self.prms.args["verbose"]:
                    print("○ Warning message: Category colours legend track cannot be created since there is no "
                          "categories.", file=sys.stdout)
                return None
        except Exception as error:
            raise lovis4u.Manager.lovis4uError("Unable to create a categories colour legend track track object.") \
                from error


class Canvas:
    """An Image object holds canvas;.

    Attributes:
        canvas (reportlab.pdfgen.canvas.Canvas): canvas object of the reportlab library.

    """

    def __init__(self, filename: str, width: float, height: float, parameters):
        """Create a Canvas object.

        Arguments:
            filename (str): path and name of the output pdf.
            width (float): width of the canvas.
            height (float): height of the pdf.

        """
        self.prms = parameters
        self.canvas = reportlab.pdfgen.canvas.Canvas(filename, pagesize=(width, height))
        self.canvas.setTitle("lovis4u output")
        self.canvas.setSubject("🎨")
        self.canvas.setCreator("lovis4u | The Atkinson Lab 4U")

    def save(self) -> None:
        """Save canvas as a pdf file.

        Returns:
            None

        """
        if not os.path.exists(self.prms.args["output_dir"]):
            os.mkdir(self.prms.args["output_dir"])
        self.canvas.save()
        return None
