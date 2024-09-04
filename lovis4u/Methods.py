"""
This module provides some methods (e.g. colors tranformation, data copying) used by the tool.
"""
import os
import re
import shutil
import sys

import Bio.SeqIO
import colorsys

from reportlab.lib.units import cm, mm
import reportlab.pdfbase.pdfmetrics
import reportlab.pdfbase.ttfonts

import matplotlib.colors

import lovis4u.Manager


def adjust_paths(to: str) -> None:
    """Change paths in the internal config files for linux or mac.

    Arguments:
        to (str): mac | linux

    Returns:
        None

    """
    internal_dir = os.path.join(os.path.dirname(__file__), "lovis4u_data")
    config_files = ["standard.cfg"]
    for config_file in config_files:
        config_file_path = os.path.join(internal_dir, config_file)
        with open(config_file_path, "r+") as config:
            if to == "linux":
                if not os.path.exists(os.path.join(internal_dir, "bin/mmseqs_linux")):
                    os.system(f"unzip -q -d {os.path.join(internal_dir, 'bin/')} "
                              f"{os.path.join(internal_dir, 'bin/mmseqs_linux.zip')}")
                config_txt = re.sub(r"mmseqs_mac/bin/mmseqs", "mmseqs_linux/bin/mmseqs", config.read())
            else:
                config_txt = re.sub(r"mmseqs_linux/bin/mmseqs", "mmseqs_mac/bin/mmseqs", config.read())
            config.seek(0)
            config.truncate()
            config.write(config_txt)
    print(f"⚙ mmseqs path was adjusted to {to}")
    return None


def copy_package_data() -> None:
    """Copy the lovis4u package data folder to your current dir.

    Returns:
        None

    """
    try:
        users_dir = os.path.join(os.getcwd(), "lovis4u_data")
        internal_dir = os.path.join(os.path.dirname(__file__), "lovis4u_data")
        if os.path.exists(users_dir):
            raise lovis4u.Manager.lovis4uError("lovis4u_data folder already exists.")
        shutil.copytree(internal_dir, users_dir, ignore=shutil.ignore_patterns("help*", ".*"))
        print("⚙ lovis4u_data folder was copied to the current working directory", file=sys.stdout)
        return None
    except Exception as error:
        raise lovis4u.Manager.lovis4uError(f"Unable to copy lovis4u folder in your working dir.") from error


def scale_lightness(hex_c: str, scale_l: float) -> str:
    """Helper function to get darker version of input colour

    Arguments:
        hex_c (str): input HEX colour
        scale_l (float): scale of lightness

    Returns:
        str: new HEX colour

    """
    rgb = matplotlib.colors.hex2color(hex_c)
    h, l, s = colorsys.rgb_to_hls(*rgb)
    hex_c = matplotlib.colors.rgb2hex(colorsys.hls_to_rgb(h, min(1, l * scale_l), s=s))
    return hex_c


def get_colour(name: str, parameters: dict) -> str:
    """Get HEX colour by its name

    Arguments:
        name (str): name of a colour.
        parameters (dict): Parameters' object dict.

    Returns:
        str: HEX colour.

    """
    hex_c = parameters.args["palette"][parameters.args[name]]
    return hex_c


def get_colour_rgba(name: str, parameters: dict) -> tuple:
    """Get rgba colour by its name

    Arguments:
        name (str): name of a colour.
        parameters (dict): Parameters' object dict.

    Returns:
        tuple: RGBA colour

    """
    return *matplotlib.colors.hex2color(get_colour(name, parameters)), parameters.args[f"{name}_alpha"]


def str_height_to_size(height: float, font_type: str) -> float:
    """Transform string height to the font size.

    Arguments:
        height (float): available height of the string.
        font_type (str): font type (see config file; at this moment only regular is available).

    Returns:
        float: font size defined by height.

    """

    face = reportlab.pdfbase.pdfmetrics.getFont(font_type).face
    font_size = (1000 * 1.38 * height) / (face.ascent - face.descent)
    return font_size


def update_path_extension(path: str, new_extension: str) -> str:
    """Get path basename and replace its extension

    Arguments:
        path (str): path to a file
        new_extension (str): new extension

    Returns:
        str: basename of a file with new extension.

    """
    updated_filename = f"{os.path.splitext(os.path.basename(path))[0]}.{new_extension}"
    return updated_filename


def nt_to_x_transform(nt: int, locus, layout: dict, mode: str) -> float:
    """Transform nucleotide coordinate to x page coordinate.

    Arguments:
        nt (int): Nucleotide coordinate.
        locus (lovis4u.DataProcessing.Locus): Corresponding locus object.
        layout (dict): Layout of the canvas.
        mode (str): Mode whether coordinate should be centered to the nt or on the side.

    Returns:
        float: Corresponding page coordinate.

    """
    passed_x = layout["loci_tracks_left_border"]
    for c_i in range(len(locus.coordinates)):
        coordinate = locus.coordinates[c_i]
        coordinate_region_width = (coordinate["end"] - coordinate["start"] + 1) * layout["width_per_nt"]
        if coordinate["start"] <= nt <= coordinate["end"]:
            if coordinate["strand"] == 1:
                relative_nt = nt - coordinate["start"]
            else:
                relative_nt = coordinate["end"] - nt
            if mode == "start":
                pass
            elif mode == "center":
                relative_nt += 0.5
            elif mode == "end":
                relative_nt += 1
            relative_x = relative_nt * layout["width_per_nt"]
            global_x = passed_x + relative_x
            break
        passed_x += coordinate_region_width + layout["x_gap_between_regions"]
        if c_i < len(locus.coordinates) - 1:
            if locus.circular:
                if coordinate["strand"] == locus.coordinates[c_i + 1]["strand"] == 1:
                    if coordinate["end"] == locus.length and locus.coordinates[c_i + 1]["start"] == 1:
                        passed_x -= layout["x_gap_between_regions"]
                elif coordinate["strand"] == locus.coordinates[c_i + 1]["strand"] == -1:
                    if coordinate["start"] == 1 and locus.coordinates[c_i + 1]["end"] == locus.length:
                        passed_x -= layout["x_gap_between_regions"]
    return global_x


def region_nt_to_x_transform(nt_start: int, nt_end: int, locus, layout: dict) -> dict:
    """Transform region coordinates to x page coordinates.

    Arguments:
        nt_start (int): 1-based nucleotide start coordinate.
        nt_end (int): 1-based nucleotide end coordinate.
        locus (lovis4u.DataProcessing.Locus): Corresponding locus object.
        layout (dict): Layout of the canvas.

    Returns:
        dict: Dictionary with corresponding start and end page coordinates.

    """
    for coordinate in locus.coordinates:
        if coordinate["start"] <= nt_start <= coordinate["end"] and coordinate["start"] <= nt_end <= coordinate["end"]:
            if coordinate["strand"] == -1:
                nt_start, nt_end = nt_end, nt_start
            x_start = nt_to_x_transform(nt_start, locus, layout, "start")
            x_end = nt_to_x_transform(nt_end, locus, layout, "end")
            break
    return dict(start=x_start, end=x_end)


def feature_nt_to_x_transform(nt_start: int, nt_end: int, feature_strand: int, locus, layout: dict) -> dict:
    """Transform feature coordinates to x page coordinates.

    Arguments:
        nt_start (int): 1-based nucleotide start coordinate.
        nt_end (int): 1-based nucleotide end coordinate.
        feature_strand (int): 1 | -1 corresponding to plus or minus strand, respectively.
        locus (lovis4u.DataProcessing.Locus): Corresponding locus object.
        layout (dict): Layout of the canvas.

    Returns:
        dict: Dictionary with corresponding start and end page coordinates.

    """
    left_out, right_out = False, False
    for coordinate in locus.coordinates:
        if coordinate["start"] <= nt_start <= coordinate["end"] or coordinate["start"] <= nt_end <= coordinate["end"]:
            if nt_start < coordinate["start"]:
                if coordinate["strand"] == 1:
                    left_out = True
                else:
                    right_out = True
                nt_start = coordinate["start"]
            if nt_end > coordinate["end"]:
                if coordinate["strand"] == 1:
                    right_out = True
                else:
                    left_out = True
                nt_end = coordinate["end"]
            x_coordinates = region_nt_to_x_transform(nt_start, nt_end, locus, layout)
            x_coordinates["center"] = (x_coordinates["start"] + x_coordinates["end"]) / 2
            x_coordinates["orient"] = feature_strand * coordinate["strand"]
            x_coordinates["lout"] = left_out
            x_coordinates["rout"] = right_out
            break
    return x_coordinates
