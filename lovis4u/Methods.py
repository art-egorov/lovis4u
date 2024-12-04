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

import pandas as pd
import pyhmmer.plan7
import pyhmmer.easel
import pyhmmer

import progress.bar

import requests
import tarfile
import lovis4u.Manager


def adjust_paths(to: str) -> None:
    """Change paths in the internal config files for linux or mac.

    Arguments:
        to (str): mac | linux

    Returns:
        None

    """
    internal_dir = os.path.join(os.path.dirname(__file__), "lovis4u_data")
    config_files = ["standard.cfg", "A4L.cfg", "A4p1.cfg", "A4p2.cfg"]
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


def str_font_size_to_height(font_size: float, font_type: str) -> float:
    """Transform string font size to height.

    Arguments:
        font_size (float): font_size.
        font_type (str): font type (see config file; at this moment only regular is available).

    Returns:
        float:  height of the string.

    """

    face = reportlab.pdfbase.pdfmetrics.getFont(font_type).face
    height = font_size * (face.ascent - face.descent) / (1000 * 1.38)
    return height


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
            try:
                x_coordinates["orient"] = feature_strand * coordinate["strand"]
            except:
                x_coordinates["orient"] = 1
            x_coordinates["lout"] = left_out
            x_coordinates["rout"] = right_out
            break
    return x_coordinates


def run_pyhmmer(query_fasta: str, query_size: int, prms: lovis4u.Manager.Parameters) -> pd.DataFrame:
    """Run pyhmmer hmmscan for a set of query proteins

    Arguments:
        query_fasta (str): Path to a query fasta file.
        query_size (int): Number of query proteins.
        prms (LoVis4u.Manager.Parameters): Parameters class object that holds all arguments.

    Returns:
        pd.DataFrame: Table with hmmscan search results.

    """
    with pyhmmer.easel.SequenceFile(query_fasta, digital=True) as seqs_file:
        query_proteins = seqs_file.read_block()
    num_of_query_proteins = query_size

    databases_cname = prms.args["hmm_config_names"]
    databases_short_names = prms.args["database_names"]

    databases_names = {"hmm_defence_df": "DefenceFinder and CasFinder databases",
                       "hmm_defence_padloc": "PADLOC database", "hmm_virulence": "virulence factor database (VFDB)",
                       "hmm_anti_defence": "anti-prokaryotic immune systems database (dbAPIS)",
                       "hmm_amr": "AMRFinderPlus Database"}
    hmmscan_output_folder = os.path.join(prms.args["output_dir"], "hmmscan")
    alignment_table_rows = []
    for db_ind, db_name in enumerate(databases_cname):
        if prms.args["defence_models"] == "DefenseFinder" and db_name == "hmm_defence_padloc":
            continue
        if prms.args["defence_models"] == "PADLOC" and db_name == "hmm_defence_df":
            continue
        db_alignment_table_rows = []
        db_shortname = databases_short_names[db_ind]
        db_path = prms.args[db_name]
        if databases_cname[db_ind] in databases_names.keys():
            db_full_name = databases_names[databases_cname[db_ind]]
        else:
            db_full_name = db_shortname
        if not os.path.exists(db_path):
            print(db_path)
            print(f"  ⦿ Database {db_full_name} was not found.", file=sys.stdout)
            continue
        hmm_files = [fp for fp in os.listdir(db_path) if os.path.splitext(fp)[1].lower() == ".hmm" and fp[0] != "."]
        hmms = []
        for hmm_file in hmm_files:
            hmms.append(pyhmmer.plan7.HMMFile(os.path.join(db_path, hmm_file)).read())
        if prms.args["verbose"]:
            print(f"  ⦿ Running pyhmmer hmmscan versus {db_full_name}...", file=sys.stdout)
            bar = progress.bar.FillingCirclesBar("   ", max=num_of_query_proteins, suffix="%(index)d/%(max)d")
        for hits in pyhmmer.hmmscan(query_proteins, hmms, E=1e-3, cpus=0):
            if prms.args["verbose"]:
                bar.next()
            for hit in hits:
                if hit.included:
                    for domain in hit.domains.reported:
                        if domain.i_evalue < prms.args["hmmscan_evalue"]:
                            alignment = domain.alignment
                            hit_name = hit.name.decode()
                            hit_description = hit.description
                            if hit.description:
                                hit_description = hit_description.decode()
                                if hit_description == "NA":
                                    hit_description = ""
                            else:
                                hit_description = ""

                            if hit_name != hit_description and hit_name not in hit_description and hit_description:
                                hname = f"{hit_name} {hit_description}"
                            elif hit_description:
                                hname = hit_description
                            else:
                                hname = hit_name
                            alignment_row = dict(query=alignment.target_name.decode(),
                                                 target_db=db_shortname, target=hname, t_name=hit_name,
                                                 t_description=hit_description,
                                                 hit_evalue=hit.evalue, di_evalue=domain.i_evalue,
                                                 q_from=alignment.target_from, q_to=alignment.target_to,
                                                 qlen=alignment.target_length, t_from=alignment.hmm_from,
                                                 t_to=alignment.hmm_to, tlen=alignment.hmm_length)
                            alignment_row["q_cov"] = round((alignment_row["q_to"] - alignment_row["q_from"]) / \
                                                           alignment_row["qlen"], 2)
                            alignment_row["t_cov"] = round((alignment_row["t_to"] - alignment_row["t_from"]) / \
                                                           alignment_row["tlen"], 2)
                            if alignment_row["q_cov"] >= prms.args["hmmscan_query_coverage_cutoff"] and \
                                    alignment_row["t_cov"] >= prms.args["hmmscan_hmm_coverage_cutoff"]:
                                db_alignment_table_rows.append(alignment_row)
        if prms.args["verbose"]:
            bar.finish()
        alignment_table_rows += db_alignment_table_rows
        db_alignment_table = pd.DataFrame(db_alignment_table_rows)
        if not db_alignment_table.empty:
            db_alignment_table = db_alignment_table.sort_values(by="hit_evalue", ascending=True)
            db_alignment_table = db_alignment_table.drop_duplicates(subset="query", keep="first").set_index("query")
            db_alignment_table.to_csv(os.path.join(hmmscan_output_folder, f"{db_shortname}.tsv"), sep="\t",
                                      index_label="query")
        n_hits = len(db_alignment_table.index)
        if prms.args["verbose"]:
            print(f"    Number of hits: {n_hits}", file=sys.stdout)
    alignment_table = pd.DataFrame(alignment_table_rows)
    if not alignment_table.empty:
        alignment_table = alignment_table.sort_values(by="hit_evalue", ascending=True)
        alignment_table = alignment_table.drop_duplicates(subset="query", keep="first").set_index("query")

    return alignment_table

def download_file_with_progress(url: str, local_folder: str) -> None:
    """Function for downloading a particular file from a web server.

    Arguments:
        url (str): Link to the file.
        local_folder (str): Path to a folder where file will be saved.

    Returns:
        None

    """
    try:
        response = requests.head(url)
        file_size = int(response.headers.get('content-length', 0))
        # Extract the original file name from the URL
        file_name = os.path.basename(url)
        local_path = os.path.join(local_folder, file_name)
        # Stream the file download and show progress bar
        with requests.get(url, stream=True) as r, open(local_path, 'wb') as f:
            bar = progress.bar.FillingCirclesBar(" ", max=file_size // 8192, suffix='%(percent)d%%')
            downloaded_size = 0
            for chunk in r.iter_content(chunk_size=8192):
                if chunk:
                    downloaded_size += len(chunk)
                    f.write(chunk)
                    bar.next()
            bar.finish()
        # Verify that the file was fully downloaded
        if downloaded_size != file_size:
            raise lovis4u.Manager.lovis4uError(f"Downloaded file size ({downloaded_size} bytes) does not match "
                                               f"expected size ({file_size} bytes).")
        print(f"⦿ File was saved to {local_path}")
        if file_name.endswith('.tar.gz'):
            with tarfile.open(local_path, 'r:gz') as tar:
                tar.extractall(path=local_folder)
            print(f"⦿ Folder was successfully unarchived")
            os.remove(local_path)
    except Exception as error:
        raise lovis4u.Manager.lovis4uError(f"Unable to get file from the {url}.") from error


def get_HMM_models(parameters) -> None:
    """Download HMM models

    Returns:
        None

    """
    try:
        url =  parameters["hmm_models_url"]
        internal_dir = os.path.join(os.path.dirname(__file__), "lovis4u_data")
        if os.path.exists(os.path.join(internal_dir, "HMMs")):
            print(f"○ HMMs folder already exists and will be rewritten...", file=sys.stdout)
        # Add checking if it's already downloaded
        print(f"○ Downloading HMM models...\n"
              f"  Source: {url}", file=sys.stdout)
        download_file_with_progress(url, internal_dir)
        return None
    except Exception as error:
        raise lovis4u.Manager.lovis4uError(f"Unable to download HMM models.") from error
