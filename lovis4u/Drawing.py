"""
This module provides visualisation of loci annotation.
"""
from reportlab.lib.units import cm, mm
import reportlab.pdfbase.pdfmetrics
import reportlab.pdfbase.ttfonts
import reportlab.pdfgen.canvas
import reportlab.rl_config
import reportlab.pdfgen

reportlab.rl_config.warnOnMissingFontGlyphs = 0

import matplotlib.colors
import numpy as np
import copy
import lovis4u.Methods


class Track:
    """Parent class for visualisation Tracks.

    Attributes:
        layout (dict): Layout built by CanvasManager's define_layout() method.
        track_data (dict): a dictionary with prepared track specific data.
        prms (lovis4u.Manager.Parameters): Parameters' class object.

    """

    def __init__(self, layout: dict, track_data: dict, parameters: lovis4u.Manager.Parameters):
        """Parent's constructor for creating a Track object.

        Arguments:
            layout (dict): Layout built by CanvasManager's define_layout() method.
            track_data (dict): a dictionary with prepared track specific data.
            parameters (lovis4u.Manager.Parameters): Parameters' class object.

        """
        self.layout = layout
        self.track_data = track_data
        self.prms = parameters

    def draw(self, canvas: reportlab.pdfgen.canvas.Canvas) -> None:
        """Empy parent's method for track drawing.

        Arguments:
            canvas (reportlab.pdfgen.canvas.Canvas): a canvas object.

        Returns:
            None

        """
        pass


class CrossTrack:
    """Parent class for Cross-Tracks visualisation.

    Attributes:
        layout (dict): Layout built by CanvasManager's define_layout() method.
        tracks (dict): List with track objects participated in CrossTrack.
        prms (lovis4u.Manager.Parameters): Parameters' class object.

    """

    def __init__(self, layout, tracks, parameters):
        """Parent's constructor for creating a CrossTrack object.

        Attributes:
            layout (dict): Layout built by CanvasManager's define_layout() method.
            tracks (dict): List with track objects participated in CrossTrack.
            parameters (lovis4u.Manager.Parameters): Parameters' class object.

        """
        self.layout = layout
        self.tracks = tracks
        self.prms = parameters

    def draw(self, canvas: reportlab.pdfgen.canvas.Canvas):
        """Empy parent's method for cross track drawing.

        Arguments:
            canvas (reportlab.pdfgen.canvas.Canvas): a canvas object.

        Returns:
            None

        """
        pass


class HomologyTrack(CrossTrack):
    """Track that handle visualisation of homology lines between homologous features on neighbours' loci.

    Attributes:
        layout (dict): Layout built by CanvasManager's define_layout() method.
        tracks (dict): List with track objects participated in CrossTrack.
        prms (lovis4u.Manager.Parameters): Parameters' class object.

    """

    def __init__(self, layout, tracks, parameters):
        """Create a HomologyTrack.

        Attributes:
            layout (dict): Layout built by CanvasManager's define_layout() method.
            tracks (dict): List with track objects participated in CrossTrack.
            parameters (lovis4u.Manager.Parameters): Parameters' class object.

        """
        super().__init__(layout, tracks, parameters)

    def draw(self, canvas: reportlab.pdfgen.canvas.Canvas) -> None:
        """Draw HomologyTrack on canvas.

        Arguments:
            canvas (reportlab.pdfgen.canvas.Canvas): a canvas object.

        Returns:
            None

        """
        try:
            for track in self.tracks:
                track.track_data["y_top"] = self.layout["figure_height"] - track.layout["inverse_y_coordinate"]
                track.track_data["feature_upper"] = track.track_data["y_top"] - \
                                                    (track.track_data["n_label_rows"] * track.track_data[
                                                        "f_label_height"] *
                                                     (1 + track.prms.args["feature_label_gap"]))
            num_of_loci_tracks = len(self.tracks)
            for ti in range(num_of_loci_tracks - 1):
                current_track = self.tracks[ti]
                current_track_features = current_track.track_data["features"]
                next_track = self.tracks[ti + 1]
                next_track_features = next_track.track_data["features"]
                for ctf in current_track_features:
                    ctf_group = ctf["group"]
                    next_track_same_group_features = [i for i in next_track_features if i["group"] == ctf_group]
                    for ntf in next_track_same_group_features:
                        cty_u = current_track.track_data["feature_upper"]
                        cty_c = current_track.track_data["feature_upper"] - 0.5 * self.prms.args["feature_height"] * mm
                        cty_b = current_track.track_data["feature_upper"] - self.prms.args["feature_height"] * mm
                        nty_u = next_track.track_data["feature_upper"]
                        nty_c = next_track.track_data["feature_upper"] - 0.5 * self.prms.args["feature_height"] * mm
                        nty_b = next_track.track_data["feature_upper"] - self.prms.args["feature_height"] * mm
                        ct_arrow_len = min(
                            self.prms.args["feature_height"] * mm * self.prms.args["feature_arrow_length"],
                            (ctf["coordinates"]["end"] - ctf["coordinates"]["start"]))
                        nt_arrow_len = min(
                            self.prms.args["feature_height"] * mm * self.prms.args["feature_arrow_length"],
                            (ntf["coordinates"]["end"] - ntf["coordinates"]["start"]))
                        canvas.setLineCap(0)
                        canvas.setLineJoin(1)
                        canvas.setLineWidth(self.prms.args["homology_line_width"])
                        canvas.setStrokeColorRGB(*lovis4u.Methods.get_colour_rgba("homology_stroke_colour", self.prms))
                        canvas.setFillColorRGB(*lovis4u.Methods.get_colour_rgba("homology_fill_colour", self.prms))
                        p = canvas.beginPath()

                        if ctf["coordinates"]["orient"] == 1:
                            cts, cte = ctf["coordinates"]["start"], ctf["coordinates"]["end"]
                            p.moveTo(ctf["coordinates"]["start"], cty_b)
                            if not ctf["coordinates"]["rout"]:
                                p.lineTo(ctf["coordinates"]["end"] - ct_arrow_len, cty_b)
                                p.lineTo(ctf["coordinates"]["end"], cty_c)
                            else:
                                p.lineTo(ctf["coordinates"]["end"], cty_b)
                        elif ctf["coordinates"]["orient"] == -1:
                            cts, cte = ctf["coordinates"]["end"], ctf["coordinates"]["start"]
                            p.moveTo(ctf["coordinates"]["end"], cty_b)
                            if not ctf["coordinates"]["lout"]:
                                p.lineTo(ctf["coordinates"]["start"] + ct_arrow_len, cty_b)
                                p.lineTo(ctf["coordinates"]["start"], cty_c)
                            else:
                                p.lineTo(ctf["coordinates"]["start"], cty_b)
                        if ntf["coordinates"]["orient"] == 1:
                            nts, nte = ntf["coordinates"]["end"], ntf["coordinates"]["start"]
                            if not ntf["coordinates"]["rout"]:
                                p.lineTo(ntf["coordinates"]["end"], nty_c)
                                p.lineTo(ntf["coordinates"]["end"] - nt_arrow_len, nty_u)
                            else:
                                p.lineTo(ntf["coordinates"]["end"], nty_u)
                            p.lineTo(ntf["coordinates"]["start"], nty_u)
                            if nte >= cts:
                                p.lineTo(ntf["coordinates"]["start"], nty_c)
                        elif ntf["coordinates"]["orient"] == -1:
                            nts, nte = ntf["coordinates"]["start"], ntf["coordinates"]["end"]
                            if not ntf["coordinates"]["lout"]:
                                p.lineTo(ntf["coordinates"]["start"], nty_c)
                                p.lineTo(ntf["coordinates"]["start"] + nt_arrow_len, nty_u)
                            else:
                                p.lineTo(ntf["coordinates"]["start"], nty_u)
                            p.lineTo(ntf["coordinates"]["end"], nty_u)
                            if nte <= cts:
                                p.lineTo(ntf["coordinates"]["end"], nty_c)
                        if (nte <= cts and ctf["coordinates"]["orient"] == 1) or (
                                nte >= cts and ctf["coordinates"]["orient"] == -1):
                            p.lineTo(cts, cty_c)
                        p.lineTo(cts, cty_b)
                        canvas.drawPath(p, stroke=1, fill=1)
        except Exception as error:
            raise lovis4u.Manager.lovis4uError("Unable to draw homology line track.") from error


class LocusVis(Track):
    """LocusVis track object that handles each locus visualisation.

    Attributes:
        layout (dict): Layout built by CanvasManager's define_layout() method.
        track_data (dict): a dictionary with prepared track specific data.
        prms (lovis4u.Manager.Parameters): Parameters' class object.

    """

    def __init__(self, layout: dict, track_data: dict, parameters):
        """Create a LocusVis object.

        Arguments:
            layout (dict): Layout built by CanvasManager's define_layout() method.
            track_data (dict): a dictionary with prepared track specific data.
            parameters (lovis4u.Manager.Parameters): Parameters' class object.


        """
        super().__init__(layout, track_data, parameters)

    def draw(self, canvas: reportlab.pdfgen.canvas.Canvas) -> None:
        """Draw a LocusVis track.

        Arguments:
            canvas (reportlab.pdfgen.canvas.Canvas): a canvas object.

        Returns:
            None

        """
        try:
            y_track_bottom = self.layout["current_y_coordinate"]  - self.track_data["track_height"]
            feature_height = self.prms.args["feature_height"] * mm
            y_feature_upper = self.layout["current_y_coordinate"] - (self.track_data["n_label_rows"] *
                                                                     self.track_data["f_label_height"] *
                                                                     (1 + self.prms.args["feature_label_gap"]))
            y_feature_bottom = y_feature_upper - feature_height
            y_feature_center = y_feature_upper - feature_height * 0.5

            # Sequence label
            canvas.setFillColorRGB(*lovis4u.Methods.get_colour_rgba("locus_label_colour", self.prms))
            canvas.setFont(self.prms.args["locus_label_description_font_face"], self.prms.args["locus_label_font_size"])
            if self.prms.args["locus_label_position"] == "left":
                if self.prms.args["locus_label_style"] == "full" and self.track_data["locus_description"]:
                    label_bottom = y_feature_upper - self.prms.args["locus_label_height"]
                    canvas.drawRightString(self.layout["locus_label_right_border"], label_bottom,
                                           self.track_data["locus_description"])
                    label_bottom = y_feature_bottom
                else:
                    label_bottom = y_feature_bottom + (feature_height - self.prms.args["locus_label_height"]) * 0.5
                    label = self.track_data["locus_description"]
                if self.prms.args["locus_label_style"] != "description":
                    label = self.track_data["locus_id"]
                    canvas.setFont(self.prms.args["locus_label_id_font_face"], self.prms.args["locus_label_font_size"])
                canvas.drawRightString(self.layout["locus_label_right_border"], label_bottom, label)
            elif self.prms.args["locus_label_position"] == "bottom":
                label_bottom = y_track_bottom
                current_left = self.layout["locus_label_left_border"]
                if self.prms.args["locus_label_style"] == "full" and self.track_data["locus_description"]:
                    canvas.drawString(current_left, label_bottom, self.track_data["locus_description"])
                    current_left += self.track_data["locus_description_width"] + self.track_data["two_space_width"]
                canvas.setFont(self.prms.args["locus_label_id_font_face"], self.prms.args["locus_label_font_size"])
                canvas.drawString(current_left, label_bottom, self.track_data["locus_id"])
                current_left += self.track_data["locus_id_width"] + self.track_data["two_space_width"]
                canvas.drawString(current_left, label_bottom, self.track_data["text_coordinates"])

            # Middle line
            if self.prms.args["draw_middle_line"]:
                canvas.setLineWidth(self.prms.args["x_axis_line_width"])
                canvas.setStrokeColorRGB(*lovis4u.Methods.get_colour_rgba("x_axis_line_colour", self.prms))
                canvas.setLineCap(0)
                canvas.setFillColorRGB(*lovis4u.Methods.get_colour_rgba("x_axis_line_colour", self.prms))
                p = canvas.beginPath()
                for md_line_coordinates in self.track_data["middle_line_coordinates"]:
                    p.moveTo(md_line_coordinates["start"], y_feature_center)
                    p.lineTo(md_line_coordinates["end"], y_feature_center)
                canvas.drawPath(p, stroke=1, fill=0)

            # Category annotation
            if self.track_data["functions_coordinates"]:
                for feature_function, ff_region in self.track_data["functions_coordinates"].items():
                    feature_colour = self.track_data["category_colours"][feature_function]
                    canvas.setFillColorRGB(*matplotlib.colors.hex2color(feature_colour),
                                           self.prms.args["category_annotation_alpha"])
                    canvas.setLineJoin(1)
                    y_upper_sausage = y_feature_bottom - self.prms.args["feature_bottom_gap"] * mm
                    y_bottom_sausage = y_upper_sausage - self.prms.args["category_annotation_line_width"] * mm
                    for ffr in ff_region:
                        p = canvas.beginPath()
                        p.moveTo(ffr[0], y_bottom_sausage)
                        p.lineTo(ffr[0], y_upper_sausage)
                        p.lineTo(ffr[1], y_upper_sausage)
                        p.lineTo(ffr[1], y_bottom_sausage)
                        p.lineTo(ffr[0], y_bottom_sausage)
                        p.close()
                        canvas.drawPath(p, stroke=0, fill=1)
            # Features
            canvas.setLineCap(0)
            canvas.setLineJoin(1)
            canvas.setLineWidth(self.prms.args["feature_stroke_width"])
            if self.track_data["clean_features_coordinates"]:
                for f_data in self.track_data["features"]:
                    f_data_copy = copy.deepcopy(f_data)
                    f_data_copy["stroke_colour"] = None
                    f_data_copy["fill_colour"] = (*matplotlib.colors.hex2color(self.prms.args["palette"]["white"]), 1)
                    self.__plot_feature(canvas, f_data_copy, y_center=y_feature_center, height=feature_height)
            for f_data in self.track_data["features"]:
                f_data["stroke_colour"] = *matplotlib.colors.hex2color(f_data["stroke_colour"]), self.prms.args[
                    "feature_stroke_colour_alpha"]
                if f_data["fill_colour"]:
                    f_data["fill_colour"] = *matplotlib.colors.hex2color(f_data["fill_colour"]), self.prms.args[
                        "feature_fill_colour_alpha"]
                else:
                    f_data["fill_colour"] = None
                self.__plot_feature(canvas, f_data, y_center=y_feature_center, height=feature_height)
                if f_data["label_width"]:
                    canvas.setFillColorRGB(*f_data["stroke_colour"])
                    canvas.setFont(self.prms.args["feature_label_font_face"],
                                   self.track_data["f_label_font_size"])
                    fx_center = f_data["coordinates"]["center"]
                    canvas.drawString(f_data["label_position"][0], f_data["label_y_bottom"] + y_feature_upper,
                                      f_data["label"])
                    canvas.setLineWidth(self.prms.args["feature_stroke_width"])
                    underline_colour = f_data["stroke_colour"]
                    canvas.setStrokeColorRGB(*underline_colour)
                    canvas.setLineCap(1)
                    if f_data["label_row"] > 0:
                        p = canvas.beginPath()
                        for ls, le in f_data["label_line_coordinates"]:
                            p.moveTo(fx_center, ls + y_feature_upper)
                            p.lineTo(fx_center, le + y_feature_upper)
                        canvas.drawPath(p, stroke=1, fill=0)
                        l_start = f_data["coordinates"]["start"]
                        l_end = min(f_data["coordinates"]["end"], fx_center + self.track_data["feature_label_gap"])
                        l_end = f_data["coordinates"]["end"]
                        ly = y_feature_upper + f_data["label_y_bottom"] - self.track_data["feature_label_gap"] * 0.5
                        canvas.line(l_start, ly, l_end, ly)
                    else:
                        overlapping = min(f_data["coordinates"]["end"], f_data["label_position"][1]) - (
                            max(f_data["coordinates"]["start"], f_data["label_position"][0]))
                        if overlapping / (f_data["label_position"][1] - f_data["label_position"][0]) < 1:
                            l_start = max(f_data["coordinates"]["start"], fx_center -
                                          self.track_data["feature_label_gap"])
                            l_start = f_data["coordinates"]["start"]
                            l_end = min(f_data["coordinates"]["end"], fx_center + self.track_data["feature_label_gap"])
                            l_end = f_data["coordinates"]["end"]
                            ly = y_feature_upper + f_data["label_y_bottom"] - self.track_data["feature_label_gap"] * 0.5
                            canvas.line(l_start, ly, l_end, ly)
            # Axis ticks
            if self.prms.args["draw_individual_x_axis"]:
                canvas.setLineWidth(self.prms.args["x_axis_line_width"])
                canvas.setStrokeColorRGB(*lovis4u.Methods.get_colour_rgba("x_axis_line_colour", self.prms))
                canvas.setLineCap(1)
                canvas.setFillColorRGB(*lovis4u.Methods.get_colour_rgba("x_axis_line_colour", self.prms))
                canvas.setFont(self.prms.args["x_axis_ticks_labels_font_face"],
                               self.track_data["x_axis_annotation"]["label_size"])
                axis_line_y_coordinate = y_feature_bottom - self.prms.args["feature_bottom_gap"] * mm
                axis_tick_height = self.prms.args["x_axis_ticks_height"] * mm
                axis_tick_label_y_coordinate = axis_line_y_coordinate - self.prms.args["x_axis_ticks_height"] * \
                                               1.3 * mm - self.prms.args["x_axis_ticks_labels_height"] * mm
                for ati in range(len(self.track_data["x_axis_annotation"]["axis_tics_position"])):
                    tick_coordinate = self.track_data["x_axis_annotation"]["axis_tics_position"][ati]
                    tick_label_position = self.track_data["x_axis_annotation"]["tics_labels_coordinates"][ati]
                    tick_label = self.track_data["x_axis_annotation"]["axis_tics_labels"][ati]
                    canvas.drawCentredString(tick_label_position, axis_tick_label_y_coordinate, tick_label)
                    canvas.line(tick_coordinate, axis_line_y_coordinate, tick_coordinate,
                                axis_line_y_coordinate - axis_tick_height)
                for region in self.track_data["x_axis_annotation"]["axis_regions"]:
                    canvas.setLineCap(0)
                    canvas.line(region["start"], axis_line_y_coordinate, region["end"], axis_line_y_coordinate)
            return None
        except Exception as error:
            raise lovis4u.Manager.lovis4uError("Unable to draw a Locus track.") from error

    def __plot_feature(self, canvas: reportlab.pdfgen.canvas.Canvas, feature_data: dict, y_center: float,
                           height: float) -> None:
        """Helper method to plot feature polygone

        Returns:
            None

        """
        x_start = feature_data["coordinates"]["start"]
        x_end = feature_data["coordinates"]["end"]
        orientation = feature_data["coordinates"]["orient"]
        left_out = feature_data["coordinates"]["lout"]
        right_out = feature_data["coordinates"]["rout"]
        fill_colour = feature_data["fill_colour"]
        stroke_colour = feature_data["stroke_colour"]
        feature_type = feature_data["type"]
        y_center = y_center
        height = height


        canvas.setLineCap(0)
        canvas.setLineWidth(self.prms.args["feature_stroke_width"])
        arrow_length = min(height * self.prms.args["feature_arrow_length"], (x_end - x_start))
        p_noncoding = canvas.beginPath()
        if feature_type != "CDS":
            p_noncoding.moveTo(x_start, y_center + height / 2)
            p_noncoding.lineTo(x_end, y_center + height / 2)
            if not right_out:
                p_noncoding.lineTo(x_end, y_center - height / 2)
            else:
                p_noncoding.moveTo(x_end, y_center - height / 2)
            p_noncoding.lineTo(x_start, y_center - height / 2)
            if not left_out:
                p_noncoding.lineTo(x_start, y_center + height / 2)
            if not left_out and not right_out:
                p_noncoding.close()
            stroke, fill = 0, 0
            if stroke_colour:
                canvas.setStrokeColorRGB(*stroke_colour)
                stroke = 1
            if fill_colour:
                canvas.setFillColorRGB(*fill_colour)
                fill = 1
            canvas.drawPath(p_noncoding, stroke=stroke, fill=fill)
        p = canvas.beginPath()
        if orientation == 1:
            if right_out:
                p.moveTo(x_end, y_center - height / 2)
                p.lineTo(x_start, y_center - height / 2)
                p.lineTo(x_start, y_center + height / 2)
                p.lineTo(x_end, y_center + height / 2)
            else:
                p.moveTo(x_start, y_center + height / 2)
                p.lineTo(x_end - arrow_length, y_center + height / 2)
                p.lineTo(x_end, y_center)
                p.lineTo(x_end - arrow_length, y_center - height / 2)
                p.lineTo(x_start, y_center - height / 2)
                if not left_out:
                    p.lineTo(x_start, y_center + height / 2)
        elif orientation == -1:
            if left_out:
                p.moveTo(x_start, y_center - height / 2)
                p.lineTo(x_end, y_center - height / 2)
                p.lineTo(x_end, y_center + height / 2)
                p.lineTo(x_start, y_center + height / 2)
            else:
                p.moveTo(x_end, y_center + height / 2)
                p.lineTo(x_start + arrow_length, y_center + height / 2)
                p.lineTo(x_start, y_center)
                p.lineTo(x_start + arrow_length, y_center - height / 2)
                p.lineTo(x_end, y_center - height / 2)
                if not right_out:
                    p.lineTo(x_end, y_center + height / 2)
        if left_out and right_out:
            p.moveTo(x_start, y_center + height / 2)
            p.lineTo(x_end, y_center + height / 2)
            p.moveTo(x_start, y_center - height / 2)
            p.lineTo(x_end, y_center - height / 2)
        if not left_out and not right_out:
            p.close()
        stroke, fill = 0, 0
        if stroke_colour:
            canvas.setStrokeColorRGB(*stroke_colour)
            stroke = 1
        if fill_colour:
            canvas.setFillColorRGB(*fill_colour)
            fill = 1
        canvas.drawPath(p, stroke=stroke, fill=fill)
        return None


class ScaleVis(Track):
    """ScaleVis track object that handles visualisation of scale bottom line.

    Attributes:
        layout (dict): Layout built by CanvasManager's define_layout() method.
        track_data (dict): a dictionary with prepared track specific data.
        prms (lovis4u.Manager.Parameters): Parameters' class object.

    """

    def __init__(self, layout: dict, track_data: dict, parameters):
        """Create a LocusVis object.

        Arguments:
            layout (dict): Layout built by CanvasManager's define_layout() method.
            track_data (dict): a dictionary with prepared track specific data.
            parameters (lovis4u.Manager.Parameters): Parameters' class object.


        """
        super().__init__(layout, track_data, parameters)
        self.track_height = None

    def draw(self, canvas: reportlab.pdfgen.canvas.Canvas) -> None:
        """Draw a ScaleVis track.

        Arguments:
            canvas (reportlab.pdfgen.canvas.Canvas): a canvas object.

        Returns:
            None

        """
        try:
            y_upper = self.layout["current_y_coordinate"]
            y_bottom = self.layout["current_y_coordinate"] - self.track_data["track_height"]
            y_center = self.layout["current_y_coordinate"] - 0.5 * self.track_data["track_height"]

            middle_line_x_position = np.mean(self.track_data["coordinates"])

            canvas.setLineWidth(self.prms.args["scale_line_width"])
            canvas.setStrokeColorRGB(*lovis4u.Methods.get_colour_rgba("scale_line_colour", self.prms))
            canvas.setLineCap(0)
            canvas.setLineJoin(1)
            if self.track_data["style"] == "fancy":
                tick_height = self.track_data["scale_line_label_height"]

                p = canvas.beginPath()
                p.moveTo(self.track_data["coordinates"][0], y_center - 0.5 * tick_height)
                p.lineTo(self.track_data["coordinates"][0], y_center + 0.5 * tick_height)
                p.moveTo(self.track_data["coordinates"][0], y_center)
                p.lineTo(middle_line_x_position - 0.5 * self.track_data["scale_line_label_width"] -
                         self.track_data["space_width"], y_center)
                p.moveTo(self.track_data["coordinates"][1], y_center - 0.5 * tick_height)
                p.lineTo(self.track_data["coordinates"][1], y_center + 0.5 * tick_height)
                p.moveTo(self.track_data["coordinates"][1], y_center)
                p.lineTo(middle_line_x_position + 0.5 * self.track_data["scale_line_label_width"] +
                         self.track_data["space_width"], y_center)
            else:
                tick_height = self.prms.args["scale_line_tics_height"] * mm

                p = canvas.beginPath()
                p.moveTo(self.track_data["coordinates"][0], y_upper - tick_height)
                p.lineTo(self.track_data["coordinates"][0], y_upper)
                p.moveTo(self.track_data["coordinates"][0], y_upper - tick_height * 0.5)
                p.lineTo(self.track_data["coordinates"][1], y_upper - tick_height * 0.5)
                p.lineTo(self.track_data["coordinates"][1], y_upper)
                p.lineTo(self.track_data["coordinates"][1], y_upper - tick_height)
            canvas.drawPath(p, stroke=1, fill=0)
            canvas.setFillColorRGB(*lovis4u.Methods.get_colour_rgba("scale_line_colour", self.prms))
            canvas.setFont(self.prms.args["scale_line_label_font_face"],
                           self.track_data["scale_line_label_font_size"])
            canvas.drawCentredString(middle_line_x_position, y_bottom, self.track_data["scale_label"])
            return None
        except Exception as error:
            raise lovis4u.Manager.lovis4uError("Unable to draw a scale track.") from error


class ColorLegendVis(Track):
    """ColorLegend track object that handles visualisation of legend to feature's category colours.

    Attributes:
        layout (dict): Layout built by CanvasManager's define_layout() method.
        track_data (dict): a dictionary with prepared track specific data.
        prms (lovis4u.Manager.Parameters): Parameters' class object.

    """

    def __init__(self, layout: dict, track_data: dict, parameters):
        """Create a ColorLegend object.

        Arguments:
            layout (dict): Layout built by CanvasManager's define_layout() method.
            track_data (dict): a dictionary with prepared track specific data.
            parameters (lovis4u.Manager.Parameters): Parameters' class object.

        """
        super().__init__(layout, track_data, parameters)
        self.track_height = None

    def draw(self, canvas: reportlab.pdfgen.canvas.Canvas) -> None:
        """Draw a ColorLegend track.

        Arguments:
            canvas (reportlab.pdfgen.canvas.Canvas): a canvas object.

        Returns:
            None

        """
        try:
            y_upper = self.layout["current_y_coordinate"]
            canvas.setLineWidth(self.prms.args["scale_line_width"])
            canvas.setFont(self.prms.args["colour_legend_font_face"],
                           self.track_data["colour_legend_label_size"])
            for label_dict in self.track_data["labels"]:
                yl = y_upper + label_dict["relative_y"]
                yt = y_upper + label_dict["relative_y_text"]
                canvas.setFillColorRGB(*matplotlib.colors.hex2color(label_dict["colour"]),
                                       self.prms.args["category_annotation_alpha"])
                canvas.rect(label_dict["label_x"], yl, label_dict["label_width"],
                            self.track_data["line_width"], fill=1, stroke=0)
                canvas.setFillColorRGB(*lovis4u.Methods.get_colour_rgba("colour_legend_label_colour", self.prms))
                canvas.drawString(label_dict["label_x"], yt, label_dict["label"])
            return None
        except Exception as error:
            raise lovis4u.Manager.lovis4uError("Unable to draw a colour legend track.") from error
