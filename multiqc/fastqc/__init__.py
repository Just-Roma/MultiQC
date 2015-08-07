#!/usr/bin/env python

""" MultiQC module to parse output from FastQC
"""

import collections
import json
import logging
import os
import re
import shutil
import zipfile

import multiqc

class MultiqcModule(multiqc.BaseMultiqcModule):

    def __init__(self, analysis_dir, output_dir):

        # Initialise the parent object
        super(MultiqcModule, self).__init__()

        # Static variables
        self.name = "FastQC"
        self.analysis_dir = analysis_dir
        self.output_dir = output_dir
        self.data_dir = os.path.join(output_dir, 'report_data')

        # Find and load any FastQC reports
        fastqc_raw_data = {}
        for root, dirnames, filenames in os.walk(analysis_dir):
            # Extracted FastQC directory
            if root[-7:] == '_fastqc' and 'fastqc_data.txt' in filenames:
                s_name = os.path.basename(root)
                s_name = s_name[:-7]
                d_path = os.path.join(root, 'fastqc_data.txt')
                with open (d_path, "r") as f:
                    fastqc_raw_data[s_name] = f.read()

            # Zipped FastQC report
            for f in filenames:
                if f[-11:] == '_fastqc.zip':
                    s_name = f[:-11]
                    d_name = f[:-4]
                    fqc_zip = zipfile.ZipFile(os.path.join(root, f))
                    try:
                        with fqc_zip.open(os.path.join(d_name, 'fastqc_data.txt')) as f:
                            fastqc_raw_data[s_name] = f.read()
                    except KeyError:
                        pass # Can't find fastqc_raw_data.txt in the zip file

        if len(fastqc_raw_data) == 0:
            logging.debug("Could not find any FastQC reports in {}".format(analysis_dir))
            return

        logging.info("Found {} FastQC reports".format(len(fastqc_raw_data)))

        self.sections = list()

        # Section 1 - Basic Stats
        parsed_stats = self.fastqc_basic_stats(fastqc_raw_data)
        stats_table = self.fastqc_stats_table(parsed_stats)
        self.sections.append({
            'name': 'Basic Stats',
            'content': stats_table
        })

        # Section 2 - Quality Histograms
        histogram_data = self.fastqc_seq_quality(fastqc_raw_data)
        self.sections.append({
            'name': 'Sequence Quality Histograms',
            'content': self.fastqc_quality_overlay_plot(histogram_data)
        })

        # Section 3 - Per-base sequence content
        seq_content = self.fastqc_seq_content(fastqc_raw_data)
        self.sections.append({
            'name': 'Per Base Sequence Content',
            'content': self.fastqc_seq_heatmap(seq_content)
        })

        # Section 4 - GC Content
        gc_data = self.fastqc_gc_content(fastqc_raw_data)
        self.sections.append({
            'name': 'Per Sequence GC Content',
            'content': self.fastqc_gc_overlay_plot(gc_data)
        })

        # Section 5 - Adapter Content
        adapter_data = self.fastqc_adapter_content(fastqc_raw_data)
        self.sections.append({
            'name': 'Adapter Content',
            'content': self.fastqc_adapter_overlay_plot(adapter_data)
        })


        # Copy across the required module files (CSS / Javascript etc)
        self.init_modfiles()

    def fastqc_basic_stats(self, fastqc_raw_data):
        """ Parse fastqc_data.txt for basic stats.
        Returns a 2D dict with basic stats, sample names as first keys,
        then statistic type as second key. """
        parsed_data = {}
        for s, data in fastqc_raw_data.iteritems():
            parsed_data[s] = {}

            dups = re.search("#Total Deduplicated Percentage\s+([\d\.]+)", data)
            if dups:
                parsed_data[s]['percent_duplicates'] = "{0:.2f}%".format(100 - float(dups.group(1)))

            seqlen = re.search("Sequence length\s+([\d-]+)", data)
            if seqlen:
                parsed_data[s]['sequence_length'] = seqlen.group(1)

            gc = re.search("%GC\s+(\d+)", data)
            if gc:
                parsed_data[s]['percent_gc'] = "{}%".format(gc.group(1))

            numseq = re.search("Total Sequences\s+(\d+)", data)
            if numseq:
                parsed_data[s]['total_sequences_m'] = "{0:.1f}".format(float(numseq.group(1))/1000000)

            # Work out the average sequence length as the range is a bit useless
            if '-' in parsed_data[s]['sequence_length']:
                in_module = False
                total_count = 0
                total_bp = 0
                for l in data.splitlines():
                    if l[:30] == ">>Sequence Length Distribution":
                        in_module = True
                    elif l == ">>END_MODULE":
                        in_module = False
                    elif in_module is True:
                        len_matches = re.search("([\d-]+)\s+([\d\.E]+)", l)
                        try:
                            avg_len = len_matches.group(1)
                            if '-' in avg_len:
                                maxlen = int(avg_len.split("-",1)[1])
                                minlen = int(avg_len.split("-",1)[0])
                                avg_len = ((maxlen - minlen)/2) + minlen
                            total_bp += float(len_matches.group(2)) * avg_len
                            total_count += float(len_matches.group(2))
                        except AttributeError:
                            pass
                if total_count > 0:
                    parsed_data[s]['sequence_length'] = '{:.0f}'.format(total_bp / total_count)

        return parsed_data

    def fastqc_stats_table(self, parsed_stats):
        """ Take the parsed stats from the FastQC report and turn it into a
        nice attractive HTML table. """

        # Order the table by the sample names
        parsed_stats = collections.OrderedDict(sorted(parsed_stats.items()))

        stats_table_headers = {
            'percent_duplicates': '<span data-toggle="tooltip" title="% Duplicate Reads">%&nbsp; Dups</span>',
            'sequence_length': '<span data-toggle="tooltip" title="Average Sequence Length (bp)">Length</span>',
            'percent_gc': '<span data-toggle="tooltip" title="Average % GC Content">%&nbsp;GC</span>',
            'total_sequences_m': '<span data-toggle="tooltip" title="Total Sequences (millions)">M Seqs</span>'
        }
        header_attrs = {
            'percent_duplicates': 'class="chroma-col" data-chroma-scale="RdYlGn-rev" data-chroma-max="100" data-chroma-min="0"',
            'sequence_length': 'class="chroma-col" data-chroma-scale="RdYlGn" data-chroma-min="0"',
            'percent_gc': 'class="chroma-col"  data-chroma-scale="PRGn" data-chroma-max="100" data-chroma-min="0"',
            'total_sequences_m': 'class="chroma-col" data-chroma-scale="Blues" data-chroma-min="0"'
        }
        cell_attrs = {
            'percent_duplicates': 'class="text-right"',
            'sequence_length': 'class="text-right"',
            'percent_gc': 'class="text-right"',
            'total_sequences_m': 'class="text-right"'
        }

        # Use the base class dict_to_table to make a HTML table
        parsed_stats = self.dict_to_table(parsed_stats,
            colheaders=stats_table_headers,
            table_attrs='class="table table-hover table-bordered table-responsive table-condensed table-nostretch"',
            header_attrs=header_attrs,
            cell_attrs=cell_attrs)

        return parsed_stats

    def fastqc_seq_quality(self, fastqc_raw_data):
        """ Parse the 'Per base sequence quality' data from fastqc_data.txt
        Returns a 2D dict, sample names as first keys, then a dict of lists
        containing base, mean, median, lower_quart, upper_quart, 10_percentile
        and 90_percentile. """

        parsed_data = {}
        for s, data in fastqc_raw_data.iteritems():
            parsed_data[s] = {}
            parsed_data[s]['base'] = list()
            parsed_data[s]['mean'] = list()
            parsed_data[s]['median'] = list()
            parsed_data[s]['lower_quart'] = list()
            parsed_data[s]['upper_quart'] = list()
            parsed_data[s]['10_percentile'] = list()
            parsed_data[s]['90_percentile'] = list()
            in_module = False
            for l in data.splitlines():
                if l[:27] == ">>Per base sequence quality":
                    in_module = True
                elif l == ">>END_MODULE":
                    in_module = False
                elif in_module is True:
                    quals = re.search("([\d-]+)\s+([\d\.]+)\s+([\d\.]+)\s+([\d\.]+)\s+([\d\.]+)\s+([\d\.]+)\s+([\d\.]+)", l)
                    try:
                        parsed_data[s]['base'].append(quals.group(1))
                        parsed_data[s]['mean'].append(float(quals.group(2)))
                        parsed_data[s]['median'].append(float(quals.group(3)))
                        parsed_data[s]['lower_quart'].append(float(quals.group(4)))
                        parsed_data[s]['upper_quart'].append(float(quals.group(5)))
                        parsed_data[s]['10_percentile'].append(float(quals.group(6)))
                        parsed_data[s]['90_percentile'].append(float(quals.group(7)))
                    except AttributeError:
                        pass

        return parsed_data

    def fastqc_quality_overlay_plot (self, parsed_data):

        data = list()
        for s in sorted(parsed_data):
            pairs = list()
            for k, p in enumerate(parsed_data[s]['base']):
                pairs.append([int(p.split('-', 1)[0]), parsed_data[s]['mean'][k]])
            data.append({
                'name': s,
                'data': pairs
            })

        html = '<div id="fastqc_quality_overlay" style="height:500px;"></div> \
        <script type="text/javascript"> \
            fastqc_overlay_hist_data = {};\
            $(function () {{ \
                plot_xy_line_graph("#fastqc_quality_overlay", \
                    fastqc_overlay_hist_data, "Mean Quality Scores", "Phred Score", "Position (bp)", \
                    undefined, 0, undefined, undefined, "Base {{point.x}}"); \
            }}); \
        </script>'.format(json.dumps(data));

        return html


    def fastqc_seq_content(self, fastqc_raw_data):
        """ Parse the 'Per base sequence content' data from fastqc_data.txt
        Returns a 3D dict, sample names as first keys, second key the base
        position and third key with the base ([ACTG]). Values contain percentages """

        parsed_data = {}
        for s, data in fastqc_raw_data.iteritems():
            parsed_data[s] = {}
            in_module = False
            for l in data.splitlines():
                if l[:27] == ">>Per base sequence content":
                    in_module = True
                elif l == ">>END_MODULE":
                    in_module = False
                elif in_module is True:
                    seq_matches = re.search("([\d-]+)\s+([\d\.]+)\s+([\d\.]+)\s+([\d\.]+)\s+([\d\.]+)", l)
                    try:
                        bp = int(seq_matches.group(1).split('-', 1)[0])
                        parsed_data[s][bp] = {}
                        parsed_data[s][bp]['base'] = seq_matches.group(1)
                        parsed_data[s][bp]['G'] = float(seq_matches.group(2))
                        parsed_data[s][bp]['A'] = float(seq_matches.group(3))
                        parsed_data[s][bp]['T'] = float(seq_matches.group(4))
                        parsed_data[s][bp]['C'] = float(seq_matches.group(5))
                    except AttributeError:
                        if l[:1] != '#':
                            raise

        return parsed_data

    def fastqc_seq_heatmap (self, parsed_data):

        # Order the table by the sample names
        parsed_data = collections.OrderedDict(sorted(parsed_data.items()))

        html = '<canvas id="fastqc_seq_heatmap" height="300px" width="800px" style="width:100%;"></canvas> \n\
        <ul id="fastqc_seq_heatmap_key">\n\
            <li>%G: <span id="fastqc_seq_heatmap_key_g">?</span> <span id="fastqc_seq_heatmap_key_colourbar_g" class="heatmap_colourbar"><span></span></span></li>\n\
            <li>%A: <span id="fastqc_seq_heatmap_key_a">?</span> <span id="fastqc_seq_heatmap_key_colourbar_a" class="heatmap_colourbar"><span></span></span></li>\n\
            <li>%T: <span id="fastqc_seq_heatmap_key_t">?</span> <span id="fastqc_seq_heatmap_key_colourbar_t" class="heatmap_colourbar"><span></span></span></li>\n\
            <li>%C: <span id="fastqc_seq_heatmap_key_c">?</span> <span id="fastqc_seq_heatmap_key_colourbar_c" class="heatmap_colourbar"><span></span></span></li>\n\
        </ul>\n\
        <p id="fastqc_seq_heatmap_footer"><small class="text-muted">Fuzzy after resizing the window? Hit your browser\'s refresh button..</small></p>\n\
        <script type="text/javascript"> \n\
            fastqc_seq_content_data = {};\n\
            $(function () {{ \n\
                var num_samples = Object.keys(fastqc_seq_content_data).length; \n\
                var num_data = undefined; \n\
                var ypos = 0; \n\
                var xpos_start = 0; \n\
                var labels_width = 0; \n\
                var max_bp = 0; \n\
                // Convert the CSS percentage size into pixels \n\
                var c_width = $("#fastqc_seq_heatmap").width(); \n\
                var c_height = $("#fastqc_seq_heatmap").height(); \n\
                if( c_height / num_samples > 14){{ \n\
                    c_height = num_samples * 14; \n\
                }} \n\
                $("#fastqc_seq_heatmap").prop("width", c_width); \n\
                $("#fastqc_seq_heatmap").prop("height", c_height); \n\
                var s_height = c_height / num_samples; \n\
                var canvas = document.getElementById("fastqc_seq_heatmap"); \n\
                if (canvas.getContext) {{ \n\
                    var ctx = canvas.getContext("2d"); \n\
                    ctx.font = (parseInt(s_height)-2)+"px Arial"; \n\
                    // First, do labels and get max base pairs \n\
                    $.each(fastqc_seq_content_data, function(name, s){{ \n\
                        var txt_w = ctx.measureText(name).width; \n\
                        if(txt_w > labels_width) {{ \n\
                            labels_width = txt_w; \n\
                        }} \n\
                        ctx.fillText(name, 0, ypos+parseInt(s_height)-2); \n\
                        ypos += s_height; \n\
                        $.each(s, function(bp, v){{ \n\
                            bp = parseInt(bp); \n\
                            if(bp > max_bp){{ \n\
                                max_bp = bp; \n\
                            }} \n\
                        }}); \n\
                    }}); \n\
                    labels_width += 5; \n\
                    ypos = 0; \n\
                    xpos_start = labels_width; \n\
                    c_width -= labels_width; \n\
                    $.each(fastqc_seq_content_data, function(name, s){{ \n\
                        if (num_data == undefined){{ \n\
                            var s_width = c_width / Object.keys(s).length; \n\
                        }} \n\
                        var xpos = xpos_start; \n\
                        var last_bp = 0; \n\
                        $.each(s, function(bp, v){{ \n\
                            bp = parseInt(bp); \n\
                            var this_width = (bp - last_bp) * (c_width / max_bp); \n\
                            last_bp = bp; \n\
                            var r = (v["G"] / 100)*255; \n\
                            var g = (v["A"] / 100)*255; \n\
                            var b = (v["T"] / 100)*255; \n\
                            ctx.fillStyle = chroma(r,g,b).css(); \n\
                            ctx.fillRect (xpos, ypos, this_width, s_height); \n\
                            xpos += this_width; \n\
                        }}); \n\
                        ypos += s_height; \n\
                    }}); \n\
                }} \n\
                \n\
                // http://stackoverflow.com/questions/6735470/get-pixel-color-from-canvas-on-mouseover \n\
                $("#fastqc_seq_heatmap").mousemove(function(e) {{ \n\
                    var pos = findPos(this); \n\
                    var x = e.pageX - pos.x; \n\
                    var y = e.pageY - pos.y; \n\
                    var coord = "x=" + x + ", y=" + y; \n\
                    var c = this.getContext("2d"); \n\
                    var p = c.getImageData(x, y, 1, 1).data;  \n\
                    var seq_g = Math.round((p[0]/255)*100); \n\
                    var seq_a = Math.round((p[1]/255)*100); \n\
                    var seq_t = Math.round((p[2]/255)*100); \n\
                    var seq_c = 100 - (seq_g + seq_a + seq_t); \n\
                    $("#fastqc_seq_heatmap_key_g").text(seq_g+"%"); \n\
                    $("#fastqc_seq_heatmap_key_a").text(seq_a+"%"); \n\
                    $("#fastqc_seq_heatmap_key_t").text(seq_t+"%"); \n\
                    $("#fastqc_seq_heatmap_key_c").text(seq_c+"%"); \n\
                    $("#fastqc_seq_heatmap_key_colourbar_g span").css("margin-left", seq_g); \n\
                    $("#fastqc_seq_heatmap_key_colourbar_a span").css("margin-left", seq_a); \n\
                    $("#fastqc_seq_heatmap_key_colourbar_t span").css("margin-left", seq_t); \n\
                    $("#fastqc_seq_heatmap_key_colourbar_c span").css("margin-left", seq_c); \n\
                }}); \n\
            }}); \n\
        </script>'.format(json.dumps(parsed_data))

        return html


    def fastqc_gc_content(self, fastqc_raw_data):
        """ Parse the 'Per sequence GC content' data from fastqc_data.txt
        Returns a 2D dict, sample names as first keys, then a dict of with keys
        containing percentage and values containing counts """

        parsed_data = {}
        for s, data in fastqc_raw_data.iteritems():
            parsed_data[s] = {}
            in_module = False
            for l in data.splitlines():
                if l[:25] == ">>Per sequence GC content":
                    in_module = True
                elif l == ">>END_MODULE":
                    in_module = False
                elif in_module is True:
                    gc_matches = re.search("([\d]+)\s+([\d\.E]+)", l)
                    try:
                        parsed_data[s][int(gc_matches.group(1))] = float(gc_matches.group(2))
                    except AttributeError:
                        pass

        return parsed_data

    def fastqc_gc_overlay_plot (self, parsed_data):
        data = list()
        for s in sorted(parsed_data):
            pairs = list()
            for k, p in iter(sorted(parsed_data[s].iteritems())):
                pairs.append([int(k), p])
            data.append({
                'name': s,
                'data': pairs
            })

        html = '<div id="fastqc_gc_overlay" style="height:500px;"></div> \
        <script type="text/javascript"> \
            fastqc_overlay_gc_data = {};\
            $(function () {{ \
                plot_xy_line_graph("#fastqc_gc_overlay", \
                    fastqc_overlay_gc_data, "Per Sequence GC Content", "Count", "%GC", \
                    undefined, 0, 100, 0, "{{point.x}}% GC"); \
            }}); \
        </script>'.format(json.dumps(data));

        return html


    def fastqc_adapter_content(self, fastqc_raw_data):
        """ Parse the 'Adapter Content' data from fastqc_data.txt
        Returns a 2D dict, sample names as first keys, then a dict of with keys
        containing adapter type and values containing counts """

        parsed_data = {}
        for s, data in fastqc_raw_data.iteritems():
            parsed_data[s] = dict()
            in_module = False
            adapter_types = []
            for l in data.splitlines():
                if l[:17] == ">>Adapter Content":
                    in_module = True
                elif l == ">>END_MODULE":
                    in_module = False
                elif in_module is True:
                    if l[:1] == '#':
                        adapter_types = l[1:].split("\t")[1:]
                        for a in adapter_types:
                            parsed_data[s][a] = dict()
                    else:
                        cols = l.split("\t")
                        pos = int(cols[0].split('-', 1)[0])
                        for idx, val in enumerate(cols[1:]):
                            parsed_data[s][adapter_types[idx]][pos] = float(val)
        return parsed_data

    def fastqc_adapter_overlay_plot (self, parsed_data):

        data = list()
        for s in sorted(parsed_data):
            for a, d in parsed_data[s].iteritems():
                pairs = list()
                for base, count in iter(sorted(d.iteritems())):
                    pairs.append([base, count])
                data.append({
                    'name': '{} - {}'.format(s, a),
                    'data': pairs
                })

        html = '<div id="fastqc_adapter_overlay" style="height:500px;"></div> \
        <script type="text/javascript"> \
            fastqc_adapter_data = {};\
            $(function () {{ \
                plot_xy_line_graph("#fastqc_adapter_overlay", \
                    fastqc_adapter_data, "Adapter Content", "% of Sequences", "Position", \
                    100, 0, undefined, undefined, "Base {{point.x}}"); \
            }}); \
        </script>'.format(json.dumps(data));

        return html


    def init_modfiles(self):
        """ Copy the required assets into the output directory.
        Need to do this manually as we don't want to over-write
        the existing directory tree."""

        self.css = [ os.path.join('assets', 'css', 'multiqc_fastqc.css') ]
        self.js = [ os.path.join('assets', 'js', 'multiqc_fastqc.js') ]

        for f in self.css + self.js:
            d = os.path.join(self.output_dir, os.path.dirname(f))
            if not os.path.exists(d):
                os.makedirs(d)
            if not os.path.exists(os.path.join(self.output_dir, f)):
                shutil.copy(os.path.join(os.path.dirname(__file__), f), os.path.join(self.output_dir, f))
