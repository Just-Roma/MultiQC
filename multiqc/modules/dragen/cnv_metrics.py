'''""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
This module gathers CNV metrics data and prepares it for the output report.
It relies on the following official sources:
https://support-docs.illumina.com/SW/DRAGEN_v41/Content/SW/DRAGEN/CNVMetrics.htm
https://support-docs.illumina.com/SW/DRAGEN_v41/Content/SW/DRAGEN/CNVSampleCorrelationSexGenotyper.htm
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""'''

import logging
import re
from collections import defaultdict

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.utils.util_functions import write_data_file

from .utils import (check_duplicate_samples, clean_headers, make_general_stats,
                    make_own_table_plot, make_parsing_log_report,
                    order_headers)

# Initialise the logger.
log = logging.getLogger(__name__)


NAMESPACE = "DRAGEN CNV"


'''"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The SINGLE_HEADER is used to define the common settings for all headers:
https://github.com/ewels/MultiQC/blob/master/docs/plots.md#creating-a-table
'''
SINGLE_HEADER = {
    "namespace": NAMESPACE,  # Name for grouping. Prepends desc and is in Config Columns modal
    "title": None,  # Short title, table column title
    "description": None,  # Longer description, goes in mouse hover text
    "max": None,  # Minimum value in range, for bar / colour coding
    "min": 0,  # Maximum value in range, for bar / colour coding
    "ceiling": None,  # Maximum value for automatic bar limit
    "floor": None,  # Minimum value for automatic bar limit
    "minRange": None,  # Minimum range for automatic bar
    "scale": "GnBu",  # Colour scale for colour coding. False to disable.
    "bgcols": None,  # Dict with values: background colours for categorical data.
    "colour": "255, 150, 30",  # Colour for column grouping
    "suffix": None,  # Suffix for value (eg. "%")
    "format": "{:,.1f}",  # Value format string - default 1 decimal place
    "cond_formatting_rules": None,  # Rules for conditional formatting table cell values.
    "cond_formatting_colours": None,  # Styles for conditional formatting of table cell values
    "shared_key": None,  # See the link for description
    "modify": None,  # Lambda function to modify values, special case, see below
    "hidden": True,  # Set to True to hide the column on page load
    "hidden_own": False,  # Set to True to hide all columns in own plot.
    "exclude": True,  # True to exclude all metrics from the general html table.
}
'''""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The EXTRA_HEADER contains the logical keys of the SINGLE_HEADER. These are just extensions,
but can be used as if they were a part of the SINGLE_HEADER. Extra configurations are not
set by default. Only the SINGLE_HEADER is used as a basis for a header. But you can always
append an extra config(eg hidden_own, exclude) to the SINGLE_HEADER to overcome this issue.

Hints and rules:
- "hidden" is used for own CNV section if "hidden_own" is not provided.

- "exclude" and "exclude_own" may be useful in cases where html size matters.
  Eg. many samples but only few metrics of interest.

- "order_priority" can be used to arrange table's columns. Only int and float are valid.
'''
EXTRA_HEADER = {
    "hidden_own": False,  # For non-general plot in the CNV section.
    "exclude": False,  # True to exclude metric from the general html table.
    "exclude_own": False,  # True to exclude from own CNV section html table.
    "order_priority": None,  # Used to specify columns' order in all tables.
}
'''"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The METRICS is the container for "high level" representations of the standard metrics.
They are designed to be seen as simple identifiers of the well known metrics.
You can set any possible real/virtual configuration defined in the SINGLE_HEADER.

Some rules:

- "modify" must conform to the following general structure:
  lambda x: (arbitrary expression) if isinstance(x, str) else (math expression)
  The x is either number or string. The latter is some value, which could not be
  converted to int or float. NA for instance. Such values can be modified in the
  (arbitrary expression) if needed. If not then just return the input.

- The metrics are sensitive and must conform to the output format in csv.

- If you wish to add new metrics for the 'CNV Summary' section, please take a look
  at the "make_metric_id".
'''

V2 = "(second value)"  # Used to create a unique ID if a metric has two values.

SEX_GENOTYPER_SECTION = "SEX GENOTYPER"
CNV_SUMMARY_SECTION = "CNV SUMMARY"

# Used to order panel of samples. Please do not exceed it.
RESERVED_START = 100

# Used to make a number "narrower" to not overlap with adjacent right columns.
MODIFY = lambda x: x if isinstance(x, str) else x * 0.000001
PREFIX = "M"
DESCR_UNIT = "millions"

METRICS = {
    # The SEX_GENOTYPER_SECTION defines only the headers for the case sample.
    # The creation of headers for panel of samples is done by the make_panel_headers.
    SEX_GENOTYPER_SECTION: {
        "karyotype": {
            "order_priority": 0,
            "exclude": False,
            "title": "Kar",
            "colour": "0, 255, 0",
            "hidden": False,
            "description": "Estimated sex karyotype for the case sample.",
            "cond_formatting_rules": {
                "red": [{"s_contains": ""}],  # Each value is red by default.
                "green": [{"s_eq": "XX"}, {"s_eq": "YX"}, {"s_eq": "XY"}],
            },
            "cond_formatting_colours": [
                {"red": "#FF0000"},
                {"green": "#00FF00"},
            ],
        },
        "confidence": {
            "order_priority": 1,
            "exclude": False,
            "title": "KarConf",
            "colour": "0, 255, 0",
            "hidden": False,
            "suffix": " %",
            "format": "{:,.4f}",
            "max": 1.0,
            "cond_formatting_rules": {
                "red": [{"lt": 0}, {"gt": 1}],
            },
            "cond_formatting_colours": [
                {"red": "#FF0000"},
            ],
            "description": "Confidence metric ranging from 0.0 to 1.0. If the case sample sex is specified, this metric is 0.0.",
        },
    },
    CNV_SUMMARY_SECTION: {
        "bases in reference genome": {
            "title": config.base_count_prefix + " bases",
            "hidden_own": True,
            "scale": "RdYlGn",
            "colour": "0, 0, 255",
            "format": "{:,.0f}",
            "description": "Bases ({}) in reference genome in use.".format(config.base_count_desc),
            "modify": lambda x: x if isinstance(x, str) else x * config.base_count_multiplier,
        },
        "average alignment coverage over genome": {
            "order_priority": 2,
            "title": "AvgCov",
            "scale": "YlOrRd",
            "colour": "150, 50, 0",
            "suffix": " %",
            "description": "Average alignment coverage over genome.",
        },
        "number of alignment records": {
            "order_priority": 3,
            "title": PREFIX + " Aln records",
            "scale": "BrBG",
            "colour": "0, 55, 100",
            "format": "{:,.0f}",
            "modify": MODIFY,
            "description": "Number ({}) of alignment records processed.".format(DESCR_UNIT),
        },
        "coverage uniformity": {
            "order_priority": 4,
            "title": "CovUnif",
            "scale": "Reds",
            "colour": "0, 0, 255",
            "suffix": " %",
            "format": "{:,.2f}",
            "description": "Coverage uniformity.",
        },
        "number of target intervals": {
            "order_priority": 5,
            "title": "Trg intervals",
            "scale": "Greens",
            "colour": "0, 255, 200",
            "format": "{:,.0f}",
            "description": "Number of target intervals.",
        },
        "number of segments": {
            "order_priority": 6,
            "title": "Segments",
            "scale": "Blues",
            "colour": "250, 0, 0",
            "format": "{:,.0f}",
            "description": "Number of segments.",
        },
        "number of filtered records (total)": {
            "order_priority": 7.0,
            "title": PREFIX + " Total",
            "scale": "Purples",
            "colour": "255, 0, 255",
            "format": "{:,.2f}",
            "modify": MODIFY,
            "description": "Number ({}) of filtered records (total).".format(DESCR_UNIT),
        },
        "number of filtered records (total)"
        + V2: {
            "order_priority": 7.1,
            "title": "Total%",
            "scale": "Oranges",
            "colour": "255, 0, 255",
            "suffix": " %",
            "description": "Percentage of filtered records (total).",
        },
        "number of filtered records (duplicates)": {
            "order_priority": 7.2,
            "title": "Duplicates",
            "colour": "255, 0, 255",
            "scale": "Greens",
            "format": "{:,.2f}",
            "description": "Number of filtered records (due to duplicates).",
        },
        "number of filtered records (duplicates)"
        + V2: {
            "order_priority": 7.3,
            "title": "Duplicates%",
            "colour": "255, 0, 255",
            "suffix": " %",
            "description": "Percentage of filtered records (due to duplicates).",
        },
        "number of filtered records (mapq)": {
            "order_priority": 7.4,
            "title": PREFIX + " MAPQ",
            "colour": "255, 0, 255",
            "scale": "Reds",
            "format": "{:,.2f}",
            "modify": MODIFY,
            "description": "Number ({}) of filtered records (due to MAPQ).".format(DESCR_UNIT),
        },
        "number of filtered records (mapq)"
        + V2: {
            "order_priority": 7.5,
            "title": "MAPQ%",
            "colour": "255, 0, 255",
            "suffix": " %",
            "description": "Percentage of filtered records (due to MAPQ).",
        },
        "number of filtered records (unmapped)": {
            "order_priority": 7.6,
            "title": PREFIX + " Unmapped",
            "scale": "PiYG",
            "colour": "255, 0, 255",
            "format": "{:,.2f}",
            "modify": MODIFY,
            "description": "Number ({}) of filtered records (due to being unmapped).".format(DESCR_UNIT),
        },
        "number of filtered records (unmapped)"
        + V2: {
            "order_priority": 7.7,
            "title": "Unmapped%",
            "colour": "255, 0, 255",
            "suffix": " %",
            "description": "Percentage of filtered records (due to being unmapped).",
        },
        "number of amplifications": {
            "order_priority": 8.0,
            "title": "Amplifs",
            "scale": "RdYlGn",
            "colour": "255, 255, 0",
            "format": "{:,.0f}",
            "description": "Number of amplifications.",
        },
        "number of passing amplifications": {
            "order_priority": 8.1,
            "title": "Amplifs PASS",
            "scale": "PiYG",
            "colour": "255, 255, 0",
            "format": "{:,.0f}",
            "description": "Number of PASS amplifications.",
        },
        "number of passing amplifications"
        + V2: {
            "order_priority": 8.2,
            "title": "Amplifs PASS%",
            "colour": "255, 255, 0",
            "suffix": " %",
            "description": "Percentage of PASS amplifications.",
        },
        "number of deletions": {
            "order_priority": 9.0,
            "title": "Deletions",
            "scale": "Reds",
            "colour": "0, 0, 255",
            "format": "{:,.0f}",
            "description": "Number of deletions.",
        },
        "number of passing deletions": {
            "order_priority": 9.1,
            "title": "Deletions PASS",
            "scale": "BrBG",
            "colour": "0, 0, 255",
            "format": "{:,.0f}",
            "description": "Number of PASS deletions.",
        },
        "number of passing deletions"
        + V2: {
            "order_priority": 9.2,
            "title": "Deletions PASS%",
            "scale": "Oranges",
            "colour": "0, 0, 255",
            "suffix": " %",
            "description": "Percentage of PASS deletions.",
        },
        "number of de novo calls": {
            "order_priority": 10.0,
            "title": "de Novo",
            "scale": "Greens",
            "colour": "255, 0, 0",
            "format": "{:,.0f}",
        },
        "number of passing de novo calls": {
            "order_priority": 10.1,
            "title": "de Novo PASS",
            "colour": "255, 0, 0",
            "format": "{:,.0f}",
            "description": "Number of passing de novo calls.",
        },
        "number of passing de novo calls"
        + V2: {
            "order_priority": 10.2,
            "title": "de Novo PASS%",
            "scale": "Purples",
            "colour": "255, 0, 0",
            "suffix": " %",
            "description": "Percentage of passing de novo calls.",
        },
    },
}


# The TABLE_CONFIG defines configs for the whole own table.
TABLE_CONFIG = {
    "namespace": NAMESPACE,  # Name for grouping. Prepends desc and is in Config Columns modal
    "id": None,  # ID used for the table
    "table_title": "CNV Metrics",  # Title of the table. Used in the column config modal
    "save_file": False,  # Whether to save the table data to a file
    "raw_data_fn": None,  # File basename to use for raw data file
    "sortRows": True,  # Whether to sort rows alphabetically
    "only_defined_headers": True,  # Only show columns that are defined in the headers config
    "col1_header": None,  # The header used for the first column with sample names.
    "no_beeswarm": False,  # Force a table to always be plotted (beeswarm by default if many rows)
}


class DragenCNVMetrics(BaseMultiqcModule):
    """Public members of the DragenCNVMetrics module are available to other dragen modules.
    To avoid cluttering up the global dragen space, only one method is defined."""

    def add_cnv_metrics(self):
        """The main function of the dragen CNV metrics module.
        The public members of the BaseMultiqcModule and dragen modules defined in
        MultiqcModule are available within it. Returns either a view object containing
        a list of sample names, or an empty set if data could not be collected or
        became empty after ignoring samples."""

        cnv_data = {}
        cnv_headers = {}

        # Stores cleaned sample names with references to file handlers.
        # Used later to check for duplicates and inform about found ones.
        all_samples = defaultdict(list)

        for file in self.find_log_files("dragen/cnv_metrics"):
            out = cnv_parser(file, cnv_headers)
            if out["success"]:
                self.add_data_source(file, section="stats")
                cleaned_sample = self.clean_s_name(out["sample_name"], file)
                all_samples[cleaned_sample].append(file)
                # Add/overwrite the sample.
                cnv_data[cleaned_sample] = out["data"]

        cnv_data = self.ignore_samples(cnv_data)

        if not cnv_data:
            return set()

        check_duplicate_samples(all_samples, log)

        # Write data into the general table.
        gen_data, gen_headers = make_general_stats(cnv_data, cnv_headers)
        self.general_stats_addcols(gen_data, gen_headers, namespace=NAMESPACE)

        # Write data to file.
        self.write_data_file(cnv_data, "dragen_cnv_metrics")

        # Reporting found info/warnings/errors, which were collected while
        # calling the cnv_parser. You can disable it anytime, if it is not wanted.
        make_parsing_log_report("cnv_metrics", log_data, log)

        self.add_section(
            name="CNV Metrics",
            anchor="dragen-cnv-metrics-own-section",
            description="""
            DRAGEN CNV outputs metrics in CSV format.
            The output follows the general convention for QC metrics reporting in DRAGEN.
            The CNV metrics are output to a file with a *.cnv_metrics.csv file extension.
            The following list summarizes the metrics that are output from a CNV run.
            Press the `Help` button for details.
            """,
            helptext="""
            Sex Genotyper Metrics:

            * Estimated sex karyotype for the sample, along with a confidence metric ranging from 0.0 to 1.0.
            If the sample sex is specified, this metric is 0.0.
            * If using a panel of normals, all panel samples are also reported.

            CNV Summary Metrics:

            * Bases in reference genome in use.
            * Average alignment coverage over genome.
            * Number of alignment records processed.
            * Number of filtered records (total).
            * Number of filtered records (due to duplicates).
            * Number of filtered records (due to MAPQ).
            * Number of filtered records (due to being unmapped).
            * Number of target intervals.
            * Number of normal samples.
            * Number of segments.
            * Number of amplifications.
            * Number of deletions.
            * Number of PASS amplifications.
            * Number of PASS deletions.
            """,
            plot=make_own_table_plot(
                cnv_data,
                cnv_headers,
                {config: val for config, val in TABLE_CONFIG.items() if val is not None},
            ),
        )
        return cnv_data.keys()


def construct_coverage_parser():
    """Isolation for the parsing codeblock.
    Returns:
    * cnv_metrics_parser: parser for CNV data stored in *.cnv_metrics.csv files.
    * data with found info/warnings/errors."""

    '''"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    The log_data is the container for found info, warnings and errors.
    Collected data is stored in groups to make the output log file more readable.

    warnings:
    - invalid_file_names, which do not conform to the:
        <output prefix>.cnv_metrics.csv

    - invalid_file_lines, which do not conform to:
        SEX GENOTYPER,,<sample>,<value1>,<value2> or
        CNV SUMMARY,,<metric>,<value1> or
        CNV SUMMARY,,<metric>,<value1>,<value2>

    debug/info:
    - unusual_values are those except for int and float in the CNV SUMMARY.
      And those except for float in SEX GENOTYPER's second value.
    '''
    log_data = {
        "invalid_file_names": defaultdict(list),
        "invalid_file_lines": defaultdict(lambda: defaultdict(list)),
        "unusual_values": defaultdict(lambda: defaultdict(dict)),
    }

    def make_metric_id(metric):
        """Tries to guarantee consistency to a certain degree.
        Metric-specific peculiarities can be handled here."""

        # Single backslashes are escaped if presented.
        metric_id = re.sub("\s+", " ", metric).strip().lower()

        return metric_id

    def make_header(metric, section):
        """Creates a single header for the given metric with section."""
        configs = {
            config: val for config, val in SINGLE_HEADER.items() if val is not None
        }

        if metric in METRICS[section]:
            configs.update(METRICS[section][metric])

        # If not, set the bare minimum.
        else:
            configs["title"] = metric
            configs["description"] = metric

        return configs

    # Used to order found panel samples.
    order_counter = RESERVED_START

    def make_panel_headers(sample):
        """Creates headers for both values for the given panel sample."""

        nonlocal order_counter

        configs = {
            config: val for config, val in SINGLE_HEADER.items() if val is not None
        }

        configs["order_priority"] = order_counter
        configs["title"] = "Kar " + sample
        configs["colour"] = "0, 255, 255"
        configs["description"] = (
            "Estimated sex karyotype for the panel sample " + sample
        )
        configs["cond_formatting_rules"] = {
            "red": [{"s_contains": ""}],  # Each value is red by default.
            "green": [{"s_eq": "XX"}, {"s_eq": "YX"}, {"s_eq": "XY"}],
        }
        configs["cond_formatting_colours"] = [
            {"red": "#FF0000"},
            {"green": "#00FF00"},
        ]

        order_counter += 1

        configs2 = {
            config: val for config, val in SINGLE_HEADER.items() if val is not None
        }
        configs2["order_priority"] = order_counter
        configs2["title"] = "KarConf " + sample
        configs2["colour"] = "0, 255, 255"
        configs2["description"] = (
            "Confidence metric for the panel sample "
            + sample
            + " ranging from 0.0 to 1.0. If the sample sex is specified, this metric is 0.0."
        )
        configs2["suffix"] = " %"
        configs2["format"] = "{:,.4f}"
        configs2["max"] = 1.0
        configs2["cond_formatting_rules"] = {
            "red": [{"lt": 0}, {"gt": 1}],
        }
        configs2["cond_formatting_colours"] = [
            {"red": "#FF0000"},
        ]

        order_counter += 1

        return configs, configs2

    # Precompile file and line regexes to speed up execution.
    FILE_RGX = re.compile(r"^([^.]+)\.cnv_metrics\.csv$")
    LINE_RGX = re.compile("^(SEX GENOTYPER|CNV SUMMARY),,([^,]+),([^,]+),?([^,]*)$")

    def cnv_metrics_parser(file_handler, cnv_headers):
        """Parser for *.cnv_metrics.csv files.
        Input:  file_handler with necessary info - file name/content/root.
                cnv_headers will be filled with newfound metrics.
        Output: {"success": 0} if file name test failed. Otherwise:
                {"success": 0/1,
                 "sample_name": <output prefix>,
                 "data": extracted_metrics_with_values
                }, with success 0 if all file's lines are invalid."""

        root_name = file_handler["root"]
        file_name = file_handler["fn"]

        # Official structure of files: *.cnv_metrics.csv
        # Accepted structure: <output prefix>.cnv_metrics.csv
        file_match = FILE_RGX.search(file_name)

        if not file_match:
            log_data["invalid_file_names"][root_name].append(file_name)
            return {"success": 0}

        sample = file_match.group(1)
        success = 0
        data = {}

        for line in file_handler["f"].splitlines():
            # Check the general structure.
            line_match = LINE_RGX.search(line)

            # If line is fine then extract the necessary fields.
            if line_match:
                success = 1
                section, metric, value1, value2 = (
                    line_match.group(1),
                    line_match.group(2),
                    line_match.group(3),
                    line_match.group(4),
                )

                if section == "CNV SUMMARY":
                    metric_id = make_metric_id(metric)

                    """ Check the extracted values. These are int/float.
                    Type conversion is only performed to int and float.
                    Other types of data are left unchanged as strings."""
                    try:
                        value1 = int(value1)
                    except ValueError:
                        try:
                            value1 = float(value1)
                        except ValueError:
                            log_data["unusual_values"][root_name][file_name][metric] = value1

                    data[metric_id] = value1

                    if metric_id not in cnv_headers:
                        cnv_headers[metric_id] = make_header(metric_id, CNV_SUMMARY_SECTION)

                    if value2:
                        try:
                            value2 = float(value2)
                        except ValueError:
                            try:
                                value2 = int(value2)
                            except ValueError:
                                log_data["unusual_values"][root_name][file_name][metric] = value2

                        metric_id += V2
                        data[metric_id] = value2
                        if metric_id not in cnv_headers:
                            cnv_headers[metric_id] = make_header(metric_id, CNV_SUMMARY_SECTION)

                # Else "SEX GENOTYPER"
                else:
                    """
                    SEX GENOTYPER,,SAMPLE_1,XX,0.9842

                    Acc. to the Illumina docs this section currently has only one metric.
                    The metric itself is presumably just an (im)proper subset of a sample name.
                    Besides the case sample a panel of normals samples can be also included.
                    Only the case samples can be stored in one column. Panel of normals are
                    just some random strings, hence can not be combined in the same column.
                    """

                    # Check if case sample.
                    if metric in sample:
                        metric_id = "Karyotype case sample"
                        if metric_id not in cnv_headers:
                            cnv_headers[metric_id] = make_header("karyotype", SEX_GENOTYPER_SECTION)

                        if metric_id + V2 not in cnv_headers:
                            cnv_headers[metric_id + V2] = make_header("confidence", SEX_GENOTYPER_SECTION)

                    # Panel of normals.
                    else:
                        metric_id = "Karyotype panel sample" + metric
                        configs, configs2 = make_panel_headers(metric)
                        if metric_id not in cnv_headers:
                            cnv_headers[metric_id] = configs

                        if metric_id + V2 not in cnv_headers:
                            cnv_headers[metric_id + V2] = configs2

                    data[metric_id] = value1

                    # The confidence shall be float.
                    if value2:
                        try:
                            value2 = float(value2)
                        except ValueError:
                            try:
                                value2 = int(value2)
                            except ValueError:
                                log_data["unusual_values"][root_name][file_name][metric + V2] = value2

                        data[metric_id + V2] = value2

            # Otherwise check if line is empty. If not then report it and go to the next line.
            else:
                if not re.search("^\s*$", line):
                    log_data["invalid_file_lines"][root_name][file_name].append(line)
                continue

        return {"success": success, "sample_name": sample, "data": data}

    return cnv_metrics_parser, log_data


cnv_parser, log_data = construct_coverage_parser()
