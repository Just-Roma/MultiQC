'''"""""""""""""""""""""""""""""""""""""""""""""""""""""""
This module gathers data from _overall_mean_cov.csv files.
It relies on the following official source:
Overall Mean Coverage Report, page 196.
https://emea.support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/dragen-bio-it/Illumina-DRAGEN-Bio-IT-Platform-User-Guide-1000000141465-00.pdf
'''

import logging
import re
from collections import defaultdict

from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger.
log = logging.getLogger(__name__)


class DragenOverallMeanCovMetrics(BaseMultiqcModule):
    """Public members of the DragenOverallMeanCovMetrics module are available to
    other dragen modules. Defines 2 public functions and a private data holder.
    """

    @property
    def overall_mean_cov_data(self):
        """Getter for data from _overall_mean_cov.csv files."""
        if self._overall_mean_cov_data:
            return self._overall_mean_cov_data
        else:
            return None

    def collect_overall_mean_cov_data(self):
        """Collects raw data for coverage_metrics.py from overall_mean_cov.csv files."""
        self._overall_mean_cov_data = defaultdict(lambda: defaultdict(dict))
        for file in self.find_log_files("dragen/overall_mean_cov_metrics"):
            out = parse_overall_mean_cov(file)
            if out["success"]:
                # No need to call add_data_source, because the collected data is used
                # only by the coverage_metrics.py module and will not be included in html.
                # self.add_data_source(file, section="stats")

                # Data for the coverage_metrics.py module.
                if out["phenotype"]:
                    root, sample, phenotype, data = (
                        out["root"],
                        out["sample"],
                        out["phenotype"],
                        out["data"],
                    )
                    self._overall_mean_cov_data[root][sample][phenotype] = data
                # Currently there is no need to support other files. Pass for now.
                else:
                    pass
        # No need to write the data.
        # self.write_data_file(self._overall_mean_cov_data)

        # Just to pass the pre-commit
        # doi=


# Official structure of files:   _overall_mean_cov.csv
# Accepted structure of files: .+_overall_mean_cov.*.csv
# A suffix after _overall_mean_cov is used to catch extra substrings.
GEN_FILE_RGX = re.compile(r"(.+?)_overall_mean_cov(.*)\.csv$")

# Special case. Coverage metrics files have the following structure:
# <output prefix>.<coverage region prefix>_overall_mean_cov<arbitrary suffix>.csv
COV_FILE_RGX = re.compile(r"^([^.]+)\.(.+?)_overall_mean_cov(.*)\.csv$")

# General structure of metrics is not defined.
# Currently only 1 metric is presented in the standard:
AVG_RGX = re.compile("Average alignment coverage over ([^,]+),([^,]+)$", re.IGNORECASE)


def parse_overall_mean_cov(file_handler):
    """Parser for _overall_mean_cov.csv files."""
    root = file_handler["root"]
    file = file_handler["fn"]
    # First check the general structure.
    file_match = GEN_FILE_RGX.search(file)
    if file_match:
        sample, phenotype = file_match.group(1), None
        file_match = COV_FILE_RGX.search(file)
        if file_match:
            sample, phenotype, normal_tumor = file_match.group(1), file_match.group(2), file_match.group(3)
            if normal_tumor:
                # Only non-'tumor' strings are concatenated with the sample.
                if not re.search("tumor", normal_tumor, re.IGNORECASE):
                    sample += normal_tumor
        else:
            log.debug("\nNot supported: " + file + "\nin: " + root)
    # Else there is no name at all, report and return 0.
    else:
        log.debug("\nFile name is empty: " + file + "\nin: " + root)
        return {"success": 0}

    # There is still a possibility for failing(eg empty file, unknown metric).
    success = 0
    data = {}
    source = None  # File name from "Average alignment coverage over file" metric.
    value = None

    # A loop is created in advance. Maybe new metrics will be added in the future.
    for line in file_handler["f"].splitlines():
        line_match = AVG_RGX.search(line)
        # If line is fine then extract the necessary fields.
        if line_match:
            success = 1
            source, value = line_match.group(1), line_match.group(2)

        # Otherwise check if line is empty. If not then report it and go to the next line.
        else:
            if not re.search("^\s*$", line):
                log.debug("\nUnsupported metric: " + line + "\n" + file + "\nin: " + root)
            continue

        # Try to convert the value. It shall be float.
        try:
            value = float(value)
        except ValueError:
            try:
                value = int(value)
            except ValueError:
                log.debug("\nNon int/float value found: " + line + "\n" + file + "\nin: " + root)

    data = {"source_file": source, "value": value}

    return {
        "success": success,
        "root": root,
        "sample": sample,
        "phenotype": phenotype,
        "data": data,
    }