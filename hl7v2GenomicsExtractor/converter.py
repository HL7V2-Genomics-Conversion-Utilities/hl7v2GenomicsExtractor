import xml.etree.ElementTree as ET
import xmltodict
import pyranges
from .hl7v2_message_generator import _get_hl7v2_message
import logging
import sys
import pandas as pd
from .common import *
general_logger = logging.getLogger('hl7v2GenomicsExtractor.general')
"""Hl7 V2 Genomics Extractor"""


class Converter(object):
    """Creates a new Converter Object.

    Parameters
    ----------

    **filename** (required): Path to a file or qci xml file.
    Valid path and filename without whitespace must be provided.
    
    **ref_build** (required): Genome Reference Consortium genome
    assembly. Must be one of 'GRCh37' or 'GRCh38'.

    **patient_id** (optional): Supplied patient ID is inserted into
    generated HL7V2 output. Alphanumeric string without whitespace. if \
    not provided, header of first sample column is used.

    **Conversion region** (optional): Subset of the file to be
    converted into HL7V2. If absent, the entire file is converted. Can
    be supplied as either a parameter (conv_region_dict) or as a BED
    file (conv_region_filename):

       **conv_region_dict** : Array of regions (e.g. '{"Chromosome":
       ["X", "X", "M"],"Start": [50000, 55000, 50000],"End": [52000,
       60600, 60025]}'). Ranges must be `0-based \
       <https://www.biostars.org/p/84686/>`_
       (or 0-start, half-open) and based on GRCh37 or GRCh38 \
       reference sequences.

       **conv_region_filename**: Valid path and filename without
       whitespace must be provided. Must be a valid BED file with first 3
       columns: <chr> <start> <stop>. Values in <chr> field must align
       with values in #CHROM field. Ranges must be based on GRCh37 or
       GRCh38 reference sequences.

    **region_studied_filename** (optional): Subset of patient's genome
    that was studied in the generation of the file. Valid path and
    filename without whitespace must be provided. Must be a valid BED
    file, with first 3 columns: <chr> <start> <stop>. Values in <chr>
    field must align with values in #CHROM field. Ranges must be
    based on GRCh37 or GRCh38 reference sequences.

    **source_class** (optional)(default value = germline): An \
    assertion as to whether variants in the file are in the \
    germline (i.e. inherited), somatic (e.g. arose spontaneously \
    as cancer mutations), or mixed (i.e. may be a combination of \
    germline and/or somatic).

    **seed** (optional)(default value = 1000): Used to set the \
    starting integer count for OBX-1 (sequence number)

    **variant_analysis_method** (optional)(default value = Sequencing): Used \
    to set the variant analysis method.

    **report_filename** (optional): Path to a text-based report file.

    Returns

    -------

    Object

    An Instance of Converter.

    """

    def __init__(self, filename=None, ref_build=None, patient_id=None,
                 conv_region_filename=None, conv_region_dict=None,
                 region_studied_filename=None,
                 source_class=None, seed=1000,
                 variant_analysis_method="Sequencing", report_filename=None):

        super(Converter, self).__init__()
        if not filename:
            raise Exception('You must provide a file')

        if not ref_build or ref_build not in ["GRCh37", "GRCh38"]:
            raise Exception(
                'You must provide build number ("GRCh37" or "GRCh38")')

        if not validate_filename(filename):
            raise Exception('Either filename or extension is not correct')

        if is_xml_file(filename):
            try:
                with open(filename, encoding='utf-8') as fd:
                    xml_dict = xmltodict.parse(fd.read())
                    report = xml_dict['report']
                self._xml_reader = report
            except FileNotFoundError:
                raise
            except BaseException:
                self._generate_exception("Please provide valid 'filename'")

        self.report = None
        if report_filename and is_txt_file(report_filename):
            try:
                with open(report_filename, encoding='utf-8') as fd:
                    self.report = fd.readlines()
            except FileNotFoundError:
                raise
            except BaseException:
                self._generate_exception(
                    "Please provide valid 'report_filename'")

        if not patient_id and self._xml_reader is not None:
            patient_id = 'patient_id'

        if conv_region_filename:
            try:
                self.conversion_region = pyranges.read_bed(
                    conv_region_filename)
            except FileNotFoundError:
                raise
            except BaseException:
                self._generate_exception(
                    "Please provide valid 'conv_region_filename'")
        elif conv_region_dict:
            try:
                self.conversion_region = pyranges.from_dict(conv_region_dict)
            except FileNotFoundError:
                raise
            except BaseException:
                self._generate_exception(
                    "Please provide valid 'conv_region_dict'")
        else:
            self.conversion_region = None

        if region_studied_filename:
            try:
                self.region_studied = pyranges.read_bed(
                    region_studied_filename)
            except FileNotFoundError:
                raise
            except BaseException:
                self._generate_exception(
                    "Please provide valid 'region_studied_filename'")
        else:
            self.region_studied = None

        if not validate_seed(seed):
            raise Exception("Please provide a valid seed")

        if(source_class is not None and
           source_class.title() not in Genomic_Source_Class.set_()):
            raise Exception(
                'Please provide a valid Source Class ("germline" or "somatic")'
            )

        if variant_analysis_method.lower() not in SEQUENCING_TO_CODE.keys():
            raise Exception(
                'Please provide a valid variant analysis method'
            )

        self.patient_id = patient_id
        self.ref_build = ref_build
        self.conv_region_filename = conv_region_filename
        self.source_class = source_class
        self.seed = seed
        self.variant_analysis_method = variant_analysis_method
        general_logger.info("Converter class instantiated successfully")

    def convert(self, output_filename='hl7v2.txt'):
        """ Generates HL7 V2 format data as output_filename \
        or hl7v2.txt if it is not provided

        Parameters
        ----------
        output_filename:
            Path to output hl7v2 text file.

        """
        general_logger.info("Starting HL7V2 Conversion")
        hl7v2_oru_message =\
            _get_hl7v2_message(
                self._xml_reader, self.ref_build,
                self.patient_id, self.conversion_region,
                self.source_class, self.region_studied,
                self.seed,
                self.variant_analysis_method, self.report, output_filename)
        general_logger.info("Completed HL7V2 Conversion")
        return hl7v2_oru_message

    def _generate_exception(self, msg):
        general_logger.error(msg, exc_info=True)
        raise Exception(msg, sys.exc_info)
