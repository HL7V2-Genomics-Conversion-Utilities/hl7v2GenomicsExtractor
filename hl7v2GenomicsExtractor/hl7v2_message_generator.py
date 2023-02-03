import logging
from .gene_ref_seq import _get_ref_seq_by_chrom
from .hl7v2_helper import _HL7V2_Helper
from .common import *

invalid_record_logger = logging.getLogger("hl7v2GenomicsExtractor.invalidrecord")
general_logger = logging.getLogger("hl7v2GenomicsExtractor.general")


def _get_chrom(chrom_index):
    switcher = {
        23: 'X',
        24: 'Y',
        25: 'M'
    }
    return switcher.get(chrom_index, str(chrom_index))


def _fix_regions_chrom(region):
    if region:
        region.Chromosome = region.Chromosome.apply(
            extract_chrom_identifier)


def _add_record_variants(
        record, ref_seq, patientID,
        hl7v2_helper, source_class,
        ref_build, variant_analysis_method):
    spdi_representation = get_spdi_representation(record, ref_seq)
    annotation_record =\
        get_annotations(record, spdi_representation)

    if annotation_record is not None:
        hl7v2_helper.add_variant_obv(record, ref_seq,
                                     source_class, annotation_record,
                                     spdi_representation, ref_build,
                                     variant_analysis_method)


def _add_region_studied(
        region_studied, hl7v2_helper, chrom, ref_seq, patientID):
    if(region_studied and
       not region_studied[chrom].empty):
        general_logger.info("Adding region Studied OBX for %s", chrom)
        general_logger.debug("Region Examined %s", region_studied[chrom])

        hl7v2_helper.add_regionstudied_obv(
            ref_seq, region_studied[chrom])


def create_region_studied_obxs(
        xml_reader, ref_build, patientID,
        conversion_region, region_studied, source_class,
        report, output_filename, hl7v2_helper):
    _fix_regions_chrom(conversion_region)
    _fix_regions_chrom(region_studied)

    if conversion_region:
        if region_studied:
            region_studied = region_studied.intersect(conversion_region)
            general_logger.debug(
                "Final Conmputed Reportable Query Regions: %s", region_studied)
    general_logger.info("Start adding the Region studied OBXs")

    if xml_reader.get("variant") is not None:
        if not isinstance(xml_reader.get("variant"), list):
            xml_reader["variant"] = [xml_reader.get("variant")]
        variants = xml_reader.get("variant")
    else:
        return

    chrom_index = 1
    prev_add_chrom = ""
    for variant in variants:
        CHROM = extract_chrom_identifier(get_chromosome(variant))
        if(not (conversion_region and
           conversion_region[CHROM].empty)):
            if prev_add_chrom != CHROM and (
                    region_studied):
                chrom = _get_chrom(chrom_index)
                while prev_add_chrom != CHROM:
                    current_ref_seq = _get_ref_seq_by_chrom(
                        ref_build, chrom)
                    _add_region_studied(
                        region_studied, hl7v2_helper, chrom,
                        current_ref_seq, patientID)
                    prev_add_chrom = chrom
                    chrom_index += 1
                    chrom = _get_chrom(chrom_index)
    hl7v2_helper.add_final_region_studied_obx_segment(report)


def create_variant_obxs(
        xml_reader, ref_build, patientID,
        conversion_region, source_class,
        variant_analysis_method, output_filename, hl7v2_helper):
    _fix_regions_chrom(conversion_region)

    general_logger.info("Start adding the Variant OBXs")

    if xml_reader.get("variant") is not None:
        if not isinstance(xml_reader.get("variant"), list):
            xml_reader["variant"] = [xml_reader.get("variant")]
            print(xml_reader["variant"])
        variants = xml_reader.get("variant")
    else:
        return

    for variant in variants:
        CHROM = extract_chrom_identifier(get_chromosome(variant))
        if(not (conversion_region and
           conversion_region[CHROM].empty)):
            ref_seq = _get_ref_seq_by_chrom(ref_build, CHROM)
            POS = get_position(variant)
            REF = get_reference(variant)
            end = POS + len(REF) - 1
            if(not conversion_region or
               conversion_region[CHROM, POS - 1: end].empty is False):
                _add_record_variants(
                    variant, ref_seq, patientID, hl7v2_helper,
                    source_class,
                    ref_build,
                    variant_analysis_method)


def _get_hl7v2_message(
        xml_reader, ref_build, patientID,
        conversion_region, source_class, region_studied,
        seed,
        variant_analysis_method, report, output_filename):

    hl7v2_helper = _HL7V2_Helper(patientID, seed)
    general_logger.debug("Finished Initializing empty HL7V2 message")

    create_region_studied_obxs(
        xml_reader, ref_build, patientID,
        conversion_region, region_studied, source_class,
        report, output_filename, hl7v2_helper)

    hl7v2_helper.add_section1_components(xml_reader)

    create_variant_obxs(
        xml_reader, ref_build, patientID,
        conversion_region, source_class,
        variant_analysis_method, output_filename,
        hl7v2_helper)

    if hl7v2_helper.message.obx[-1].obx_4[0].value[0] == "2":
        hl7v2_helper.message.obx[hl7v2_helper.final_rs_index].obx_5 =\
            'Positive'

    general_logger.info(
        f'Export the HL7V2 message object to the file {output_filename}')
    hl7v2_helper.export_hl7v2_message(output_filename)
    general_logger.info("Completed conversion")
    return hl7v2_helper.message
