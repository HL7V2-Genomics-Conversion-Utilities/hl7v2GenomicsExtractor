import pandas as pd
import logging
import re
from enum import Enum

general_logger = logging.getLogger("hl7v2GenomicsExtractor.general")

GERMLINE = 'Germline'
SOMATIC = 'Somatic'

CLIN_SIG_TO_CODE = {
    'pathogenic': 'LA6668-3',
    'likely pathogenic': 'LA26332-9',
    'uncertain significance': 'LA26333-7',
    'likely benign': 'LA26334-5',
    'benign': 'LA6675-8'
}

GENOTYPE_TO_ALLELIC_STATE = {
    "Het": "LA6706-1^Heterozygous^LN",
    "Hom": "LA6705-3^Homozygous^LN"
}

# Translation Impact to Molecular Consequence
TI_TO_MC_CODE = {
    "wild type": "LA9658-1",
    "deletion": "LA6692-3",
    "duplication": "LA6686-5",
    "frameshift": "LA6694-9",
    "initiating methionine": "LA6695-6",
    "insertion": "LA6687-3",
    "insertion and deletion": "LA9659-9",
    "missense": "LA6698-0",
    "nonsense": "LA6699-8",
    "silent": "LA6700-4",
    "stop codon mutation": "LA6701-2"
}

SEQUENCING_TO_CODE = {
    "sequencing": "LA26398-0",
    "oligo acgh": "LA26399-8",
    "snp array": "LA26400-4",
    "bac acgh": "LA26401-2",
    "curated": "LA26402-0",
    "digital array": "LA26403-8",
    "fish": "LA26404-6",
    "gene expression array": "LA26405-3",
    "karyotyping": "LA26406-1",
    "maph": "LA26407-9",
    "maldi-tof": "LA26408-7",
    "merging": "LA26808-8",
    "multiple complete digestion": "LA26414-5",
    "mlpa": "LA26415-2",
    "optical mapping": "LA26417-8",
    "pcr": "LA26418-6",
    "qpcr (real-time PCR)": "LA26419-4",
    "roma": "LA26420-2",
    "denaturing high-pressure liquid chromatography (dhplc)": "LA26809-6",
    "dna hybridization": "LA26810-4",
    "computational analysis": "LA26811-2",
    "single-stranded conformational polymorphism (sscp)": "LA26812-0",
    "restriction fragment length polymorphism (rflp)": "LA26813-8"
}


class Genomic_Source_Class(Enum):

    @classmethod
    def set_(cls):
        return set(map(lambda c: c.value, cls))

    GERMLINE = GERMLINE
    SOMATIC = SOMATIC


def get_variant_name_xml(variant):
    if(is_present_xml(variant, 'transcriptchange') and
       is_present_xml(variant.get('transcriptchange'), 'transcript') and
       is_present_xml(variant.get('transcriptchange'), 'change') and
       is_present_xml(variant, 'proteinchange') and
       is_present_xml(variant.get('proteinchange'), 'change')):
        display_name =\
            (f'{variant.get("transcriptchange").get("transcript")}' +
             f'({variant.get("gene")}):' +
             f'{variant.get("transcriptchange").get("change")} ' +
             f'({variant.get("proteinchange").get("change")})')
    elif(is_present_xml(variant, 'transcriptchange') and
         is_present_xml(variant.get('transcriptchange'), 'transcript') and
         is_present_xml(variant.get('transcriptchange'), 'change') and
         not is_present_xml(variant, 'proteinchange') and
         not is_present_xml(variant.get('proteinchange'), 'change')):
        display_name =\
            (f'{variant.get("transcriptchange").get("transcript")}' +
             f'({variant.get("gene")}):' +
             f'{variant.get("transcriptchange").get("chnage")}')
    else:
        display_name =\
            (f'{variant.get("chromosome")}:' +
             f'{variant.get("genomicchange").get("change")}')
    return display_name


def extract_chrom_identifier(chrom):
    chrom = chrom.upper().replace("CHR", "")
    if chrom == "MT":
        chrom = "M"
    return chrom


def get_spdi_representation(record, ref_seq):
    return (f'{ref_seq}:{int(record.get("position")) - 1}:' +
            f'{record.get("reference")}:{record.get("alternate")}')


def get_chromosome(variant):
    return str(variant.get("chromosome"))


def get_position(variant):
    return int(variant.get("position"))


def get_reference(variant):
    return str(variant.get("reference"))


def get_alternate(variant):
    return str(variant.get("alternate"))


def get_variant_display_name(record):
    return get_variant_name_xml(record)


def get_as_af(record):
    allelic_state =\
        GENOTYPE_TO_ALLELIC_STATE.get(str(record.get("genotype")))
    allelic_frequency =\
        (float(record.get("allelefraction")) / 100)

    return allelic_state, allelic_frequency


def validate_filename(filename):
    pattern = r'^[\w,\s\-,/,\.]+\.(xml)$'
    result = re.match(pattern, filename)
    return bool(result)


def is_xml_file(filename):
    extension = filename.split('.')[-1]
    if extension == 'xml':
        return True
    else:
        return False


def is_txt_file(filename):
    extension = filename.split('.')[-1]
    if extension == 'txt':
        return True
    else:
        return False


def validate_chrom_identifier(chrom):
    chrom = extract_chrom_identifier(chrom)
    pattern = '^[1-9]$|^1[0-9]$|^2[0-2]$|^[XYM]$'
    result = re.match(pattern, chrom)
    return bool(result)


def validate_seed(seed):
    if isinstance(seed, int) and not isinstance(seed, bool):
        return True
    return False


def get_alphabet_index(n):
    a_index = ""
    while(n != 0):
        remainder = n % 26
        n = n // 26
        if remainder == 0:
            remainder = 26
            n = n - 1
        a_index = chr(96 + remainder) + a_index
    return a_index


def is_present(annotation, component):
    if(not annotation[component].empty and
       not pd.isna(annotation.iloc[0][component])):
        return True
    else:
        return False


def get_annotations(record, spdi_representation):
    return get_annotations_from_xml(record, spdi_representation)
    


def is_present_xml(variant, component):
    if variant is not None and variant.get(component) is not None:
        return True
    else:
        return False


def get_annotations_from_xml(variant, spdi_representation):
    if(is_present_xml(variant, 'transcriptchange') and
       is_present_xml(variant.get('transcriptchange'), 'change') and
       is_present_xml(variant.get('transcriptchange'), 'transcript')):
        dna_change = (f'{variant.get("transcriptchange").get("transcript")}:' +
                      f'{variant.get("transcriptchange").get("change")}')
    elif(is_present_xml(variant, 'transcriptchange') and
         is_present_xml(variant.get('transcriptchange'), 'change') and
         not is_present_xml(variant.get('transcriptchange'), 'transcript')):
        dna_change = f'{variant.get("transcriptchange").get("change")}'
    else:
        dna_change = spdi_representation

    if(is_present_xml(variant, 'proteinchange') and
       is_present_xml(variant.get('proteinchange'), 'change') and
       is_present_xml(variant.get('proteinchange'), 'protein')):
        amino_acid_change =\
            (f'{variant.get("proteinchange").get("protein")}:' +
             f'{variant.get("proteinchange").get("change")}')
    elif(is_present_xml(variant, 'proteinchange') and
         is_present_xml(variant.get('proteinchange'), 'change') and
         not is_present_xml(variant.get('proteinchange'), 'protein')):
        amino_acid_change = f'{variant.get("proteinchange").get("change")}'
    else:
        amino_acid_change = None

    if is_present_xml(variant, 'assessment'):
        assessment = str(variant.get("assessment"))
        code = CLIN_SIG_TO_CODE.get(assessment.lower())
        clin_sig = f'{code}^{assessment}^LN'
    else:
        clin_sig = "not specified"

    phenotype = None
    phenotype_description = None

    if is_present_xml(variant, 'gene'):
        gene_studied = f'^{variant.get("gene")}^HGNC'
    else:
        gene_studied = 'HGNC:0000^NoGene^HGNC'

    if(is_present_xml(variant, 'transcriptchange') and
       is_present_xml(variant.get('transcriptchange'), 'transcript')):
        transcript_ref_seq =\
            (f'{variant.get("transcriptchange").get("transcript")}^' +
             f'{variant.get("transcriptchange").get("transcript")}^RefSeq-T')
    else:
        transcript_ref_seq = None

    if(is_present_xml(variant, "genomicchange") and
       is_present_xml(variant.get("genomicchange"), "change")):
        genomic_dna_change = f'{variant.get("genomicchange").get("change")}'
    else:
        genomic_dna_change = None

    if(is_present_xml(variant, "proteinchange") and
       is_present_xml(variant.get("proteinchange"), "translationimpact")):
        translation_impact =\
            str(variant.get("proteinchange").get("translationimpact"))
        code = TI_TO_MC_CODE.get(translation_impact.lower())
        translation_value = translation_impact

        if translation_impact.lower() == "in-frame":
            translation_impact_updated = "".join(translation_impact.split("-"))
            variation = str(variant.get("variation")).lower()
            if variation == "substitution":
                variation = "insertion and deletion"
            translation_code = " ".join(
                [translation_impact_updated, variation])
            translation_value = translation_code
        molecular_consequence = f"{code}^{translation_value}^LN"
    else:
        molecular_consequence = None

    dna_region = None

    if(is_present_xml(variant, "proteinchange") and
       is_present_xml(variant.get("proteinchange"), "protein")):
        protein_ref_seq = f'{variant.get("proteinchange").get("protein")}'
    else:
        protein_ref_seq = None

    if is_present_xml(variant, "dbsnp"):
        db_snp_id = f'{variant.get("dbsnp")}^{variant.get("dbsnp")}^dbSNP'
    else:
        db_snp_id = None

    if is_present_xml(variant, "readDepth"):
        read_depth = f'{int(variant.get("readDepth"))}'
    else:
        read_depth = None

    return {'dna_change': dna_change,
            'amino_acid_change': amino_acid_change,
            'clin_sig': clin_sig, 'phenotype': phenotype,
            'gene_studied': gene_studied,
            'transcript_ref_seq': transcript_ref_seq,
            'molecular_consequence': molecular_consequence,
            'genomic_dna_change': genomic_dna_change,
            'dna_region': dna_region,
            'protein_ref_seq': protein_ref_seq, 'db_snp_id': db_snp_id,
            'phenotype_description': phenotype_description,
            'read_depth': read_depth}


def _error_log_allelicstate(record):
    general_logger.error(
        "Cannot Determine AllelicState for: %s , considered sample: %s",
        record,
        record.samples[0].data)
