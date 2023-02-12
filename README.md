## **hl7v2GenomicsExtractor**

### Introduction

HL7 Version 2 messaging format is the predominant means by which labs send structured results to Electronic Health Records (EHRs). EHRs can directly import HL7 V2 formatted results that conform to the [HL7 Lab Results Interface Implementation Guide](https://www.hl7.org/documentcenter/public/standards/dstu/V251_IG_LRI_R1_STU3_2018JUN.pdf). Here, we provide an open source utility for converting variants from VCF/XML format into HL7 V2 format.

### Install
Before installing hl7v2GenomicsExtractor you need to install cython and wheel.
```
pip install cython wheel
```
Now, install hl7v2GenomicsExtractor binary from pip.
```
pip install hl7v2GenomicsExtractor
```

```python
>>> import hl7v2GenomicsExtractor
>>> hl7v2GenomicsExtractor_converter = hl7v2GenomicsExtractor.Converter(filename="sample.xml", ref_build="GRCh37", patient_id=1234, seed=1, source_class="somatic", variant_analysis_method="sequencing")
>>> hl7v2GenomicsExtractor_converter.convert()
```

### License and Limitations

Software is available for use under an [Apache 2.0 license](https://opensource.org/licenses/Apache-2.0), and is intended solely for experimental use, to help further Genomics-EHR integration exploration. Software is expressly not ready to be used with identifiable patient data or in delivering care to patients. Code issues should be tracked here.
