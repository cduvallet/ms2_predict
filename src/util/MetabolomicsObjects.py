#!/usr/bin/env python
"""
This file contains class definitions for various objects used in metabolomics
data processing.
"""

class MS2:
    '''
    A container for all MS2 metadata. Final MS2 objects may have empty
    or additional fields, depending on any anomalies in the data.
    '''
    def __init__(self):
        self.peak_counter = None
        self.mono_mass = None
        self.spectra_type = None
        self.inchi_key = None
        self.frequency = None
        self.instrument_type = None
        self.energy_field = None
        self.base_peak = None
        self.sample_mass_units = None
        self.chromatography_type = None
        self.searchable = None
        self.sample_mass = None
        self.derivative_mw = None
        self.retention_time = None
        self.updated_at = None
        self.sample_assessment = None
        self.derivative_formula = None
        self.derivative_type = None
        self.database_id = None
        self.ref_text = None
        self.mass_charge = None
        self.collision_energy_voltage = None
        self.sample_concentration = None
        self.spectra_assessment = None
        self.solvent = None
        self.nucleus_y = None
        self.ionization_mode = None
        self.collection_date = None
        self.nucleus = None
        self.sample_temperature_units = None
        self.sample_ph = None
        self.spectra_id = None
        self.c_ms_id = None
        self.sample_concentration_units = None
        self.sample_source = None
        self.nil_classes = None
        self.nucleus_x = None
        self.database = None
        self.notes = None
        self.created_at = None
        self.sample_temperature = None
        self.ri_type = None
        self.pubmed_id = None
        self.molecule_id = None
        self.column_type = None
        self.retention_index = None
        self.collision_energy_level = None
        self.references = None
        self.name = None
        self.accession = None
        self.chemical_formula = None
        self.monoisotopic_molecular_weight = None
        self.iupac_name = None
        self.traditional_iupac = None
        self.cas_registry = None
        self.smiles = None
        self.inchi = None
        self.inchikey = None
        self.taxonomy = None
        self.biofluid_locations = None
        self.id_dict = None
        self.peaks = None

    def __str__(self):
        return 'MS2'


class MSPeak:
    '''
    Container for data about an ms-ms peak. Final MSPeak objects may have
    empty or additional fields, depending on any anomalies in the data.
    '''
    def __init__(self):
        self.id = None
        self.ms_ms_id = None
        self.intensity = None
        self.mass_charge = None

    def __str__(self):
        return 'MSPeak'


class Metabolite:
    '''
    Container for metadata about a metabolite. Final Metabolite objects may have
    empty or additional fields, depending on any anomalies in the data.
    '''
    def __init__(self):
        self.accession = None
        self.secondary_accessions = None
        self.name = None
        self.chemical_formula = None
        self.monisotopic_molecular_weight = None
        self.iupac_name = None
        self.traditional_iupac = None
        self.cas_registry = None
        self.smiles = None
        self.inchi = None
        self.inchikey = None
        self.biofluid_locations = None
        self.id_dict = None
        self.taxonomy_dict = None
        self.MS2 = []

    def __str__(self):
        return 'Metabolite'
