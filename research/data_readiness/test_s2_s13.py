"""Actual S2/S13 preservation checks, plus a synthetic rounding edge case."""
import hashlib
import json
from pathlib import Path
import shutil
import tempfile
import unittest

from research.data_readiness import extract_s2_s13 as extract
from research.data_readiness.prepare_observations import read_csv


class SupplementTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory()
        cls.addClassCleanup(cls.temp.cleanup)
        cls.base = Path(cls.temp.name)
        cls.out = cls.base / 'extracted'
        cls.before = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in extract.RAW.iterdir() if p.is_file()}
        cls.report = extract.build(output=cls.out)
        cls.books, _ = extract.verified_books(extract.RAW)

    def rows(self, name):
        return read_csv(self.out / (name + '.csv'))

    def test_all_s2_cells_and_source_rows_are_retained(self):
        sheet = self.books['pnas.2521522122.sd02.xlsx']['Dataset S2']
        rows = self.rows('s2_genomes')
        self.assertEqual(len(rows), 2628)
        self.assertEqual(len({r['source_row'] for r in rows}), 2628)
        schema = json.loads((self.out / 's2_columns.json').read_text())
        self.assertEqual(schema['amino_acid_identity_to_bg0019_percent']['group_header_cell'], 'P3')
        self.assertEqual(schema['gtdb_species_second']['group_header_cell'], 'J3')
        for row in rows:
            for c, field in enumerate(extract.S2_FIELDS, 2):
                value = sheet.cell(int(row['source_row']), c).value
                self.assertEqual(row[field], '' if value is None else str(value), (row['source_row'], field))

    def test_whitespace_collisions_are_reported_without_merging_records(self):
        self.assertEqual(self.report['s2']['raw_unique_identifiers'], 2628)
        self.assertEqual(self.report['s2']['whitespace_modified_search_labels'], 1611)
        self.assertEqual(self.report['s2']['stripped_collision_groups'], 6)
        rows = self.rows('s2_whitespace_collisions')
        self.assertEqual(json.loads(rows[0]['source_rows']), [8, 1212])
        self.assertNotEqual(*json.loads(rows[0]['raw_identifiers']))

    def test_all_six_target_records_keep_literal_identity_and_type(self):
        rows = self.rows('s2_target_identities')
        self.assertEqual([int(r['source_row']) for r in rows], [49, 271, 544, 561, 562, 2455])
        self.assertEqual(rows[-1]['genome'], 'Segatella_copri_BgF5_2 ')
        self.assertEqual([r['genome_type_or_purpose'] for r in rows], ['reference_genome', 'reference_genome', 'isolate', 'mag', 'mag', 'isolate'])
        self.assertEqual(self.report['s2']['genome_types_literal'], {'isolate': 74, 'reference_genome': 2043, 'mag': 484, 'MAG': 27})

    def test_f5_annotations_and_gtdb_versions_are_not_reconciled(self):
        rows = {r['source_row']: r for r in self.rows('s2_target_identities')}
        self.assertEqual((rows['544']['bifunctional_presence_first_pass'], rows['2455']['bifunctional_presence_first_pass']), ('1', '0'))
        self.assertEqual(rows['2455']['bifunctional_presence_second_pass'], '')
        self.assertEqual(rows['2455']['amino_acid_identity_to_bg0019_percent'], '')
        self.assertEqual(rows['544']['amino_acid_identity_to_bg0019_percent'], '92.884')
        self.assertEqual((rows['561']['gtdb_species_first'], rows['561']['gtdb_species_second']), ('copri_C', 'copri_N'))
        self.assertEqual((rows['544']['gtdb_version_first'], rows['544']['gtdb_version_second']), ('r220', 'r226'))

    def test_all_s13_measurements_labels_and_units_match_source_cells(self):
        sheet = self.books['pnas.2521522122.sd13.xlsx']['Dataset S13']
        rows = self.rows('s13_measurements')
        self.assertEqual(len(rows), 648)
        self.assertEqual(len({r['value_cell'] for r in rows}), 648)
        for row in rows:
            for field, anchor in [('value', 'value_cell'), ('analyte_raw', 'analyte_cell'), ('units_raw', 'units_cell'), ('food_input', 'diet_cell'), ('replicate_label', 'replicate_cell')]:
                self.assertEqual(row[field], str(sheet[row[anchor]].value))
        self.assertEqual(next(r for r in rows if r['value_cell'] == 'D56')['analyte_raw'], 't-Glucose ')

    def test_diet_grids_are_complete_without_invented_mouse_or_batch_join(self):
        for section, n in [('total_monosaccharides', 126), ('free_monosaccharides', 126), ('linkages', 396)]:
            report = self.report['s13'][section]
            self.assertEqual(report['records'], n)
            self.assertEqual(report['finite'], n)
            self.assertEqual((report['missing'], report['zero'], report['negative']), (0, 0, 0))
            self.assertEqual(report['replicate_labels'], [1, 2, 3])
        self.assertFalse({'mouse_id', 'experiment', 'dose', 'batch'} & self.rows('s13_measurements')[0].keys())

    def test_all_published_s13_summaries_reproduce_mean_and_sample_sd(self):
        report = self.report['s13_summary_checks']
        self.assertEqual((report['groups'], report['mean_matches'], report['sample_sd_matches']), (216, 216, 216))
        self.assertEqual({r['n_reported_replicates'] for r in self.rows('s13_summary_diagnostics')}, {'3'})

    def test_published_pvalues_are_retained_without_classification(self):
        sheet = self.books['pnas.2521522122.sd13.xlsx']['Dataset S13']
        rows = self.rows('s13_published_pvalues')
        self.assertEqual(len(rows), 72)
        for row in rows:
            self.assertEqual(row['reported_value'], str(sheet[row['value_cell']].value))
        self.assertFalse(any('significant' in k or 'threshold' in k for k in rows[0]))

    def test_scientific_notation_rounding_does_not_use_large_absolute_tolerance(self):
        values, tolerances = extract.parse_summary('2.38E-01 ± 6.38E-02')
        self.assertEqual(values, [0.238, 0.0638])
        self.assertEqual(tolerances, [0.0005, 0.00005])
        vals = [0.9e-12, 1e-12, 1.1e-12]
        self.assertTrue(extract.summarize_source_values(vals, '1.00E-12 ± 1.00E-13')['calculated_mean'] > 0)
        self.assertTrue(extract.summarize_source_values(vals, '1.00E-12 ± 1.00E-13')['mean_matches_printed_precision'])
        self.assertFalse(extract.summarize_source_values(vals, '1.05E-12 ± 1.00E-13')['mean_matches_printed_precision'])
        with self.assertRaises(ValueError):
            extract.summarize_source_values([1, None, 3], '2 ± 1')

    def test_tampered_original_is_rejected_before_extraction(self):
        altered = self.base / 'altered'
        shutil.copytree(extract.RAW, altered)
        path = altered / 'pnas.2521522122.sd13.xlsx'
        path.write_bytes(path.read_bytes() + b'changed')
        with self.assertRaisesRegex(ValueError, 'source integrity failure'):
            extract.build(raw=altered, output=self.base / 'rejected')

    def test_repeat_outputs_are_identical_and_originals_unchanged(self):
        second = self.base / 'second'
        extract.build(output=second)
        self.assertEqual({p.name: p.read_bytes() for p in self.out.iterdir()}, {p.name: p.read_bytes() for p in second.iterdir()})
        self.assertEqual(self.before, {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in extract.RAW.iterdir() if p.is_file()})


if __name__ == '__main__':
    unittest.main()
