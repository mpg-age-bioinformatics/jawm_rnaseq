"""Run the checksum script against small files without JAWM or reference data."""
import ast
import copy
import gzip
import hashlib
import os
import re
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[1]
TREE = ast.parse((ROOT / 'submodules/_rnaseq.py').read_text())
PROCESS = next(node.value for node in TREE.body if isinstance(node, ast.Assign)
               and any(isinstance(target, ast.Name) and target.id == 'raw_md5sums'
                       for target in node.targets))
SCRIPT = ast.literal_eval(next(item.value for item in PROCESS.keywords if item.arg == 'script'))
HELPER = next(node for node in TREE.body if isinstance(node, ast.FunctionDef)
              and node.name == 'configure_raw_md5sums')


class ProcessStub:
    def __init__(self):
        self.var = {}
        self.script = SCRIPT

    def clone(self):
        return copy.deepcopy(self)


class RawChecksumTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.output = self.root / 'metadata' / "checksums' report.tsv"

    def configure(self, files, md5sums_file=None):
        namespace = {'os': os, 're': re, 'raw_md5sums': ProcessStub()}
        exec(compile(ast.Module(body=[HELPER], type_ignores=[]), '<helper>', 'exec'), namespace)
        return namespace['configure_raw_md5sums'](files, self.output, md5sums_file)

    def execute(self, files, md5sums_file=None):
        process = self.configure(files, md5sums_file)
        script = process.script
        for key, value in process.var.items():
            script = script.replace('{{' + key + '}}', value)
        return subprocess.run([sys.executable, '-c', script], capture_output=True, text=True)

    def test_paths_and_exact_compressed_bytes(self):
        files = [self.root / "one's reads.fastq.gz", self.root / 'two' / 'reads.fastq.gz']
        for index, path in enumerate(files):
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(gzip.compress(b'ACGT' * (index + 1), mtime=0))
        result = self.execute(files + files[:1])
        self.assertEqual(result.returncode, 0, result.stderr)
        rows = [line.split('\t') for line in self.output.read_text().splitlines()]
        self.assertEqual(rows, [[str(path), hashlib.md5(path.read_bytes()).hexdigest()]
                                for path in files])

    def test_missing_input_does_not_publish_partial_report(self):
        first = self.root / 'exists.fastq.gz'
        first.write_bytes(b'first')
        result = self.execute([first, self.root / 'missing.fastq.gz'])
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse(self.output.exists())

    def test_external_checksums_match_basename_in_both_formats(self):
        raw = self.root / "reads with spaces.fastq.gz"
        raw.write_bytes(gzip.compress(b'ACGT', mtime=0))
        checksum = hashlib.md5(raw.read_bytes()).hexdigest()
        reference = self.root / "external checksums.txt"
        for entry in [f"{checksum.upper()}  elsewhere/{raw.name}\n",
                      f"{checksum} *elsewhere/{raw.name}\n",
                      f"elsewhere/{raw.name}\t{checksum}\n"]:
            with self.subTest(entry=entry):
                reference.write_text(entry)
                result = self.execute([raw], reference)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertEqual(self.output.read_text(), f"{raw}\t{checksum}\n")

    def test_missing_mismatched_or_conflicting_external_checksums_fail(self):
        raw = self.root / 'reads.fastq.gz'
        raw.write_bytes(b'ACGT')
        checksum = hashlib.md5(raw.read_bytes()).hexdigest()
        reference = self.root / 'external.txt'
        for entry in [f"{checksum}  other.fastq.gz\n",
                      f"{'0' * 32}  {raw.name}\n",
                      f"{checksum}  a/{raw.name}\n{'0' * 32}  b/{raw.name}\n"]:
            with self.subTest(entry=entry):
                reference.write_text(entry)
                result = self.execute([raw], reference)
                self.assertNotEqual(result.returncode, 0)
                self.assertFalse(self.output.exists())

    def test_inline_checksums_support_both_formats_and_literal_blocks(self):
        raw = self.root / "reads' with spaces.fastq.gz"
        raw.write_bytes(b'ACGT')
        checksum = hashlib.md5(raw.read_bytes()).hexdigest()
        entries = [f"{checksum.upper()}  elsewhere/{raw.name}",
                   f"{checksum} *elsewhere/{raw.name}",
                   f"elsewhere/{raw.name}\t{checksum}"]
        for entry in entries:
            for content in [entry, entry + '\n', '# recorded checksums\n\n' + entry + '\n']:
                with self.subTest(content=content):
                    process = self.configure([raw], content)
                    self.assertEqual(process.var['md5sums_file'], '')
                    self.assertEqual(process.var['md5sums_content'], content)
                    mounts = [value for key, value in process.var.items()
                              if key.startswith('map.input_folder_')]
                    self.assertEqual(set(mounts), {str(self.root), str(self.root.resolve())})
                    result = self.execute([raw], content)
                    self.assertEqual(result.returncode, 0, result.stderr)
                    self.assertEqual(self.output.read_text(), f"{raw}\t{checksum}\n")

    def test_inline_block_verifies_multiple_files(self):
        files = [self.root / 'first.fastq.gz', self.root / 'second.fastq.gz']
        for index, path in enumerate(files):
            path.write_bytes(f'reads {index}'.encode())
        content = ''.join(f"{hashlib.md5(path.read_bytes()).hexdigest()}  {path.name}\n"
                          for path in files)
        result = self.execute(files, content)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(len(self.output.read_text().splitlines()), 2)

    def test_inline_verification_fails_without_replacing_existing_report(self):
        raw = self.root / 'reads.fastq.gz'
        raw.write_bytes(b'ACGT')
        checksum = hashlib.md5(raw.read_bytes()).hexdigest()
        self.output.parent.mkdir()
        self.output.write_text('previous report\n')
        for content in [f"{checksum}  other.fastq.gz\n",
                        f"{'0' * 32}  {raw.name}\n",
                        f"{checksum}  a/{raw.name}\n{'0' * 32}  b/{raw.name}\n",
                        '# no entries\n', 'invalid entry\n']:
            with self.subTest(content=content):
                process = self.configure([raw], content)
                condition = next(item.value for item in PROCESS.keywords if item.arg == 'when')
                should_run = eval(compile(ast.Expression(condition), '<when>', 'eval'), {'os': os})
                self.assertTrue(should_run(process))
                result = self.execute([raw], content)
                self.assertNotEqual(result.returncode, 0)
                self.assertEqual(self.output.read_text(), 'previous report\n')


    def test_provider_download_captures_checksums_before_preprocessing(self):
        download = next(node.value for node in TREE.body if isinstance(node, ast.Assign)
                        and any(isinstance(target, ast.Name) and target.id == 'download_link'
                                for target in node.targets))
        script = ast.literal_eval(next(item.value for item in download.keywords
                                       if item.arg == 'script'))
        for valid in [True, False]:
            with self.subTest(valid=valid):
                folder = self.root / str(valid)
                nested = folder / 'delivery'
                nested.mkdir(parents=True)
                raw = nested / 'reads with spaces.fastq.gz'
                raw.write_bytes(b'provider data')
                checksum = hashlib.md5(raw.read_bytes()).hexdigest() if valid else '0' * 32
                (nested / 'provider.md5').write_text(f"{checksum}  delivery/{raw.name}\n")
                rendered = script.replace('{{source_folder}}', str(folder))
                rendered = rendered.replace('{{download_link}}', ':')
                rendered = rendered.replace('{{md5sums_file}}', 'delivery/provider.md5')
                result = subprocess.run(['bash', '-c', rendered], capture_output=True, text=True)
                self.assertEqual(result.returncode == 0, valid, result.stderr)
                self.assertEqual((folder / 'fastq.md5sums.check.log').exists(), valid)
                self.assertEqual((folder / 'fastq.md5sums.txt').read_text(),
                                 f"{checksum}  {raw.name}\n")


    def test_workflow_selects_original_names_and_individual_lanes(self):
        source = self.root / 'originals'
        source.mkdir()
        files = [source / name for name in ['lane1_R1.gz', 'lane2_R1.gz', 'lane1_R2.gz']]
        for index, path in enumerate(files):
            path.write_bytes(f'lane {index}'.encode())
        groups = 'lane1_R1.gz,lane2_R1.gz;lane1_R2.gz;control\\nlane1_R1.gz;;treated'
        workflow = ast.parse((ROOT / 'rnaseq.py').read_text())
        assignment = next(node for node in ast.walk(workflow) if isinstance(node, ast.Assign)
                          and any(isinstance(target, ast.Name) and target.id == 'raw_files'
                                  for target in node.targets))
        namespace = {'os': os, 'var': {'groups': groups},
                     '_rnaseq': SimpleNamespace(preprocess=SimpleNamespace(
                         var={'source_folder': str(source)}))}
        exec(compile(ast.Module(body=[assignment], type_ignores=[]), '<workflow>', 'exec'), namespace)
        result = self.execute(namespace['raw_files'])
        self.assertEqual(result.returncode, 0, result.stderr)
        rows = [line.split('\t') for line in self.output.read_text().splitlines()]
        self.assertEqual(rows, [[str(path), hashlib.md5(path.read_bytes()).hexdigest()]
                                for path in files])



if __name__ == '__main__':
    unittest.main()
