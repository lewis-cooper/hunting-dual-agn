"""Check that concurrent control workers use complete, separate dust-map configs."""

from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
import unittest


ROOT = Path(__file__).resolve().parents[1]


class DustmapWorkerTests(unittest.TestCase):
    def test_parallel_workers_do_not_share_their_configuration_file(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            program = '''
import json, sys
from pathlib import Path
from types import SimpleNamespace
sys.path.insert(0, str(Path(sys.argv[1])/'chance_alignments/engine'))
from score_controls import configure_production
root = Path(sys.argv[2])
configure_production(SimpleNamespace(repo=Path(sys.argv[1]), output=root/'output',
    models=root/'models/multisurvey_psf', dustmaps=root/'dustmaps'))
from dustmaps.config import config
print(json.dumps(dict(path=config.fname, data_dir=config.get('data_dir'))))
'''
            command = [sys.executable, '-B', '-c', program, str(ROOT), str(root)]
            def run_worker(_):
                return json.loads(subprocess.check_output(command, text=True, stderr=subprocess.PIPE))
            with ThreadPoolExecutor(max_workers=4) as pool:
                outputs = list(pool.map(run_worker, range(4)))
            self.assertEqual(len({row['path'] for row in outputs}), 4)
            for row in outputs:
                self.assertEqual(row['data_dir'], str((root/'dustmaps').resolve()))
                self.assertEqual(json.loads(Path(row['path']).read_text()),
                                 {'data_dir': row['data_dir']})


if __name__ == '__main__':
    unittest.main()
