"""Verify visible heartbeats and clean completion without network access."""

from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
from tempfile import TemporaryDirectory
from threading import Event, enumerate as threads
import unittest

from dual_agn.progress import report_progress


class HeartbeatOutput(StringIO):
    def __init__(self):
        super().__init__()
        self.heartbeat_received = Event()

    def write(self, text):
        result = super().write(text)
        if 'WAITING ' in text:
            self.heartbeat_received.set()
        return result


class ProgressTests(unittest.TestCase):
    def test_long_stage_prints_a_heartbeat_and_received_bytes(self):
        output = HeartbeatOutput()
        with TemporaryDirectory() as directory, redirect_stdout(output):
            path = Path(directory)/'export.partial'
            path.write_bytes(b'x' * 1024**2)
            with report_progress('test export', interval_seconds=.01, output_path=path):
                self.assertTrue(output.heartbeat_received.wait(1.))
                self.assertNotIn('DONE test export', output.getvalue())
        text = output.getvalue()
        self.assertIn('START test export', text)
        self.assertIn('received 1.00 MiB', text)
        self.assertIn('database progress unknown', text)
        self.assertIn('DONE test export', text)
        self.assertFalse(any(thread.name == 'pipeline-progress' for thread in threads()))

    def test_failure_reports_type_without_exception_details(self):
        output = StringIO()
        with redirect_stdout(output), self.assertRaises(ValueError):
            with report_progress('test failure', interval_seconds=1.):
                raise ValueError('private connection detail')
        self.assertIn('FAILED (ValueError) test failure', output.getvalue())
        self.assertNotIn('private connection detail', output.getvalue())
        self.assertNotIn('DONE test failure', output.getvalue())
        self.assertFalse(any(thread.name == 'pipeline-progress' for thread in threads()))

    def test_keyboard_interrupt_reports_cancellation(self):
        output = StringIO()
        with redirect_stdout(output), self.assertRaises(KeyboardInterrupt):
            with report_progress('test cancellation', interval_seconds=1.):
                raise KeyboardInterrupt
        self.assertIn('CANCELLED test cancellation', output.getvalue())
        self.assertNotIn('DONE test cancellation', output.getvalue())

    def test_invalid_interval_fails_before_starting_a_thread(self):
        for interval in (0., -1., float('nan'), float('inf')):
            with self.subTest(interval=interval), self.assertRaisesRegex(ValueError, 'finite and positive'):
                with report_progress('invalid', interval_seconds=interval):
                    self.fail('Invalid interval was accepted')


if __name__ == '__main__':
    unittest.main()
