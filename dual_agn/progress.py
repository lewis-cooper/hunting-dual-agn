"""Small terminal progress messages; a heartbeat is not a completion estimate."""

from contextlib import contextmanager
from datetime import datetime
import math
import os
from pathlib import Path
from threading import Event, Thread
from time import monotonic


def print_status(message):
    print(f'[{datetime.now().astimezone():%H:%M:%S %Z}] {message}', flush=True)


def command_label(command):
    if '-m' in command:
        return command[command.index('-m') + 1].removeprefix('scripts.').replace('_', ' ')
    script = next((Path(value).stem for value in command if value.endswith('.py')), 'command')
    return script.replace('_', ' ')


@contextmanager
def report_progress(label, interval_seconds=None, output_path=None):
    interval = float(interval_seconds if interval_seconds is not None
                     else os.environ.get('PAPER1_PROGRESS_INTERVAL_SECONDS', '30'))
    if not math.isfinite(interval) or interval <= 0:
        raise ValueError('Progress interval must be finite and positive')
    started = monotonic()
    stopped = Event()

    def elapsed():
        seconds = int(monotonic() - started)
        minutes, seconds = divmod(seconds, 60)
        hours, minutes = divmod(minutes, 60)
        return f'{hours:02d}:{minutes:02d}:{seconds:02d}'

    def heartbeat():
        while not stopped.wait(interval):
            detail = ''
            if output_path is not None:
                try:
                    size = Path(output_path).stat().st_size
                except OSError:
                    size = 0
                detail = f'; received {size / 1024**2:.2f} MiB; database progress unknown'
            print_status(f'WAITING {label}; elapsed {elapsed()}{detail} (heartbeat)')

    print_status(f'START {label}')
    thread = Thread(target=heartbeat, name='pipeline-progress', daemon=True)
    thread.start()
    outcome = 'DONE'
    try:
        yield
    except BaseException as error:
        outcome = 'CANCELLED' if isinstance(error, KeyboardInterrupt) else f'FAILED ({type(error).__name__})'
        raise
    finally:
        stopped.set()
        thread.join()
        print_status(f'{outcome} {label}; elapsed {elapsed()}')
