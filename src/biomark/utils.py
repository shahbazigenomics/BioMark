# BioMark - Utility Functions
# Author: Amir Shahbazi
# GitHub: shahbazigenomics/BioMark

import time
import os
import threading
import json
import numpy as np
import psutil


class BioMarkEncoder(json.JSONEncoder):
    """Custom JSON encoder to handle numpy types"""
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.bool_):
            return bool(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super().default(obj)


def measure_cpu_time(fn, *args, **kwargs):
    """Execute a callable and record wall-clock time plus peak RAM"""
    peak_ram   = [0.0]
    stop_event = threading.Event()

    def _ram_monitor():
        proc = psutil.Process(os.getpid())
        while not stop_event.is_set():
            try:
                mem_gb = proc.memory_info().rss / (1024**3)
                if mem_gb > peak_ram[0]:
                    peak_ram[0] = mem_gb
            except psutil.NoSuchProcess:
                break
            time.sleep(0.05)

    monitor = threading.Thread(target=_ram_monitor, daemon=True)
    monitor.start()

    t0           = time.perf_counter()
    error        = None
    return_value = None
    success      = False

    try:
        return_value = fn(*args, **kwargs)
        success      = True
    except Exception as exc:
        error = str(exc)
    finally:
        stop_event.set()
        monitor.join(timeout=1)

    wall_time = round(time.perf_counter() - t0, 4)
    return {
        "success"          : success,
        "wall_time_seconds": wall_time,
        "peak_ram_gb"      : round(peak_ram[0], 3),
        "return_value"     : return_value,
        "error"            : error,
    }


def measure_io_speed(path=".", size_mb=64):
    """Measure sequential write and read speed"""
    tmp  = os.path.join(path, "_biomark_io_test.tmp")
    data = os.urandom(size_mb * 1024 * 1024)

    write_speed = 0.0
    read_speed  = 0.0

    try:
        t0 = time.perf_counter()
        with open(tmp, "wb") as fh:
            fh.write(data)
        write_speed = round(size_mb / max(time.perf_counter() - t0, 1e-6), 1)

        t0 = time.perf_counter()
        with open(tmp, "rb") as fh:
            _ = fh.read()
        read_speed = round(size_mb / max(time.perf_counter() - t0, 1e-6), 1)

    except OSError:
        pass
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)

    return {
        "write_speed_mbs": write_speed,
        "read_speed_mbs" : read_speed,
    }


def calculate_score(wall_time, peak_ram_gb, read_speed, write_speed,
                    baseline_time=10.0, baseline_ram=8.0,
                    baseline_read=500.0, baseline_write=300.0):
    """Compute a 0-100 composite performance score"""
    cpu_score   = min(1.0, baseline_time / max(wall_time,   0.001)) * 100
    ram_score   = min(1.0, baseline_ram  / max(peak_ram_gb, 0.001)) * 100
    read_score  = min(1.0, read_speed    / baseline_read)           * 100
    write_score = min(1.0, write_speed   / baseline_write)          * 100

    composite = (
        0.40 * cpu_score  +
        0.30 * ram_score  +
        0.20 * read_score +
        0.10 * write_score
    )
    return int(round(composite))
