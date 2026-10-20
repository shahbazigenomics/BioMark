# BioMark - System Hardware Detection
# Author: Amir Shahbazi
# GitHub: shahbazigenomics

import platform
import subprocess
import psutil
import os
from datetime import datetime

# ─── Get CPU Info ──────────────────────────────────────────
def get_cpu_info():
    cpu = {}
    cpu["architecture"] = platform.machine()
    cpu["python_processor"] = platform.processor()
    cpu["physical_cores"] = psutil.cpu_count(logical=False)
    cpu["logical_cores"] = psutil.cpu_count(logical=True)

    # Apple Silicon detection
    if platform.system() == "Darwin":
        try:
            result = subprocess.run(
                ["sysctl", "-n", "machdep.cpu.brand_string"],
                capture_output=True, text=True
            )
            cpu["brand"] = result.stdout.strip()

            # Detect Apple chip model (M1/M2/M3/M4/M5)
            result2 = subprocess.run(
                ["sysctl", "-n", "hw.model"],
                capture_output=True, text=True
            )
            cpu["mac_model"] = result2.stdout.strip()

            # Get chip name via system_profiler
            result3 = subprocess.run(
                ["system_profiler", "SPHardwareDataType"],
                capture_output=True, text=True
            )
            for line in result3.stdout.split("\n"):
                if "Chip" in line or "Processor" in line:
                    cpu["apple_chip"] = line.strip()
                    break

        except Exception as e:
            cpu["brand"] = "Unknown"
            cpu["error"] = str(e)

    return cpu


# ─── Get RAM Info ──────────────────────────────────────────
def get_ram_info():
    ram = {}
    mem = psutil.virtual_memory()

    ram["total_gb"] = round(mem.total / (1024 ** 3), 1)
    ram["available_gb"] = round(mem.available / (1024 ** 3), 1)
    ram["used_gb"] = round(mem.used / (1024 ** 3), 1)
    ram["percent_used"] = mem.percent

    # Apple unified memory bandwidth
    if platform.system() == "Darwin":
        try:
            result = subprocess.run(
                ["sysctl", "-n", "hw.memsize"],
                capture_output=True, text=True
            )
            ram["hw_memsize_gb"] = round(
                int(result.stdout.strip()) / (1024 ** 3), 1
            )
        except:
            pass

    return ram


# ─── Get SSD Info ──────────────────────────────────────────
def get_storage_info():
    storage = {}
    disk = psutil.disk_usage("/")

    storage["total_gb"] = round(disk.total / (1024 ** 3), 1)
    storage["used_gb"] = round(disk.used / (1024 ** 3), 1)
    storage["free_gb"] = round(disk.free / (1024 ** 3), 1)
    storage["percent_used"] = disk.percent

    return storage


# ─── Get GPU Info ──────────────────────────────────────────
def get_gpu_info():
    gpu = {}

    if platform.system() == "Darwin":
        try:
            result = subprocess.run(
                ["system_profiler", "SPDisplaysDataType"],
                capture_output=True, text=True
            )
            lines = result.stdout.split("\n")
            for i, line in enumerate(lines):
                if "Chipset Model" in line:
                    gpu["model"] = line.split(":")[1].strip()
                if "VRAM" in line:
                    gpu["vram"] = line.split(":")[1].strip()
                if "Metal" in line:
                    gpu["metal"] = line.split(":")[1].strip()
        except Exception as e:
            gpu["error"] = str(e)

    return gpu


# ─── Get OS Info ───────────────────────────────────────────
def get_os_info():
    os_info = {}
    os_info["system"] = platform.system()
    os_info["release"] = platform.release()
    os_info["version"] = platform.version()
    os_info["python_version"] = platform.python_version()

    if platform.system() == "Darwin":
        os_info["macos_version"] = platform.mac_ver()[0]

    return os_info


# ─── Get Full System Profile ───────────────────────────────
def get_full_system_profile():
    print("\n🖥️  Detecting system hardware...")

    profile = {
        "timestamp": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "cpu": get_cpu_info(),
        "ram": get_ram_info(),
        "storage": get_storage_info(),
        "gpu": get_gpu_info(),
        "os": get_os_info(),
    }

    return profile


# ─── Print System Profile ──────────────────────────────────
def print_system_profile(profile):
    print("\n" + "=" * 55)
    print("🖥️   BioMark System Profile")
    print("=" * 55)

    # OS
    os_i = profile["os"]
    print(f"\n  💻 Operating System")
    print(f"     System    : {os_i.get('system', 'N/A')}")
    print(f"     macOS     : {os_i.get('macos_version', 'N/A')}")
    print(f"     Python    : {os_i.get('python_version', 'N/A')}")

    # CPU
    cpu = profile["cpu"]
    print(f"\n  ⚡ Processor")
    print(f"     Chip      : {cpu.get('apple_chip', cpu.get('brand', 'N/A'))}")
    print(f"     Model     : {cpu.get('mac_model', 'N/A')}")
    print(f"     Arch      : {cpu.get('architecture', 'N/A')}")
    print(f"     P-Cores   : {cpu.get('physical_cores', 'N/A')}")
    print(f"     L-Cores   : {cpu.get('logical_cores', 'N/A')}")

    # RAM
    ram = profile["ram"]
    print(f"\n  🧠 Memory (RAM)")
    print(f"     Total     : {ram.get('total_gb', 'N/A')} GB")
    print(f"     Available : {ram.get('available_gb', 'N/A')} GB")
    print(f"     Used      : {ram.get('used_gb', 'N/A')} GB")
    print(f"     Usage     : {ram.get('percent_used', 'N/A')}%")

    # Storage
    st = profile["storage"]
    print(f"\n  💾 Storage (SSD)")
    print(f"     Total     : {st.get('total_gb', 'N/A')} GB")
    print(f"     Used      : {st.get('used_gb', 'N/A')} GB")
    print(f"     Free      : {st.get('free_gb', 'N/A')} GB")
    print(f"     Usage     : {st.get('percent_used', 'N/A')}%")

    # GPU
    gpu = profile["gpu"]
    print(f"\n  🎮 GPU")
    print(f"     Model     : {gpu.get('model', 'N/A')}")
    print(f"     VRAM      : {gpu.get('vram', 'Unified (shared RAM)')}")
    print(f"     Metal     : {gpu.get('metal', 'N/A')}")

    print("\n" + "=" * 55)
