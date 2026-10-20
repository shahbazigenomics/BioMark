"""Tests for hardware detection in biomark.system_info (no real Mac required)."""
from types import SimpleNamespace
from unittest.mock import patch

from biomark import system_info

SYSTEM_PROFILER_OUTPUT = (
    "Hardware:\n\n    Hardware Overview:\n\n"
    "      Model Name: MacBook Pro\n"
    "      Model Identifier: Mac17,9\n"
    "      Chip: Apple M5 Pro\n"
)


def _result(stdout="", returncode=0):
    return SimpleNamespace(stdout=stdout, returncode=returncode)


def _fake_run(perflevels):
    """Build a fake subprocess.run. perflevels = {0: ("Super", "5"), ...} (test fixture values)."""
    def fake(cmd, **kwargs):
        if cmd[0] == "system_profiler":
            return _result(SYSTEM_PROFILER_OUTPUT)
        key = cmd[-1]
        if key == "machdep.cpu.brand_string":
            return _result("Apple M5 Pro")
        if key == "hw.model":
            return _result("Mac17,9")
        for level, (name, count) in perflevels.items():
            if key == f"hw.perflevel{level}.name":
                return _result(name)
            if key == f"hw.perflevel{level}.physicalcpu":
                return _result(count)
        return _result("", returncode=1)
    return fake


def _darwin_cpu(perflevels):
    with patch.object(system_info.platform, "system", return_value="Darwin"), \
         patch.object(system_info.subprocess, "run", side_effect=_fake_run(perflevels)):
        return system_info.get_cpu_info()


def test_chip_name_has_no_label_prefix():
    cpu = _darwin_cpu({0: ("Super", "5"), 1: ("Performance", "10")})
    assert cpu["apple_chip"] == "Apple M5 Pro"


def test_core_mix_read_from_perflevels():
    cpu = _darwin_cpu({0: ("Super", "5"), 1: ("Performance", "10")})
    assert cpu["core_mix"] == [
        {"name": "Super", "cores": 5},
        {"name": "Performance", "cores": 10},
    ]


def test_core_mix_empty_when_sysctl_has_no_perflevels():
    cpu = _darwin_cpu({})
    assert cpu["core_mix"] == []


def test_non_mac_has_core_counts_but_no_core_mix():
    with patch.object(system_info.platform, "system", return_value="Linux"):
        cpu = system_info.get_cpu_info()
    assert "core_mix" not in cpu
    assert cpu["physical_cores"] is None or cpu["physical_cores"] >= 1


def test_linux_pretty_name_reads_os_release(tmp_path):
    f = tmp_path / "os-release"
    f.write_text('NAME="Ubuntu"\nPRETTY_NAME="Ubuntu 24.04 LTS"\n')
    assert system_info._linux_pretty_name(str(f)) == "Ubuntu 24.04 LTS"


def test_linux_pretty_name_missing_file_returns_none(tmp_path):
    assert system_info._linux_pretty_name(str(tmp_path / "nope")) is None


def test_os_version_on_linux_uses_distro_name():
    with patch.object(system_info.platform, "system", return_value="Linux"), \
         patch.object(system_info, "_linux_pretty_name", return_value="Ubuntu 24.04 LTS"):
        info = system_info.get_os_info()
    assert info["os_version"] == "Ubuntu 24.04 LTS"
    assert "macos_version" not in info


def test_os_version_on_linux_falls_back_to_kernel_release():
    with patch.object(system_info.platform, "system", return_value="Linux"), \
         patch.object(system_info, "_linux_pretty_name", return_value=None), \
         patch.object(system_info.platform, "release", return_value="6.8.0"):
        info = system_info.get_os_info()
    assert info["os_version"] == "6.8.0"


def test_os_version_on_mac_matches_macos_version():
    with patch.object(system_info.platform, "system", return_value="Darwin"), \
         patch.object(system_info.platform, "mac_ver", return_value=("27.0.1", ("", "", ""), "arm64")):
        info = system_info.get_os_info()
    assert info["os_version"] == info["macos_version"] == "27.0.1"
