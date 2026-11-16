#!/usr/bin/env python3
"""Create a release notice containing Cargo dependency license texts."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LICENSE_PREFIXES = ("license", "copying", "unlicense", "notice")


def license_files(package: dict) -> list[Path]:
    if package.get("license_file"):
        candidates = [Path(package["license_file"])]
    else:
        package_dir = Path(package["manifest_path"]).parent
        candidates = [
            path
            for path in package_dir.iterdir()
            if path.is_file() and path.name.lower().startswith(LICENSE_PREFIXES)
        ]
    return sorted(
        (path for path in candidates if path.is_file()),
        key=lambda path: path.name.lower(),
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    command = [
        "cargo",
        "metadata",
        "--locked",
        "--format-version",
        "1",
        "--manifest-path",
        str(ROOT / "rust/Cargo.toml"),
    ]
    metadata = json.loads(
        subprocess.check_output(command, cwd=ROOT, text=True)
    )
    packages = sorted(
        (package for package in metadata["packages"] if package.get("source")),
        key=lambda package: (package["name"].lower(), package["version"]),
    )
    sections: list[str] = [
        "Third-party dependency license notices",
        "========================================",
        "",
        "Generated from rust/Cargo.lock and the license files shipped with each",
        "resolved Cargo package. The DockQ project license is in the separate",
        "LICENSE file.",
        "",
    ]

    missing: list[str] = []
    for package in packages:
        files = license_files(package)
        if not files:
            missing.append(f"{package['name']} {package['version']}")
            continue
        sections.extend(
            [
                f"{package['name']} {package['version']}",
                f"Declared license: {package.get('license') or 'unspecified'}",
                f"Package source: {package.get('repository') or package.get('source')}",
                "",
            ]
        )
        for path in files:
            sections.extend(
                [
                    f"----- {package['name']} / {path.name} -----",
                    path.read_text(encoding="utf-8", errors="replace").rstrip(),
                    "",
                ]
            )

    if missing:
        print(
            "Missing license files for: " + ", ".join(missing),
            file=sys.stderr,
        )
        return 1

    output = args.output if args.output.is_absolute() else ROOT / args.output
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(sections) + "\n", encoding="utf-8")
    print(f"Wrote licenses for {len(packages)} Cargo packages to {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
