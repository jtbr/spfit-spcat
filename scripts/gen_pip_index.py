#!/usr/bin/env python3
"""Generate a PEP 503 pip index from GitHub release wheel assets.

Writes pip/index.html and pip/<package>/index.html, collecting .whl assets
from all releases so historical versions remain installable.
Deploy pip/ to the gh-pages branch; pip users then install with:
  pip install pickett --extra-index-url https://jtbr.github.io/spfit-spcat/pip/
"""

import json
import os
import sys
from pathlib import Path
from urllib.request import Request, urlopen

REPO = os.environ.get("GITHUB_REPOSITORY", "jtbr/spfit-spcat")
PACKAGE = "pickett"


def gh_get(path):
    url = f"https://api.github.com/{path}"
    headers = {
        "Accept": "application/vnd.github+json",
        "X-GitHub-Api-Version": "2022-11-28",
    }
    token = os.environ.get("GH_TOKEN")
    if token:
        headers["Authorization"] = f"Bearer {token}"
    with urlopen(Request(url, headers=headers)) as resp:
        return json.loads(resp.read())


def fetch_wheels():
    wheels = []
    page = 1
    while True:
        releases = gh_get(f"repos/{REPO}/releases?per_page=100&page={page}")
        if not releases:
            break
        for release in releases:
            for asset in release.get("assets", []):
                if asset["name"].endswith(".whl"):
                    wheels.append((asset["name"], asset["browser_download_url"]))
        page += 1
    return wheels


def write_index(out_dir: Path, wheels: list):
    out_dir.mkdir(parents=True, exist_ok=True)
    pkg_dir = out_dir / PACKAGE
    pkg_dir.mkdir(exist_ok=True)

    (out_dir / "index.html").write_text(
        "<!DOCTYPE html><html><head><title>Simple Index</title></head><body>\n"
        f'<a href="{PACKAGE}/">{PACKAGE}</a>\n'
        "</body></html>\n"
    )

    links = "\n".join(f'<a href="{url}">{name}</a>' for name, url in sorted(wheels))
    (pkg_dir / "index.html").write_text(
        f"<!DOCTYPE html><html><head><title>Links for {PACKAGE}</title></head><body>\n"
        f"<h1>Links for {PACKAGE}</h1>\n"
        f"{links}\n"
        "</body></html>\n"
    )

    print(f"Wrote index with {len(wheels)} wheels across all releases")


def main():
    wheels = fetch_wheels()
    if not wheels:
        print("No wheels found in any release.", file=sys.stderr)
        sys.exit(1)
    write_index(Path("pip"), wheels)


if __name__ == "__main__":
    main()
