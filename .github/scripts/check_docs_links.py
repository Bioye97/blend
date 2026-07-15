#!/usr/bin/env python3
"""Check generated Sphinx HTML links without requiring network access."""

from pathlib import Path
from urllib.parse import unquote, urlsplit
import re
import sys


LINK_RE = re.compile(r"\b(?:href|src)=[\"']([^\"']+)[\"']")
ID_RE = re.compile(r"\b(?:id|name)=[\"']([^\"']+)[\"']")
GITHUB_BLOB_RE = re.compile(r"^/([^/]+)/([^/]+)/blob/([^/]+)/(.+)$")


def fail(message):
    print(message, file=sys.stderr)
    return 1


def main():
    if len(sys.argv) != 2:
        return fail("usage: check_docs_links.py HTML_DIR")

    html_root = Path(sys.argv[1]).resolve()
    repo_root = Path.cwd().resolve()
    if not html_root.is_dir():
        return fail(f"{html_root} is not a directory")

    html_files = list(html_root.rglob("*.html"))
    ids_by_page = {}
    links_by_page = {}

    for path in html_files:
        text = path.read_text(encoding="utf-8", errors="replace")
        resolved = path.resolve()
        ids_by_page[resolved] = set(ID_RE.findall(text))
        links_by_page[resolved] = LINK_RE.findall(text)

    errors = []
    for path, links in links_by_page.items():
        for raw in links:
            if not raw:
                continue

            split = urlsplit(raw)
            if split.scheme in ("http", "https"):
                if split.netloc == "github.com":
                    match = GITHUB_BLOB_RE.match(split.path)
                    if match:
                        target = repo_root / unquote(match.group(4).split("#", 1)[0])
                        if not target.exists():
                            errors.append((path, raw, "missing GitHub source target"))
                continue
            if raw.startswith(("mailto:", "javascript:", "data:")) or split.scheme or split.netloc:
                continue

            link_path = split.path
            fragment = unquote(split.fragment)
            if not link_path:
                target = path
            else:
                target = (path.parent / unquote(link_path)).resolve()
                try:
                    target.relative_to(html_root)
                except ValueError:
                    errors.append((path, raw, "points outside documentation tree"))
                    continue
                if target.is_dir():
                    target = target / "index.html"

            if not target.exists():
                errors.append((path, raw, "missing target"))
                continue
            if target.suffix == ".html" and fragment and fragment not in ids_by_page.get(target, set()):
                errors.append((path, raw, f"missing anchor #{fragment}"))

    if errors:
        for source, link, message in errors[:200]:
            print(f"{source.relative_to(html_root)}: {link}: {message}", file=sys.stderr)
        if len(errors) > 200:
            print(f"... {len(errors) - 200} more", file=sys.stderr)
        return fail(f"{len(errors)} broken documentation links")

    print(f"checked {len(html_files)} HTML pages; documentation links ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
