#!/usr/bin/env python3
"""
nb2md.py - render an executed Jupyter notebook as an MkDocs markdown page.

The published page is generated from the notebook itself, so the site and
the notebook can never drift apart: markdown cells are copied through
verbatim (keeping the numbered equations), code cells become fenced
`python` blocks, and the stored outputs - text and figures - are emitted
underneath each one.

Why a converter instead of the mkdocs-jupyter plugin: the output is plain
markdown, so the site keeps working even on a MkDocs build without plugin
support, and it needs no dependency beyond the Python standard library.
See README.md for the mkdocs-jupyter alternative.

Usage:
    python3 tools/nb2md.py EQE/eqe_analysis.ipynb docs/eqe.md \\
        --assets docs/assets/nb --assets-url ../assets/nb \\
        --rewrite ./figures/=../assets/ \\
        --repo-url https://github.com/USER/REPO

Cells tagged `hide-in-docs` in their notebook metadata are skipped.
"""
import argparse
import base64
import json
import os
import re
import sys

COLAB_BASE = "https://colab.research.google.com/github"


def cell_source(cell):
    src = cell.get("source", "")
    return src if isinstance(src, str) else "".join(src)


def out_text(output):
    """Plain text carried by a stream / text-bearing output, or ''."""
    if output.get("output_type") == "stream":
        t = output.get("text", "")
        return t if isinstance(t, str) else "".join(t)
    data = output.get("data", {})
    t = data.get("text/plain", "")
    if isinstance(t, list):
        t = "".join(t)
    # A bare "<Figure>" repr adds nothing next to the image itself.
    return "" if t.strip() in ("<Figure>", "") else t


def out_png(output):
    """Base64 PNG carried by an output, or None."""
    data = output.get("data", {})
    png = data.get("image/png")
    if png is None:
        return None
    return png if isinstance(png, str) else "".join(png)


def apply_rewrites(text, rewrites):
    for old, new in rewrites:
        text = text.replace(old, new)
    return text


def convert(nb_path, md_path, assets_dir, assets_url, rewrites, repo_url,
            branch="main"):
    with open(nb_path) as f:
        nb = json.load(f)

    os.makedirs(assets_dir, exist_ok=True)
    os.makedirs(os.path.dirname(md_path) or ".", exist_ok=True)

    stem = os.path.splitext(os.path.basename(nb_path))[0]
    # Clear this notebook's previously generated images. Stale files are
    # only a tidiness problem - every image still in use is rewritten
    # below - so a filesystem that refuses deletes must not fail the build.
    for old in os.listdir(assets_dir):
        if old.startswith(stem + "_") and old.endswith(".png"):
            try:
                os.remove(os.path.join(assets_dir, old))
            except OSError as exc:
                print(f"  note: could not remove stale {old} ({exc.strerror})")

    parts = [
        f"<!-- GENERATED FILE - do not edit. Produced from {nb_path} by "
        "tools/nb2md.py (see tools/build_docs.sh). -->"
    ]

    # Header: where this page comes from and how to run it.
    if repo_url:
        repo_path = repo_url.rstrip("/").split("github.com/")[-1]
        nb_url = f"{repo_url.rstrip('/')}/blob/{branch}/{nb_path}"
        colab_url = f"{COLAB_BASE}/{repo_path}/blob/{branch}/{nb_path}"
        parts.append(
            '!!! info "Generated from a Jupyter notebook"\n'
            f"    This page is `{nb_path}`, rendered with its stored outputs.\n"
            f"    [Run it in Google Colab]({colab_url}) or\n"
            f"    [view the notebook on GitHub]({nb_url}).\n"
        )

    n_img = 0
    for cell in nb.get("cells", []):
        tags = (cell.get("metadata") or {}).get("tags") or []
        if "hide-in-docs" in tags:
            continue

        src = cell_source(cell).rstrip("\n")

        if cell.get("cell_type") == "markdown":
            if src.strip():
                parts.append(apply_rewrites(src, rewrites))
            continue

        if cell.get("cell_type") != "code":
            continue

        if src.strip():
            parts.append("```python\n" + src + "\n```")

        for output in cell.get("outputs", []):
            png = out_png(output)
            if png is not None:
                n_img += 1
                name = f"{stem}_{n_img:02d}.png"
                with open(os.path.join(assets_dir, name), "wb") as fh:
                    fh.write(base64.b64decode(png))
                parts.append(f"![Output {n_img}]({assets_url}/{name})")
                continue

            text = out_text(output).rstrip("\n")
            if text:
                parts.append("```text\n" + text + "\n```")

            if output.get("output_type") == "error":
                tb = "\n".join(output.get("traceback", []))
                tb = re.sub(r"\x1b\[[0-9;]*m", "", tb)  # strip ANSI colours
                parts.append("```text\n" + tb.rstrip("\n") + "\n```")

    body = "\n\n".join(parts).rstrip("\n") + "\n"

    with open(md_path, "w") as f:
        f.write(body)

    print(f"{nb_path} -> {md_path}  ({n_img} figures -> {assets_dir}/)")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("notebook")
    ap.add_argument("output")
    ap.add_argument("--assets", default="docs/assets/nb",
                    help="directory to write extracted output images into")
    ap.add_argument("--assets-url", default="../assets/nb",
                    help="URL prefix for those images, as seen from the page")
    ap.add_argument("--rewrite", action="append", default=[], metavar="OLD=NEW",
                    help="literal path rewrite applied to markdown cells "
                         "(repeatable), e.g. ./figures/=../assets/")
    ap.add_argument("--repo-url", default="",
                    help="GitHub repo URL, for the Colab / source links")
    ap.add_argument("--branch", default="main")
    args = ap.parse_args(argv)

    rewrites = []
    for r in args.rewrite:
        if "=" not in r:
            ap.error(f"--rewrite expects OLD=NEW, got {r!r}")
        old, new = r.split("=", 1)
        rewrites.append((old, new))

    convert(args.notebook, args.output, args.assets, args.assets_url,
            rewrites, args.repo_url, args.branch)


if __name__ == "__main__":
    sys.exit(main())
