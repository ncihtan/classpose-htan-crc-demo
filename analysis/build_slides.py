#!/usr/bin/env python3
"""
Build a Google Slides presentation from classpose_meeting_deck.md.
Uses the gws CLI to create the presentation and upload images to Drive.
"""

import json
import re
import subprocess
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).parent
DECK_PATH = PROJECT_ROOT / "classpose_meeting_deck.md"
FIGURES_DIR = PROJECT_ROOT / "figures"

# Slide dimensions (EMU) — 16:9 widescreen
SLIDE_W = 9144000
SLIDE_H = 5143500

# Layout areas (in EMU)
TITLE_TOP    = 274638
TITLE_LEFT   = 457200
TITLE_W      = 8229600
TITLE_H      = 900000

# Content area when there's an image (left half)
BODY_LEFT_IMG   = 457200
BODY_TOP_IMG    = 1200000
BODY_W_IMG      = 3900000
BODY_H_IMG      = 3700000

# Image area (right side)
IMG_LEFT  = 4700000
IMG_TOP   = 800000
IMG_W     = 4200000
IMG_H     = 4100000

# Content area when no image (full width)
BODY_LEFT_FULL = 457200
BODY_TOP_FULL  = 1200000
BODY_W_FULL    = 8229600
BODY_H_FULL    = 3700000


def run_gws(*args) -> dict:
    cmd = ["gws"] + list(args)
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"ERROR: {' '.join(cmd)}\n{result.stderr}", file=sys.stderr)
        sys.exit(1)
    return json.loads(result.stdout)


def parse_deck(path: Path) -> list[dict]:
    """Parse markdown into list of slide dicts with title, bullets, image."""
    text = path.read_text()
    raw_slides = re.split(r"\n---\n", text)
    slides = []

    for raw in raw_slides:
        raw = raw.strip()
        if not raw:
            continue

        slide = {"title": "", "bullets": [], "image": None, "is_title_slide": False}

        # Title slide (starts with # not ##)
        if raw.startswith("# "):
            lines = raw.splitlines()
            slide["is_title_slide"] = True
            slide["title"] = lines[0].lstrip("# ").strip()
            slide["bullets"] = [l.lstrip("#").strip() for l in lines[1:] if l.strip()]
            slides.append(slide)
            continue

        # Extract ## title
        title_match = re.search(r"^## (.+)$", raw, re.MULTILINE)
        if title_match:
            slide["title"] = title_match.group(1).strip()

        # Extract image path
        img_match = re.search(r"!\[.*?\]\((figures/[^\)]+)\)", raw)
        if img_match:
            img_path = PROJECT_ROOT / img_match.group(1)
            if img_path.exists():
                slide["image"] = img_path
            else:
                print(f"  WARNING: image not found: {img_path.name}")

        # Extract bullets and body text (lines starting with - or plain paragraphs)
        bullets = []
        in_bullets = False
        for line in raw.splitlines():
            if line.startswith("## ") or line.startswith("!["):
                continue
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith("- "):
                bullets.append(stripped[2:])
                in_bullets = True
            elif stripped.startswith("**") and stripped.endswith("**") and not in_bullets:
                # Section header within slide — treat as bold bullet
                bullets.append(stripped)
            elif stripped.startswith("|"):
                # Table row — flatten to text
                cells = [c.strip() for c in stripped.strip("|").split("|") if c.strip() and c.strip() != "---"]
                if cells:
                    bullets.append("  |  ".join(cells))
            elif in_bullets and not stripped.startswith("#"):
                # Continuation paragraph after bullets
                bullets.append(stripped)

        slide["bullets"] = bullets
        slides.append(slide)

    return slides


def upload_image(image_path: Path, folder_id: str) -> str:
    """Upload image to Drive and return a public URL."""
    print(f"  Uploading {image_path.name}...")
    result = run_gws(
        "drive", "files", "create",
        "--json", json.dumps({"name": image_path.name, "parents": [folder_id]}),
        "--upload", str(image_path),
    )
    file_id = result["id"]

    # Make publicly readable
    run_gws(
        "drive", "permissions", "create",
        "--params", json.dumps({"fileId": file_id}),
        "--json", json.dumps({"role": "reader", "type": "anyone"}),
    )

    return f"https://drive.google.com/uc?id={file_id}"


def pt(pixels: int) -> int:
    """Already in EMU, pass through."""
    return pixels


def make_text_box_request(slide_id: str, box_id: str, text: str,
                           left: int, top: int, width: int, height: int,
                           font_size: int = 18, bold: bool = False) -> list[dict]:
    requests = [
        {"createShape": {
            "objectId": box_id,
            "shapeType": "TEXT_BOX",
            "elementProperties": {
                "pageObjectId": slide_id,
                "size": {"width": {"magnitude": width, "unit": "EMU"},
                         "height": {"magnitude": height, "unit": "EMU"}},
                "transform": {"scaleX": 1, "scaleY": 1,
                               "translateX": left, "translateY": top, "unit": "EMU"},
            }
        }},
        {"insertText": {
            "objectId": box_id,
            "text": text,
        }},
        {"updateTextStyle": {
            "objectId": box_id,
            "style": {
                "fontSize": {"magnitude": font_size, "unit": "PT"},
                "bold": bold,
                "foregroundColor": {"opaqueColor": {"rgbColor": {"red": 0.13, "green": 0.13, "blue": 0.13}}},
            },
            "fields": "fontSize,bold,foregroundColor",
        }},
        {"updateParagraphStyle": {
            "objectId": box_id,
            "style": {"spaceAbove": {"magnitude": 4, "unit": "PT"},
                      "spaceBelow": {"magnitude": 4, "unit": "PT"}},
            "fields": "spaceAbove,spaceBelow",
        }},
    ]
    return requests


def strip_markdown(text: str) -> str:
    """Remove bold/italic markdown markers for Slides text."""
    text = re.sub(r"\*\*(.+?)\*\*", r"\1", text)
    text = re.sub(r"\*(.+?)\*", r"\1", text)
    text = re.sub(r"`(.+?)`", r"\1", text)
    return text


def build_slide_requests(slide: dict, slide_index: int, image_url: str | None) -> list[dict]:
    slide_id = f"slide_{slide_index}"
    requests = []

    # Create slide
    requests.append({"createSlide": {
        "objectId": slide_id,
        "insertionIndex": slide_index,
        "slideLayoutReference": {"predefinedLayout": "BLANK"},
    }})

    # Background: white
    requests.append({"updatePageProperties": {
        "objectId": slide_id,
        "pageProperties": {
            "pageBackgroundFill": {
                "solidFill": {"color": {"rgbColor": {"red": 1.0, "green": 1.0, "blue": 1.0}}}
            }
        },
        "fields": "pageBackgroundFill",
    }})

    title = strip_markdown(slide["title"])
    bullets = [strip_markdown(b) for b in slide["bullets"]]

    if slide.get("is_title_slide"):
        # Centred title slide
        requests += make_text_box_request(
            slide_id, f"{slide_id}_title", title,
            left=TITLE_LEFT, top=1800000, width=TITLE_W, height=1200000,
            font_size=36, bold=True,
        )
        if bullets:
            subtitle = "\n".join(bullets)
            requests += make_text_box_request(
                slide_id, f"{slide_id}_sub", subtitle,
                left=TITLE_LEFT, top=3200000, width=TITLE_W, height=600000,
                font_size=20, bold=False,
            )
        return requests

    # Title bar
    requests += make_text_box_request(
        slide_id, f"{slide_id}_title", title,
        left=TITLE_LEFT, top=TITLE_TOP, width=TITLE_W, height=TITLE_H,
        font_size=24, bold=True,
    )

    # Divider line under title
    requests.append({"createLine": {
        "objectId": f"{slide_id}_line",
        "lineCategory": "STRAIGHT",
        "elementProperties": {
            "pageObjectId": slide_id,
            "size": {"width": {"magnitude": TITLE_W, "unit": "EMU"},
                     "height": {"magnitude": 0, "unit": "EMU"}},
            "transform": {"scaleX": 1, "scaleY": 1,
                           "translateX": TITLE_LEFT, "translateY": 1190000, "unit": "EMU"},
        }
    }})
    requests.append({"updateLineProperties": {
        "objectId": f"{slide_id}_line",
        "lineProperties": {
            "lineFill": {"solidFill": {"color": {"rgbColor": {"red": 0.2, "green": 0.4, "blue": 0.7}}}},
            "weight": {"magnitude": 2, "unit": "PT"},
        },
        "fields": "lineFill,weight",
    }})

    # Body text
    if bullets:
        bullet_text = "\n".join(f"• {b}" if not b.startswith("•") else b for b in bullets)
        if image_url:
            requests += make_text_box_request(
                slide_id, f"{slide_id}_body", bullet_text,
                left=BODY_LEFT_IMG, top=BODY_TOP_IMG,
                width=BODY_W_IMG, height=BODY_H_IMG,
                font_size=14,
            )
        else:
            requests += make_text_box_request(
                slide_id, f"{slide_id}_body", bullet_text,
                left=BODY_LEFT_FULL, top=BODY_TOP_FULL,
                width=BODY_W_FULL, height=BODY_H_FULL,
                font_size=15,
            )

    # Image
    if image_url:
        requests.append({"createImage": {
            "objectId": f"{slide_id}_img",
            "url": image_url,
            "elementProperties": {
                "pageObjectId": slide_id,
                "size": {"width": {"magnitude": IMG_W, "unit": "EMU"},
                         "height": {"magnitude": IMG_H, "unit": "EMU"}},
                "transform": {"scaleX": 1, "scaleY": 1,
                               "translateX": IMG_LEFT, "translateY": IMG_TOP, "unit": "EMU"},
            }
        }})

    return requests


def main():
    print("Parsing deck...")
    slides = parse_deck(DECK_PATH)
    print(f"  Found {len(slides)} slides")

    # Create Drive folder for images
    print("\nCreating Drive folder for images...")
    folder = run_gws(
        "drive", "files", "create",
        "--json", json.dumps({
            "name": "ClassPose Meeting Deck Images",
            "mimeType": "application/vnd.google-apps.folder",
        }),
    )
    folder_id = folder["id"]
    print(f"  Folder ID: {folder_id}")

    # Upload images
    print("\nUploading images...")
    image_urls: dict[str, str] = {}
    for slide in slides:
        if slide["image"] and slide["image"] not in image_urls:
            url = upload_image(slide["image"], folder_id)
            image_urls[str(slide["image"])] = url

    # Create presentation
    print("\nCreating presentation...")
    pres = run_gws(
        "slides", "presentations", "create",
        "--json", json.dumps({"title": "ClassPose × HTAN CRC: Work to Date"}),
    )
    pres_id = pres["presentationId"]
    print(f"  Presentation ID: {pres_id}")

    # Delete the default blank slide
    default_slide_id = pres["slides"][0]["objectId"]

    # Build all requests
    print("\nBuilding slide requests...")
    all_requests = [{"deleteObject": {"objectId": default_slide_id}}]

    for i, slide in enumerate(slides):
        image_url = image_urls.get(str(slide["image"])) if slide["image"] else None
        reqs = build_slide_requests(slide, i, image_url)
        all_requests.extend(reqs)
        print(f"  Slide {i+1}/{len(slides)}: {slide['title'][:50]}")

    # Send batchUpdate in chunks (API limit ~100 requests per call)
    print(f"\nSending {len(all_requests)} requests to Slides API...")
    chunk_size = 50
    for i in range(0, len(all_requests), chunk_size):
        chunk = all_requests[i:i + chunk_size]
        run_gws(
            "slides", "presentations", "batchUpdate",
            "--params", json.dumps({"presentationId": pres_id}),
            "--json", json.dumps({"requests": chunk}),
        )
        print(f"  Sent requests {i+1}–{min(i+chunk_size, len(all_requests))}")

    url = f"https://docs.google.com/presentation/d/{pres_id}/edit"
    print(f"\nDone! Presentation: {url}")
    return url


if __name__ == "__main__":
    main()
