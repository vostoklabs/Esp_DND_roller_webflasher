# ESP32 Webflasher Boilerplate

A ready-to-use, copy/paste boilerplate to publish an ESP32 browser flasher (powered by ESP Web Tools) similar to the existing Bongo Cat installer.

## Quick Start

1. Copy the entire `webflasher-boilerplate/` folder into a new repository (or make it your repo root).
2. Put your firmware binaries into the `build/` subfolders (see below).
3. Update the manifest JSONs to point to your actual file names/paths and offsets.
4. Open `index.html` and replace the title, headings, and any images in `web resources/`.
5. Test locally via a web server (HTTPS recommended) and then deploy to GitHub Pages.

## File Structure

```
webflasher-boilerplate/
├── index.html                      # Web flasher UI with multiple variants
├── manifest.json                   # Variant A manifest (default)
├── manifest-variant-b.json         # Variant B manifest
├── manifest-variant-a-inverse.json # Variant A inverse colors fix (optional)
├── manifest-variant-b-inverse.json # Variant B inverse colors fix (optional)
├── build/
│   ├── Variant A/
│   │   └── README.md               # Place bootloader.bin, partitions.bin, boot_app0.bin (optional), firmware.bin
│   ├── Variant A inverse/
│   │   └── README.md
│   ├── Variant B/
│   │   └── README.md
│   └── Variant B inverse/
│       └── README.md
└── web resources/
    └── README.md                   # Put images/screenshots used in index.html
```

## Binaries Needed

Typical ESP32 builds produce these files (names may vary by toolchain):

- bootloader.bin → offset 0x1000 (4096)
- partitions.bin → offset 0x8000 (32768)
- boot_app0.bin → offset 0xE000 (57344) (optional, depends on build)
- firmware.bin (app) → offset 0x10000 (65536)

Update the manifest `parts.path` to match your filenames, and keep the offsets correct for your build.

## Local Testing

- Python: `python -m http.server 8080` then open `http://localhost:8080/webflasher-boilerplate/`
- Node: `npx http-server -p 8080` then open `http://localhost:8080/webflasher-boilerplate/`

Chrome/Edge required (Web Serial). Some features need HTTPS; localhost usually works for testing.

## Deploy to GitHub Pages

- If this folder is your repo root, enable Pages on the `main` branch `/ (root)`.
- If this folder is within a larger repo, copy its contents into a new repo or configure Pages to serve from the folder.

## Customize

- Manifests: `manifest*.json` → set `name`, `chipFamily`, `parts` paths and offsets.
- Text/UI: edit `index.html` and replace the hero image in `web resources/`.
- Variants: rename "Variant A/B" (and inverse) to your board/display names, or remove the inverse section.

For deeper details, see `INSTRUCTIONS.md`.
