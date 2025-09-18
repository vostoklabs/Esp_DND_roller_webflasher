# Webflasher Boilerplate – Implementation Guide

This guide explains how to adapt the boilerplate to any ESP32 project and publish a browser-based flasher.

## 1) Collect Firmware Artifacts

Place the following binaries into the appropriate `build/<Variant...>/` folder(s):

- bootloader.bin → typical offset 0x1000 (4096)
- partitions.bin → typical offset 0x8000 (32768)
- boot_app0.bin → typical offset 0xE000 (57344) (optional; depends on build)
- firmware.bin (application) → typical offset 0x10000 (65536)

Notes:
- Exact offsets can vary by project. Use your build system’s output or `esptool.py` metadata.
- Arduino/PlatformIO/ESP-IDF usually produce these with the shown offsets for ESP32.

## 2) Update Manifest(s)

Edit the `manifest*.json` files at the project root of this boilerplate:

- `name`: A human-friendly label (appears in the installer dialog)
- `chipFamily`: One of `ESP32`, `ESP32S2`, `ESP32S3`, `ESP32C3` (match your target)
- `parts`: File `path` must match the location you placed binaries in `build/` and the `offset` must match your firmware layout

Example (Variant A):

```json
{
  "name": "Your Project – Variant A",
  "version": "0.0.1",
  "new_install_prompt_erase": true,
  "builds": [
    {
      "chipFamily": "ESP32",
      "parts": [
        { "path": "build/Variant A/bootloader.bin", "offset": 4096 },
        { "path": "build/Variant A/partitions.bin", "offset": 32768 },
        { "path": "build/Variant A/boot_app0.bin", "offset": 57344 },
        { "path": "build/Variant A/firmware.bin", "offset": 65536 }
      ]
    }
  ]
}
```

## 3) Customize the UI

- Open `index.html` and replace the title, headings, and texts.
- Replace `web resources/hero.png` and add any other images or GIFs.
- Rename "Variant A/B" to the actual board/display names. Remove the inverse section if not needed.

## 4) USB Device Filtering (optional)

The page pre-filters for common ESP32 USB Vendor IDs (CP210x, CH34x, FTDI, Espressif, Arduino, Adafruit, SparkFun). If your device uses different IDs, add them in `index.html` under `esp32VendorIds`.

## 5) Test Locally

Use a local server (Chrome/Edge):
- Python: `python -m http.server 8080`
- Node: `npx http-server -p 8080`

Then open `http://localhost:8080/webflasher-boilerplate/`.

If you see "HTTPS required" warnings, that refers to production hosting. `localhost` typically works for testing.

## 6) Deploy

- GitHub Pages: push the boilerplate contents to a repo, then enable Pages for the `main` branch (root). Your site will be `https://<user>.github.io/<repo>/`.
- Any static host (Netlify/Vercel) will also work.

## 7) Troubleshooting

- Stuck at "Initializing": hold BOOT while starting the install.
- Wrong colors on display: flash the matching "inverse" manifest.
- Port doesn’t appear: try different cable/port, install USB drivers (CP210x/CH34x), or reboot.
- Offsets mismatch: verify offsets produced by your build and mirror them in the manifest.

## 8) Advanced

- Multiple chips: add more entries to `builds` with different `chipFamily` values.
- Single-file images: If your build produces a single combined image, reference it once at the correct offset.
- Versioning: bump the `version` field so users can see they’re flashing a new release.

---

This boilerplate mirrors the structure and behavior of the reference implementation while keeping customization simple and minimal.
