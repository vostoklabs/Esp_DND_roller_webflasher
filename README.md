# ESP32 DnD Roller (CYD — Cheap Yellow Display)

An ESP32-based dice roller for Dungeons & Dragons running on the ESP32-2432S024R (aka Cheap Yellow Display / CYD) with a built-in 2.4" 240x320 TFT.

## Features
- On-device touchscreen UI for rolling common dice (d4, d6, d8, d10, d12, d20, d100)
- Uses `TFT_eSPI` with a custom `User_Setup.h`

## Hardware
- ESP32-2432S024R (CYD) with 2.4" TFT (240x320)
- ILI9341 driver, landscape rotation (tft.setRotation(1))

## Getting Started
1. Open `dndroller.ino` in the Arduino IDE or PlatformIO.
2. Install required libraries per includes in the sketch.
3. Flash to your ESP32 board.

## Notes
- You can flash using Arduino IDE or via a Web flasher.
- No web interface is included.

## License
MIT
