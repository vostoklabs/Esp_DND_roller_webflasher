// Types must be defined before Arduino's auto-generated prototypes
struct Point3D { float x, y, z; };
struct Point2D { int x, y; };

#include <TFT_eSPI.h>
#include <TFT_eWidget.h>
#include <LittleFS.h>
using namespace fs;
#include <stdint.h>
// Use ESP32 hardware RNG for unbiased dice rolls
#include "esp_system.h"

#define CALIBRATION_FILE "/TouchCalData1"
#define REPEAT_CAL false

TFT_eSPI tft = TFT_eSPI();
ButtonWidget* diceButtons[7];  // Increased to 7 for D100
ButtonWidget* quantityUpBtn;
ButtonWidget* quantityDownBtn;
ButtonWidget* rollBtn;

char diceLabels[7][6] = {"D4", "D6", "D8", "D10", "D12", "D20", "D100"};
int diceSides[7] = {4, 6, 8, 10, 12, 20, 100};
uint8_t buttonCount = 7;

// Game state variables
int selectedDiceIndex = -1;  // Track which dice is selected
int diceQuantity = 1;        // Number of dice to roll (1-10)
int rollResults[10];         // Store individual dice results
int totalResult = 0;         // Sum of all dice

// Button position storage for redrawing
struct ButtonPos {
  int x, y, w, h;
};
ButtonPos diceButtonPos[7];

// Unbiased uniform integer generation using ESP32 hardware RNG
// Returns number in [minInclusive, maxInclusive]
static inline int uniformIntInclusive(int minInclusive, int maxInclusive) {
  if (maxInclusive <= minInclusive) return minInclusive;
  uint32_t span = (uint32_t)(maxInclusive - minInclusive + 1);
  // Rejection sampling to avoid modulo bias
  uint32_t limit = UINT32_MAX - (UINT32_MAX % span);
  uint32_t r;
  do {
    r = esp_random();
  } while (r >= limit);
  return (int)(minInclusive + (r % span));
}

static inline int rollUnbiasedDie(int sides) {
  if (sides <= 1) return 1;
  return uniformIntInclusive(1, sides);
}

// Helper to draw perfectly centered button labels (vertical + horizontal)
void drawCenteredLabel(int x, int y, int w, int h, const char* label, uint8_t textSize, uint16_t color, uint16_t bg) {
  tft.setTextSize(textSize);
  tft.setTextColor(color, bg);
  tft.setTextDatum(MC_DATUM);
  tft.drawString(label, x + w / 2, y + h / 2);
  tft.setTextDatum(TL_DATUM);
}

// Fixed 3D Animation System with correct polyhedra coordinates

// D4 - Regular tetrahedron vertices (centered and symmetric)
// Using coordinates of a regular tetra at ±1, which keeps shape rigid
Point3D tetraVertices[4] = {
  { 1.0f,  1.0f,  1.0f},
  {-1.0f, -1.0f,  1.0f},
  {-1.0f,  1.0f, -1.0f},
  { 1.0f, -1.0f, -1.0f}
};

// Tetrahedron edges (6 edges) - FIXED
int tetraEdges[6][2] = {
  {0,1}, {0,2}, {0,3},  // From top to other vertices
  {1,2}, {1,3}, {2,3}   // Connect all base vertices
};

// Tetrahedron triangular faces (for solid rendering)
// Winding chosen consistently so normals point outward
int tetraFaces[4][3] = {
  {0,1,2},
  {0,3,1},
  {0,2,3},
  {1,3,2}
};

// D6 - Cube vertices (correctly centered)
Point3D cubeVertices[8] = {
  {-0.5f, -0.5f, -0.5f}, { 0.5f, -0.5f, -0.5f}, 
  { 0.5f,  0.5f, -0.5f}, {-0.5f,  0.5f, -0.5f},  // back face
  {-0.5f, -0.5f,  0.5f}, { 0.5f, -0.5f,  0.5f}, 
  { 0.5f,  0.5f,  0.5f}, {-0.5f,  0.5f,  0.5f}   // front face
};

// Cube edges (12 edges)
int cubeEdges[12][2] = {
  {0,1}, {1,2}, {2,3}, {3,0},  // back face
  {4,5}, {5,6}, {6,7}, {7,4},  // front face
  {0,4}, {1,5}, {2,6}, {3,7}   // connecting edges
};

// Cube faces split into triangles (12 triangles total), outward winding
int cubeTriFaces[12][3] = {
  // Front z=+0.5
  {4,5,6}, {4,6,7},
  // Back z=-0.5 (reverse to point -Z)
  {0,2,1}, {0,3,2},
  // Right x=+0.5
  {1,2,6}, {1,6,5},
  // Left x=-0.5 (reverse)
  {0,7,3}, {0,4,7},
  // Top y=+0.5
  {3,6,2}, {3,7,6},
  // Bottom y=-0.5
  {0,1,5}, {0,5,4}
};

// D8 - Octahedron vertices (FIXED - proper regular octahedron)
Point3D octaVertices[6] = {
  { 0.0f,  0.0f,  0.707f},  // Top
  { 0.0f,  0.0f, -0.707f},  // Bottom
  { 0.707f,  0.0f,  0.0f},  // Right
  {-0.707f,  0.0f,  0.0f},  // Left
  { 0.0f,  0.707f,  0.0f},  // Front
  { 0.0f, -0.707f,  0.0f}   // Back
};

// Octahedron edges (12 edges) - FIXED
int octaEdges[12][2] = {
  {0,2}, {0,3}, {0,4}, {0,5},  // Top to sides
  {1,2}, {1,3}, {1,4}, {1,5},  // Bottom to sides  
  {2,4}, {4,3}, {3,5}, {5,2}   // Side square connections
};

// Octahedron faces (8 triangles), outward winding
int octaTriFaces[8][3] = {
  {0,2,4}, {0,4,3}, {0,3,5}, {0,5,2},
  {1,4,2}, {1,3,4}, {1,5,3}, {1,2,5}
};

// D10 - Pentagonal trapezohedron geometry (12 vertices)
Point3D d10Vertices[12]; // filled in initD10Geometry()

// D10 edges (30 edges)
int d10Edges[30][2];     // filled in initD10Geometry()

void initD10Geometry() {
  // Parameters (tweak for look)
  const float poleZ = 0.85f;
  const float ringZ = 0.22f;
  const float ringR = 0.62f;
  const float deg2rad = 0.01745329252f;

  // Poles
  d10Vertices[0] = {0.0f, 0.0f,  poleZ};
  d10Vertices[1] = {0.0f, 0.0f, -poleZ};

  // Upper ring (5 vertices)
  for (int i = 0; i < 5; i++) {
    float a = (72.0f * i) * deg2rad;
    d10Vertices[2 + i] = { ringR * cos(a), ringR * sin(a), ringZ };
  }

  // Lower ring (5 vertices), rotated by 36°
  for (int i = 0; i < 5; i++) {
    float a = (72.0f * i + 36.0f) * deg2rad;
    d10Vertices[7 + i] = { ringR * cos(a), ringR * sin(a), -ringZ };
  }

  // Build edges
  int k = 0;
  for (int i = 0; i < 5; i++) {
    int ui = 2 + i;
    int ui1 = 2 + ((i + 1) % 5);
    int li = 7 + i;
    int li1 = 7 + ((i + 1) % 5);
    // Pole connections
    d10Edges[k][0] = 0; d10Edges[k++][1] = ui;
    d10Edges[k][0] = 1; d10Edges[k++][1] = li;
    // Ring edges
    d10Edges[k][0] = ui;  d10Edges[k++][1] = ui1;
    d10Edges[k][0] = li;  d10Edges[k++][1] = li1;
    // Cross edges
    d10Edges[k][0] = ui;  d10Edges[k++][1] = li;
    d10Edges[k][0] = ui1; d10Edges[k++][1] = li;
  }
}

// D10 triangular faces (built in setup)
int d10TriFaces[20][3];
int d10TriCount = 0;

// D12 - Dodecahedron vertices (PROPER coordinates with normalization)
const float phi = 1.618034f;  // Golden ratio
const float norm = 0.525731f;  // Normalization factor

Point3D dodecaVertices[20] = {
  // 8 cube vertices (±1, ±1, ±1) normalized
  { 0.577f,  0.577f,  0.577f}, {-0.577f,  0.577f,  0.577f}, 
  { 0.577f, -0.577f,  0.577f}, {-0.577f, -0.577f,  0.577f},
  { 0.577f,  0.577f, -0.577f}, {-0.577f,  0.577f, -0.577f}, 
  { 0.577f, -0.577f, -0.577f}, {-0.577f, -0.577f, -0.577f},
  // 12 vertices on coordinate planes
  { 0.000f,  0.935f,  0.357f}, { 0.000f, -0.935f,  0.357f},  // XY plane
  { 0.000f,  0.935f, -0.357f}, { 0.000f, -0.935f, -0.357f},
  { 0.357f,  0.000f,  0.935f}, {-0.357f,  0.000f,  0.935f},  // XZ plane
  { 0.357f,  0.000f, -0.935f}, {-0.357f,  0.000f, -0.935f},
  { 0.935f,  0.357f,  0.000f}, {-0.935f,  0.357f,  0.000f},  // YZ plane
  { 0.935f, -0.357f,  0.000f}, {-0.935f, -0.357f,  0.000f}
};

// Dodecahedron edges (30 edges) - Proper connectivity
int dodecaEdges[30][2] = {
  // Each vertex connects to exactly 3 others
  {0,8}, {0,12}, {0,16}, {1,8}, {1,13}, {1,17}, {2,9}, {2,12}, 
  {2,18}, {3,9}, {3,13}, {3,19}, {4,10}, {4,14}, {4,16}, {5,10}, 
  {5,15}, {5,17}, {6,11}, {6,14}, {6,18}, {7,11}, {7,15}, {7,19},
  {8,10}, {9,11}, {12,13}, {14,15}, {16,18}, {17,19}
};

// Dodecahedron faces triangulated (built in setup)
int dodecaTriFaces[108][3]; // 36 faces but we may add reinforcement tris
int dodecaTriCount = 0;

// D20 - Icosahedron vertices (FIXED - proper icosahedron)
Point3D icosaVertices[12] = {
  // Standard icosahedron coordinates
  {-0.525f,  0.0f,    0.850f}, { 0.525f,  0.0f,    0.850f}, 
  {-0.525f,  0.0f,   -0.850f}, { 0.525f,  0.0f,   -0.850f},
  { 0.0f,    0.850f,  0.525f}, { 0.0f,    0.850f, -0.525f}, 
  { 0.0f,   -0.850f,  0.525f}, { 0.0f,   -0.850f, -0.525f},
  { 0.850f,  0.525f,  0.0f},   {-0.850f,  0.525f,  0.0f}, 
  { 0.850f, -0.525f,  0.0f},   {-0.850f, -0.525f,  0.0f}
};

// Icosahedron edges (30 edges) - FIXED
int icosaEdges[30][2] = {
  {0,1}, {0,4}, {0,6}, {0,9}, {0,11}, 
  {1,4}, {1,6}, {1,8}, {1,10},
  {2,3}, {2,5}, {2,7}, {2,9}, {2,11},
  {3,5}, {3,7}, {3,8}, {3,10},
  {4,5}, {4,8}, {4,9},
  {5,8}, {5,9},
  {6,7}, {6,10}, {6,11},
  {7,10}, {7,11},
  {8,10}, {9,11}
};

// Icosahedron faces (built by enumerating triangles from edges)
int icosaTriFaces[20][3];
int icosaTriCount = 0;

// Build icosahedron triangular faces from vertices and edge list
void buildIcosaFaces() {
  icosaTriCount = 0;
  const float PLANE_EPSILON = 1e-4f;
  for (int i = 0; i < 12; i++) {
    for (int j = i + 1; j < 12; j++) {
      if (!hasEdge(i, j, icosaEdges, 30)) continue;
      for (int k = j + 1; k < 12; k++) {
        if (!hasEdge(j, k, icosaEdges, 30)) continue;
        if (!hasEdge(k, i, icosaEdges, 30)) continue;

        Point3D vi = icosaVertices[i];
        Point3D vj = icosaVertices[j];
        Point3D vk = icosaVertices[k];
        float ux = vj.x - vi.x, uy = vj.y - vi.y, uz = vj.z - vi.z;
        float vx = vk.x - vi.x, vy = vk.y - vi.y, vz = vk.z - vi.z;
        float nx = uy * vz - uz * vy;
        float ny = uz * vx - ux * vz;
        float nz = ux * vy - uy * vx;

        // Skip nearly-degenerate triangles
        float nlen2 = nx*nx + ny*ny + nz*nz;
        if (nlen2 < 1e-6f) continue;

        bool allBack = true;
        bool allFront = true;
        for (int m = 0; m < 12; m++) {
          if (m == i || m == j || m == k) continue;
          Point3D vm = icosaVertices[m];
          float dx = vm.x - vi.x;
          float dy = vm.y - vi.y;
          float dz = vm.z - vi.z;
          float d = nx * dx + ny * dy + nz * dz;
          if (d > PLANE_EPSILON) allBack = false;   // point in front of plane
          if (d < -PLANE_EPSILON) allFront = false; // point behind plane
          if (!allBack && !allFront) break; // points on both sides, not a face
        }

        if (allBack || allFront) {
          // Ensure outward orientation relative to origin
          Point3D c = { (vi.x + vj.x + vk.x) / 3.0f,
                        (vi.y + vj.y + vk.y) / 3.0f,
                        (vi.z + vj.z + vk.z) / 3.0f };
          float outDot = nx * c.x + ny * c.y + nz * c.z;
          int a = i, b = j, cidx = k;
          if (outDot < 0.0f) { int tmp = b; b = cidx; cidx = tmp; }

          if (icosaTriCount < 20) {
            icosaTriFaces[icosaTriCount][0] = a;
            icosaTriFaces[icosaTriCount][1] = b;
            icosaTriFaces[icosaTriCount][2] = cidx;
            icosaTriCount++;
          }
        }
      }
    }
  }
}

// Animation state
float angleX = 0;
float angleY = 0;
float angleZ = 0;
bool animationActive = false;
bool isRolling = false;  // Track if we're in rolling animation
unsigned long lastAnimationTime = 0;
unsigned long animationStartTime = 0;

// Global orthographic scale so all dice share the same on-screen size
const float ORTHO_SCALE = 40.0f;  // tuned to roughly match D20 apparent size
// Per-die scale adjustment to visually match D20 footprint
const float D4_SCALE_FACTOR = 0.60f;  // shrink D4 further

void displayResults() {
  // Clear results area (top right, where ROLL button used to be)
  tft.fillRect(125, 5, 190, 35, TFT_BLACK);
  
  if (selectedDiceIndex == -1) return;
  
  tft.setTextSize(2);
  tft.setTextColor(TFT_YELLOW, TFT_BLACK);
  
  // Display results at top right
  if (diceQuantity > 1) {
    // Show individual results and total in compact format
    tft.setCursor(130, 8);
    tft.print("Rolls: ");
    for (int i = 0; i < diceQuantity && i < 6; i++) { // Fewer to fit with bigger text
      tft.printf("%d ", rollResults[i]);
    }
    if (diceQuantity > 6) tft.print("...");
    
    // Show total on second line
    tft.setCursor(130, 24);
    tft.setTextSize(2);
    tft.setTextColor(TFT_GREEN, TFT_BLACK);
    tft.printf("Total: %d", totalResult);
  } else {
    // Single die result - show prominently
    tft.setCursor(130, 8);
    tft.setTextSize(2);
    tft.print("Result: ");
    tft.setTextSize(3);
    tft.setTextColor(TFT_GREEN, TFT_BLACK);
    tft.printf("%d", rollResults[0]);
  }
}

// Fixed rotation with proper matrix multiplication
Point3D rotatePoint(Point3D point, float rx, float ry, float rz) {
  // Pre-calculate trigonometric values
  float cosX = cos(rx), sinX = sin(rx);
  float cosY = cos(ry), sinY = sin(ry);
  float cosZ = cos(rz), sinZ = sin(rz);
  
  Point3D result;
  
  // Apply rotations in order: X, then Y, then Z for stability
  // Rotate around X axis
  float y1 = point.y * cosX - point.z * sinX;
  float z1 = point.y * sinX + point.z * cosX;
  
  // Rotate around Y axis
  float x2 = point.x * cosY + z1 * sinY;
  float z2 = -point.x * sinY + z1 * cosY;
  
  // Rotate around Z axis
  result.x = x2 * cosZ - y1 * sinZ;
  result.y = x2 * sinZ + y1 * cosZ;
  result.z = z2;
  
  return result;
}

// Improved 3D to 2D projection
Point2D project3D(Point3D point) {
  int centerX = 220;  
  int centerY = 110;
  float scale = 140.0;  // DOUBLED from 70 to 140 for 2x size
  float distance = 3.0;  // Camera distance
  
  // Perspective projection
  float z = point.z + distance;
  if (z < 0.1f) z = 0.1f;  // Prevent division by zero
  
  float projX = (point.x * scale) / z;
  float projY = (point.y * scale) / z;
  
  Point2D result;
  result.x = centerX + (int)projX;
  result.y = centerY - (int)projY;
  return result;
}

// Orthographic projection (no perspective) — preserves shape
Point2D projectOrtho(Point3D point, float scale) {
  int centerX = 220;
  int centerY = 110;
  Point2D result;
  result.x = centerX + (int)(point.x * scale);
  result.y = centerY - (int)(point.y * scale);
  return result;
}

// Simple 16-bit (RGB565) color shading helper (0..1 factor)
static inline uint16_t shadeColor(uint16_t color, float factor) {
  if (factor < 0.0f) factor = 0.0f;
  if (factor > 1.0f) factor = 1.0f;
  uint8_t r = (color >> 11) & 0x1F;
  uint8_t g = (color >> 5) & 0x3F;
  uint8_t b = color & 0x1F;
  r = (uint8_t)(r * factor);
  g = (uint8_t)(g * factor);
  b = (uint8_t)(b * factor);
  return (uint16_t)((r << 11) | (g << 5) | b);
}

// Draw a solid D4 (tetrahedron) with painter's algorithm (back-to-front)
void drawSolidD4(Point3D rotated[4], Point2D projected[4], uint16_t baseColor) {
  struct FaceInfo { int a, b, c; float avgZ; float shade; } faces[4];

  // Precompute a simple directional light
  const float lx = 0.3f, ly = 0.6f, lz = 1.0f;
  const float llen = sqrtf(lx*lx + ly*ly + lz*lz);
  const float nxL = lx / llen, nyL = ly / llen, nzL = lz / llen;

  for (int i = 0; i < 4; i++) {
    int a = tetraFaces[i][0];
    int b = tetraFaces[i][1];
    int c = tetraFaces[i][2];

    Point3D v0 = rotated[a];
    Point3D v1 = rotated[b];
    Point3D v2 = rotated[c];

    // Face normal (right-hand rule)
    float ux = v1.x - v0.x, uy = v1.y - v0.y, uz = v1.z - v0.z;
    float vx = v2.x - v0.x, vy = v2.y - v0.y, vz = v2.z - v0.z;
    float nx = uy * vz - uz * vy;
    float ny = uz * vx - ux * vz;
    float nz = ux * vy - uy * vx;
    float nlen = sqrtf(nx*nx + ny*ny + nz*nz);
    if (nlen < 1e-6f) nlen = 1.0f;

    // Simple Lambert shading; flip normal toward viewer if needed
    float ndotl = (nx/ nlen) * nxL + (ny/ nlen) * nyL + (nz/ nlen) * nzL;
    if (ndotl < 0.0f) ndotl = 0.0f;  // no negative light
    float shade = 0.35f + 0.65f * ndotl; // keep minimum brightness

    float avgZ = (v0.z + v1.z + v2.z) / 3.0f;
    faces[i] = {a, b, c, avgZ, shade};
  }

  // Sort faces back-to-front (largest z first)
  for (int i = 0; i < 3; i++) {
    for (int j = i + 1; j < 4; j++) {
      if (faces[j].avgZ > faces[i].avgZ) {
        FaceInfo tmp = faces[i];
        faces[i] = faces[j];
        faces[j] = tmp;
      }
    }
  }

  // Draw filled triangles
  for (int i = 0; i < 4; i++) {
    Point2D p0 = projected[faces[i].a];
    Point2D p1 = projected[faces[i].b];
    Point2D p2 = projected[faces[i].c];
    uint16_t color = shadeColor(baseColor, faces[i].shade);
    tft.fillTriangle(p0.x, p0.y, p1.x, p1.y, p2.x, p2.y, color);
    // Optional edge overlay for crispness
    tft.drawLine(p0.x, p0.y, p1.x, p1.y, TFT_BLACK);
    tft.drawLine(p1.x, p1.y, p2.x, p2.y, TFT_BLACK);
    tft.drawLine(p2.x, p2.y, p0.x, p0.y, TFT_BLACK);
  }
}

// Simple 2px-thick line using offsets
static inline void drawThickLine2(int x0, int y0, int x1, int y1, uint16_t color) {
  tft.drawLine(x0, y0, x1, y1, color);
  tft.drawLine(x0 + 1, y0, x1 + 1, y1, color);
  tft.drawLine(x0, y0 + 1, x1, y1 + 1, color);
}

// Utility: check whether an undirected edge exists in an edge list
static inline bool hasEdge(int a, int b, int (*edges)[2], int numEdges) {
  for (int i = 0; i < numEdges; i++) {
    int e0 = edges[i][0];
    int e1 = edges[i][1];
    if ((e0 == a && e1 == b) || (e0 == b && e1 == a)) return true;
  }
  return false;
}

// Generic front-face-only wireframe using triangle faces
void drawWireframeFrontFaces(
  Point3D rotated[], Point2D projected[],
  int (*edges)[2], int numEdges,
  int (*triFaces)[3], int numTriFaces,
  uint16_t color
) {
  bool edgeVisible[64];
  for (int i = 0; i < 64; i++) edgeVisible[i] = false;
  const float FACE_EPSILON = 0.02f;
  const float viewX = 0.0f, viewY = 0.0f, viewZ = -1.0f; // viewer looking along -Z

  for (int f = 0; f < numTriFaces; f++) {
    int a = triFaces[f][0];
    int b = triFaces[f][1];
    int c = triFaces[f][2];

    Point3D v0 = rotated[a];
    Point3D v1 = rotated[b];
    Point3D v2 = rotated[c];

    float ux = v1.x - v0.x, uy = v1.y - v0.y, uz = v1.z - v0.z;
    float vx = v2.x - v0.x, vy = v2.y - v0.y, vz = v2.z - v0.z;
    float nx = uy * vz - uz * vy;
    float ny = uz * vx - ux * vz;
    float nz = ux * vy - uy * vx;
    // Ensure outward normal regardless of triangle winding by checking center
    Point3D ctri = { (v0.x + v1.x + v2.x) / 3.0f,
                     (v0.y + v1.y + v2.y) / 3.0f,
                     (v0.z + v1.z + v2.z) / 3.0f };
    float outDot = nx * ctri.x + ny * ctri.y + nz * ctri.z;
    if (outDot < 0.0f) { nx = -nx; ny = -ny; nz = -nz; }
    float facing = nx * viewX + ny * viewY + nz * viewZ;
    if (facing > FACE_EPSILON) {
      int faceEdges[3][2] = { {a,b}, {b,c}, {c,a} };
      for (int e = 0; e < 3; e++) {
        int u = faceEdges[e][0];
        int v = faceEdges[e][1];
        // Mark if the physical edge exists; otherwise draw the triangle edge anyway
        for (int i = 0; i < numEdges; i++) {
          int e0 = edges[i][0];
          int e1 = edges[i][1];
          if ((e0 == u && e1 == v) || (e0 == v && e1 == u)) {
            edgeVisible[i] = true;
            break;
          }
        }
      }
    }
  }

  for (int i = 0; i < numEdges; i++) {
    if (!edgeVisible[i]) continue;
    int a = edges[i][0];
    int b = edges[i][1];
    drawThickLine2(projected[a].x, projected[a].y, projected[b].x, projected[b].y, color);
  }
}

// Wireframe D4 showing only front-facing edges (orthographic)
void drawWireframeD4Ortho(Point3D rotated[4], uint16_t color) {
  // Project using orthographic with per-die scale
  Point2D proj[4];
  float scale = ORTHO_SCALE * D4_SCALE_FACTOR;
  for (int i = 0; i < 4; i++) proj[i] = projectOrtho(rotated[i], scale);

  // Mark edges that belong to front-facing faces only
  bool edgeVisible[6] = {false,false,false,false,false,false};
  const float FACE_EPSILON = 0.02f; // threshold to reduce flicker near silhouette
  const float viewX = 0.0f, viewY = 0.0f, viewZ = 1.0f; // looking along +Z

  for (int f = 0; f < 4; f++) {
    int a = tetraFaces[f][0];
    int b = tetraFaces[f][1];
    int c = tetraFaces[f][2];

    Point3D v0 = rotated[a];
    Point3D v1 = rotated[b];
    Point3D v2 = rotated[c];

    // Face normal
    float ux = v1.x - v0.x, uy = v1.y - v0.y, uz = v1.z - v0.z;
    float vx = v2.x - v0.x, vy = v2.y - v0.y, vz = v2.z - v0.z;
    float nx = uy * vz - uz * vy;
    float ny = uz * vx - ux * vz;
    float nz = ux * vy - uy * vx;

    // Facing if normal has positive dot with view direction
    float facing = nx * viewX + ny * viewY + nz * viewZ;
    if (facing > FACE_EPSILON) {
      int faceEdges[3][2] = { {a,b}, {b,c}, {c,a} };
      for (int e = 0; e < 3; e++) {
        int u = faceEdges[e][0];
        int v = faceEdges[e][1];
        for (int i = 0; i < 6; i++) {
          int e0 = tetraEdges[i][0];
          int e1 = tetraEdges[i][1];
          if ((e0 == u && e1 == v) || (e0 == v && e1 == u)) {
            edgeVisible[i] = true;
            break;
          }
        }
      }
    }
  }

  // Draw only visible edges (front-facing)
  for (int i = 0; i < 6; i++) {
    if (!edgeVisible[i]) continue;
    int a = tetraEdges[i][0];
    int b = tetraEdges[i][1];
    drawThickLine2(proj[a].x, proj[a].y, proj[b].x, proj[b].y, color);
  }
}

void animateDice() {
  if (selectedDiceIndex == -1) return;
  
  // Animation timing
  if (millis() - lastAnimationTime < 50) return;  // 20 FPS
  lastAnimationTime = millis();
  
  // Clear animation area
  tft.fillRect(125, 45, 190, 130, TFT_BLACK);
  
  // Clip drawing to the dice area to prevent edge popping/flicker
  // Use absolute coordinates (vpDatum=false) so our projected points remain screen-based
  tft.setViewport(125, 45, 190, 130, false);
  
  // Variable rotation speed based on state
  if (isRolling) {
    // Fast spin when rolling
    angleX += 0.3;
    angleY += 0.25;
    angleZ += 0.2;
  } else if (animationActive) {
    // Slow down after roll
    float slowdownFactor = 1.0 - ((millis() - animationStartTime - 2000) / 1000.0);
    if (slowdownFactor < 0.1) slowdownFactor = 0.1;
    angleX += 0.3 * slowdownFactor;
    angleY += 0.25 * slowdownFactor;
    angleZ += 0.2 * slowdownFactor;
  } else {
    // Slow idle spin
    angleX += 0.02;
    angleY += 0.015;
    angleZ += 0.01;
  }
  
  // Choose geometry based on dice type
  uint16_t color;
  Point3D* vertices;
  int numVertices;
  int (*edges)[2];
  int numEdges;
  
  switch(selectedDiceIndex) {
    case 0: // D4
      color = TFT_RED;
      vertices = tetraVertices;
      numVertices = 4;
      edges = tetraEdges;
      numEdges = 6;
      break;
    case 1: // D6
      color = TFT_GREEN;
      vertices = cubeVertices;
      numVertices = 8;
      edges = cubeEdges;
      numEdges = 12;
      break;
    case 2: // D8
      color = TFT_BLUE;
      vertices = octaVertices;
      numVertices = 6;
      edges = octaEdges;
      numEdges = 12;
      break;
    case 3: // D10
      color = TFT_CYAN;
      vertices = d10Vertices;
      numVertices = 12;
      edges = d10Edges;
      numEdges = 30;
      break;
    case 4: // D12
      color = TFT_MAGENTA;
      vertices = dodecaVertices;
      numVertices = 20;
      edges = dodecaEdges;
      numEdges = 30;
      break;
    case 5: // D20
      color = TFT_YELLOW;
      vertices = icosaVertices;
      numVertices = 12;
      edges = icosaEdges;
      numEdges = 30;
      break;
    default: // D100 (same as D20)
      color = TFT_WHITE;
      vertices = icosaVertices;
      numVertices = 12;
      edges = icosaEdges;
      numEdges = 30;
      break;
  }
  
  // Project all vertices (and keep rotated positions for face normals)
  Point3D rotatedVertices[20];
  Point2D projectedVertices[20];  // Max needed for dodecahedron
  for (int i = 0; i < numVertices; i++) {
    rotatedVertices[i] = rotatePoint(vertices[i], angleX, angleY, angleZ);
    // Use orthographic for consistent size
    projectedVertices[i] = projectOrtho(rotatedVertices[i], ORTHO_SCALE);
  }

  // D4: orthographic wireframe with hidden edges to preserve shape
  if (selectedDiceIndex == 0) {
    drawWireframeD4Ortho(rotatedVertices, color);
    // Restore full-screen drawing
    tft.resetViewport();
    return;
  }
  
  // Draw wireframe for other dice: front-facing only
  if (selectedDiceIndex == 1) { // D6
    drawWireframeFrontFaces(rotatedVertices, projectedVertices, cubeEdges, 12, cubeTriFaces, 12, color);
  } else if (selectedDiceIndex == 2) { // D8
    drawWireframeFrontFaces(rotatedVertices, projectedVertices, octaEdges, 12, octaTriFaces, 8, color);
  } else if (selectedDiceIndex == 3) { // D10
    drawWireframeFrontFaces(rotatedVertices, projectedVertices, d10Edges, 30, d10TriFaces, d10TriCount, color);
  } else if (selectedDiceIndex == 4) { // D12
    drawWireframeFrontFaces(rotatedVertices, projectedVertices, dodecaEdges, 30, dodecaTriFaces, dodecaTriCount, color);
  } else { // D20 and D100
    drawWireframeFrontFaces(rotatedVertices, projectedVertices, icosaEdges, 30, icosaTriFaces, icosaTriCount, color);
  }

  // Restore full-screen drawing
  tft.resetViewport();
}

bool resultsShown = false;

void rollDice() {
  if (selectedDiceIndex == -1) return;
  
  // Start 3D animation
  animationActive = true;
  isRolling = true;  // Start fast spin
  animationStartTime = millis();
  resultsShown = false;
  lastAnimationTime = millis();
  
  // Reset rotation angles for fresh animation
  angleX = 0;
  angleY = 0;
  angleZ = 0;
  
  totalResult = 0;
  int sides = diceSides[selectedDiceIndex];
  
  // Roll multiple dice
  for (int i = 0; i < diceQuantity; i++) {
    rollResults[i] = rollUnbiasedDie(sides);
    totalResult += rollResults[i];
  }
  
  // Animation will play, then results will show
}

void updateQuantityDisplay() {
  // Clear quantity display area (between - and + buttons)
  tft.fillRect(252, 187, 16, 28, TFT_BLACK);
  tft.setTextSize(2);
  tft.setTextColor(TFT_WHITE, TFT_BLACK);
  tft.setTextDatum(MC_DATUM);
  tft.drawNumber(diceQuantity, 260, 201); // center of 16x28 area
  tft.setTextDatum(TL_DATUM);
}

void updateDiceSelection() {
  // Update all dice button colors to show selection
  for (int i = 0; i < buttonCount; i++) {
    uint16_t fill = (i == selectedDiceIndex) ? TFT_GREEN : TFT_BLUE;
    diceButtons[i]->initButtonUL(diceButtonPos[i].x, diceButtonPos[i].y, 
                                 diceButtonPos[i].w, diceButtonPos[i].h,
                                 TFT_WHITE, fill, TFT_WHITE, "", 2);
    diceButtons[i]->drawSmoothButton(false, 3, TFT_BLACK);
    drawCenteredLabel(diceButtonPos[i].x, diceButtonPos[i].y,
                      diceButtonPos[i].w, diceButtonPos[i].h,
                      diceLabels[i], 2, TFT_WHITE, fill);
  }
}

// Dice selection actions
void btn0_action() { selectedDiceIndex = 0; updateDiceSelection(); }
void btn1_action() { selectedDiceIndex = 1; updateDiceSelection(); }
void btn2_action() { selectedDiceIndex = 2; updateDiceSelection(); }
void btn3_action() { selectedDiceIndex = 3; updateDiceSelection(); }
void btn4_action() { selectedDiceIndex = 4; updateDiceSelection(); }
void btn5_action() { selectedDiceIndex = 5; updateDiceSelection(); }
void btn6_action() { selectedDiceIndex = 6; updateDiceSelection(); }

void quantityUp_action() {
  if (diceQuantity < 10) {
    diceQuantity++;
    updateQuantityDisplay();
  }
}

void quantityDown_action() {
  if (diceQuantity > 1) {
    diceQuantity--;
    updateQuantityDisplay();
  }
}

void roll_action() {
  rollDice();
}

void (*btnActions[])() = {
  btn0_action, btn1_action, btn2_action, btn3_action, 
  btn4_action, btn5_action, btn6_action
};

void setupDiceButtons() {
  // Left column - Rectangular dice buttons
  int btnWidth = 110;
  int btnHeight = 30;
  int spacing = 3;
  int leftX = 5;
  int startY = 5;

  // Create dice buttons in left column
  for (int i = 0; i < buttonCount; i++) {
    int x = leftX;
    int y = startY + i * (btnHeight + spacing);

    // Store button position for later reference
    diceButtonPos[i] = {x, y, btnWidth, btnHeight};

    diceButtons[i] = new ButtonWidget(&tft);
    diceButtons[i]->initButtonUL(x, y, btnWidth, btnHeight, TFT_WHITE, TFT_BLUE, TFT_WHITE, "", 2);
    diceButtons[i]->setPressAction(btnActions[i]);
    diceButtons[i]->drawSmoothButton(false, 3, TFT_BLACK);
    drawCenteredLabel(x, y, btnWidth, btnHeight, diceLabels[i], 2, TFT_WHITE, TFT_BLUE);
  }
  
  // Right column - ROLL button at bottom
  int rollX = 125;
  int rollY = 180;
  rollBtn = new ButtonWidget(&tft);
  rollBtn->initButtonUL(rollX, rollY, 90, 35, TFT_WHITE, TFT_RED, TFT_WHITE, "", 3);
  rollBtn->setPressAction(roll_action);
  rollBtn->drawSmoothButton(false, 3, TFT_BLACK);
  drawCenteredLabel(rollX, rollY, 90, 35, "ROLL", 3, TFT_WHITE, TFT_RED);
  
  // Quantity controls next to ROLL button
  int qtyY = 180;
  quantityDownBtn = new ButtonWidget(&tft);
  quantityDownBtn->initButtonUL(220, qtyY, 30, 35, TFT_WHITE, TFT_CYAN, TFT_BLACK, "", 3);
  quantityDownBtn->setPressAction(quantityDown_action);
  quantityDownBtn->drawSmoothButton(false, 2, TFT_BLACK);
  drawCenteredLabel(220, qtyY, 30, 35, "-", 3, TFT_BLACK, TFT_CYAN);
  
  quantityUpBtn = new ButtonWidget(&tft);
  quantityUpBtn->initButtonUL(270, qtyY, 30, 35, TFT_WHITE, TFT_CYAN, TFT_BLACK, "", 3);
  quantityUpBtn->setPressAction(quantityUp_action);
  quantityUpBtn->drawSmoothButton(false, 2, TFT_BLACK);
  drawCenteredLabel(270, qtyY, 30, 35, "+", 3, TFT_BLACK, TFT_CYAN);
  
  // Initial quantity display
  updateQuantityDisplay();
}

void touch_calibrate() {
  uint16_t calData[5];
  bool calDataOK = false;

  if (!LittleFS.begin()) {
    LittleFS.format();
    LittleFS.begin();
  }

  if (LittleFS.exists(CALIBRATION_FILE)) {
    if (!REPEAT_CAL) {
      File f = LittleFS.open(CALIBRATION_FILE, "r");
      if (f && f.readBytes((char *)calData, 14) == 14)
        calDataOK = true;
      f.close();
    } else {
      LittleFS.remove(CALIBRATION_FILE);
    }
  }

  if (calDataOK) {
    tft.setTouch(calData);
  } else {
    tft.fillScreen(TFT_BLACK);
    tft.setCursor(20, 0);
    tft.setTextFont(2);
    tft.setTextSize(1);
    tft.setTextColor(TFT_WHITE, TFT_BLACK);
    tft.println("Touch corners as indicated");
    tft.calibrateTouch(calData, TFT_MAGENTA, TFT_BLACK, 15);
    File f = LittleFS.open(CALIBRATION_FILE, "w");
    if (f) {
      f.write((const unsigned char *)calData, 14);
      f.close();
    }
    tft.setTouch(calData);
  }
}

void setup() {
  Serial.begin(115200);
  randomSeed(analogRead(0));
  tft.begin();
  tft.setRotation(1);  // Landscape 90° right
  tft.fillScreen(TFT_BLACK);
  tft.setTextFont(1);
  tft.setTextSize(1);

  // Initialize generated geometry/face tables
  initD10Geometry();
  // Build D10 faces (two triangles per kite face)
  d10TriCount = 0;
  for (int i = 0; i < 5; i++) {
    int ui  = 2 + i;
    int ui1 = 2 + ((i + 1) % 5);
    int li  = 7 + i;
    // upper kite split
    d10TriFaces[d10TriCount][0] = 0; d10TriFaces[d10TriCount][1] = ui;  d10TriFaces[d10TriCount++][2] = li;
    d10TriFaces[d10TriCount][0] = 0; d10TriFaces[d10TriCount][1] = li;  d10TriFaces[d10TriCount++][2] = ui1;
    // lower kite split (toward bottom pole)
    d10TriFaces[d10TriCount][0] = 1; d10TriFaces[d10TriCount][1] = li;  d10TriFaces[d10TriCount++][2] = ui;
    d10TriFaces[d10TriCount][0] = 1; d10TriFaces[d10TriCount][1] = ui;  d10TriFaces[d10TriCount++][2] = ((7 + ((i + 4) % 5))); // previous lower ring vertex
  }

  // Dodeca triangulation from known pentagon cycles
  const int pentagons[12][5] = {
    {0,8,10,4,16}, {1,8,10,5,17}, {2,9,11,6,18}, {3,9,11,7,19},
    {0,12,13,1,8}, {2,12,13,3,9}, {4,14,15,5,10}, {6,14,15,7,11},
    {0,16,18,2,12}, {4,16,18,6,14}, {1,17,19,3,13}, {5,17,19,7,15}
  };
  dodecaTriCount = 0;
  for (int p = 0; p < 12; p++) {
    int v0 = pentagons[p][0];
    for (int i = 1; i < 4; i++) {
      int v1 = pentagons[p][i];
      int v2 = pentagons[p][i+1];
      dodecaTriFaces[dodecaTriCount][0] = v0;
      dodecaTriFaces[dodecaTriCount][1] = v1;
      dodecaTriFaces[dodecaTriCount++][2] = v2;
    }
    // Add reinforcement triangles to improve culling stability near pentagon center
    int v1 = pentagons[p][1], v2 = pentagons[p][2], v3 = pentagons[p][3], v4 = pentagons[p][4];
    dodecaTriFaces[dodecaTriCount][0] = v1; dodecaTriFaces[dodecaTriCount][1] = v2; dodecaTriFaces[dodecaTriCount++][2] = v3;
    dodecaTriFaces[dodecaTriCount][0] = v1; dodecaTriFaces[dodecaTriCount][1] = v3; dodecaTriFaces[dodecaTriCount++][2] = v4;
  }

  // Build robust icosa triangular faces from edges (avoids indexing mismatch)
  buildIcosaFaces();

  touch_calibrate();
  setupDiceButtons();
}

void loop() {
  // Run dice animation continuously when a dice is selected
  if (selectedDiceIndex != -1) {
    animateDice();
    
    // Check if rolling animation should transition to slowdown
    if (isRolling && (millis() - animationStartTime > 2000)) {
      isRolling = false;  // Stop fast spin after 2 seconds
    }
    
    // Check if animation should end and show results
    if (!resultsShown && animationActive && (millis() - animationStartTime > 3000)) {
      animationActive = false;
      resultsShown = true;
      displayResults();
    }
  }
  
  static uint32_t lastScan = 0;
  if (millis() - lastScan >= 50) {
    uint16_t x, y;
    bool touched = tft.getTouch(&x, &y);
    
    // Check all buttons
    ButtonWidget* allButtons[] = {
      diceButtons[0], diceButtons[1], diceButtons[2], diceButtons[3],
      diceButtons[4], diceButtons[5], diceButtons[6],
      quantityUpBtn, quantityDownBtn, rollBtn
    };
    int totalButtons = 10;
    
    if (touched) {
      for (int i = 0; i < totalButtons; i++) {
        if (allButtons[i]->contains(x, y)) {
          allButtons[i]->press(true);
          allButtons[i]->pressAction();
          break;
        }
      }
    } else {
      // Release all buttons
      for (int i = 0; i < totalButtons; i++) {
        allButtons[i]->press(false);
      }
    }
    
    lastScan = millis();
  }
}