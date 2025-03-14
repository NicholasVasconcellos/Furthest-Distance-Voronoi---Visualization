# Farthest Neighbor Voronoi

Finds regions where each point is farthest from a specific point than any other point in a 2D dataset.

## How It Works

1. **Convex Hull**: First identifies the outermost points (only these have farthest Voronoi regions).
2. **Region Calculation**: For each hull point, creates its region by:
   - Starting with a large bounding box
   - Clipping away areas where other points would be farther

## Implementation Details

- Uses Andrew's monotone chain algorithm (O(n log n)) for the convex hull
- Creates regions through half-plane intersection
- Displays results with matplotlib

## Usage

You can modify the `sample_points` list to test different inputs.
