import math
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon


def orientation(p, q, r):
    # Determine the orientation of the triplet (p, q, r)
    return (q[0] - p[0]) * (r[1] - p[1]) - (q[1] - p[1]) * (r[0] - p[0])


def convex_hull(points):
    # Find the convex hull of a set of 2D points
    pts = sorted(points)  # Sort points by x, then by y
    lower = []
    for p in pts:
        while len(lower) >= 2 and orientation(lower[-2], lower[-1], p) <= 0:
            lower.pop()  # Remove last point if it doesn't make a left turn
        lower.append(p)
    upper = []
    for p in reversed(pts):
        while len(upper) >= 2 and orientation(upper[-2], upper[-1], p) <= 0:
            upper.pop()  # Remove last point if it doesn't make a left turn
        upper.append(p)
    lower.pop()
    upper.pop()
    return lower + upper  # Combine lower and upper hull


def line_from_two_points(a, b):
    # Get the line equation (A, B, C) from two points a and b
    (x1, y1), (x2, y2) = a, b
    dx, dy = (x2 - x1), (y2 - y1)
    length = math.hypot(dx, dy)
    A = dy / length
    B = -dx / length
    C = A * x1 + B * y1
    return (A, B, C)


def signed_distance(line, p):
    # Get the signed distance from point p to the line (A, B, C)
    (A, B, C) = line
    return A * p[0] + B * p[1] - C


def intersect_lines(L1, L2):
    # Find the intersection point of two lines L1 and L2
    (A1, B1, C1) = L1
    (A2, B2, C2) = L2
    det = A1 * B2 - B1 * A2  # Calculate the determinant
    if abs(det) < 1e-14:
        return None  # Lines are parallel or coincident
    x = (B2 * C1 - B1 * C2) / det
    y = (A1 * C2 - A2 * C1) / det
    return (x, y)


def halfplane_intersect(region, hline):
    # Get the intersection of a polygon region with a half-plane defined by hline
    outpoly = []
    n = len(region)
    if n == 0:
        return []
    for i in range(n):
        cur = region[i]
        nxt = region[(i + 1) % n]
        dc = signed_distance(hline, cur)
        dn = signed_distance(hline, nxt)
        if dc >= 0:
            outpoly.append(cur)  # Add current point if it's inside the half-plane
        if dc * dn < 0:
            inter_pt = intersect_lines(
                (hline[0], hline[1], hline[2]), line_from_two_points(cur, nxt)
            )
            if inter_pt is not None:
                outpoly.append(inter_pt)  # Add intersection point if it exists
    return outpoly


def farthest_voronoi(points):
    # Compute the farthest-Voronoi diagram for a set of points
    hull = convex_hull(points)  # Get the convex hull of the points
    minx = min(p[0] for p in points)
    maxx = max(p[0] for p in points)
    miny = min(p[1] for p in points)
    maxy = max(p[1] for p in points)
    M = 10.0  # Margin for bounding box
    minx, maxx = minx - M, maxx + M
    miny, maxy = miny - M, maxy + M
    bounding_poly = [
        (minx, miny),
        (minx, maxy),
        (maxx, maxy),
        (maxx, miny),
    ]  # Define bounding polygon

    def sqnorm(p):
        return p[0] * p[0] + p[1] * p[1]  # Compute squared norm of a point

    result = {}
    for q in hull:
        region = bounding_poly[:]  # Start with the bounding polygon
        q2 = sqnorm(q)
        for r in hull:
            if r == q:
                continue
            r2 = sqnorm(r)
            A = r[0] - q[0]
            B = r[1] - q[1]
            C = (r2 - q2) / 2.0
            hline = (A, B, C)
            region = halfplane_intersect(region, hline)  # Intersect with half-plane
            if len(region) == 0:
                break
        result[q] = region  # Store the region for point q
    return result


if __name__ == "__main__":
    # Example sample points
    sample_points = [
        (0, 0),
        (1, 0),
        (2, 0),
        (2.5, 0.5),
        (1, 2),
        (0.5, 1.5),
        (-1, 1.2),
        (-2, 0.5),
        (0, 3),
        (1.5, 3),
        (2, 2.8),
    ]

    # Compute farthest-Voronoi cells
    fv = farthest_voronoi(sample_points)

    # Plot the points
    xs = [p[0] for p in sample_points]
    ys = [p[1] for p in sample_points]
    plt.plot(xs, ys, "o")  # Plot points with default marker style

    # Plot each hull cell as polygon patches
    for site, poly in fv.items():
        if len(poly) < 3:
            continue  # Skip degenerate or empty cells
        # Add a Polygon patch with default style
        patch = Polygon(poly, closed=True, fill=False)
        plt.gca().add_patch(patch)

    # highlight hull vertices
    hull = convex_hull(sample_points)
    hx = [p[0] for p in hull] + [hull[0][0]]
    hy = [p[1] for p in hull] + [hull[0][1]]
    plt.plot(hx, hy)  # Highlight hull boundary with default style

    plt.title("Farthest-Voronoi diagram") 
    # Set aspect ratio
    plt.gca().set_aspect("equal", adjustable="datalim")  
    plt.show()
