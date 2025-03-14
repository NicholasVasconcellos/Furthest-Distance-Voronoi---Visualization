import math
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

def orientation(p, q, r):
    return (q[0] - p[0])*(r[1] - p[1]) - (q[1] - p[1])*(r[0] - p[0])

def convex_hull(points):
    pts = sorted(points)
    lower = []
    for p in pts:
        while len(lower) >= 2 and orientation(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)
    upper = []
    for p in reversed(pts):
        while len(upper) >= 2 and orientation(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)
    lower.pop()
    upper.pop()
    return lower + upper

def line_from_two_points(a, b):
    (x1, y1), (x2, y2) = a, b
    dx, dy = (x2 - x1), (y2 - y1)
    length = math.hypot(dx, dy)
    A =  dy / length
    B = -dx / length
    C = A*x1 + B*y1
    return (A, B, C)

def signed_distance(line, p):
    (A, B, C) = line
    return A*p[0] + B*p[1] - C

def intersect_lines(L1, L2):
    (A1,B1,C1) = L1
    (A2,B2,C2) = L2
    det = A1*B2 - B1*A2
    if abs(det) < 1e-14:
        return None
    x = (B2*C1 - B1*C2)/det
    y = (A1*C2 - A2*C1)/det
    return (x,y)

def halfplane_intersect(region, hline):
    outpoly = []
    n = len(region)
    if n == 0:
        return []
    for i in range(n):
        cur = region[i]
        nxt = region[(i+1)%n]
        dc = signed_distance(hline, cur)
        dn = signed_distance(hline, nxt)
        if dc >= 0:
            outpoly.append(cur)
        if dc*dn < 0:  
            inter_pt = intersect_lines((hline[0],hline[1],hline[2]),
                                       line_from_two_points(cur, nxt))
            if inter_pt is not None:
                outpoly.append(inter_pt)
    return outpoly

def farthest_voronoi(points):
    hull = convex_hull(points)
    minx = min(p[0] for p in points); maxx = max(p[0] for p in points)
    miny = min(p[1] for p in points); maxy = max(p[1] for p in points)
    M = 10.0
    minx, maxx = minx - M, maxx + M
    miny, maxy = miny - M, maxy + M
    bounding_poly = [(minx,miny), (minx,maxy), (maxx,maxy), (maxx,miny)]

    def sqnorm(p): return p[0]*p[0] + p[1]*p[1]

    result = {}
    for q in hull:
        region = bounding_poly[:]
        q2 = sqnorm(q)
        for r in hull:
            if r == q:
                continue
            r2 = sqnorm(r)
            A = (r[0] - q[0])
            B = (r[1] - q[1])
            C = (r2 - q2)/2.0
            hline = (A,B,C)
            region = halfplane_intersect(region, hline)
            if len(region) == 0:
                break
        result[q] = region
    return result

if __name__ == "__main__":
    # Example sample points
    sample_points = [
        (0,0), (1,0), (2,0), (2.5,0.5),
        (1,2), (0.5,1.5), (-1,1.2), (-2,0.5),
        (0,3), (1.5,3), (2,2.8)
    ]

    # Compute farthest-Voronoi cells
    fv = farthest_voronoi(sample_points)

    # Plot the points
    xs = [p[0] for p in sample_points]
    ys = [p[1] for p in sample_points]
    plt.plot(xs, ys, 'o')  # default marker style

    # Plot each hull cell as polygon patches
    for site, poly in fv.items():
        if len(poly) < 3:
            continue  # degenerate or empty cell
        # Add a Polygon patch with default style
        patch = Polygon(poly, closed=True, fill=False)
        plt.gca().add_patch(patch)

    # Optionally, highlight hull vertices
    hull = convex_hull(sample_points)
    hx = [p[0] for p in hull] + [hull[0][0]]
    hy = [p[1] for p in hull] + [hull[0][1]]
    plt.plot(hx, hy)  # highlight hull boundary in default style

    plt.title("Farthest-Voronoi diagram (clipped) for sample points")
    plt.gca().set_aspect('equal', adjustable='datalim')
    plt.show()
