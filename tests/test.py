import numpy as np
def rotate_around_origin(x, y, angle):
    s = np.sin(angle)
    c = np.cos(angle)
    return x * c - y * s, x * s + y * c

def area(a, b, c):
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])

def convex_hull(points):
    P = sorted(points)
    upper, lower = [P[0]], [P[0]]
    for i in range(len(P)):
        if i == len(P) - 1 or area(P[0], P[i], P[-1]) > 0:
            while len(upper) >= 2 and area(upper[-2], upper[-1], P[i]) < 0:
                del upper[-1]
            upper.append(P[i])

        if i == len(P) - 1 or area(P[0], P[i], P[-1]) < 0:
            while len(lower) >= 2 and area(lower[-2], lower[-1], P[i]) > 0:
                del lower[-1]
            lower.append(P[i])
    lower = lower[::-1]
    upper.extend(lower[1:-1])
    return upper

def remove_surface_points(pts_cloud):
    d = {}
    for x in pts_cloud:
        a = int(x[1] * 100.)
        if a in d:
            d[a] += 1
        else:
            d[a] = 1

    maxi, maxiZ = 0, 0
    for a, b in d.items():
        if b > maxi:
            maxi = b
            maxiZ = a

    s1 = 0; s2 = 0; s3 = 0;
    x1 = 0; x2 = 0; x3 = 0;
    r = 2
    for i in range(len(pts_cloud)):
        a = int(pts_cloud[i][1] * 100.)
        if a < maxiZ - r or a > maxiZ + r:
            s1 += 1
            if a < maxiZ - r:
                s2 += 1
        else:
            s3 += 1
            x1 += pts_cloud[i][0]
            x2 += pts_cloud[i][1]
            x3 += pts_cloud[i][2]

    if s2 * 3 < s1:
        for i in range(len(pts_cloud)):
            a = int(pts_cloud[i][1] * 100.)
            if a <= maxiZ + r:
                x1 += pts_cloud[i][0]
                x2 += pts_cloud[i][1]
                x3 += pts_cloud[i][2]

    return x1, x2, x3, s1, s2, s3


def smallest_rectangle(P):
    N_ROTATIONS = 1000
    PI = 3.14159265
    min_area = float('inf')
    for i in range(N_ROTATIONS):
        X, Y = [], []
        angle = (2. * PI / N_ROTATIONS) * i
        for p in P:
            x, y = rotate_around_origin(p[0], p[1], angle)
            X.append(x)
            Y.append(y)
        minx, maxx = np.min(X), np.max(X)
        miny, maxy = np.min(Y), np.max(Y)
        area = (maxx - minx) * (maxy - miny)
        if area < min_area:
            min_area = area
            tmp = [(minx, maxy), (maxx, maxy), (minx, miny), (maxx, miny)]
            best_angle = angle
            res = [rotate_around_origin(x[0], x[1], -angle) for x in tmp]

    return res[0], res[-1], best_angle


if __name__ == '__main__':

    P = []
    for line in open('input.txt'):
        x = line.split(' ')
        if len(x) == 1:
            N = int(x[0])
        else:
            P.append((float(x[0]), float(x[1]), float(x[2])))

    res = smallest_rectangle(P)
    print('%.5f %.5f %.5f %.5f %.5f' % (res[0][0], res[0][1], res[1][0], res[1][1], res[2]))

    res = convex_hull(P)
    for x, y, _ in res:
        print('%.5f %.5f' % (x, y))

    res = remove_surface_points(P)
    for x in res:
        print('%.5f ' % x, end="")
