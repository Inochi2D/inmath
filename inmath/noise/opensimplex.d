/**
    OpenSimplex2 implementation, based on open_simplex_2.
    This implementation is based on K.jpeg's faster OpenSimplex2 variant.

    See: https://github.com/KdotJPG/OpenSimplex2
*/
module inmath.noise.opensimplex;
import inmath.math;

private {
    struct LatticePoint2D {
        int xsv, ysv;
        double dx, dy;
        this(int xsv, int ysv) {
            this.xsv = xsv; this.ysv = ysv;
            double ssv = (xsv + ysv) * -0.211324865405187;
            this.dx = -xsv - ssv;
            this.dy = -ysv - ssv;
        }
    }

    struct LatticePoint3D {
        double dxr, dyr, dzr;
        int xrv, yrv, zrv;
        LatticePoint3D* nextOnFailure, nextOnSuccess;

        this(int xrv, int yrv, int zrv, int lattice) {
            this.dxr = -xrv + lattice * 0.5; this.dyr = -yrv + lattice * 0.5; this.dzr = -zrv + lattice * 0.5;
            this.xrv = xrv + lattice * 1024; this.yrv = yrv + lattice * 1024; this.zrv = zrv + lattice * 1024;
        }
    }

    struct LatticePoint4D {
        int xsv, ysv, zsv, wsv;
        double dx, dy, dz, dw;
        double xsi, ysi, zsi, wsi;
        double ssiDelta;
        this(int xsv, int ysv, int zsv, int wsv) {
            this.xsv = xsv + 409; this.ysv = ysv + 409; this.zsv = zsv + 409; this.wsv = wsv + 409;
            double ssv = (xsv + ysv + zsv + wsv) * 0.309016994374947;
            this.dx = -xsv - ssv;
            this.dy = -ysv - ssv;
            this.dz = -zsv - ssv;
            this.dw = -wsv - ssv;
            this.xsi = xsi = 0.2 - xsv;
            this.ysi = ysi = 0.2 - ysv;
            this.zsi = zsi = 0.2 - zsv;
            this.wsi = wsi = 0.2 - wsv;
            this.ssiDelta = (0.8 - xsv - ysv - zsv - wsv) * 0.309016994374947;
        }
    }

    /*
     * Gradients
     */
    struct Grad2 {
        double dx, dy;
        this(double dx, double dy) {
            this.dx = dx; this.dy = dy;
        }
    }

    struct Grad3 {
        double dx, dy, dz;
        this(double dx, double dy, double dz) {
            this.dx = dx; this.dy = dy; this.dz = dz;
        }
    }

    struct Grad4 {
        double dx, dy, dz, dw;
        this(double dx, double dy, double dz,  double dw) {
            this.dx = dx; this.dy = dy; this.dz = dz; this.dw = dw;
        }
    }


    enum int PSIZE = 2048;
    int PMASK = 2047;

    __gshared short[] perm;
    __gshared Grad2[] permGrad2;
    __gshared Grad3[] permGrad3;
    __gshared Grad4[] permGrad4;
    

    __gshared LatticePoint2D*[] LOOKUP_2D;
    __gshared LatticePoint3D*[] LOOKUP_3D;
    __gshared LatticePoint4D*[] VERTICES_4D;
    __gshared Grad2[] GRADIENTS_2D;
    __gshared Grad3[] GRADIENTS_3D;
    __gshared Grad4[] GRADIENTS_4D;
    
    enum double N2 = 0.01001634121365712;
    enum double N3 = 0.030485933181293584;
    enum double N4 = 0.009202377986303158;

    double _osnoise2(double xs, double ys) {
        double value = 0;

        // Get base points and offsets
        int xsb = cast(int)floor(xs), ysb = cast(int)floor(ys);
        double xsi = xs - xsb, ysi = ys - ysb;

        // Index to point list
        int index = cast(int)((ysi - xsi) / 2 + 1);

        double ssi = (xsi + ysi) * -0.211324865405187;
        double xi = xsi + ssi, yi = ysi + ssi;

        // Point contributions
        for (int i = 0; i < 3; i++) {
            LatticePoint2D* c = LOOKUP_2D[index + i];

            double dx = xi + c.dx, dy = yi + c.dy;
            double attn = 0.5 - dx * dx - dy * dy;
            if (attn <= 0) continue;

            int pxm = (xsb + c.xsv) & PMASK, pym = (ysb + c.ysv) & PMASK;
            Grad2 grad = permGrad2[perm[pxm] ^ pym];
            double extrapolation = grad.dx * dx + grad.dy * dy;

            attn *= attn;
            value += attn * attn * extrapolation;
        }

        return value;
    }

    double _osnoise3(double xr, double yr, double zr) {

        // Get base and offsets inside cube of first lattice.
        int xrb = cast(int)floor(xr), yrb = cast(int)floor(yr), zrb = cast(int)floor(zr);
        double xri = xr - xrb, yri = yr - yrb, zri = zr - zrb;

        // Identify which octant of the cube we're in. This determines which cell
        // in the other cubic lattice we're in, and also narrows down one point on each.
        int xht = cast(int)(xri + 0.5), yht = cast(int)(yri + 0.5), zht = cast(int)(zri + 0.5);
        int index = (xht << 0) | (yht << 1) | (zht << 2);

        // Point contributions
        double value = 0;
        LatticePoint3D* c = LOOKUP_3D[index];
        while (c !is null) {
            double dxr = xri + c.dxr, dyr = yri + c.dyr, dzr = zri + c.dzr;
            double attn = 0.5 - dxr * dxr - dyr * dyr - dzr * dzr;
            if (attn < 0) {
                c = c.nextOnFailure;
            } else {
                int pxm = (xrb + c.xrv) & PMASK, pym = (yrb + c.yrv) & PMASK, pzm = (zrb + c.zrv) & PMASK;
                Grad3 grad = permGrad3[perm[perm[pxm] ^ pym] ^ pzm];
                double extrapolation = grad.dx * dxr + grad.dy * dyr + grad.dz * dzr;

                attn *= attn;
                value += attn * attn * extrapolation;
                c = c.nextOnSuccess;
            }
        }
        return value;
    }

    double _osnoise4(double xs, double ys, double zs, double ws) {
        double value = 0;

        // Get base points and offsets
        int xsb = cast(int)floor(xs), ysb = cast(int)floor(ys), zsb = cast(int)floor(zs), wsb = cast(int)floor(ws);
        double xsi = xs - xsb, ysi = ys - ysb, zsi = zs - zsb, wsi = ws - wsb;

        // If we're in the lower half, flip so we can repeat the code for the upper half. We'll flip back later.
        double siSum = xsi + ysi + zsi + wsi;
        double ssi = siSum * 0.309016994374947; // Prep for vertex contributions.
        bool inLowerHalf = (siSum < 2);
        if (inLowerHalf) {
            xsi = 1 - xsi; ysi = 1 - ysi; zsi = 1 - zsi; wsi = 1 - wsi;
            siSum = 4 - siSum;
        }

        // Consider opposing vertex pairs of the octahedron formed by the central cross-section of the stretched tesseract
        double aabb = xsi + ysi - zsi - wsi, abab = xsi - ysi + zsi - wsi, abba = xsi - ysi - zsi + wsi;
        double aabbScore = abs(aabb), ababScore = abs(abab), abbaScore = abs(abba);

        // Find the closest point on the stretched tesseract as if it were the upper half
        int vertexIndex, via, vib;
        double asi, bsi;
        if (aabbScore > ababScore && aabbScore > abbaScore) {
            if (aabb > 0) {
                asi = zsi; bsi = wsi; vertexIndex = 0b0011; via = 0b0111; vib = 0b1011;
            } else {
                asi = xsi; bsi = ysi; vertexIndex = 0b1100; via = 0b1101; vib = 0b1110;
            }
        } else if (ababScore > abbaScore) {
            if (abab > 0) {
                asi = ysi; bsi = wsi; vertexIndex = 0b0101; via = 0b0111; vib = 0b1101;
            } else {
                asi = xsi; bsi = zsi; vertexIndex = 0b1010; via = 0b1011; vib = 0b1110;
            }
        } else {
            if (abba > 0) {
                asi = ysi; bsi = zsi; vertexIndex = 0b1001; via = 0b1011; vib = 0b1101;
            } else {
                asi = xsi; bsi = wsi; vertexIndex = 0b0110; via = 0b0111; vib = 0b1110;
            }
        }
        if (bsi > asi) {
            via = vib;
            double temp = bsi;
            bsi = asi;
            asi = temp;
        }
        if (siSum + asi > 3) {
            vertexIndex = via;
            if (siSum + bsi > 4) {
                vertexIndex = 0b1111;
            }
        }

        // Now flip back if we're actually in the lower half.
        if (inLowerHalf) {
            xsi = 1 - xsi; ysi = 1 - ysi; zsi = 1 - zsi; wsi = 1 - wsi;
            vertexIndex ^= 0b1111;
        }

        // Five points to add, total, from five copies of the A4 lattice.
        for (int i = 0; i < 5; i++) {

            // Update xsb/etc. and add the lattice point's contribution.
            LatticePoint4D* c = VERTICES_4D[vertexIndex];
            xsb += c.xsv; ysb += c.ysv; zsb += c.zsv; wsb += c.wsv;
            double xi = xsi + ssi, yi = ysi + ssi, zi = zsi + ssi, wi = wsi + ssi;
            double dx = xi + c.dx, dy = yi + c.dy, dz = zi + c.dz, dw = wi + c.dw;
            double attn = 0.5 - dx * dx - dy * dy - dz * dz - dw * dw;
            if (attn > 0) {
                int pxm = xsb & PMASK, pym = ysb & PMASK, pzm = zsb & PMASK, pwm = wsb & PMASK;
                Grad4 grad = permGrad4[perm[perm[perm[pxm] ^ pym] ^ pzm] ^ pwm];
                double ramped = grad.dx * dx + grad.dy * dy + grad.dz * dz + grad.dw * dw;

                attn *= attn;
                value += attn * attn * ramped;
            }

            // Maybe this helps the compiler/JVM/LLVM/etc. know we can end the loop here. Maybe not.
            if (i == 4) break;

            // Update the relative skewed coordinates to reference the vertex we just added.
            // Rather, reference its counterpart on the lattice copy that is shifted down by
            // the vector <-0.2, -0.2, -0.2, -0.2>
            xsi += c.xsi; ysi += c.ysi; zsi += c.zsi; wsi += c.wsi;
            ssi += c.ssiDelta;

            // Next point is the closest vertex on the 4-simplex whose base vertex is the aforementioned vertex.
            double score0 = 1.0 + ssi * (-1.0 / 0.309016994374947); // Seems slightly faster than 1.0-xsi-ysi-zsi-wsi
            vertexIndex = 0b0000;
            if (xsi >= ysi && xsi >= zsi && xsi >= wsi && xsi >= score0) {
                vertexIndex = 0b0001;
            }
            else if (ysi > xsi && ysi >= zsi && ysi >= wsi && ysi >= score0) {
                vertexIndex = 0b0010;
            }
            else if (zsi > xsi && zsi > ysi && zsi >= wsi && zsi >= score0) {
                vertexIndex = 0b0100;
            }
            else if (wsi > xsi && wsi > ysi && wsi > zsi && wsi >= score0) {
                vertexIndex = 0b1000;
            }
        }

        return value;
    }
}

/**
    Seeds the OpenSimplex algorithm
*/
void osseed(ulong seed) {
    perm = new short[PSIZE];
    permGrad2 = new Grad2[PSIZE];
    permGrad3 = new Grad3[PSIZE];
    permGrad4 = new Grad4[PSIZE];
    short[] source = new short[PSIZE];
    for (short i = 0; i < PSIZE; i++)
        source[i] = i;
    for (int i = PSIZE - 1; i >= 0; i--) {
        seed = seed * 6_364_136_223_846_793_005L + 1_442_695_040_888_963_407L;
        int r = cast(int)((seed + 31) % (i + 1));
        if (r < 0)
            r += (i + 1);
        perm[i] = source[r];
        permGrad2[i] = GRADIENTS_2D[perm[i]];
        permGrad3[i] = GRADIENTS_3D[perm[i]];
        permGrad4[i] = GRADIENTS_4D[perm[i]];
        source[r] = source[i];
    }
}

/**
    2D Simplex noise, standard lattice orientation.
*/
double osnoise2(double x, double y) {

    // Get points for A2* lattice
    double s = 0.366025403784439 * (x + y);
    double xs = x + s, ys = y + s;

    return _osnoise2(xs, ys);
}

/**
    2D Simplex noise, with Y pointing down the main diagonal.
    Might be better for a 2D sandbox style game, where Y is vertical.
    Probably slightly less optimal for heightmaps or continent maps.
*/
double osnoise2XY(double x, double y) {

    // Skew transform and rotation baked into one.
    double xx = x * 0.7071067811865476;
    double yy = y * 1.224744871380249;

    return _osnoise2(yy + xx, yy - xx);
}

/**
    3D Re-oriented 4-point BCC noise, classic orientation.
    Proper substitute for 3D Simplex in light of Forbidden Formulae.
    Use noise3_XYBeforeZ or noise3_XZBeforeY instead, wherever appropriate.
*/
double osnoise3(double x, double y, double z) {

    // Re-orient the cubic lattices via rotation, to produce the expected look on cardinal planar slices.
    // If texturing objects that don't tend to have cardinal plane faces, you could even remove this.
    // Orthonormal rotation. Not a skew transform.
    double r = (2.0 / 3.0) * (x + y + z);
    double xr = r - x, yr = r - y, zr = r - z;

    // Evaluate both lattices to form a BCC lattice.
    return _osnoise3(xr, yr, zr);
}

/**
    3D Re-oriented 4-point BCC noise, with better visual isotropy in (X, Y).
    Recommended for 3D terrain and time-varied animations.
    The Z coordinate should always be the "different" coordinate in your use case.
    If Y is vertical in world coordinates, call noise3_XYBeforeZ(x, z, Y) or use noise3_XZBeforeY.
    If Z is vertical in world coordinates, call noise3_XYBeforeZ(x, y, Z).
    For a time varied animation, call noise3_XYBeforeZ(x, y, T).
*/
double osnoise3XYZ(double x, double y, double z) {

    // Re-orient the cubic lattices without skewing, to make X and Y triangular like 2D.
    // Orthonormal rotation. Not a skew transform.
    double xy = x + y;
    double s2 = xy * -0.211324865405187;
    double zz = z * 0.577350269189626;
    double xr = x + s2 - zz, yr = y + s2 - zz;
    double zr = xy * 0.577350269189626 + zz;

    // Evaluate both lattices to form a BCC lattice.
    return _osnoise3(xr, yr, zr);
}

/**
    3D Re-oriented 4-point BCC noise, with better visual isotropy in (X, Z).
    Recommended for 3D terrain and time-varied animations.
    The Y coordinate should always be the "different" coordinate in your use case.
    If Y is vertical in world coordinates, call noise3_XZBeforeY(x, Y, z).
    If Z is vertical in world coordinates, call noise3_XZBeforeY(x, Z, y) or use noise3_XYBeforeZ.
    For a time varied animation, call noise3_XZBeforeY(x, T, y) or use noise3_XYBeforeZ.
*/
double osnoise3XZY(double x, double y, double z) {

    // Re-orient the cubic lattices without skewing, to make X and Z triangular like 2D.
    // Orthonormal rotation. Not a skew transform.
    double xz = x + z;
    double s2 = xz * -0.211324865405187;
    double yy = y * 0.577350269189626;
    double xr = x + s2 - yy; double zr = z + s2 - yy;
    double yr = xz * 0.577350269189626 + yy;

    // Evaluate both lattices to form a BCC lattice.
    return _osnoise3(xr, yr, zr);
}

/**
* 4D OpenSimplex2F noise, classic lattice orientation.
*/
double osnoise4(double x, double y, double z, double w) {

    // Get points for A4 lattice
    double s = -0.138196601125011 * (x + y + z + w);
    double xs = x + s, ys = y + s, zs = z + s, ws = w + s;

    return _osnoise4(xs, ys, zs, ws);
}

/**
    4D OpenSimplex2F noise, with XY and ZW forming orthogonal triangular-based planes.
    Recommended for 3D terrain, where X and Y (or Z and W) are horizontal.
    Recommended for noise(x, y, sin(time), cos(time)) trick.
*/
double osnoise4OrthoXYZW(double x, double y, double z, double w) {

    double s2 = (x + y) * -0.178275657951399372 + (z + w) * 0.215623393288842828;
    double t2 = (z + w) * -0.403949762580207112 + (x + y) * -0.375199083010075342;
    double xs = x + s2, ys = y + s2, zs = z + t2, ws = w + t2;

    return _osnoise4(xs, ys, zs, ws);
}

/**
    4D OpenSimplex2F noise, with XZ and YW forming orthogonal triangular-based planes.
    Recommended for 3D terrain, where X and Z (or Y and W) are horizontal.
*/
double osnoise4OrthoXZYW(double x, double y, double z, double w) {

    double s2 = (x + z) * -0.178275657951399372 + (y + w) * 0.215623393288842828;
    double t2 = (y + w) * -0.403949762580207112 + (x + z) * -0.375199083010075342;
    double xs = x + s2, ys = y + t2, zs = z + s2, ws = w + t2;

    return _osnoise4(xs, ys, zs, ws);
}

/**
    4D OpenSimplex2F noise, with XYZ oriented like noise3_Classic,
    and W for an extra degree of freedom. W repeats eventually.
    Recommended for time-varied animations which texture a 3D object (W=time)
*/
double osnoise4XYZW(double x, double y, double z, double w) {

    double xyz = x + y + z;
    double ww = w * 0.2236067977499788;
    double s2 = xyz * -0.16666666666666666 + ww;
    double xs = x + s2, ys = y + s2, zs = z + s2, ws = -0.5 * xyz + ww;

    return _osnoise4(xs, ys, zs, ws);
}

shared static this() {
    LOOKUP_2D = new LatticePoint2D*[4];
    LOOKUP_3D = new LatticePoint3D*[8];
    VERTICES_4D = new LatticePoint4D*[16];

    LOOKUP_2D[0] = new LatticePoint2D(1, 0);
    LOOKUP_2D[1] = new LatticePoint2D(0, 0);
    LOOKUP_2D[2] = new LatticePoint2D(1, 1);
    LOOKUP_2D[3] = new LatticePoint2D(0, 1);

    for (int i = 0; i < 8; i++) {
        int i1, j1, k1, i2, j2, k2;
        i1 = (i >> 0) & 1; j1 = (i >> 1) & 1; k1 = (i >> 2) & 1;
        i2 = i1 ^ 1; j2 = j1 ^ 1; k2 = k1 ^ 1;

        // The two points within this octant, one from each of the two cubic half-lattices.
        LatticePoint3D* c0 = new LatticePoint3D(i1, j1, k1, 0);
        LatticePoint3D* c1 = new LatticePoint3D(i1 + i2, j1 + j2, k1 + k2, 1);

        // Each single step away on the first half-lattice.
        LatticePoint3D* c2 = new LatticePoint3D(i1 ^ 1, j1, k1, 0);
        LatticePoint3D* c3 = new LatticePoint3D(i1, j1 ^ 1, k1, 0);
        LatticePoint3D* c4 = new LatticePoint3D(i1, j1, k1 ^ 1, 0);

        // Each single step away on the second half-lattice.
        LatticePoint3D* c5 = new LatticePoint3D(i1 + (i2 ^ 1), j1 + j2, k1 + k2, 1);
        LatticePoint3D* c6 = new LatticePoint3D(i1 + i2, j1 + (j2 ^ 1), k1 + k2, 1);
        LatticePoint3D* c7 = new LatticePoint3D(i1 + i2, j1 + j2, k1 + (k2 ^ 1), 1);

        // First two are guaranteed.
        c0.nextOnFailure = c0.nextOnSuccess = c1;
        c1.nextOnFailure = c1.nextOnSuccess = c2;

        // Once we find one on the first half-lattice, the rest are out.
        // In addition, knowing c2 rules out c5.
        c2.nextOnFailure = c3; c2.nextOnSuccess = c6;
        c3.nextOnFailure = c4; c3.nextOnSuccess = c5;
        c4.nextOnFailure = c4.nextOnSuccess = c5;

        // Once we find one on the second half-lattice, the rest are out.
        c5.nextOnFailure = c6; c5.nextOnSuccess = null;
        c6.nextOnFailure = c7; c6.nextOnSuccess = null;
        c7.nextOnFailure = c7.nextOnSuccess = null;

        LOOKUP_3D[i] = c0;
    }

    for (int i = 0; i < 16; i++) {
        VERTICES_4D[i] = new LatticePoint4D((i >> 0) & 1, (i >> 1) & 1, (i >> 2) & 1, (i >> 3) & 1);
    }

    
    GRADIENTS_2D = new Grad2[PSIZE];
    Grad2[] grad2 = [
        Grad2( 0.130526192220052,  0.99144486137381),
        Grad2( 0.38268343236509,   0.923879532511287),
        Grad2( 0.608761429008721,  0.793353340291235),
        Grad2( 0.793353340291235,  0.608761429008721),
        Grad2( 0.923879532511287,  0.38268343236509),
        Grad2( 0.99144486137381,   0.130526192220051),
        Grad2( 0.99144486137381,  -0.130526192220051),
        Grad2( 0.923879532511287, -0.38268343236509),
        Grad2( 0.793353340291235, -0.60876142900872),
        Grad2( 0.608761429008721, -0.793353340291235),
        Grad2( 0.38268343236509,  -0.923879532511287),
        Grad2( 0.130526192220052, -0.99144486137381),
        Grad2(-0.130526192220052, -0.99144486137381),
        Grad2(-0.38268343236509,  -0.923879532511287),
        Grad2(-0.608761429008721, -0.793353340291235),
        Grad2(-0.793353340291235, -0.608761429008721),
        Grad2(-0.923879532511287, -0.38268343236509),
        Grad2(-0.99144486137381,  -0.130526192220052),
        Grad2(-0.99144486137381,   0.130526192220051),
        Grad2(-0.923879532511287,  0.38268343236509),
        Grad2(-0.793353340291235,  0.608761429008721),
        Grad2(-0.608761429008721,  0.793353340291235),
        Grad2(-0.38268343236509,   0.923879532511287),
        Grad2(-0.130526192220052,  0.99144486137381)
    ];

    for (int i = 0; i < grad2.length; i++) {
        grad2[i].dx /= N2; grad2[i].dy /= N2;
    }
    for (int i = 0; i < PSIZE; i++) {
        GRADIENTS_2D[i] = grad2[i % grad2.length];
    }

    GRADIENTS_3D = new Grad3[PSIZE];
    Grad3[] grad3 = [
        Grad3(-2.22474487139,      -2.22474487139,      -1.0),
        Grad3(-2.22474487139,      -2.22474487139,       1.0),
        Grad3(-3.0862664687972017, -1.1721513422464978,  0.0),
        Grad3(-1.1721513422464978, -3.0862664687972017,  0.0),
        Grad3(-2.22474487139,      -1.0,                -2.22474487139),
        Grad3(-2.22474487139,       1.0,                -2.22474487139),
        Grad3(-1.1721513422464978,  0.0,                -3.0862664687972017),
        Grad3(-3.0862664687972017,  0.0,                -1.1721513422464978),
        Grad3(-2.22474487139,      -1.0,                 2.22474487139),
        Grad3(-2.22474487139,       1.0,                 2.22474487139),
        Grad3(-3.0862664687972017,  0.0,                 1.1721513422464978),
        Grad3(-1.1721513422464978,  0.0,                 3.0862664687972017),
        Grad3(-2.22474487139,       2.22474487139,      -1.0),
        Grad3(-2.22474487139,       2.22474487139,       1.0),
        Grad3(-1.1721513422464978,  3.0862664687972017,  0.0),
        Grad3(-3.0862664687972017,  1.1721513422464978,  0.0),
        Grad3(-1.0,                -2.22474487139,      -2.22474487139),
        Grad3( 1.0,                -2.22474487139,      -2.22474487139),
        Grad3( 0.0,                -3.0862664687972017, -1.1721513422464978),
        Grad3( 0.0,                -1.1721513422464978, -3.0862664687972017),
        Grad3(-1.0,                -2.22474487139,       2.22474487139),
        Grad3( 1.0,                -2.22474487139,       2.22474487139),
        Grad3( 0.0,                -1.1721513422464978,  3.0862664687972017),
        Grad3( 0.0,                -3.0862664687972017,  1.1721513422464978),
        Grad3(-1.0,                 2.22474487139,      -2.22474487139),
        Grad3( 1.0,                 2.22474487139,      -2.22474487139),
        Grad3( 0.0,                 1.1721513422464978, -3.0862664687972017),
        Grad3( 0.0,                 3.0862664687972017, -1.1721513422464978),
        Grad3(-1.0,                 2.22474487139,       2.22474487139),
        Grad3( 1.0,                 2.22474487139,       2.22474487139),
        Grad3( 0.0,                 3.0862664687972017,  1.1721513422464978),
        Grad3( 0.0,                 1.1721513422464978,  3.0862664687972017),
        Grad3( 2.22474487139,      -2.22474487139,      -1.0),
        Grad3( 2.22474487139,      -2.22474487139,       1.0),
        Grad3( 1.1721513422464978, -3.0862664687972017,  0.0),
        Grad3( 3.0862664687972017, -1.1721513422464978,  0.0),
        Grad3( 2.22474487139,      -1.0,                -2.22474487139),
        Grad3( 2.22474487139,       1.0,                -2.22474487139),
        Grad3( 3.0862664687972017,  0.0,                -1.1721513422464978),
        Grad3( 1.1721513422464978,  0.0,                -3.0862664687972017),
        Grad3( 2.22474487139,      -1.0,                 2.22474487139),
        Grad3( 2.22474487139,       1.0,                 2.22474487139),
        Grad3( 1.1721513422464978,  0.0,                 3.0862664687972017),
        Grad3( 3.0862664687972017,  0.0,                 1.1721513422464978),
        Grad3( 2.22474487139,       2.22474487139,      -1.0),
        Grad3( 2.22474487139,       2.22474487139,       1.0),
        Grad3( 3.0862664687972017,  1.1721513422464978,  0.0),
        Grad3( 1.1721513422464978,  3.0862664687972017,  0.0)
    ];
    for (int i = 0; i < grad3.length; i++) {
        grad3[i].dx /= N3; grad3[i].dy /= N3; grad3[i].dz /= N3;
    }
    for (int i = 0; i < PSIZE; i++) {
        GRADIENTS_3D[i] = grad3[i % grad3.length];
    }

    GRADIENTS_4D = new Grad4[PSIZE];
    Grad4[] grad4 = [
        Grad4(-0.753341017856078,    -0.37968289875261624,  -0.37968289875261624,  -0.37968289875261624),
        Grad4(-0.7821684431180708,   -0.4321472685365301,   -0.4321472685365301,    0.12128480194602098),
        Grad4(-0.7821684431180708,   -0.4321472685365301,    0.12128480194602098,  -0.4321472685365301),
        Grad4(-0.7821684431180708,    0.12128480194602098,  -0.4321472685365301,   -0.4321472685365301),
        Grad4(-0.8586508742123365,   -0.508629699630796,     0.044802370851755174,  0.044802370851755174),
        Grad4(-0.8586508742123365,    0.044802370851755174, -0.508629699630796,     0.044802370851755174),
        Grad4(-0.8586508742123365,    0.044802370851755174,  0.044802370851755174, -0.508629699630796),
        Grad4(-0.9982828964265062,   -0.03381941603233842,  -0.03381941603233842,  -0.03381941603233842),
        Grad4(-0.37968289875261624,  -0.753341017856078,    -0.37968289875261624,  -0.37968289875261624),
        Grad4(-0.4321472685365301,   -0.7821684431180708,   -0.4321472685365301,    0.12128480194602098),
        Grad4(-0.4321472685365301,   -0.7821684431180708,    0.12128480194602098,  -0.4321472685365301),
        Grad4( 0.12128480194602098,  -0.7821684431180708,   -0.4321472685365301,   -0.4321472685365301),
        Grad4(-0.508629699630796,    -0.8586508742123365,    0.044802370851755174,  0.044802370851755174),
        Grad4( 0.044802370851755174, -0.8586508742123365,   -0.508629699630796,     0.044802370851755174),
        Grad4( 0.044802370851755174, -0.8586508742123365,    0.044802370851755174, -0.508629699630796),
        Grad4(-0.03381941603233842,  -0.9982828964265062,   -0.03381941603233842,  -0.03381941603233842),
        Grad4(-0.37968289875261624,  -0.37968289875261624,  -0.753341017856078,    -0.37968289875261624),
        Grad4(-0.4321472685365301,   -0.4321472685365301,   -0.7821684431180708,    0.12128480194602098),
        Grad4(-0.4321472685365301,    0.12128480194602098,  -0.7821684431180708,   -0.4321472685365301),
        Grad4( 0.12128480194602098,  -0.4321472685365301,   -0.7821684431180708,   -0.4321472685365301),
        Grad4(-0.508629699630796,     0.044802370851755174, -0.8586508742123365,    0.044802370851755174),
        Grad4( 0.044802370851755174, -0.508629699630796,    -0.8586508742123365,    0.044802370851755174),
        Grad4( 0.044802370851755174,  0.044802370851755174, -0.8586508742123365,   -0.508629699630796),
        Grad4(-0.03381941603233842,  -0.03381941603233842,  -0.9982828964265062,   -0.03381941603233842),
        Grad4(-0.37968289875261624,  -0.37968289875261624,  -0.37968289875261624,  -0.753341017856078),
        Grad4(-0.4321472685365301,   -0.4321472685365301,    0.12128480194602098,  -0.7821684431180708),
        Grad4(-0.4321472685365301,    0.12128480194602098,  -0.4321472685365301,   -0.7821684431180708),
        Grad4( 0.12128480194602098,  -0.4321472685365301,   -0.4321472685365301,   -0.7821684431180708),
        Grad4(-0.508629699630796,     0.044802370851755174,  0.044802370851755174, -0.8586508742123365),
        Grad4( 0.044802370851755174, -0.508629699630796,     0.044802370851755174, -0.8586508742123365),
        Grad4( 0.044802370851755174,  0.044802370851755174, -0.508629699630796,    -0.8586508742123365),
        Grad4(-0.03381941603233842,  -0.03381941603233842,  -0.03381941603233842,  -0.9982828964265062),
        Grad4(-0.6740059517812944,   -0.3239847771997537,   -0.3239847771997537,    0.5794684678643381),
        Grad4(-0.7504883828755602,   -0.4004672082940195,    0.15296486218853164,   0.5029860367700724),
        Grad4(-0.7504883828755602,    0.15296486218853164,  -0.4004672082940195,    0.5029860367700724),
        Grad4(-0.8828161875373585,    0.08164729285680945,   0.08164729285680945,   0.4553054119602712),
        Grad4(-0.4553054119602712,   -0.08164729285680945,  -0.08164729285680945,   0.8828161875373585),
        Grad4(-0.5029860367700724,   -0.15296486218853164,   0.4004672082940195,    0.7504883828755602),
        Grad4(-0.5029860367700724,    0.4004672082940195,   -0.15296486218853164,   0.7504883828755602),
        Grad4(-0.5794684678643381,    0.3239847771997537,    0.3239847771997537,    0.6740059517812944),
        Grad4(-0.3239847771997537,   -0.6740059517812944,   -0.3239847771997537,    0.5794684678643381),
        Grad4(-0.4004672082940195,   -0.7504883828755602,    0.15296486218853164,   0.5029860367700724),
        Grad4( 0.15296486218853164,  -0.7504883828755602,   -0.4004672082940195,    0.5029860367700724),
        Grad4( 0.08164729285680945,  -0.8828161875373585,    0.08164729285680945,   0.4553054119602712),
        Grad4(-0.08164729285680945,  -0.4553054119602712,   -0.08164729285680945,   0.8828161875373585),
        Grad4(-0.15296486218853164,  -0.5029860367700724,    0.4004672082940195,    0.7504883828755602),
        Grad4( 0.4004672082940195,   -0.5029860367700724,   -0.15296486218853164,   0.7504883828755602),
        Grad4( 0.3239847771997537,   -0.5794684678643381,    0.3239847771997537,    0.6740059517812944),
        Grad4(-0.3239847771997537,   -0.3239847771997537,   -0.6740059517812944,    0.5794684678643381),
        Grad4(-0.4004672082940195,    0.15296486218853164,  -0.7504883828755602,    0.5029860367700724),
        Grad4( 0.15296486218853164,  -0.4004672082940195,   -0.7504883828755602,    0.5029860367700724),
        Grad4( 0.08164729285680945,   0.08164729285680945,  -0.8828161875373585,    0.4553054119602712),
        Grad4(-0.08164729285680945,  -0.08164729285680945,  -0.4553054119602712,    0.8828161875373585),
        Grad4(-0.15296486218853164,   0.4004672082940195,   -0.5029860367700724,    0.7504883828755602),
        Grad4( 0.4004672082940195,   -0.15296486218853164,  -0.5029860367700724,    0.7504883828755602),
        Grad4( 0.3239847771997537,    0.3239847771997537,   -0.5794684678643381,    0.6740059517812944),
        Grad4(-0.6740059517812944,   -0.3239847771997537,    0.5794684678643381,   -0.3239847771997537),
        Grad4(-0.7504883828755602,   -0.4004672082940195,    0.5029860367700724,    0.15296486218853164),
        Grad4(-0.7504883828755602,    0.15296486218853164,   0.5029860367700724,   -0.4004672082940195),
        Grad4(-0.8828161875373585,    0.08164729285680945,   0.4553054119602712,    0.08164729285680945),
        Grad4(-0.4553054119602712,   -0.08164729285680945,   0.8828161875373585,   -0.08164729285680945),
        Grad4(-0.5029860367700724,   -0.15296486218853164,   0.7504883828755602,    0.4004672082940195),
        Grad4(-0.5029860367700724,    0.4004672082940195,    0.7504883828755602,   -0.15296486218853164),
        Grad4(-0.5794684678643381,    0.3239847771997537,    0.6740059517812944,    0.3239847771997537),
        Grad4(-0.3239847771997537,   -0.6740059517812944,    0.5794684678643381,   -0.3239847771997537),
        Grad4(-0.4004672082940195,   -0.7504883828755602,    0.5029860367700724,    0.15296486218853164),
        Grad4( 0.15296486218853164,  -0.7504883828755602,    0.5029860367700724,   -0.4004672082940195),
        Grad4( 0.08164729285680945,  -0.8828161875373585,    0.4553054119602712,    0.08164729285680945),
        Grad4(-0.08164729285680945,  -0.4553054119602712,    0.8828161875373585,   -0.08164729285680945),
        Grad4(-0.15296486218853164,  -0.5029860367700724,    0.7504883828755602,    0.4004672082940195),
        Grad4( 0.4004672082940195,   -0.5029860367700724,    0.7504883828755602,   -0.15296486218853164),
        Grad4( 0.3239847771997537,   -0.5794684678643381,    0.6740059517812944,    0.3239847771997537),
        Grad4(-0.3239847771997537,   -0.3239847771997537,    0.5794684678643381,   -0.6740059517812944),
        Grad4(-0.4004672082940195,    0.15296486218853164,   0.5029860367700724,   -0.7504883828755602),
        Grad4( 0.15296486218853164,  -0.4004672082940195,    0.5029860367700724,   -0.7504883828755602),
        Grad4( 0.08164729285680945,   0.08164729285680945,   0.4553054119602712,   -0.8828161875373585),
        Grad4(-0.08164729285680945,  -0.08164729285680945,   0.8828161875373585,   -0.4553054119602712),
        Grad4(-0.15296486218853164,   0.4004672082940195,    0.7504883828755602,   -0.5029860367700724),
        Grad4( 0.4004672082940195,   -0.15296486218853164,   0.7504883828755602,   -0.5029860367700724),
        Grad4( 0.3239847771997537,    0.3239847771997537,    0.6740059517812944,   -0.5794684678643381),
        Grad4(-0.6740059517812944,    0.5794684678643381,   -0.3239847771997537,   -0.3239847771997537),
        Grad4(-0.7504883828755602,    0.5029860367700724,   -0.4004672082940195,    0.15296486218853164),
        Grad4(-0.7504883828755602,    0.5029860367700724,    0.15296486218853164,  -0.4004672082940195),
        Grad4(-0.8828161875373585,    0.4553054119602712,    0.08164729285680945,   0.08164729285680945),
        Grad4(-0.4553054119602712,    0.8828161875373585,   -0.08164729285680945,  -0.08164729285680945),
        Grad4(-0.5029860367700724,    0.7504883828755602,   -0.15296486218853164,   0.4004672082940195),
        Grad4(-0.5029860367700724,    0.7504883828755602,    0.4004672082940195,   -0.15296486218853164),
        Grad4(-0.5794684678643381,    0.6740059517812944,    0.3239847771997537,    0.3239847771997537),
        Grad4(-0.3239847771997537,    0.5794684678643381,   -0.6740059517812944,   -0.3239847771997537),
        Grad4(-0.4004672082940195,    0.5029860367700724,   -0.7504883828755602,    0.15296486218853164),
        Grad4( 0.15296486218853164,   0.5029860367700724,   -0.7504883828755602,   -0.4004672082940195),
        Grad4( 0.08164729285680945,   0.4553054119602712,   -0.8828161875373585,    0.08164729285680945),
        Grad4(-0.08164729285680945,   0.8828161875373585,   -0.4553054119602712,   -0.08164729285680945),
        Grad4(-0.15296486218853164,   0.7504883828755602,   -0.5029860367700724,    0.4004672082940195),
        Grad4( 0.4004672082940195,    0.7504883828755602,   -0.5029860367700724,   -0.15296486218853164),
        Grad4( 0.3239847771997537,    0.6740059517812944,   -0.5794684678643381,    0.3239847771997537),
        Grad4(-0.3239847771997537,    0.5794684678643381,   -0.3239847771997537,   -0.6740059517812944),
        Grad4(-0.4004672082940195,    0.5029860367700724,    0.15296486218853164,  -0.7504883828755602),
        Grad4( 0.15296486218853164,   0.5029860367700724,   -0.4004672082940195,   -0.7504883828755602),
        Grad4( 0.08164729285680945,   0.4553054119602712,    0.08164729285680945,  -0.8828161875373585),
        Grad4(-0.08164729285680945,   0.8828161875373585,   -0.08164729285680945,  -0.4553054119602712),
        Grad4(-0.15296486218853164,   0.7504883828755602,    0.4004672082940195,   -0.5029860367700724),
        Grad4( 0.4004672082940195,    0.7504883828755602,   -0.15296486218853164,  -0.5029860367700724),
        Grad4( 0.3239847771997537,    0.6740059517812944,    0.3239847771997537,   -0.5794684678643381),
        Grad4( 0.5794684678643381,   -0.6740059517812944,   -0.3239847771997537,   -0.3239847771997537),
        Grad4( 0.5029860367700724,   -0.7504883828755602,   -0.4004672082940195,    0.15296486218853164),
        Grad4( 0.5029860367700724,   -0.7504883828755602,    0.15296486218853164,  -0.4004672082940195),
        Grad4( 0.4553054119602712,   -0.8828161875373585,    0.08164729285680945,   0.08164729285680945),
        Grad4( 0.8828161875373585,   -0.4553054119602712,   -0.08164729285680945,  -0.08164729285680945),
        Grad4( 0.7504883828755602,   -0.5029860367700724,   -0.15296486218853164,   0.4004672082940195),
        Grad4( 0.7504883828755602,   -0.5029860367700724,    0.4004672082940195,   -0.15296486218853164),
        Grad4( 0.6740059517812944,   -0.5794684678643381,    0.3239847771997537,    0.3239847771997537),
        Grad4( 0.5794684678643381,   -0.3239847771997537,   -0.6740059517812944,   -0.3239847771997537),
        Grad4( 0.5029860367700724,   -0.4004672082940195,   -0.7504883828755602,    0.15296486218853164),
        Grad4( 0.5029860367700724,    0.15296486218853164,  -0.7504883828755602,   -0.4004672082940195),
        Grad4( 0.4553054119602712,    0.08164729285680945,  -0.8828161875373585,    0.08164729285680945),
        Grad4( 0.8828161875373585,   -0.08164729285680945,  -0.4553054119602712,   -0.08164729285680945),
        Grad4( 0.7504883828755602,   -0.15296486218853164,  -0.5029860367700724,    0.4004672082940195),
        Grad4( 0.7504883828755602,    0.4004672082940195,   -0.5029860367700724,   -0.15296486218853164),
        Grad4( 0.6740059517812944,    0.3239847771997537,   -0.5794684678643381,    0.3239847771997537),
        Grad4( 0.5794684678643381,   -0.3239847771997537,   -0.3239847771997537,   -0.6740059517812944),
        Grad4( 0.5029860367700724,   -0.4004672082940195,    0.15296486218853164,  -0.7504883828755602),
        Grad4( 0.5029860367700724,    0.15296486218853164,  -0.4004672082940195,   -0.7504883828755602),
        Grad4( 0.4553054119602712,    0.08164729285680945,   0.08164729285680945,  -0.8828161875373585),
        Grad4( 0.8828161875373585,   -0.08164729285680945,  -0.08164729285680945,  -0.4553054119602712),
        Grad4( 0.7504883828755602,   -0.15296486218853164,   0.4004672082940195,   -0.5029860367700724),
        Grad4( 0.7504883828755602,    0.4004672082940195,   -0.15296486218853164,  -0.5029860367700724),
        Grad4( 0.6740059517812944,    0.3239847771997537,    0.3239847771997537,   -0.5794684678643381),
        Grad4( 0.03381941603233842,   0.03381941603233842,   0.03381941603233842,   0.9982828964265062),
        Grad4(-0.044802370851755174, -0.044802370851755174,  0.508629699630796,     0.8586508742123365),
        Grad4(-0.044802370851755174,  0.508629699630796,    -0.044802370851755174,  0.8586508742123365),
        Grad4(-0.12128480194602098,   0.4321472685365301,    0.4321472685365301,    0.7821684431180708),
        Grad4( 0.508629699630796,    -0.044802370851755174, -0.044802370851755174,  0.8586508742123365),
        Grad4( 0.4321472685365301,   -0.12128480194602098,   0.4321472685365301,    0.7821684431180708),
        Grad4( 0.4321472685365301,    0.4321472685365301,   -0.12128480194602098,   0.7821684431180708),
        Grad4( 0.37968289875261624,   0.37968289875261624,   0.37968289875261624,   0.753341017856078),
        Grad4( 0.03381941603233842,   0.03381941603233842,   0.9982828964265062,    0.03381941603233842),
        Grad4(-0.044802370851755174,  0.044802370851755174,  0.8586508742123365,    0.508629699630796),
        Grad4(-0.044802370851755174,  0.508629699630796,     0.8586508742123365,   -0.044802370851755174),
        Grad4(-0.12128480194602098,   0.4321472685365301,    0.7821684431180708,    0.4321472685365301),
        Grad4( 0.508629699630796,    -0.044802370851755174,  0.8586508742123365,   -0.044802370851755174),
        Grad4( 0.4321472685365301,   -0.12128480194602098,   0.7821684431180708,    0.4321472685365301),
        Grad4( 0.4321472685365301,    0.4321472685365301,    0.7821684431180708,   -0.12128480194602098),
        Grad4( 0.37968289875261624,   0.37968289875261624,   0.753341017856078,     0.37968289875261624),
        Grad4( 0.03381941603233842,   0.9982828964265062,    0.03381941603233842,   0.03381941603233842),
        Grad4(-0.044802370851755174,  0.8586508742123365,   -0.044802370851755174,  0.508629699630796),
        Grad4(-0.044802370851755174,  0.8586508742123365,    0.508629699630796,    -0.044802370851755174),
        Grad4(-0.12128480194602098,   0.7821684431180708,    0.4321472685365301,    0.4321472685365301),
        Grad4( 0.508629699630796,     0.8586508742123365,   -0.044802370851755174, -0.044802370851755174),
        Grad4( 0.4321472685365301,    0.7821684431180708,   -0.12128480194602098,   0.4321472685365301),
        Grad4( 0.4321472685365301,    0.7821684431180708,    0.4321472685365301,   -0.12128480194602098),
        Grad4( 0.37968289875261624,   0.753341017856078,     0.37968289875261624,   0.37968289875261624),
        Grad4( 0.9982828964265062,    0.03381941603233842,   0.03381941603233842,   0.03381941603233842),
        Grad4( 0.8586508742123365,   -0.044802370851755174, -0.044802370851755174,  0.508629699630796),
        Grad4( 0.8586508742123365,   -0.044802370851755174,  0.508629699630796,    -0.044802370851755174),
        Grad4( 0.7821684431180708,   -0.12128480194602098,   0.4321472685365301,    0.4321472685365301),
        Grad4( 0.8586508742123365,    0.508629699630796,    -0.044802370851755174, -0.044802370851755174),
        Grad4( 0.7821684431180708,    0.4321472685365301,   -0.12128480194602098,   0.4321472685365301),
        Grad4( 0.7821684431180708,    0.4321472685365301,    0.4321472685365301,   -0.12128480194602098),
        Grad4( 0.753341017856078,     0.37968289875261624,   0.37968289875261624,   0.37968289875261624)
    ];

    for (int i = 0; i < grad4.length; i++) {
        grad4[i].dx /= N4; grad4[i].dy /= N4; grad4[i].dz /= N4; grad4[i].dw /= N4;
    }
    for (int i = 0; i < PSIZE; i++) {
        GRADIENTS_4D[i] = grad4[i % grad4.length];
    }
}