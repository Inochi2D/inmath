/**
    inmath.lerpolate

    Authors: David Herberth, Inochi2D Project
    License: MIT
*/
module inmath.interpolate;

private {
    import inmath.linalg : Vector, dot, vec2, vec3, vec4, quat;
    import inmath.util : isVector, isQuaternion;
    import inmath.math : almostEqual, acos, sin, sqrt, clamp, PI, exp, log;
    import std.conv : to;
}

@safe pure nothrow @nogc:

/// Linearly interpolates between `a` and `b`, using scalar `t` from 0..1 
T lerp(T)(T a, T b, float t) {
    return a * (1 - t) + b * t;
}

/// Logarithmically interpolates between `a` and `b`, using scalar `t` from 0..1
/// Also known as eerp
T lorp(T)(T a, T b, float t) {
    return a * exp( t * log(b / a));
}

/// Linearly interpolates between to vectors or quaternions.
T slerp(T)(T a, T b, float t) if(isVector!T || isQuaternion!T) {
    static if(isVector!T) {
        real theta = acos(dot(a, b));
    } else {
        real theta = acos(
            // this is a workaround, acos returning -nan on certain values near +/-1
            clamp(a.w * b.w + a.x * b.x + a.y * b.y + a.z * b.z, -1, 1)
        );
    }
    
    if(almostEqual(theta, 0)) {
        return a;
    } else if(almostEqual(theta, PI)) { // 180°?
        return lerp(a, b, t);
    } else { // slerp
        real sintheta = sin(theta);
        return (sin((1.0-t)*theta)/sintheta)*a + (sin(t*theta)/sintheta)*b;
    }
}


/// Normalized quaternion linear interpolation.
quat nlerp(quat a, quat b, float t) {
    // TODO: tests
    float dot = a.w * b.w + a.x * b.x + a.y * b.y + a.z * b.z;

    quat result;
    if(dot < 0) { // Determine the "shortest route"...
        result = a - (b + a) * t; // use -b instead of b
    } else {
        result = a + (b - a) * t;
    }
    result.normalize();

    return result;
}

@("vector lerp, slerp")
unittest {
    vec2 v2_1 = vec2(1.0f);
    vec2 v2_2 = vec2(0.0f);
    vec3 v3_1 = vec3(1.0f);
    vec3 v3_2 = vec3(0.0f);
    vec4 v4_1 = vec4(1.0f);
    vec4 v4_2 = vec4(0.0f);
    
    assert(lerp(v2_1, v2_2, 0.5f).vector == [0.5f, 0.5f]);
    assert(lerp(v2_1, v2_2, 0.0f) == v2_1);
    assert(lerp(v2_1, v2_2, 1.0f) == v2_2);
    assert(lerp(v3_1, v3_2, 0.5f).vector == [0.5f, 0.5f, 0.5f]);
    assert(lerp(v3_1, v3_2, 0.0f) == v3_1);
    assert(lerp(v3_1, v3_2, 1.0f) == v3_2);
    assert(lerp(v4_1, v4_2, 0.5f).vector == [0.5f, 0.5f, 0.5f, 0.5f]);
    assert(lerp(v4_1, v4_2, 0.0f) == v4_1);
    assert(lerp(v4_1, v4_2, 1.0f) == v4_2);
    
    assert(slerp(v2_1, v2_2, 0.0).vector == v2_1.vector);
    assert(slerp(v2_1, v2_2, 1.0).vector == v2_2.vector);
    assert(slerp(v3_1, v3_2, 0.0).vector == v3_1.vector);
    assert(slerp(v3_1, v3_2, 1.0).vector == v3_2.vector);
    assert(slerp(v4_1, v4_2, 0.0).vector == v4_1.vector);
    assert(slerp(v4_1, v4_2, 1.0).vector == v4_2.vector);
}

@("real lerp")
unittest {
    real r1 = 0.0;
    real r2 = 1.0;
    assert(lerp(r1, r2, 0.5f) == 0.5);
    assert(lerp(r1, r2, 0.0f) == r1);
    assert(lerp(r1, r2, 1.0f) == r2);
    
    assert(lerp(0.0, 1.0, 0.5f) == 0.5);
    assert(lerp(0.0, 1.0, 0.0f) == 0.0);
    assert(lerp(0.0, 1.0, 1.0f) == 1.0);
    
    assert(lerp(0.0f, 1.0f, 0.5f) == 0.5f);
    assert(lerp(0.0f, 1.0f, 0.0f) == 0.0f);
    assert(lerp(0.0f, 1.0f, 1.0f) == 1.0f);
}

@("quat lerp")
unittest {
    quat q1 = quat(1.0f, 1.0f, 1.0f, 1.0f);
    quat q2 = quat(0.0f, 0.0f, 0.0f, 0.0f);

    assert(lerp(q1, q2, 0.0f).quaternion == q1.quaternion);
    assert(lerp(q1, q2, 0.5f).quaternion == [0.5f, 0.5f, 0.5f, 0.5f]);
    assert(lerp(q1, q2, 1.0f).quaternion == q2.quaternion);
    
    assert(slerp(q1, q2, 0.0f).quaternion == q1.quaternion);
    assert(slerp(q1, q2, 1.0f).quaternion == q2.quaternion);
}

/// Nearest interpolation of two points.
T nearest(T)(T x, T y, float t) {
    if(t < 0.5f) { return x; }
    else { return y; } 
}

@("nearest lerp")
unittest {
    assert(nearest(0.0, 1.0, 0.5f) == 1.0);
    assert(nearest(0.0, 1.0, 0.4f) == 0.0);
    assert(nearest(0.0, 1.0, 0.6f) == 1.0);
}

/// Catmull-rom interpolation between four points.
T catmullrom(T)(T p0, T p1, T p2, T p3, float t) {
    return 0.5f * ((2 * p1) + 
                   (-p0 + p2) * t +
                   (2 * p0 - 5 * p1 + 4 * p2 - p3) * t^^2 +
                   (-p0 + 3 * p1 - 3 * p2 + p3) * t^^3);
}

/// Catmull-derivatives of the interpolation between four points.
T dcatmullrom(T)(T p0, T p1, T p2, T p3, float t) {
    return 0.5f * ((2 * p1) +
                   (-p0 + p2) +
                   2 * (2 * p0 - 5 * p1 + 4 * p2 - p3) * t +
                   3 * (-p0 + 3 * p1 - 3 * p2 + p3) * t^^2);
}

/// Hermite interpolation (cubic hermite spline).
T hermite(T)(T x, T tx, T y, T ty, float t) {
    float h1 = 2 * t^^3 - 3 * t^^2 + 1;
    float h2 = -2* t^^3 + 3 * t^^2;
    float h3 = t^^3 - 2 * t^^2 + t;
    float h4 = t^^3 - t^^2;
    return h1 * x + h3 * tx + h2 * y + h4 * ty;
}

/// Cubic interpolation
T cubic(T)(T p0, T p1, T p2, T p3, float t) {
    T a = -0.5 * p0 + 1.5 * p1 - 1.5 * p2 + 0.5 * p3;
    T b = p0 - 2.5 * p1 + 2 * p2 - 0.5 * p3;
    T c = -0.5 * p0 + 0.5 * p2;
    T d = p1;
    
    return a * (t ^^ 3) + b * (t ^^ 2) + c * t + d;
}