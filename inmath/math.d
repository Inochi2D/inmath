/**
    inmath.math

    Provides nearly all GLSL functions, according to spec 4.1,
    it also publically imports other useful functions (from std.math, core.stdc.math, std.alogrithm) 
    so you only have to import this file to get all mathematical functions you need.

    Publically imports: PI, sin, cos, tan, asin, acos, atan, atan2, sinh, cosh, tanh, 
    asinh, acosh, atanh, pow, exp, log, exp2, log2, sqrt, abs, floor, trunc, round, ceil, modf,
    fmodf, min, max.

    Authors: David Herberth, Inochi2D Project
    License: MIT
*/

module inmath.math;

public {
    import std.math : PI, sin, cos, tan, asin, acos, atan, atan2,
                      sinh, cosh, tanh, asinh, acosh, atanh,
                      pow, exp, log, exp2, log2, sqrt,
                      floor, trunc, round, ceil, modf;
    alias roundEven = round;
    alias fract = floor;
    //import core.stdc.math : fmodf;
    import std.algorithm : min, max;
}

private {
    import std.conv : to;
    import std.algorithm : all;
    import std.range : zip;
    import std.traits : CommonType, Unqual;
    import std.range : ElementType;
    import smath = std.math;
    
    import inmath.util : isVector, isQuaternion, isMatrix;

    version(unittest) {
        import inmath.linalg : vec2, vec2i, vec3, vec3i, quat;
    }
}

@nogc nothrow:

/// PI / 180 at compiletime, used for degrees/radians conversion.
public enum real PI_180 = PI / 180;

/// 180 / PI at compiletime, used for degrees/radians conversion.
public enum real _180_PI = 180 / PI;

/// Modulus. Returns x - y * floor(x/y).
T mod(T)(T x, T y) { // std.math.floor is not pure
    return x - y * floor(x/y);
}

@safe pure:

extern (C) { float fmodf(float x, float y); }

/// Calculates the absolute value.
T abs(T)(T t) if(!isVector!T && !isQuaternion!T && !isMatrix!T) {
    return smath.abs(t);
}

/// Calculates the absolute value per component.
T abs(T)(T vec) if(isVector!T) {
    Unqual!T ret;

    foreach(i, element; vec.vector) {
        ret.vector[i] = abs(element);
    }
    
    return ret;
}

/// ditto
T abs(T)(T quat) if(isQuaternion!T) {
    Unqual!T ret;

    ret.quaternion[0] = abs(quat.quaternion[0]);
    ret.quaternion[1] = abs(quat.quaternion[1]);
    ret.quaternion[2] = abs(quat.quaternion[2]);
    ret.quaternion[3] = abs(quat.quaternion[3]);

    return ret;
}

unittest {
    assert(abs(0) == 0);
    assert(abs(-1) == 1);
    assert(abs(1) == 1);
    assert(abs(0.0) == 0.0);
    assert(abs(-1.0) == 1.0);
    assert(abs(1.0) == 1.0);
    
    assert(abs(vec3i(-1, 0, -12)) == vec3(1, 0, 12));
    assert(abs(vec3(-1, 0, -12)) == vec3(1, 0, 12));
    assert(abs(vec3i(12, 12, 12)) == vec3(12, 12, 12));

    assert(abs(quat(-1.0f, 0.0f, 1.0f, -12.0f)) == quat(1.0f, 0.0f, 1.0f, 12.0f));
}

/// Returns 1/sqrt(x), results are undefined if x <= 0.
real inversesqrt(real x) {
    return 1 / sqrt(x);
}

/// Returns 1.0 if x > 0, 0.0 if x = 0, or -1.0 if x < 0.
float sign(T)(T x) {
    if(x > 0) {
        return 1.0f;
    } else if(x == 0) {
        return 0.0f;
    } else { // if x < 0
        return -1.0f;
    }
}

unittest {
    assert(almostEqual(inversesqrt(1.0f), 1.0f));
    assert(almostEqual(inversesqrt(10.0f), (1/sqrt(10.0f))));
    assert(almostEqual(inversesqrt(2_342_342.0f), (1/sqrt(2_342_342.0f))));
    
    assert(sign(-1) == -1.0f);
    assert(sign(0) == 0.0f);
    assert(sign(1) == 1.0f);
    assert(sign(0.5) == 1.0f);
    assert(sign(-0.5) == -1.0f);
    
    assert(mod(12.0, 27.5) == 12.0);
    assert(mod(-12.0, 27.5) == 15.5);
    assert(mod(12.0, -27.5) == -15.5);
}

/// Compares to values and returns true if the difference is epsilon or smaller.
bool almostEqual(T, S)(T a, S b, float epsilon = 0.000001f) if(!isVector!T && !isQuaternion!T) {
    if(abs(a-b) <= epsilon) {
        return true;
    }
    return abs(a-b) <= epsilon * abs(b);
}

/// ditto
bool almostEqual(T, S)(T a, S b, float epsilon = 0.000001f) if(isVector!T && isVector!S && T.dimension == S.dimension) {
    foreach(i; 0..T.dimension) {
        if(!almostEqual(a.vector[i], b.vector[i], epsilon)) {
            return false;
        }
    }
    return true;
}

bool almostEqual(T)(T a, T b, float epsilon = 0.000001f) if(isQuaternion!T) {
    foreach(i; 0..4) {
        if(!almostEqual(a.quaternion[i], b.quaternion[i], epsilon)) {
            return false;
        }
    }
    return true;
}

unittest {
    assert(almostEqual(0, 0));
    assert(almostEqual(1, 1));
    assert(almostEqual(-1, -1));    
    assert(almostEqual(0f, 0.000001f, 0.000001f));
    assert(almostEqual(1f, 1.1f, 0.1f));
    assert(!almostEqual(1f, 1.1f, 0.01f));

    assert(almostEqual(vec2i(0, 0), vec2(0.0f, 0.0f)));
    assert(almostEqual(vec2(0.0f, 0.0f), vec2(0.000001f, 0.000001f)));
    assert(almostEqual(vec3(0.0f, 1.0f, 2.0f), vec3i(0, 1, 2)));

    assert(almostEqual(quat(0.0f, 0.0f, 0.0f, 0.0f), quat(0.0f, 0.0f, 0.0f, 0.0f)));
    assert(almostEqual(quat(0.0f, 0.0f, 0.0f, 0.0f), quat(0.000001f, 0.000001f, 0.000001f, 0.000001f)));
}

/// Converts degrees to radians.
real radians(real degrees) {
    return PI_180 * degrees;
}

/// Compiletime version of $(I radians).
real cradians(real degrees)() {
    return radians(degrees);
}

/// Converts radians to degrees.
real degrees(real radians) {
    return _180_PI * radians;
}

/// Compiletime version of $(I degrees).
real cdegrees(real radians)() {
    return degrees(radians);
}

unittest {
    assert(almostEqual(radians(to!(real)(0)), 0));
    assert(almostEqual(radians(to!(real)(90)), PI/2));
    assert(almostEqual(radians(to!(real)(180)), PI));
    assert(almostEqual(radians(to!(real)(360)), 2*PI));
    
    assert(almostEqual(degrees(to!(real)(0)), 0));
    assert(almostEqual(degrees(to!(real)(PI/2)), 90));
    assert(almostEqual(degrees(to!(real)(PI)), 180));
    assert(almostEqual(degrees(to!(real)(2*PI)), 360));

    assert(almostEqual(degrees(radians(to!(real)(12))), 12));
    assert(almostEqual(degrees(radians(to!(real)(100))), 100));
    assert(almostEqual(degrees(radians(to!(real)(213))), 213));
    assert(almostEqual(degrees(radians(to!(real)(399))), 399));
    
    /+static+/ assert(almostEqual(cdegrees!PI, 180));
    /+static+/ assert(almostEqual(cradians!180, PI));
}

/// Returns min(max(x, min_val), max_val), Results are undefined if min_val > max_val.
CommonType!(T1, T2, T3) clamp(T1, T2, T3)(T1 x, T2 min_val, T3 max_val) {
    return min(max(x, min_val), max_val);
}

unittest {
    assert(clamp(-1, 0, 2) == 0);
    assert(clamp(0, 0, 2) == 0);
    assert(clamp(1, 0, 2) == 1);
    assert(clamp(2, 0, 2) == 2);
    assert(clamp(3, 0, 2) == 2);
}

/// Returns 0.0 if x < edge, otherwise it returns 1.0.
float step(T1, T2)(T1 edge, T2 x) {
    return x < edge ? 0.0f:1.0f;
}

/// Returns 0.0 if x <= edge0 and 1.0 if x >= edge1 and performs smooth 
/// hermite interpolation between 0 and 1 when edge0 < x < edge1. 
/// This is useful in cases where you would want a threshold function with a smooth transition.
CommonType!(T1, T2, T3) smoothstep(T1, T2, T3)(T1 edge0, T2 edge1, T3 x) {
    auto t = clamp((x - edge0) / (edge1 - edge0), 0, 1);
    return t * t * (3 - 2 * t);
}

unittest {
    assert(step(0, 1) == 1.0f);
    assert(step(0, 10) == 1.0f);
    assert(step(1, 0) == 0.0f);
    assert(step(10, 0) == 0.0f);
    assert(step(1, 1) == 1.0f);
    
    assert(smoothstep(1, 0, 2) == 0);
    assert(smoothstep(1.0, 0.0, 2.0) == 0);
    assert(smoothstep(1.0, 0.0, 0.5) == 0.5);
    assert(almostEqual(smoothstep(0.0, 2.0, 0.5), 0.15625, 0.00001));
}
