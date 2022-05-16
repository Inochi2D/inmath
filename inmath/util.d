/**
    inmath.util

    Authors: David Herberth, Inochi2D Project
    License: MIT
*/
module inmath.util;

private {
    import inmath.linalg : Vector, Matrix, Quaternion;
    import inmath.plane : PlaneT;

    static import std.compiler;

}

static if (std.compiler.version_major > 2 ||
           std.compiler.version_minor > 68)
{
    private import std.meta : AliasSeq;
    public alias TypeTuple = AliasSeq;
} else {
    public import std.typetuple : TypeTuple;
}

private void isVector_impl(T, int d)(Vector!(T, d) vec) {}

/// If T is a vector, this evaluates to true, otherwise false.
template isVector(T) {
    enum isVector = is(typeof(isVector_impl(T.init)));
}

private void isMatrix_impl(T, int r, int c)(Matrix!(T, r, c) mat) {}

/// If T is a matrix, this evaluates to true, otherwise false.
template isMatrix(T) {
    enum isMatrix = is(typeof(isMatrix_impl(T.init)));
}

private void isQuaternion_impl(T)(Quaternion!(T) qu) {}

/// If T is a quaternion, this evaluates to true, otherwise false.
template isQuaternion(T) {
    enum isQuaternion = is(typeof(isQuaternion_impl(T.init)));
}

private void isPlane_impl(T)(PlaneT!(T) p) {}

/// If T is a plane, this evaluates to true, otherwise false.
template isPlane(T) {
    enum isPlane = is(typeof(isPlane_impl(T.init)));
}


unittest {
    // I need to import it here like this, otherwise you'll get a compiler
    // or a linker error depending where inmath.util gets imported
    import inmath.linalg;
    import inmath.plane;

    assert(isVector!vec2);
    assert(isVector!vec3);
    assert(isVector!vec3d);
    assert(isVector!vec4i);
    assert(!isVector!int);
    assert(!isVector!mat34);
    assert(!isVector!quat);

    assert(isMatrix!mat2);
    assert(isMatrix!mat34);
    assert(isMatrix!mat4);
    assert(!isMatrix!float);
    assert(!isMatrix!vec3);
    assert(!isMatrix!quat);

    assert(isQuaternion!quat);
    assert(!isQuaternion!vec2);
    assert(!isQuaternion!vec4i);
    assert(!isQuaternion!mat2);
    assert(!isQuaternion!mat34);
    assert(!isQuaternion!float);

    assert(isPlane!Plane);
    assert(!isPlane!vec2);
    assert(!isPlane!quat);
    assert(!isPlane!mat4);
    assert(!isPlane!float);
}

template TupleRange(int from, int to) if (from <= to) {
    static if (from >= to) {
        alias TupleRange = TypeTuple!();
    } else {
        alias TupleRange = TypeTuple!(from, TupleRange!(from + 1, to));
    }
}
