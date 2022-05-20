/**
    inmath.util

    Authors: David Herberth, Inochi2D Project
    License: MIT
*/
module inmath.util;

private {
    import inmath.linalg : Vector, Matrix, Quaternion, Rect;
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

/// If T is a vector, this evaluates to true, otherwise false.
enum isVector(T) = is(T : Vector!U, U...);

/// If T is a matrix, this evaluates to true, otherwise false.
enum isMatrix(T) = is(T : Matrix!U, U...);

/// If T is a quaternion, this evaluates to true, otherwise false.
enum isQuaternion(T) = is(T : Quaternion!U, U...);

/// If T is a rect, this evaluates to true, otherwise false.
enum isRect(T) = is(T : Rect!U, U...);

/// If T is a plane, this evaluates to true, otherwise false.
enum isPlane(T) = is(T : PlaneT!U, U...);


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
