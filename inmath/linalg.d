/**
    inmath.linalg

    Special thanks to:
    $(UL
        $(LI Tomasz Stachowiak (h3r3tic): allowed me to use parts of $(LINK2 https://bitbucket.org/h3r3tic/boxen/src/default/src/xf/omg, omg).)
        $(LI Jakob Øvrum (jA_cOp): improved the code a lot!)
        $(LI Florian Boesch (___doc__): helps me to understand opengl/complex maths better, see: $(LINK http://codeflow.org/).)
        $(LI #D on freenode: answered general questions about D.)
    )

    Authors: David Herberth, Inochi2D Project
    License: MIT

    Note: All methods marked with pure are weakly pure since, they all access an instance member.
    All static methods are strongly pure.
*/


module inmath.linalg;

private {
    import std.math : isNaN, isInfinity;
    import std.conv : to;
    import std.traits : isIntegral, isFloatingPoint, isStaticArray, isDynamicArray, isImplicitlyConvertible, isArray;
    import std.string : format, rightJustify;
    import std.array : join;
    import std.algorithm : max, min, reduce;
    import inmath.math : clamp, PI, sqrt, sin, cos, acos, tan, asin, atan2, almostEqual;
    import inmath.util : isVector, isMatrix, isQuaternion, isRect, TupleRange;
}

version(NoReciprocalMul) {
    private enum rmul = false;
} else {
    private enum rmul = true;
}

/// Base template for all vector-types.
/// Params:
/// type = all values get stored as this type
/// dimension = specifies the dimension of the vector, can be 1, 2, 3 or 4
/// Examples:
/// ---
/// alias Vector!(int, 3) vec3i;
/// alias Vector!(float, 4) vec4;
/// alias Vector!(real, 2) vec2r;
/// ---
struct Vector(type, int dimension_) {
static assert(dimension > 0, "0 dimensional vectors don't exist.");

    alias vt = type; /// Holds the internal type of the vector.
    static const int dimension = dimension_; ///Holds the dimension of the vector.

    vt[dimension] vector; /// Holds all coordinates, length conforms dimension.

    /// Returns a pointer to the coordinates.
    auto ptr() const { return vector.ptr; }

    /// Returns the current vector formatted as string, useful for printing the vector.
    const(string) toString() const {
        return format("%s", vector);
    }

    /// Gets a hash of this item
    size_t toHash() const { return typeid(this).getHash(&this); }

@safe pure nothrow:
    ///
    ref inout(vt) get_(char coord)() inout {
        return vector[coordToIndex!coord];
    }

    alias x = get_!'x'; /// static properties to access the values.
    alias u = x; /// ditto
    alias s = x; /// ditto
    alias r = x; /// ditto
    static if(dimension >= 2) {
        alias y = get_!'y'; /// ditto
        alias v = y; /// ditto
        alias t = y; /// ditto
        alias g = y; /// ditto
    }
    static if(dimension >= 3) {
        alias z = get_!'z'; /// ditto
        alias b = z; /// ditto
        alias p = z; /// ditto
    }
    static if(dimension >= 4) {
        alias w = get_!'w'; /// ditto
        alias a = w; /// ditto
        alias q = w; /// ditto
    }

    static if(dimension == 2) {
        enum Vector e1 = Vector(1.to!vt, 0.to!vt); /// canonical basis for Euclidian space
        enum Vector e2 = Vector(0.to!vt, 1.to!vt); /// ditto
    } else static if(dimension == 3) {
        enum Vector e1 = Vector(1.to!vt, 0.to!vt, 0.to!vt); /// canonical basis for Euclidian space
        enum Vector e2 = Vector(0.to!vt, 1.to!vt, 0.to!vt); /// ditto
        enum Vector e3 = Vector(0.to!vt, 0.to!vt, 1.to!vt); /// ditto
    } else static if(dimension == 4) {
        enum Vector e1 = Vector(1.to!vt, 0.to!vt, 0.to!vt, 0.to!vt); /// canonical basis for Euclidian space
        enum Vector e2 = Vector(0.to!vt, 1.to!vt, 0.to!vt, 0.to!vt); /// ditto
        enum Vector e3 = Vector(0.to!vt, 0.to!vt, 1.to!vt, 0.to!vt); /// ditto
        enum Vector e4 = Vector(0.to!vt, 0.to!vt, 0.to!vt, 1.to!vt); /// ditto
    }

    unittest {
        assert(vec2.e1.vector == [1.0, 0.0]);
        assert(vec2.e2.vector == [0.0, 1.0]);

        assert(vec3.e1.vector == [1.0, 0.0, 0.0]);
        assert(vec3.e2.vector == [0.0, 1.0, 0.0]);
        assert(vec3.e3.vector == [0.0, 0.0, 1.0]);

        assert(vec4.e1.vector == [1.0, 0.0, 0.0, 0.0]);
        assert(vec4.e2.vector == [0.0, 1.0, 0.0, 0.0]);
        assert(vec4.e3.vector == [0.0, 0.0, 1.0, 0.0]);
        assert(vec4.e4.vector == [0.0, 0.0, 0.0, 1.0]);
    }

    static void isCompatibleVectorImpl(int d)(Vector!(vt, d) vec) if(d <= dimension) {
    }

    template isCompatibleVector(T) {
        enum isCompatibleVector = is(typeof(isCompatibleVectorImpl(T.init)));
    }

    static void isCompatibleMatrixImpl(int r, int c)(Matrix!(vt, r, c) m) {
    }

    template isCompatibleMatrix(T) {
        enum isCompatibleMatrix = is(typeof(isCompatibleMatrixImpl(T.init)));
    }

    private void construct(int i, T, Tail...)(T head, Tail tail) {
        static if(i >= dimension) {
            static assert(false, "Too many arguments passed to constructor");
        } else static if(is(T : vt)) {
            vector[i] = head;
            construct!(i + 1)(tail);
        } else static if(isDynamicArray!T) {
            static assert((Tail.length == 0) && (i == 0), "dynamic array can not be passed together with other arguments");
            vector[] = head[];
        } else static if(isStaticArray!T) {
            vector[i .. i + T.length] = head[];
            construct!(i + T.length)(tail);
        } else static if(isCompatibleVector!T) {
            vector[i .. i + T.dimension] = head.vector[];
            construct!(i + T.dimension)(tail);
        } else {
            static assert(false, "Vector constructor argument must be of type " ~ vt.stringof ~ " or Vector, not " ~ T.stringof);
        }
    }

    private void construct(int i)() { // terminate
        static assert(i == dimension, "Not enough arguments passed to constructor");
    }

    /// Constructs the vector.
    /// If a single value is passed the vector, the vector will be cleared with this value.
    /// If a vector with a higher dimension is passed the vector will hold the first values up to its dimension.
    /// If mixed types are passed they will be joined together (allowed types: vector, static array, $(I vt)).
    /// Examples:
    /// ---
    /// vec4 v4 = vec4(1.0f, vec2(2.0f, 3.0f), 4.0f);
    /// vec3 v3 = vec3(v4); // v3 = vec3(1.0f, 2.0f, 3.0f);
    /// vec2 v2 = v3.xy; // swizzling returns a static array.
    /// vec3 v3_2 = vec3(1.0f); // vec3 v3_2 = vec3(1.0f, 1.0f, 1.0f);
    /// ---
    this(Args...)(Args args) {
        construct!(0)(args);
    }

    /// ditto
    this(T)(T vec) if(isVector!T && is(T.vt : vt) && (T.dimension >= dimension)) {
        foreach(i; TupleRange!(0, dimension)) {
            vector[i] = vec.vector[i];
        }
    }

    /// ditto
    this()(vt value) {
        clear(value);
    }

    /// Returns true if all values are not nan and finite, otherwise false.
    bool isFinite() const {
        static if(isIntegral!type) {
            return true;
        }
        else {
            foreach(v; vector) {
                if(isNaN(v) || isInfinity(v)) {
                    return false;
                }
            }
            return true;
        }
    }
    deprecated("Use isFinite instead of ok") alias ok = isFinite;

    /// Sets all values of the vector to value.
    void clear(vt value) {
        foreach(i; TupleRange!(0, dimension)) {
            vector[i] = value;
        }
    }

    unittest {
        vec3 vec_clear;
        assert(!vec_clear.isFinite);
        vec_clear.clear(1.0f);
        assert(vec_clear.isFinite);
        assert(vec_clear.vector == [1.0f, 1.0f, 1.0f]);
        assert(vec_clear.vector == vec3(1.0f).vector);
        vec_clear.clear(float.infinity);
        assert(!vec_clear.isFinite);
        vec_clear.clear(float.nan);
        assert(!vec_clear.isFinite);
        vec_clear.clear(1.0f);
        assert(vec_clear.isFinite);

        vec4 b = vec4(1.0f, vec_clear);
        assert(b.isFinite);
        assert(b.vector == [1.0f, 1.0f, 1.0f, 1.0f]);
        assert(b.vector == vec4(1.0f).vector);

        vec2 v2_1 = vec2(vec2(0.0f, 1.0f));
        assert(v2_1.vector == [0.0f, 1.0f]);

        vec2 v2_2 = vec2(1.0f, 1.0f);
        assert(v2_2.vector == [1.0f, 1.0f]);

        vec3 v3 = vec3(v2_1, 2.0f);
        assert(v3.vector == [0.0f, 1.0f, 2.0f]);

        vec4 v4_1 = vec4(1.0f, vec2(2.0f, 3.0f), 4.0f);
        assert(v4_1.vector == [1.0f, 2.0f, 3.0f, 4.0f]);
        assert(vec3(v4_1).vector == [1.0f, 2.0f, 3.0f]);
        assert(vec2(vec3(v4_1)).vector == [1.0f, 2.0f]);
        assert(vec2(vec3(v4_1)).vector == vec2(v4_1).vector);
        assert(v4_1.vector == vec4([1.0f, 2.0f, 3.0f, 4.0f]).vector);

        vec4 v4_2 = vec4(vec2(1.0f, 2.0f), vec2(3.0f, 4.0f));
        assert(v4_2.vector == [1.0f, 2.0f, 3.0f, 4.0f]);
        assert(vec3(v4_2).vector == [1.0f, 2.0f, 3.0f]);
        assert(vec2(vec3(v4_2)).vector == [1.0f, 2.0f]);
        assert(vec2(vec3(v4_2)).vector == vec2(v4_2).vector);
        assert(v4_2.vector == vec4([1.0f, 2.0f, 3.0f, 4.0f]).vector);

        float[2] f2 = [1.0f, 2.0f];
        float[3] f3 = [1.0f, 2.0f, 3.0f];
        float[4] f4 = [1.0f, 2.0f, 3.0f, 4.0f];
        assert(vec2(1.0f, 2.0f).vector == vec2(f2).vector);
        assert(vec3(1.0f, 2.0f, 3.0f).vector == vec3(f3).vector);
        assert(vec3(1.0f, 2.0f, 3.0f).vector == vec3(f2, 3.0f).vector);
        assert(vec4(1.0f, 2.0f, 3.0f, 4.0f).vector == vec4(f4).vector);
        assert(vec4(1.0f, 2.0f, 3.0f, 4.0f).vector == vec4(f3, 4.0f).vector);
        assert(vec4(1.0f, 2.0f, 3.0f, 4.0f).vector == vec4(f2, 3.0f, 4.0f).vector);
        // useful for: "vec4 v4 = […]" or "vec4 v4 = other_vector.rgba"

        assert(vec3(vec3i(1, 2, 3)) == vec3(1.0, 2.0, 3.0));
        assert(vec3d(vec3(1.0, 2.0, 3.0)) == vec3d(1.0, 2.0, 3.0));

        static assert(!__traits(compiles, vec3(0.0f, 0.0f)));
        static assert(!__traits(compiles, vec4(0.0f, 0.0f, 0.0f)));
        static assert(!__traits(compiles, vec4(0.0f, vec2(0.0f, 0.0f))));
        static assert(!__traits(compiles, vec4(vec3(0.0f, 0.0f, 0.0f))));
    }

    template coordToIndex(char c) {
        static if((c == 'x') || (c == 'r') || (c == 'u') || (c == 's')) {
            enum coordToIndex = 0;
        } else static if((c == 'y') || (c == 'g') || (c == 'v') || (c == 't')) {
            enum coordToIndex = 1;
        } else static if((c == 'z') || (c == 'b') || (c == 'p')) {
            static assert(dimension >= 3, "the " ~ c ~ " property is only available on vectors with a third dimension.");
            enum coordToIndex = 2;
        } else static if((c == 'w') || (c == 'a') || (c == 'q')) {
            static assert(dimension >= 4, "the " ~ c ~ " property is only available on vectors with a fourth dimension.");
            enum coordToIndex = 3;
        } else {
            static assert(false, "accepted coordinates are x, s, r, u, y, g, t, v, z, p, b, w, q and a not " ~ c ~ ".");
        }
    }

    static if(dimension == 2) { void set(vt x, vt y) { vector[0] = x; vector[1] = y; } }
    static if(dimension == 3) { void set(vt x, vt y, vt z) { vector[0] = x; vector[1] = y; vector[2] = z; } }
    static if(dimension == 4) { void set(vt x, vt y, vt z, vt w) { vector[0] = x; vector[1] = y; vector[2] = z; vector[3] = w; } }

    /// Updates the vector with the values from other.
    void update(Vector!(vt, dimension) other) {
        vector = other.vector;
    }

    unittest {
        vec2 v2 = vec2(1.0f, 2.0f);
        assert(v2.x == 1.0f);
        assert(v2.y == 2.0f);
        v2.x = 3.0f;
        v2.x += 1;
        v2.x -= 1;
        assert(v2.vector == [3.0f, 2.0f]);
        v2.y = 4.0f;
        v2.y += 1;
        v2.y -= 1;
        assert(v2.vector == [3.0f, 4.0f]);
        assert((v2.x == 3.0f) && (v2.x == v2.u) && (v2.x == v2.s) && (v2.x == v2.r));
        assert(v2.y == 4.0f);
        assert((v2.y == 4.0f) && (v2.y == v2.v) && (v2.y == v2.t) && (v2.y == v2.g));
        v2.set(0.0f, 1.0f);
        assert(v2.vector == [0.0f, 1.0f]);
        v2.update(vec2(3.0f, 4.0f));
        assert(v2.vector == [3.0f, 4.0f]);

        vec3 v3 = vec3(1.0f, 2.0f, 3.0f);
        assert(v3.x == 1.0f);
        assert(v3.y == 2.0f);
        assert(v3.z == 3.0f);
        v3.x = 3.0f;
        v3.x += 1;
        v3.x -= 1;
        assert(v3.vector == [3.0f, 2.0f, 3.0f]);
        v3.y = 4.0f;
        v3.y += 1;
        v3.y -= 1;
        assert(v3.vector == [3.0f, 4.0f, 3.0f]);
        v3.z = 5.0f;
        v3.z += 1;
        v3.z -= 1;
        assert(v3.vector == [3.0f, 4.0f, 5.0f]);
        assert((v3.x == 3.0f) && (v3.x == v3.s) && (v3.x == v3.r));
        assert((v3.y == 4.0f) && (v3.y == v3.t) && (v3.y == v3.g));
        assert((v3.z == 5.0f) && (v3.z == v3.p) && (v3.z == v3.b));
        v3.set(0.0f, 1.0f, 2.0f);
        assert(v3.vector == [0.0f, 1.0f, 2.0f]);
        v3.update(vec3(3.0f, 4.0f, 5.0f));
        assert(v3.vector == [3.0f, 4.0f, 5.0f]);

        vec4 v4 = vec4(1.0f, 2.0f, vec2(3.0f, 4.0f));
        assert(v4.x == 1.0f);
        assert(v4.y == 2.0f);
        assert(v4.z == 3.0f);
        assert(v4.w == 4.0f);
        v4.x = 3.0f;
        v4.x += 1;
        v4.x -= 1;
        assert(v4.vector == [3.0f, 2.0f, 3.0f, 4.0f]);
        v4.y = 4.0f;
        v4.y += 1;
        v4.y -= 1;
        assert(v4.vector == [3.0f, 4.0f, 3.0f, 4.0f]);
        v4.z = 5.0f;
        v4.z += 1;
        v4.z -= 1;
        assert(v4.vector == [3.0f, 4.0f, 5.0f, 4.0f]);
        v4.w = 6.0f;
        v4.w += 1;
        v4.w -= 1;
        assert(v4.vector == [3.0f, 4.0f, 5.0f, 6.0f]);
        assert((v4.x == 3.0f) && (v4.x == v4.s) && (v4.x == v4.r));
        assert((v4.y == 4.0f) && (v4.y == v4.t) && (v4.y == v4.g));
        assert((v4.z == 5.0f) && (v4.z == v4.p) && (v4.z == v4.b));
        assert((v4.w == 6.0f) && (v4.w == v4.q) && (v4.w == v4.a));
        v4.set(0.0f, 1.0f, 2.0f, 3.0f);
        assert(v4.vector == [0.0f, 1.0f, 2.0f, 3.0f]);
        v4.update(vec4(3.0f, 4.0f, 5.0f, 6.0f));
        assert(v4.vector == [3.0f, 4.0f, 5.0f, 6.0f]);
    }

    private void dispatchImpl(int i, string s, int size)(ref vt[size] result) const {
        static if(s.length > 0) {
            result[i] = vector[coordToIndex!(s[0])];
            dispatchImpl!(i + 1, s[1..$])(result);
        }
    }

    /// Implements dynamic swizzling.
    /// Returns: a Vector
    Vector!(vt, s.length) opDispatch(string s)() const {
        vt[s.length] ret;
        dispatchImpl!(0, s)(ret);
        Vector!(vt, s.length) ret_vec;
        ret_vec.vector = ret;
        return ret_vec;
    }

    unittest {
        vec2 v2 = vec2(1.0f, 2.0f);
        assert(v2.xytsy == [1.0f, 2.0f, 2.0f, 1.0f, 2.0f]);

        assert(vec3(1.0f, 2.0f, 3.0f).xybzyr == [1.0f, 2.0f, 3.0f, 3.0f, 2.0f, 1.0f]);
        assert(vec4(v2, 3.0f, 4.0f).xyzwrgbastpq == [1.0f, 2.0f, 3.0f, 4.0f,
                                                     1.0f, 2.0f, 3.0f, 4.0f,
                                                     1.0f, 2.0f, 3.0f, 4.0f]);
        assert(vec4(v2, 3.0f, 4.0f).wgyzax == [4.0f, 2.0f, 2.0f, 3.0f, 4.0f, 1.0f]);
        assert(vec4(v2.xyst).vector == [1.0f, 2.0f, 1.0f, 2.0f]);
    }

    /// Returns the squared magnitude of the vector.
    real lengthSquared() const {
        real temp = 0;

        foreach(index; TupleRange!(0, dimension)) {
            temp += vector[index]^^2;
        }

        return temp;
    }

    /// Returns the magnitude of the vector.
    real length() const {
        return sqrt(lengthSquared);
    }

    /// Normalizes the vector.
    void normalize() {
        real len = length;

        if(len != 0) {
            foreach(index; TupleRange!(0, dimension)) {
                vector[index] = cast(type)(vector[index]/len);
            }
        }
    }

    /// Returns a normalized copy of the current vector.
    Vector normalized() const {
        Vector ret;
        ret.update(this);
        ret.normalize();
        return ret;
    }

    Vector opUnary(string op : "-")() const {
        Vector ret;

        foreach(index; TupleRange!(0, dimension)) {
            ret.vector[index] = -vector[index];
        }

        return ret;
    }

    unittest {
        assert(vec2(1.0f, 1.0f) == -vec2(-1.0f, -1.0f));
        assert(vec2(-1.0f, 1.0f) == -vec2(1.0f, -1.0f));

        assert(-vec3(1.0f, 1.0f, 1.0f) == vec3(-1.0f, -1.0f, -1.0f));
        assert(-vec3(-1.0f, 1.0f, -1.0f) == vec3(1.0f, -1.0f, 1.0f));

        assert(vec4(1.0f, 1.0f, 1.0f, 1.0f) == -vec4(-1.0f, -1.0f, -1.0f, -1.0f));
        assert(vec4(-1.0f, 1.0f, -1.0f, 1.0f) == -vec4(1.0f, -1.0f, 1.0f, -1.0f));
    }

    // let the math begin!
    Vector opBinary(string op : "*")(vt r) const {
        Vector ret;

        static foreach(index; TupleRange!(0, dimension)) {
            ret.vector[index] = vector[index] * r;
        }

        return ret;
    }

    Vector opBinary(string op : "/")(vt r) const {
        Vector ret;

        static foreach(index; TupleRange!(0, dimension)) {
            ret.vector[index] = cast(vt)(vector[index] / r);
        }

        return ret;
    }

    Vector opBinary(string op)(Vector r) const if((op == "+") || (op == "-")) {
        Vector ret;

        static foreach(index; TupleRange!(0, dimension)) {
            ret.vector[index] = mixin("cast(vt)(vector[index]" ~ op ~ "r.vector[index])");
        }

        return ret;
    }

    Vector opBinary(string op : "*")(Vector r) const {
        Vector ret;

        static foreach(index; 0..dimension) {
            ret.vector[index] = vector[index] * r.vector[index];
        }

        return ret;
    }

    // vector * matrix (for matrix * vector -> struct Matrix)
    Vector!(vt, T.cols) opBinary(string op : "*", T)(T inp) const if(isCompatibleMatrix!T && (T.rows == dimension)) {
        Vector!(vt, T.cols) ret;
        ret.clear(0);

        foreach(c; TupleRange!(0, T.cols)) {
            foreach(r; TupleRange!(0, T.rows)) {
                ret.vector[c] += vector[r] * inp.matrix[r][c];
            }
        }

        return ret;
    }

    auto opBinaryRight(string op, T)(T inp) const if(!isVector!T && !isMatrix!T && !isQuaternion!T) {
        return this.opBinary!(op)(inp);
    }

    unittest {
        // vec2 v2 = vec2(1.0f, 3.0f);
        // auto v2times2 = 2 * v2;
        // assert((v2*2.5f).vector == [2.5f, 7.5f]);
        // assert((v2+vec2(3.0f, 1.0f)).vector == [4.0f, 4.0f]);
        // assert((v2-vec2(1.0f, 3.0f)).vector == [0.0f, 0.0f]);
        // assert((v2*vec2(2.0f, 2.0f)) == 8.0f);

        // vec3 v3 = vec3(1.0f, 3.0f, 5.0f);
        // assert((v3*2.5f).vector == [2.5f, 7.5f, 12.5f]);
        // assert((v3+vec3(3.0f, 1.0f, -1.0f)).vector == [4.0f, 4.0f, 4.0f]);
        // assert((v3-vec3(1.0f, 3.0f, 5.0f)).vector == [0.0f, 0.0f, 0.0f]);
        // assert((v3*vec3(2.0f, 2.0f, 2.0f)) == 18.0f);

        // vec4 v4 = vec4(1.0f, 3.0f, 5.0f, 7.0f);
        // assert((v4*2.5f).vector == [2.5f, 7.5f, 12.5f, 17.5]);
        // assert((v4+vec4(3.0f, 1.0f, -1.0f, -3.0f)).vector == [4.0f, 4.0f, 4.0f, 4.0f]);
        // assert((v4-vec4(1.0f, 3.0f, 5.0f, 7.0f)).vector == [0.0f, 0.0f, 0.0f, 0.0f]);
        // assert((v4*vec4(2.0f, 2.0f, 2.0f, 2.0f)) == 32.0f);

        // mat2 m2 = mat2(1.0f, 2.0f, 3.0f, 4.0f);
        // vec2 v2_2 = vec2(2.0f, 2.0f);
        // assert((v2_2*m2).vector == [8.0f, 12.0f]);

        // mat3 m3 = mat3(1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f, 9.0f);
        // vec3 v3_2 = vec3(2.0f, 2.0f, 2.0f);
        // assert((v3_2*m3).vector == [24.0f, 30.0f, 36.0f]);
    }

    void opOpAssign(string op : "*")(vt r) {
        foreach(index; TupleRange!(0, dimension)) {
            vector[index] *= r;
        }
    }

    void opOpAssign(string op : "/")(vt r) {
        foreach(index; TupleRange!(0, dimension)) {
            vector[index] /= r;
        }
    }

    void opOpAssign(string op)(Vector r) if((op == "+") || (op == "-")) {
        foreach(index; TupleRange!(0, dimension)) {
            mixin("vector[index]" ~ op ~ "= r.vector[index];");
        }
    }

    unittest {
        vec2 v2 = vec2(1.0f, 3.0f);
        v2 *= 2.5f;
        assert(v2.vector == [2.5f, 7.5f]);
        v2 -= vec2(2.5f, 7.5f);
        assert(v2.vector == [0.0f, 0.0f]);
        v2 += vec2(1.0f, 3.0f);
        assert(v2.vector == [1.0f, 3.0f]);
        assert(almostEqual(v2.length, sqrt(10.0f)));
        assert(v2.lengthSquared == 10.0f);
        assert((v2.length == v2.length) && (v2.lengthSquared == v2.lengthSquared));
        v2 /= 2.0f;
        assert(v2.vector == [0.5f, 1.5f]);
        assert(almostEqual(v2.normalized, vec2(1.0f/sqrt(10.0f), 3.0f/sqrt(10.0f))));

        vec3 v3 = vec3(1.0f, 3.0f, 5.0f);
        v3 *= 2.5f;
        assert(v3.vector == [2.5f, 7.5f, 12.5f]);
        v3 -= vec3(2.5f, 7.5f, 12.5f);
        assert(v3.vector == [0.0f, 0.0f, 0.0f]);
        v3 += vec3(1.0f, 3.0f, 5.0f);
        assert(v3.vector == [1.0f, 3.0f, 5.0f]);
        assert(almostEqual(v3.length, sqrt(35.0f)));
        assert(v3.lengthSquared == 35.0f);
        assert((v3.length == v3.length) && (v3.lengthSquared == v3.lengthSquared));
        v3 /= 2.0f;
        assert(v3.vector == [0.5f, 1.5f, 2.5f]);
        assert(almostEqual(v3.normalized, vec3(1.0f/sqrt(35.0f), 3.0f/sqrt(35.0f), 5.0f/sqrt(35.0f))));

        vec4 v4 = vec4(1.0f, 3.0f, 5.0f, 7.0f);
        v4 *= 2.5f;
        assert(v4.vector == [2.5f, 7.5f, 12.5f, 17.5]);
        v4 -= vec4(2.5f, 7.5f, 12.5f, 17.5f);
        assert(v4.vector == [0.0f, 0.0f, 0.0f, 0.0f]);
        v4 += vec4(1.0f, 3.0f, 5.0f, 7.0f);
        assert(v4.vector == [1.0f, 3.0f, 5.0f, 7.0f]);
        assert(almostEqual(v4.length, sqrt(84.0f)));
        assert(v4.lengthSquared == 84.0f);
        assert((v4.length == v4.length) && (v4.lengthSquared == v4.lengthSquared));
        v4 /= 2.0f;
        assert(v4.vector == [0.5f, 1.5f, 2.5f, 3.5f]);
        assert(almostEqual(v4.normalized, vec4(1.0f/sqrt(84.0f), 3.0f/sqrt(84.0f), 5.0f/sqrt(84.0f), 7.0f/sqrt(84.0f))));
    }

    int opCmp(ref const Vector vec) const {
        foreach(i, a; vector) {
            if(a < vec.vector[i]) {
                return -1;
            } else if(a > vec.vector[i]) {
                return 1;
            }
        }

        // Vectors are the same
        return 0;
    }

    bool opEquals(T)(const T vec) const if(!isArray!T && T.dimension == dimension) {
        return vector == vec.vector;
    }

    bool opEquals(T)(const(T)[] array) const if(!isArray!T && !isVector!T) {
        if(array.length != dimension) {
            return false;
        }

        foreach(index; TupleRange!(0, dimension)) {
            if(vector[index] != array[index]) {
                return false;
            }
        }

        return true;
    }

    unittest {
        assert(vec2(1.0f, 2.0f) == vec2(1.0f, 2.0f));
        assert(vec2(1.0f, 2.0f) != vec2(1.0f, 1.0f));
        assert(vec2(1.0f, 2.0f) == vec2d(1.0, 2.0));
        assert(vec2(1.0f, 2.0f) != vec2d(1.0, 1.0));
        assert(vec2(1.0f, 2.0f) == vec2(1.0f, 2.0f).vector);
        assert(vec2(1.0f, 2.0f) != vec2(1.0f, 1.0f).vector);
        assert(vec2(1.0f, 2.0f) == vec2d(1.0, 2.0).vector);
        assert(vec2(1.0f, 2.0f) != vec2d(1.0, 1.0).vector);

        assert(vec3(1.0f, 2.0f, 3.0f) == vec3(1.0f, 2.0f, 3.0f));
        assert(vec3(1.0f, 2.0f, 3.0f) != vec3(1.0f, 2.0f, 2.0f));
        assert(vec3(1.0f, 2.0f, 3.0f) == vec3d(1.0, 2.0, 3.0));
        assert(vec3(1.0f, 2.0f, 3.0f) != vec3d(1.0, 2.0, 2.0));
        assert(vec3(1.0f, 2.0f, 3.0f) == vec3(1.0f, 2.0f, 3.0f).vector);
        assert(vec3(1.0f, 2.0f, 3.0f) != vec3(1.0f, 2.0f, 2.0f).vector);
        assert(vec3(1.0f, 2.0f, 3.0f) == vec3d(1.0, 2.0, 3.0).vector);
        assert(vec3(1.0f, 2.0f, 3.0f) != vec3d(1.0, 2.0, 2.0).vector);

        assert(vec4(1.0f, 2.0f, 3.0f, 4.0f) == vec4(1.0f, 2.0f, 3.0f, 4.0f));
        assert(vec4(1.0f, 2.0f, 3.0f, 4.0f) != vec4(1.0f, 2.0f, 3.0f, 3.0f));
        assert(vec4(1.0f, 2.0f, 3.0f, 4.0f) == vec4d(1.0, 2.0, 3.0, 4.0));
        assert(vec4(1.0f, 2.0f, 3.0f, 4.0f) != vec4d(1.0, 2.0, 3.0, 3.0));
        assert(vec4(1.0f, 2.0f, 3.0f, 4.0f) == vec4(1.0f, 2.0f, 3.0f, 4.0f).vector);
        assert(vec4(1.0f, 2.0f, 3.0f, 4.0f) != vec4(1.0f, 2.0f, 3.0f, 3.0f).vector);
        assert(vec4(1.0f, 2.0f, 3.0f, 4.0f) == vec4d(1.0, 2.0, 3.0, 4.0).vector);
        assert(vec4(1.0f, 2.0f, 3.0f, 4.0f) != vec4d(1.0, 2.0, 3.0, 3.0).vector);

        assert(!vec4(float.nan).isFinite);
        if(vec4(1.0f).isFinite) { }
        else { assert(false); }
    }
}

/// Calculates the product between two vectors.
T.vt dot(T)(const T veca, const T vecb) @safe pure nothrow if(isVector!T) {
    T.vt temp = 0;

    foreach(index; TupleRange!(0, T.dimension)) {
        temp += veca.vector[index] * vecb.vector[index];
    }

    return temp;
}

/// Calculates the cross product of two 3-dimensional vectors.
T cross(T)(const T veca, const T vecb) @safe pure nothrow if(isVector!T && (T.dimension == 3)) {
   return T(veca.y * vecb.z - vecb.y * veca.z,
            veca.z * vecb.x - vecb.z * veca.x,
            veca.x * vecb.y - vecb.x * veca.y);
}

/// Calculates the distance between two vectors.
T.vt distance(T)(const T veca, const T vecb) @safe pure nothrow if(isVector!T) {
    return (veca - vecb).length;
}

@("vector dot product")
unittest {
    // dot is already tested in Vector.opBinary, so no need for testing with more vectors
    vec3 v1 = vec3(1.0f, 2.0f, -3.0f);
    vec3 v2 = vec3(1.0f, 3.0f, 2.0f);

    assert(dot(v1, v2) == 1.0f);
    assert(dot(v1, v2) == dot(v2, v1));
    assert((v1 * v2) == (v1 * v2));

    assert(cross(v1, v2).vector == [13.0f, -5.0f, 1.0f]);
    assert(cross(v2, v1).vector == [-13.0f, 5.0f, -1.0f]);

    assert(distance(vec2(0.0f, 0.0f), vec2(0.0f, 10.0f)) == 10.0);
}

/// reflect a vector using a surface normal
T reflect(T)(const T vec, const T norm) @safe pure nothrow if(isVector!T) {
    return (2 * (vec * norm) * norm) - vec;
}

@("vector reflect")
unittest
{
    assert(vec2(1,1).reflect(vec2(0,1)) == vec2(-1,1));
    assert(vec2(-1,1).reflect(vec2(0,1)) == vec2(1,1));
    assert(vec2(2,1).reflect(vec2(0,1)) == vec2(-2,1));

    assert(vec3(1,1,1).reflect(vec3(0,1,0)) == vec3(-1,1,-1));
}

/// Pre-defined vector types, the number represents the dimension and the last letter the type (none = float, d = double, i = int).
alias vec2 = Vector!(float, 2);
alias vec3 = Vector!(float, 3); /// ditto
alias vec4 = Vector!(float, 4); /// ditto

alias vec2d = Vector!(double, 2); /// ditto
alias vec3d = Vector!(double, 3); /// ditto
alias vec4d = Vector!(double, 4); /// ditto

alias vec2i = Vector!(int, 2); /// ditto
alias vec3i = Vector!(int, 3); /// ditto
alias vec4i = Vector!(int, 4); /// ditto

alias vec2u = Vector!(uint, 2); /// ditto
alias vec3u = Vector!(uint, 3); /// ditto
alias vec4u = Vector!(uint, 4); /// ditto

/*alias Vector!(ubyte, 2) vec2ub;
alias Vector!(ubyte, 3) vec3ub;
alias Vector!(ubyte, 4) vec4ub;*/


/// Base template for all matrix-types.
/// Params:
///  type = all values get stored as this type
///  rows_ = rows of the matrix
///  cols_ = columns of the matrix
/// Examples:
/// ---
/// alias Matrix!(float, 4, 4) mat4;
/// alias Matrix!(double, 3, 4) mat34d;
/// alias Matrix!(real, 2, 2) mat2r;
/// ---
struct Matrix(type, int rows_, int cols_) if((rows_ > 0) && (cols_ > 0)) {
    alias mt = type; /// Holds the internal type of the matrix;
    static const int rows = rows_; /// Holds the number of rows;
    static const int cols = cols_; /// Holds the number of columns;

    /// Holds the matrix $(RED row-major) in memory.
    mt[cols][rows] matrix; // In C it would be mt[rows][cols], D does it like this: (mt[foo])[bar]
    alias matrix this;

    unittest {
        mat2 m2 = mat2(0.0f, 1.0f, 2.0f, 3.0f);
        assert(m2[0][0] == 0.0f);
        assert(m2[0][1] == 1.0f);
        assert(m2[1][0] == 2.0f);
        assert(m2[1][1] == 3.0f);
        m2[0..1] = [2.0f, 2.0f];
        assert(m2 == [[2.0f, 2.0f], [2.0f, 3.0f]]);

        mat3 m3 = mat3(0.0f, 0.1f, 0.2f, 1.0f, 1.1f, 1.2f, 2.0f, 2.1f, 2.2f);
        assert(m3[0][1] == 0.1f);
        assert(m3[2][0] == 2.0f);
        assert(m3[1][2] == 1.2f);
        m3[0][0..$] = 0.0f;
        assert(m3 == [[0.0f, 0.0f, 0.0f],
                      [1.0f, 1.1f, 1.2f],
                      [2.0f, 2.1f, 2.2f]]);

        mat4 m4 = mat4(0.0f, 0.1f, 0.2f, 0.3f,
                       1.0f, 1.1f, 1.2f, 1.3f,
                       2.0f, 2.1f, 2.2f, 2.3f,
                       3.0f, 3.1f, 3.2f, 3.3f);
       assert(m4[0][3] == 0.3f);
       assert(m4[1][1] == 1.1f);
       assert(m4[2][0] == 2.0f);
       assert(m4[3][2] == 3.2f);
       m4[2][1..3] = [1.0f, 2.0f];
       assert(m4 == [[0.0f, 0.1f, 0.2f, 0.3f],
                     [1.0f, 1.1f, 1.2f, 1.3f],
                     [2.0f, 1.0f, 2.0f, 2.3f],
                     [3.0f, 3.1f, 3.2f, 3.3f]]);

    }

    /// Returns the pointer to the stored values as OpenGL requires it.
    /// Note this will return a pointer to a $(RED row-major) matrix,
    /// $(RED this means you've to set the transpose argument to GL_TRUE when passing it to OpenGL).
    /// Examples:
    /// ---
    /// // 3rd argument = GL_TRUE
    /// glUniformMatrix4fv(programs.main.model, 1, GL_TRUE, mat4.translation(-0.5f, -0.5f, 1.0f).ptr);
    /// ---
    auto ptr() const { return matrix[0].ptr; }

    /// Returns the current matrix formatted as flat string.
    string toString() const {
        return format("%s", matrix);
    }

    /// Returns the current matrix as pretty formatted string.
    string toPrettyString() {
        string fmtr = "%s";

        size_t rjust = max(format(fmtr, reduce!(max)(matrix[])).length,
                           format(fmtr, reduce!(min)(matrix[])).length) - 1;

        string[] outer_parts;
        foreach(mt[] row; matrix) {
            string[] inner_parts;
            foreach(mt col; row) {
                inner_parts ~= rightJustify(format(fmtr, col), rjust);
            }
            outer_parts ~= " [" ~ join(inner_parts, ", ") ~ "]";
        }

        return "[" ~ join(outer_parts, "\n")[1..$] ~ "]";
    }

@safe pure nothrow:
    static void isCompatibleMatrixImpl(int r, int c)(Matrix!(mt, r, c) m) {
    }

    template isCompatibleMatrix(T) {
        enum isCompatibleMatrix = is(typeof(isCompatibleMatrixImpl(T.init)));
    }

    static void isCompatibleVectorImpl(int d)(Vector!(mt, d) vec) {
    }

    template isCompatibleVector(T) {
        enum isCompatibleVector = is(typeof(isCompatibleVectorImpl(T.init)));
    }

    private void construct(int i, T, Tail...)(T head, Tail tail) {
        static if(i >= rows*cols) {
            static assert(false, "Too many arguments passed to constructor");
        } else static if(is(T : mt)) {
            matrix[i / cols][i % cols] = head;
            construct!(i + 1)(tail);
        } else static if(is(T == Vector!(mt, cols))) {
            static if(i % cols == 0) {
                matrix[i / cols] = head.vector;
                construct!(i + T.dimension)(tail);
            } else {
                static assert(false, "Can't convert Vector into the matrix. Maybe it doesn't align to the columns correctly or dimension doesn't fit");
            }
        } else static if(isDynamicArray!T) {
            foreach(j; 0..cols*rows)
                matrix[j / cols][j % cols] = head[j];
        } else {
            static assert(false, "Matrix constructor argument must be of type " ~ mt.stringof ~ " or Vector, not " ~ T.stringof);
        }
    }

    private void construct(int i)() { // terminate
        static assert(i == rows*cols, "Not enough arguments passed to constructor");
    }

    /// Constructs the matrix:
    /// If a single value is passed, the matrix will be cleared with this value (each column in each row will contain this value).
    /// If a matrix with more rows and columns is passed, the matrix will be the upper left nxm matrix.
    /// If a matrix with less rows and columns is passed, the passed matrix will be stored in the upper left of an identity matrix.
    /// It's also allowed to pass vectors and scalars at a time, but the vectors dimension must match the number of columns and align correctly.
    /// Examples:
    /// ---
    /// mat2 m2 = mat2(0.0f); // mat2 m2 = mat2(0.0f, 0.0f, 0.0f, 0.0f);
    /// mat3 m3 = mat3(m2); // mat3 m3 = mat3(0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f);
    /// mat3 m3_2 = mat3(vec3(1.0f, 2.0f, 3.0f), 4.0f, 5.0f, 6.0f, vec3(7.0f, 8.0f, 9.0f));
    /// mat4 m4 = mat4.identity; // just an identity matrix
    /// mat3 m3_3 = mat3(m4); // mat3 m3_3 = mat3.identity
    /// ---
    this(Args...)(Args args) {
        construct!(0)(args);
    }

    /// ditto
    this(T)(T mat) if(isMatrix!T && (T.cols >= cols) && (T.rows >= rows)) {
        foreach(r; TupleRange!(0, rows)) {
            foreach(c; TupleRange!(0, cols)) {
                matrix[r][c] = mat.matrix[r][c];
            }
        }
    }

    /// ditto
    this(T)(T mat) if(isMatrix!T && (T.cols < cols) && (T.rows < rows)) {
        makeIdentity();

        foreach(r; TupleRange!(0, T.rows)) {
            foreach(c; TupleRange!(0, T.cols)) {
                matrix[r][c] = mat.matrix[r][c];
            }
        }
    }

    /// ditto
    this()(mt value) {
        clear(value);
    }

    /// Returns true if all values are not nan and finite, otherwise false.
    bool isFinite() const {
        static if(isIntegral!type) {
            return true;
        }
        else {
            foreach(row; matrix) {
                foreach(col; row) {
                    if(isNaN(col) || isInfinity(col)) {
                        return false;
                    }
                }
            }
            return true;
        }

    }
    deprecated("Use isFinite instead of ok") alias ok = isFinite;

    /// Sets all values of the matrix to value (each column in each row will contain this value).
    void clear(mt value) {
        foreach(r; TupleRange!(0, rows)) {
            foreach(c; TupleRange!(0, cols)) {
                matrix[r][c] = value;
            }
        }
    }

    unittest {
        mat2 m2 = mat2(1.0f, 1.0f, vec2(2.0f, 2.0f));
        assert(m2.matrix == [[1.0f, 1.0f], [2.0f, 2.0f]]);
        m2.clear(3.0f);
        assert(m2.matrix == [[3.0f, 3.0f], [3.0f, 3.0f]]);
        assert(m2.isFinite);
        m2.clear(float.nan);
        assert(!m2.isFinite);
        m2.clear(float.infinity);
        assert(!m2.isFinite);
        m2.clear(0.0f);
        assert(m2.isFinite);

        mat3 m3 = mat3(1.0f);
        assert(m3.matrix == [[1.0f, 1.0f, 1.0f],
                             [1.0f, 1.0f, 1.0f],
                             [1.0f, 1.0f, 1.0f]]);

        mat4 m4 = mat4(vec4(1.0f, 1.0f, 1.0f, 1.0f),
                            2.0f, 2.0f, 2.0f, 2.0f,
                            3.0f, 3.0f, 3.0f, 3.0f,
                       vec4(4.0f, 4.0f, 4.0f, 4.0f));
        assert(m4.matrix == [[1.0f, 1.0f, 1.0f, 1.0f],
                             [2.0f, 2.0f, 2.0f, 2.0f],
                             [3.0f, 3.0f, 3.0f, 3.0f],
                             [4.0f, 4.0f, 4.0f, 4.0f]]);
        assert(mat3(m4).matrix == [[1.0f, 1.0f, 1.0f],
                                   [2.0f, 2.0f, 2.0f],
                                   [3.0f, 3.0f, 3.0f]]);
        assert(mat2(mat3(m4)).matrix == [[1.0f, 1.0f], [2.0f, 2.0f]]);
        assert(mat2(m4).matrix == mat2(mat3(m4)).matrix);
        assert(mat4(mat3(m4)).matrix == [[1.0f, 1.0f, 1.0f, 0.0f],
                                         [2.0f, 2.0f, 2.0f, 0.0f],
                                         [3.0f, 3.0f, 3.0f, 0.0f],
                                         [0.0f, 0.0f, 0.0f, 1.0f]]);

        Matrix!(float, 2, 3) mt1 = Matrix!(float, 2, 3)(1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f);
        Matrix!(float, 3, 2) mt2 = Matrix!(float, 3, 2)(6.0f, -1.0f, 3.0f, 2.0f, 0.0f, -3.0f);

        assert(mt1.matrix == [[1.0f, 2.0f, 3.0f], [4.0f, 5.0f, 6.0f]]);
        assert(mt2.matrix == [[6.0f, -1.0f], [3.0f, 2.0f], [0.0f, -3.0f]]);

        static assert(!__traits(compiles, mat2(1, 2, 1)));
        static assert(!__traits(compiles, mat3(1, 2, 3, 1, 2, 3, 1, 2)));
        static assert(!__traits(compiles, mat4(1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3)));

        auto m5 = mat2([0.0f,1,2,3]);
        assert(m5.matrix == [[0.0f, 1.0f], [2.0f, 3.0f]]);

        auto m6 = Matrix!(int, 2, 3)([0,1,2,3,4,5]);
        assert(m6.matrix == [[0, 1, 2], [3, 4, 5]]);
    }

    static if(rows == cols) {
        /// Makes the current matrix an identity matrix.
        void makeIdentity() {
            clear(0);
            foreach(r; TupleRange!(0, rows)) {
                matrix[r][r] = 1;
            }
        }

        /// Returns a identity matrix.
        static Matrix identity() {
            Matrix ret;
            ret.clear(0);

            foreach(r; TupleRange!(0, rows)) {
                ret.matrix[r][r] = 1;
            }

            return ret;
        }

        /// Transposes the current matrix;
        void transpose() {
            matrix = transposed().matrix;
        }

        unittest {
            mat2 m2 = mat2(1.0f);
            m2.transpose();
            assert(m2.matrix == mat2(1.0f).matrix);
            m2.makeIdentity();
            assert(m2.matrix == [[1.0f, 0.0f],
                                 [0.0f, 1.0f]]);
            m2.transpose();
            assert(m2.matrix == [[1.0f, 0.0f],
                                 [0.0f, 1.0f]]);
            assert(m2.matrix == m2.identity.matrix);

            mat3 m3 = mat3(1.1f, 1.2f, 1.3f,
                           2.1f, 2.2f, 2.3f,
                           3.1f, 3.2f, 3.3f);
            m3.transpose();
            assert(m3.matrix == [[1.1f, 2.1f, 3.1f],
                                 [1.2f, 2.2f, 3.2f],
                                 [1.3f, 2.3f, 3.3f]]);

            mat4 m4 = mat4(2.0f);
            m4.transpose();
            assert(m4.matrix == mat4(2.0f).matrix);
            m4.makeIdentity();
            assert(m4.matrix == [[1.0f, 0.0f, 0.0f, 0.0f],
                                 [0.0f, 1.0f, 0.0f, 0.0f],
                                 [0.0f, 0.0f, 1.0f, 0.0f],
                                 [0.0f, 0.0f, 0.0f, 1.0f]]);
            assert(m4.matrix == m4.identity.matrix);
        }

    }

    /// Returns a transposed copy of the matrix.
    Matrix!(mt, cols, rows) transposed() const {
        typeof(return) ret;

        foreach(r; TupleRange!(0, rows)) {
            foreach(c; TupleRange!(0, cols)) {
                ret.matrix[c][r] = matrix[r][c];
            }
        }

        return ret;
    }

    // transposed already tested in last unittest


    static if((rows == 2) && (cols == 2)) {
        mt det() const {
            return (matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]);
        }

        private Matrix invert(ref Matrix mat) const {
            static if(isFloatingPoint!mt && rmul) {
                mt d = 1 / det;

                mat.matrix = [[matrix[1][1]*d, -matrix[0][1]*d],
                              [-matrix[1][0]*d, matrix[0][0]*d]];
            } else {
                mt d = det;

                mat.matrix = [[matrix[1][1]/d, -matrix[0][1]/d],
                              [-matrix[1][0]/d, matrix[0][0]/d]];
            }

            return mat;
        }

        static Matrix scaling(mt x, mt y) {
            Matrix ret = Matrix.identity;

            ret.matrix[0][0] = x;
            ret.matrix[1][1] = y;

            return ret;
        }

        Matrix scale(mt x, mt y) {
            this = Matrix.scaling(x, y) * this;
            return this;
        }

        unittest {
            assert(mat2.scaling(3, 3).matrix == mat2.identity.scale(3, 3).matrix);
            assert(mat2.scaling(3, 3).matrix == [[3.0f, 0.0f], [0.0f, 3.0f]]);
        }

    } else static if((rows == 3) && (cols == 3)) {
        mt det() const {
            return (matrix[0][0] * matrix[1][1] * matrix[2][2]
                  + matrix[0][1] * matrix[1][2] * matrix[2][0]
                  + matrix[0][2] * matrix[1][0] * matrix[2][1]
                  - matrix[0][2] * matrix[1][1] * matrix[2][0]
                  - matrix[0][1] * matrix[1][0] * matrix[2][2]
                  - matrix[0][0] * matrix[1][2] * matrix[2][1]);
        }

        private Matrix invert(ref Matrix mat) const {
            static if(isFloatingPoint!mt && rmul) {
                mt d = 1 / det;
                enum op = "*";
            } else {
                mt d = det;
                enum op = "/";
            }

            mixin(`
            mat.matrix = [[(matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])`~op~`d,
                           (matrix[0][2] * matrix[2][1] - matrix[0][1] * matrix[2][2])`~op~`d,
                           (matrix[0][1] * matrix[1][2] - matrix[0][2] * matrix[1][1])`~op~`d],
                          [(matrix[1][2] * matrix[2][0] - matrix[1][0] * matrix[2][2])`~op~`d,
                           (matrix[0][0] * matrix[2][2] - matrix[0][2] * matrix[2][0])`~op~`d,
                           (matrix[0][2] * matrix[1][0] - matrix[0][0] * matrix[1][2])`~op~`d],
                          [(matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0])`~op~`d,
                           (matrix[0][1] * matrix[2][0] - matrix[0][0] * matrix[2][1])`~op~`d,
                           (matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0])`~op~`d]];
            `);

            return mat;
        }
    } else static if((rows == 4) && (cols == 4)) {
        /// Returns the determinant of the current matrix (2x2, 3x3 and 4x4 matrices).
        mt det() const {
            return (matrix[0][3] * matrix[1][2] * matrix[2][1] * matrix[3][0] - matrix[0][2] * matrix[1][3] * matrix[2][1] * matrix[3][0]
                  - matrix[0][3] * matrix[1][1] * matrix[2][2] * matrix[3][0] + matrix[0][1] * matrix[1][3] * matrix[2][2] * matrix[3][0]
                  + matrix[0][2] * matrix[1][1] * matrix[2][3] * matrix[3][0] - matrix[0][1] * matrix[1][2] * matrix[2][3] * matrix[3][0]
                  - matrix[0][3] * matrix[1][2] * matrix[2][0] * matrix[3][1] + matrix[0][2] * matrix[1][3] * matrix[2][0] * matrix[3][1]
                  + matrix[0][3] * matrix[1][0] * matrix[2][2] * matrix[3][1] - matrix[0][0] * matrix[1][3] * matrix[2][2] * matrix[3][1]
                  - matrix[0][2] * matrix[1][0] * matrix[2][3] * matrix[3][1] + matrix[0][0] * matrix[1][2] * matrix[2][3] * matrix[3][1]
                  + matrix[0][3] * matrix[1][1] * matrix[2][0] * matrix[3][2] - matrix[0][1] * matrix[1][3] * matrix[2][0] * matrix[3][2]
                  - matrix[0][3] * matrix[1][0] * matrix[2][1] * matrix[3][2] + matrix[0][0] * matrix[1][3] * matrix[2][1] * matrix[3][2]
                  + matrix[0][1] * matrix[1][0] * matrix[2][3] * matrix[3][2] - matrix[0][0] * matrix[1][1] * matrix[2][3] * matrix[3][2]
                  - matrix[0][2] * matrix[1][1] * matrix[2][0] * matrix[3][3] + matrix[0][1] * matrix[1][2] * matrix[2][0] * matrix[3][3]
                  + matrix[0][2] * matrix[1][0] * matrix[2][1] * matrix[3][3] - matrix[0][0] * matrix[1][2] * matrix[2][1] * matrix[3][3]
                  - matrix[0][1] * matrix[1][0] * matrix[2][2] * matrix[3][3] + matrix[0][0] * matrix[1][1] * matrix[2][2] * matrix[3][3]);
        }

        private Matrix invert(ref Matrix mat) const {
            static if(isFloatingPoint!mt && rmul) {
                mt d = 1 / det;
                enum op = "*";
            } else {
                mt d = det;
                enum op = "/";
            }

            mixin(`
            mat.matrix = [[(matrix[1][1] * matrix[2][2] * matrix[3][3] + matrix[1][2] * matrix[2][3] * matrix[3][1] + matrix[1][3] * matrix[2][1] * matrix[3][2]
                          - matrix[1][1] * matrix[2][3] * matrix[3][2] - matrix[1][2] * matrix[2][1] * matrix[3][3] - matrix[1][3] * matrix[2][2] * matrix[3][1])`~op~`d,
                           (matrix[0][1] * matrix[2][3] * matrix[3][2] + matrix[0][2] * matrix[2][1] * matrix[3][3] + matrix[0][3] * matrix[2][2] * matrix[3][1]
                          - matrix[0][1] * matrix[2][2] * matrix[3][3] - matrix[0][2] * matrix[2][3] * matrix[3][1] - matrix[0][3] * matrix[2][1] * matrix[3][2])`~op~`d,
                           (matrix[0][1] * matrix[1][2] * matrix[3][3] + matrix[0][2] * matrix[1][3] * matrix[3][1] + matrix[0][3] * matrix[1][1] * matrix[3][2]
                          - matrix[0][1] * matrix[1][3] * matrix[3][2] - matrix[0][2] * matrix[1][1] * matrix[3][3] - matrix[0][3] * matrix[1][2] * matrix[3][1])`~op~`d,
                           (matrix[0][1] * matrix[1][3] * matrix[2][2] + matrix[0][2] * matrix[1][1] * matrix[2][3] + matrix[0][3] * matrix[1][2] * matrix[2][1]
                          - matrix[0][1] * matrix[1][2] * matrix[2][3] - matrix[0][2] * matrix[1][3] * matrix[2][1] - matrix[0][3] * matrix[1][1] * matrix[2][2])`~op~`d],
                          [(matrix[1][0] * matrix[2][3] * matrix[3][2] + matrix[1][2] * matrix[2][0] * matrix[3][3] + matrix[1][3] * matrix[2][2] * matrix[3][0]
                          - matrix[1][0] * matrix[2][2] * matrix[3][3] - matrix[1][2] * matrix[2][3] * matrix[3][0] - matrix[1][3] * matrix[2][0] * matrix[3][2])`~op~`d,
                           (matrix[0][0] * matrix[2][2] * matrix[3][3] + matrix[0][2] * matrix[2][3] * matrix[3][0] + matrix[0][3] * matrix[2][0] * matrix[3][2]
                          - matrix[0][0] * matrix[2][3] * matrix[3][2] - matrix[0][2] * matrix[2][0] * matrix[3][3] - matrix[0][3] * matrix[2][2] * matrix[3][0])`~op~`d,
                           (matrix[0][0] * matrix[1][3] * matrix[3][2] + matrix[0][2] * matrix[1][0] * matrix[3][3] + matrix[0][3] * matrix[1][2] * matrix[3][0]
                          - matrix[0][0] * matrix[1][2] * matrix[3][3] - matrix[0][2] * matrix[1][3] * matrix[3][0] - matrix[0][3] * matrix[1][0] * matrix[3][2])`~op~`d,
                           (matrix[0][0] * matrix[1][2] * matrix[2][3] + matrix[0][2] * matrix[1][3] * matrix[2][0] + matrix[0][3] * matrix[1][0] * matrix[2][2]
                          - matrix[0][0] * matrix[1][3] * matrix[2][2] - matrix[0][2] * matrix[1][0] * matrix[2][3] - matrix[0][3] * matrix[1][2] * matrix[2][0])`~op~`d],
                          [(matrix[1][0] * matrix[2][1] * matrix[3][3] + matrix[1][1] * matrix[2][3] * matrix[3][0] + matrix[1][3] * matrix[2][0] * matrix[3][1]
                          - matrix[1][0] * matrix[2][3] * matrix[3][1] - matrix[1][1] * matrix[2][0] * matrix[3][3] - matrix[1][3] * matrix[2][1] * matrix[3][0])`~op~`d,
                           (matrix[0][0] * matrix[2][3] * matrix[3][1] + matrix[0][1] * matrix[2][0] * matrix[3][3] + matrix[0][3] * matrix[2][1] * matrix[3][0]
                          - matrix[0][0] * matrix[2][1] * matrix[3][3] - matrix[0][1] * matrix[2][3] * matrix[3][0] - matrix[0][3] * matrix[2][0] * matrix[3][1])`~op~`d,
                           (matrix[0][0] * matrix[1][1] * matrix[3][3] + matrix[0][1] * matrix[1][3] * matrix[3][0] + matrix[0][3] * matrix[1][0] * matrix[3][1]
                          - matrix[0][0] * matrix[1][3] * matrix[3][1] - matrix[0][1] * matrix[1][0] * matrix[3][3] - matrix[0][3] * matrix[1][1] * matrix[3][0])`~op~`d,
                           (matrix[0][0] * matrix[1][3] * matrix[2][1] + matrix[0][1] * matrix[1][0] * matrix[2][3] + matrix[0][3] * matrix[1][1] * matrix[2][0]
                          - matrix[0][0] * matrix[1][1] * matrix[2][3] - matrix[0][1] * matrix[1][3] * matrix[2][0] - matrix[0][3] * matrix[1][0] * matrix[2][1])`~op~`d],
                          [(matrix[1][0] * matrix[2][2] * matrix[3][1] + matrix[1][1] * matrix[2][0] * matrix[3][2] + matrix[1][2] * matrix[2][1] * matrix[3][0]
                          - matrix[1][0] * matrix[2][1] * matrix[3][2] - matrix[1][1] * matrix[2][2] * matrix[3][0] - matrix[1][2] * matrix[2][0] * matrix[3][1])`~op~`d,
                           (matrix[0][0] * matrix[2][1] * matrix[3][2] + matrix[0][1] * matrix[2][2] * matrix[3][0] + matrix[0][2] * matrix[2][0] * matrix[3][1]
                          - matrix[0][0] * matrix[2][2] * matrix[3][1] - matrix[0][1] * matrix[2][0] * matrix[3][2] - matrix[0][2] * matrix[2][1] * matrix[3][0])`~op~`d,
                           (matrix[0][0] * matrix[1][2] * matrix[3][1] + matrix[0][1] * matrix[1][0] * matrix[3][2] + matrix[0][2] * matrix[1][1] * matrix[3][0]
                          - matrix[0][0] * matrix[1][1] * matrix[3][2] - matrix[0][1] * matrix[1][2] * matrix[3][0] - matrix[0][2] * matrix[1][0] * matrix[3][1])`~op~`d,
                           (matrix[0][0] * matrix[1][1] * matrix[2][2] + matrix[0][1] * matrix[1][2] * matrix[2][0] + matrix[0][2] * matrix[1][0] * matrix[2][1]
                          - matrix[0][0] * matrix[1][2] * matrix[2][1] - matrix[0][1] * matrix[1][0] * matrix[2][2] - matrix[0][2] * matrix[1][1] * matrix[2][0])`~op~`d]];
            `);

            return mat;
        }

        // some static fun ...
        // (1) glprogramming.com/red/appendixf.html - ortographic is broken!
        // (2) http://fly.cc.fer.hr/~unreal/theredbook/appendixg.html
        // (3) http://en.wikipedia.org/wiki/Orthographic_projection_(geometry)

        static if(isFloatingPoint!mt) {
            static private mt[6] cperspective(mt width, mt height, mt fov, mt near, mt far)
                in { assert(height != 0); }
                do {
                    mt aspect = width/height;
                    mt top = near * tan(fov*(PI/360.0));
                    mt bottom = -top;
                    mt right = top * aspect;
                    mt left = -right;

                    return [left, right, bottom, top, near, far];
                }

            /// Returns a perspective matrix (4x4 and floating-point matrices only).
            static Matrix perspective(mt width, mt height, mt fov, mt near, mt far) {
                mt[6] cdata = cperspective(width, height, fov, near, far);
                return perspective(cdata[0], cdata[1], cdata[2], cdata[3], cdata[4], cdata[5]);
            }

            /// ditto
            static Matrix perspective(mt left, mt right, mt bottom, mt top, mt near, mt far)
                in {
                    assert(right-left != 0);
                    assert(top-bottom != 0);
                    assert(far-near != 0);
                }
                do {
                    Matrix ret;
                    ret.clear(0);

                    ret.matrix[0][0] = (2*near)/(right-left);
                    ret.matrix[0][2] = (right+left)/(right-left);
                    ret.matrix[1][1] = (2*near)/(top-bottom);
                    ret.matrix[1][2] = (top+bottom)/(top-bottom);
                    ret.matrix[2][2] = -(far+near)/(far-near);
                    ret.matrix[2][3] = -(2*far*near)/(far-near);
                    ret.matrix[3][2] = -1;

                    return ret;
                }

            /// Returns an inverse perspective matrix (4x4 and floating-point matrices only).
            static Matrix persperctiveInverse(mt width, mt height, mt fov, mt near, mt far) {
                mt[6] cdata = cperspective(width, height, fov, near, far);
                return persperctiveInverse(cdata[0], cdata[1], cdata[2], cdata[3], cdata[4], cdata[5]);
            }

            /// ditto
            static Matrix persperctiveInverse(mt left, mt right, mt bottom, mt top, mt near, mt far)
                in {
                    assert(near != 0);
                    assert(far != 0);
                }
                do {
                    Matrix ret;
                    ret.clear(0);

                    ret.matrix[0][0] = (right-left)/(2*near);
                    ret.matrix[0][3] = (right+left)/(2*near);
                    ret.matrix[1][1] = (top-bottom)/(2*near);
                    ret.matrix[1][3] = (top+bottom)/(2*near);
                    ret.matrix[2][3] = -1;
                    ret.matrix[3][2] = -(far-near)/(2*far*near);
                    ret.matrix[3][3] = (far+near)/(2*far*near);

                    return ret;
                }

            // (2) and (3) say this one is correct
            /// Returns an orthographic matrix (4x4 and floating-point matrices only).
            static Matrix orthographic(mt left, mt right, mt bottom, mt top, mt near, mt far)
                in {
                    assert(right-left != 0);
                    assert(top-bottom != 0);
                    assert(far-near != 0);
                }
                do {
                    Matrix ret;
                    ret.clear(0);

                    ret.matrix[0][0] = 2/(right-left);
                    ret.matrix[0][3] = -(right+left)/(right-left);
                    ret.matrix[1][1] = 2/(top-bottom);
                    ret.matrix[1][3] = -(top+bottom)/(top-bottom);
                    ret.matrix[2][2] = -2/(far-near);
                    ret.matrix[2][3] = -(far+near)/(far-near);
                    ret.matrix[3][3] = 1;

                    return ret;
                }

            // (1) and (2) say this one is correct
            /// Returns an inverse ortographic matrix (4x4 and floating-point matrices only).
            static Matrix orthographicInverse(mt left, mt right, mt bottom, mt top, mt near, mt far) {
                Matrix ret;
                ret.clear(0);

                ret.matrix[0][0] = (right-left)/2;
                ret.matrix[0][3] = (right+left)/2;
                ret.matrix[1][1] = (top-bottom)/2;
                ret.matrix[1][3] = (top+bottom)/2;
                ret.matrix[2][2] = (far-near)/-2;
                ret.matrix[2][3] = (far+near)/2;
                ret.matrix[3][3] = 1;

                return ret;
            }

            /// Returns a look at matrix (4x4 and floating-point matrices only).
            static Matrix lookAt(Vector!(mt, 3) eye, Vector!(mt, 3) target, Vector!(mt, 3) up) {
                alias vec3mt = Vector!(mt, 3);
                vec3mt look_dir = (target - eye).normalized;
                vec3mt up_dir = up.normalized;

                vec3mt right_dir = cross(look_dir, up_dir).normalized;
                vec3mt perp_up_dir = cross(right_dir, look_dir);

                Matrix ret = Matrix.identity;
                ret.matrix[0][0..3] = right_dir.vector[];
                ret.matrix[1][0..3] = perp_up_dir.vector[];
                ret.matrix[2][0..3] = (-look_dir).vector[];

                ret.matrix[0][3] = -dot(eye, right_dir);
                ret.matrix[1][3] = -dot(eye, perp_up_dir);
                ret.matrix[2][3] = dot(eye, look_dir);

                return ret;
            }

            unittest {
                mt[6] cp = cperspective(600f, 900f, 60f, 1f, 100f);
                assert(cp[4] == 1.0f);
                assert(cp[5] == 100.0f);
                assert(cp[0] == -cp[1]);
                assert((cp[0] < -0.38489f) && (cp[0] > -0.38491f));
                assert(cp[2] == -cp[3]);
                assert((cp[2] < -0.577349f) && (cp[2] > -0.577351f));

                assert(mat4.perspective(600f, 900f, 60.0, 1.0, 100.0) == mat4.perspective(cp[0], cp[1], cp[2], cp[3], cp[4], cp[5]));
                float[4][4] m4p = mat4.perspective(600f, 900f, 60.0, 1.0, 100.0).matrix;
                assert((m4p[0][0] < 2.598077f) && (m4p[0][0] > 2.598075f));
                assert(m4p[0][2] == 0.0f);
                assert((m4p[1][1] < 1.732052) && (m4p[1][1] > 1.732050));
                assert(m4p[1][2] == 0.0f);
                assert((m4p[2][2] < -1.020201) && (m4p[2][2] > -1.020203));
                assert((m4p[2][3] < -2.020201) && (m4p[2][3] > -2.020203));
                assert((m4p[3][2] < -0.9f) && (m4p[3][2] > -1.1f));

                float[4][4] m4pi = mat4.persperctiveInverse(600f, 900f, 60.0, 1.0, 100.0).matrix;
                assert((m4pi[0][0] < 0.384901) && (m4pi[0][0] > 0.384899));
                assert(m4pi[0][3] == 0.0f);
                assert((m4pi[1][1] < 0.577351) && (m4pi[1][1] > 0.577349));
                assert(m4pi[1][3] == 0.0f);
                assert(m4pi[2][3] == -1.0f);
                assert((m4pi[3][2] < -0.494999) && (m4pi[3][2] > -0.495001));
                assert((m4pi[3][3] < 0.505001) && (m4pi[3][3] > 0.504999));

                // maybe the next tests should be improved
                float[4][4] m4o = mat4.orthographic(-1.0f, 1.0f, -1.0f, 1.0f, -1.0f, 1.0f).matrix;
                assert(m4o == [[1.0f, 0.0f, 0.0f, 0.0f],
                               [0.0f, 1.0f, 0.0f, 0.0f],
                               [0.0f, 0.0f, -1.0f, 0.0f],
                               [0.0f, 0.0f, 0.0f, 1.0f]]);

                float[4][4] m4oi = mat4.orthographicInverse(-1.0f, 1.0f, -1.0f, 1.0f, -1.0f, 1.0f).matrix;
                assert(m4oi == [[1.0f, 0.0f, 0.0f, 0.0f],
                                [0.0f, 1.0f, 0.0f, 0.0f],
                                [0.0f, 0.0f, -1.0f, 0.0f],
                                [0.0f, 0.0f, 0.0f, 1.0f]]);

                //TODO: lookAt tests
            }
        }
    }

    static if((rows == cols) && (rows >= 3) && (rows <= 4)) {
        /// Returns a translation matrix (3x3 and 4x4 matrices).
        static Matrix translation(mt x, mt y, mt z) {
            Matrix ret = Matrix.identity;

            ret.matrix[0][cols-1] = x;
            ret.matrix[1][cols-1] = y;
            ret.matrix[2][cols-1] = z;

            return ret;
        }

        /// ditto
        static Matrix translation(Vector!(mt, 3) v) {
            Matrix ret = Matrix.identity;
            
            ret.matrix[0][cols-1] = v.x;
            ret.matrix[1][cols-1] = v.y;
            ret.matrix[2][cols-1] = v.z;
            
            return ret;
        }

        /// Applys a translation on the current matrix and returns $(I this) (3x3 and 4x4 matrices).
        Matrix translate(mt x, mt y, mt z) {
            this = Matrix.translation(x, y, z) * this;
            return this;
        }

        /// ditto
        Matrix translate(Vector!(mt, 3) v) {
            this = Matrix.translation(v) * this;
            return this;
        }

        /// Returns a scaling matrix (3x3 and 4x4 matrices);
        static Matrix scaling(mt x, mt y, mt z) {
            Matrix ret = Matrix.identity;

            ret.matrix[0][0] = x;
            ret.matrix[1][1] = y;
            ret.matrix[2][2] = z;

            return ret;
        }

        /// Applys a scale to the current matrix and returns $(I this) (3x3 and 4x4 matrices).
        Matrix scale(mt x, mt y, mt z) {
            this = Matrix.scaling(x, y, z) * this;
            return this;
        }

        unittest {
            mat3 m3 = mat3.identity;
            assert(m3.translate(1.0f, 2.0f, 3.0f).matrix == mat3.translation(1.0f, 2.0f, 3.0f).matrix);
            assert(mat3.translation(1.0f, 2.0f, 3.0f).matrix == [[1.0f, 0.0f, 1.0f],
                                                                 [0.0f, 1.0f, 2.0f],
                                                                 [0.0f, 0.0f, 3.0f]]);
            assert(mat3.identity.translate(0.0f, 1.0f, 2.0f).matrix == mat3.translation(0.0f, 1.0f, 2.0f).matrix);

            mat3 m31 = mat3.identity;
            assert(m31.translate(vec3(1.0f, 2.0f, 3.0f)).matrix == mat3.translation(vec3(1.0f, 2.0f, 3.0f)).matrix);
            assert(mat3.translation(vec3(1.0f, 2.0f, 3.0f)).matrix == [[1.0f, 0.0f, 1.0f],
                    [0.0f, 1.0f, 2.0f],
                    [0.0f, 0.0f, 3.0f]]);
            assert(mat3.identity.translate(vec3(0.0f, 1.0f, 2.0f)).matrix == mat3.translation(vec3(0.0f, 1.0f, 2.0f)).matrix);

            assert(m3.scaling(0.0f, 1.0f, 2.0f).matrix == mat3.scaling(0.0f, 1.0f, 2.0f).matrix);
            assert(mat3.scaling(0.0f, 1.0f, 2.0f).matrix == [[0.0f, 0.0f, 0.0f],
                                                             [0.0f, 1.0f, 0.0f],
                                                             [0.0f, 0.0f, 2.0f]]);
            assert(mat3.identity.scale(0.0f, 1.0f, 2.0f).matrix == mat3.scaling(0.0f, 1.0f, 2.0f).matrix);

            // same tests for 4x4

            mat4 m4 = mat4(1.0f);
            assert(m4.translation(1.0f, 2.0f, 3.0f).matrix == mat4.translation(1.0f, 2.0f, 3.0f).matrix);
            assert(mat4.translation(1.0f, 2.0f, 3.0f).matrix == [[1.0f, 0.0f, 0.0f, 1.0f],
                                                                 [0.0f, 1.0f, 0.0f, 2.0f],
                                                                 [0.0f, 0.0f, 1.0f, 3.0f],
                                                                 [0.0f, 0.0f, 0.0f, 1.0f]]);
            assert(mat4.identity.translate(0.0f, 1.0f, 2.0f).matrix == mat4.translation(0.0f, 1.0f, 2.0f).matrix);

            assert(m4.scaling(0.0f, 1.0f, 2.0f).matrix == mat4.scaling(0.0f, 1.0f, 2.0f).matrix);
            assert(mat4.scaling(0.0f, 1.0f, 2.0f).matrix == [[0.0f, 0.0f, 0.0f, 0.0f],
                                                             [0.0f, 1.0f, 0.0f, 0.0f],
                                                             [0.0f, 0.0f, 2.0f, 0.0f],
                                                             [0.0f, 0.0f, 0.0f, 1.0f]]);
            assert(mat4.identity.scale(0.0f, 1.0f, 2.0f).matrix == mat4.scaling(0.0f, 1.0f, 2.0f).matrix);
        }
    }


    static if((rows == cols) && (rows >= 3)) {
        static if(isFloatingPoint!mt) {
            /// Returns an identity matrix with an applied rotateAxis around an arbitrary axis (nxn matrices, n >= 3).
            static Matrix rotation(real alpha, Vector!(mt, 3) axis) {
                Matrix mult = Matrix.identity;

                if(axis.length != 1) {
                    axis.normalize();
                }

                real cosa = cos(alpha);
                real sina = sin(alpha);

                Vector!(mt, 3) temp = (1 - cosa)*axis;

                mult.matrix[0][0] = to!mt(cosa + temp.x * axis.x);
                mult.matrix[0][1] = to!mt(       temp.x * axis.y + sina * axis.z);
                mult.matrix[0][2] = to!mt(       temp.x * axis.z - sina * axis.y);
                mult.matrix[1][0] = to!mt(       temp.y * axis.x - sina * axis.z);
                mult.matrix[1][1] = to!mt(cosa + temp.y * axis.y);
                mult.matrix[1][2] = to!mt(       temp.y * axis.z + sina * axis.x);
                mult.matrix[2][0] = to!mt(       temp.z * axis.x + sina * axis.y);
                mult.matrix[2][1] = to!mt(       temp.z * axis.y - sina * axis.x);
                mult.matrix[2][2] = to!mt(cosa + temp.z * axis.z);

                return mult;
            }

            /// ditto
            static Matrix rotation(real alpha, mt x, mt y, mt z) {
                return Matrix.rotation(alpha, Vector!(mt, 3)(x, y, z));
            }

            /// Returns an identity matrix with an applied rotation around the x-axis (nxn matrices, n >= 3).
            static Matrix xRotation(real alpha) {
                Matrix mult = Matrix.identity;

                mt cosamt = to!mt(cos(alpha));
                mt sinamt = to!mt(sin(alpha));

                mult.matrix[1][1] = cosamt;
                mult.matrix[1][2] = -sinamt;
                mult.matrix[2][1] = sinamt;
                mult.matrix[2][2] = cosamt;

                return mult;
            }

            /// Returns an identity matrix with an applied rotation around the y-axis (nxn matrices, n >= 3).
            static Matrix yRotation(real alpha) {
                Matrix mult = Matrix.identity;

                mt cosamt = to!mt(cos(alpha));
                mt sinamt = to!mt(sin(alpha));

                mult.matrix[0][0] = cosamt;
                mult.matrix[0][2] = sinamt;
                mult.matrix[2][0] = -sinamt;
                mult.matrix[2][2] = cosamt;

                return mult;
            }

            /// Returns an identity matrix with an applied rotation around the z-axis (nxn matrices, n >= 3).
            static Matrix zRotation(real alpha) {
                Matrix mult = Matrix.identity;

                mt cosamt = to!mt(cos(alpha));
                mt sinamt = to!mt(sin(alpha));

                mult.matrix[0][0] = cosamt;
                mult.matrix[0][1] = -sinamt;
                mult.matrix[1][0] = sinamt;
                mult.matrix[1][1] = cosamt;

                return mult;
            }

            Matrix rotate(real alpha, Vector!(mt, 3) axis) {
                this = rotation(alpha, axis) * this;
                return this;
            }

            /// Rotates the current matrix around the x-axis and returns $(I this) (nxn matrices, n >= 3).
            Matrix rotateX(real alpha) {
                this = xRotation(alpha) * this;
                return this;
            }

            /// Rotates the current matrix around the y-axis and returns $(I this) (nxn matrices, n >= 3).
            Matrix rotateY(real alpha) {
                this = yRotation(alpha) * this;
                return this;
            }

            /// Rotates the current matrix around the z-axis and returns $(I this) (nxn matrices, n >= 3).
            Matrix rotateZ(real alpha) {
                this = zRotation(alpha) * this;
                return this;
            }

            unittest {
                assert(mat4.xRotation(0).matrix == [[1.0f, 0.0f, 0.0f, 0.0f],
                                                    [0.0f, 1.0f, -0.0f, 0.0f],
                                                    [0.0f, 0.0f, 1.0f, 0.0f],
                                                    [0.0f, 0.0f, 0.0f, 1.0f]]);
                assert(mat4.yRotation(0).matrix == [[1.0f, 0.0f, 0.0f, 0.0f],
                                                    [0.0f, 1.0f, 0.0f, 0.0f],
                                                    [0.0f, 0.0f, 1.0f, 0.0f],
                                                    [0.0f, 0.0f, 0.0f, 1.0f]]);
                assert(mat4.zRotation(0).matrix == [[1.0f, -0.0f, 0.0f, 0.0f],
                                                    [0.0f, 1.0f, 0.0f, 0.0f],
                                                    [0.0f, 0.0f, 1.0f, 0.0f],
                                                    [0.0f, 0.0f, 0.0f, 1.0f]]);
                mat4 xro = mat4.identity;
                xro.rotateX(0);
                assert(mat4.xRotation(0).matrix == xro.matrix);
                assert(xro.matrix == mat4.identity.rotateX(0).matrix);
                assert(xro.matrix == mat4.rotation(0, vec3(1.0f, 0.0f, 0.0f)).matrix);
                mat4 yro = mat4.identity;
                yro.rotateY(0);
                assert(mat4.yRotation(0).matrix == yro.matrix);
                assert(yro.matrix == mat4.identity.rotateY(0).matrix);
                assert(yro.matrix == mat4.rotation(0, vec3(0.0f, 1.0f, 0.0f)).matrix);
                mat4 zro = mat4.identity;
                xro.rotateZ(0);
                assert(mat4.zRotation(0).matrix == zro.matrix);
                assert(zro.matrix == mat4.identity.rotateZ(0).matrix);
                assert(zro.matrix == mat4.rotation(0, vec3(0.0f, 0.0f, 1.0f)).matrix);
            }
        } // isFloatingPoint


        /// Sets the translation of the matrix (nxn matrices, n >= 3).
        void translation(mt[] values...) // intended to be a property
            in { assert(values.length >= (rows-1)); }
            do {
                foreach(r; TupleRange!(0, rows-1)) {
                    matrix[r][rows-1] = values[r];
                }
            }

        /// Copyies the translation from mat to the current matrix (nxn matrices, n >= 3).
        void translation(Matrix mat) {
            foreach(r; TupleRange!(0, rows-1)) {
                matrix[r][rows-1] = mat.matrix[r][rows-1];
            }
        }

        /// Returns an identity matrix with the current translation applied (nxn matrices, n >= 3)..
        Matrix translation() {
            Matrix ret = Matrix.identity;

            foreach(r; TupleRange!(0, rows-1)) {
                ret.matrix[r][rows-1] = matrix[r][rows-1];
            }

            return ret;
        }

        unittest {
            mat3 m3 = mat3(0.0f, 1.0f, 2.0f,
                           3.0f, 4.0f, 5.0f,
                           6.0f, 7.0f, 1.0f);
            assert(m3.translation().matrix == [[1.0f, 0.0f, 2.0f], [0.0f, 1.0f, 5.0f], [0.0f, 0.0f, 1.0f]]);
            m3.translation(mat3.identity);
            assert(mat3.identity.matrix == m3.translation().matrix);
            m3.translation([2.0f, 5.0f]);
            assert(m3.translation().matrix == [[1.0f, 0.0f, 2.0f], [0.0f, 1.0f, 5.0f], [0.0f, 0.0f, 1.0f]]);
            assert(mat3.identity.matrix == mat3.identity.translation().matrix);

            mat4 m4 = mat4(0.0f, 1.0f, 2.0f, 3.0f,
                           4.0f, 5.0f, 6.0f, 7.0f,
                           8.0f, 9.0f, 10.0f, 11.0f,
                           12.0f, 13.0f, 14.0f, 1.0f);
            assert(m4.translation().matrix == [[1.0f, 0.0f, 0.0f, 3.0f],
                                       [0.0f, 1.0f, 0.0f, 7.0f],
                                       [0.0f, 0.0f, 1.0f, 11.0f],
                                       [0.0f, 0.0f, 0.0f, 1.0f]]);
            m4.translation(mat4.identity);
            assert(mat4.identity.matrix == m4.translation().matrix);
            m4.translation([3.0f, 7.0f, 11.0f]);
            assert(m4.translation().matrix == [[1.0f, 0.0f, 0.0f, 3.0f],
                                       [0.0f, 1.0f, 0.0f, 7.0f],
                                       [0.0f, 0.0f, 1.0f, 11.0f],
                                       [0.0f, 0.0f, 0.0f, 1.0f]]);
            assert(mat4.identity.matrix == mat4.identity.translation().matrix);
        }

        /// Sets the scale of the matrix (nxn matrices, n >= 3).
        void scale(mt[] values...)
            in { assert(values.length >= (rows-1)); }
            do {
                foreach(r; TupleRange!(0, rows-1)) {
                    matrix[r][r] = values[r];
                }
            }

        /// Copyies the scale from mat to the current matrix (nxn matrices, n >= 3).
        void scale(Matrix mat) {
            foreach(r; TupleRange!(0, rows-1)) {
                matrix[r][r] = mat.matrix[r][r];
            }
        }

        /// Returns an identity matrix with the current scale applied (nxn matrices, n >= 3).
        Matrix scale() {
            Matrix ret = Matrix.identity;

            foreach(r; TupleRange!(0, rows-1)) {
                ret.matrix[r][r] = matrix[r][r];
            }

            return ret;
        }

        unittest {
            mat3 m3 = mat3(0.0f, 1.0f, 2.0f,
                           3.0f, 4.0f, 5.0f,
                           6.0f, 7.0f, 1.0f);
            assert(m3.scale().matrix == [[0.0f, 0.0f, 0.0f], [0.0f, 4.0f, 0.0f], [0.0f, 0.0f, 1.0f]]);
            m3.scale(mat3.identity);
            assert(mat3.identity.matrix == m3.scale().matrix);
            m3.scale([0.0f, 4.0f]);
            assert(m3.scale().matrix == [[0.0f, 0.0f, 0.0f], [0.0f, 4.0f, 0.0f], [0.0f, 0.0f, 1.0f]]);
            assert(mat3.identity.matrix == mat3.identity.scale().matrix);

            mat4 m4 = mat4(0.0f, 1.0f, 2.0f, 3.0f,
                           4.0f, 5.0f, 6.0f, 7.0f,
                           8.0f, 9.0f, 10.0f, 11.0f,
                           12.0f, 13.0f, 14.0f, 1.0f);
            assert(m4.scale().matrix == [[0.0f, 0.0f, 0.0f, 0.0f],
                                       [0.0f, 5.0f, 0.0f, 0.0f],
                                       [0.0f, 0.0f, 10.0f, 0.0f],
                                       [0.0f, 0.0f, 0.0f, 1.0f]]);
            m4.scale(mat4.identity);
            assert(mat4.identity.matrix == m4.scale().matrix);
            m4.scale([0.0f, 5.0f, 10.0f]);
            assert(m4.scale().matrix == [[0.0f, 0.0f, 0.0f, 0.0f],
                                       [0.0f, 5.0f, 0.0f, 0.0f],
                                       [0.0f, 0.0f, 10.0f, 0.0f],
                                       [0.0f, 0.0f, 0.0f, 1.0f]]);
            assert(mat4.identity.matrix == mat4.identity.scale().matrix);
        }

        /// Copies rot into the upper left corner, the translation (nxn matrices, n >= 3).
        void rotation(Matrix!(mt, 3, 3) rot) {
            foreach(r; TupleRange!(0, 3)) {
                foreach(c; TupleRange!(0, 3)) {
                    matrix[r][c] = rot[r][c];
                }
            }
        }

        /// Returns an identity matrix with the current rotation applied (nxn matrices, n >= 3).
        Matrix!(mt, 3, 3) rotation() {
            Matrix!(mt, 3, 3) ret = Matrix!(mt, 3, 3).identity;

            foreach(r; TupleRange!(0, 3)) {
                foreach(c; TupleRange!(0, 3)) {
                    ret.matrix[r][c] = matrix[r][c];
                }
            }

            return ret;
        }

        unittest {
            mat3 m3 = mat3(0.0f, 1.0f, 2.0f,
                           3.0f, 4.0f, 5.0f,
                           6.0f, 7.0f, 1.0f);
            assert(m3.rotation().matrix == [[0.0f, 1.0f, 2.0f], [3.0f, 4.0f, 5.0f], [6.0f, 7.0f, 1.0f]]);
            m3.rotation(mat3.identity);
            assert(mat3.identity.matrix == m3.rotation().matrix);
            m3.rotation(mat3(0.0f, 1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 1.0f));
            assert(m3.rotation().matrix == [[0.0f, 1.0f, 2.0f], [3.0f, 4.0f, 5.0f], [6.0f, 7.0f, 1.0f]]);
            assert(mat3.identity.matrix == mat3.identity.rotation().matrix);

            mat4 m4 = mat4(0.0f, 1.0f, 2.0f, 3.0f,
                           4.0f, 5.0f, 6.0f, 7.0f,
                           8.0f, 9.0f, 10.0f, 11.0f,
                           12.0f, 13.0f, 14.0f, 1.0f);
            assert(m4.rotation().matrix == [[0.0f, 1.0f, 2.0f], [4.0f, 5.0f, 6.0f], [8.0f, 9.0f, 10.0f]]);
            m4.rotation(mat3.identity);
            assert(mat3.identity.matrix == m4.rotation().matrix);
            m4.rotation(mat3(0.0f, 1.0f, 2.0f, 4.0f, 5.0f, 6.0f, 8.0f, 9.0f, 10.0f));
            assert(m4.rotation().matrix == [[0.0f, 1.0f, 2.0f], [4.0f, 5.0f, 6.0f], [8.0f, 9.0f, 10.0f]]);
            assert(mat3.identity.matrix == mat4.identity.rotation().matrix);
        }

    }

    static if((rows == cols) && (rows >= 2) && (rows <= 4)) {
        /// Returns an inverted copy of the current matrix (nxn matrices, 2 >= n <= 4).
        Matrix inverse() const {
            Matrix mat;
            invert(mat);
            return mat;
        }

        /// Inverts the current matrix (nxn matrices, 2 >= n <= 4).
        void invert() {
            // workaround Issue #11238
            // uses a temporary instead of invert(this)
            Matrix temp;
            invert(temp);
            this.matrix = temp.matrix;
        }
    }

    unittest {
        mat2 m2 = mat2(1.0f, 2.0f, vec2(3.0f, 4.0f));
        assert(m2.det == -2.0f);
        assert(m2.inverse.matrix == [[-2.0f, 1.0f], [1.5f, -0.5f]]);

        mat3 m3 = mat3(1.0f, -2.0f, 3.0f,
                       7.0f, -1.0f, 0.0f,
                       3.0f, 2.0f, -4.0f);
        assert(m3.det == -1.0f);
        assert(m3.inverse.matrix == [[-4.0f, 2.0f, -3.0f],
                                     [-28.0f, 13.0f, -21.0f],
                                     [-17.0f, 8.0f, -13.0f]]);

        mat4 m4 = mat4(1.0f, 2.0f, 3.0f, 4.0f,
                       -2.0f, 1.0f, 5.0f, -2.0f,
                       2.0f, -1.0f, 7.0f, 1.0f,
                       3.0f, -3.0f, 2.0f, 0.0f);
        assert(m4.det == -8.0f);
        assert(m4.inverse.matrix == [[6.875f, 7.875f, -11.75f, 11.125f],
                                     [6.625f, 7.625f, -11.25f, 10.375f],
                                     [-0.375f, -0.375f, 0.75f, -0.625f],
                                     [-4.5f, -5.5f, 8.0f, -7.5f]]);
    }

    private void mms(mt inp, ref Matrix mat) const { // mat * scalar
        for(int r = 0; r < rows; r++) {
            for(int c = 0; c < cols; c++) {
                mat.matrix[r][c] = matrix[r][c] * inp;
            }
        }
    }

    private void masm(string op)(Matrix inp, ref Matrix mat) const { // mat + or - mat
        foreach(r; TupleRange!(0, rows)) {
            foreach(c; TupleRange!(0, cols)) {
                mat.matrix[r][c] = mixin("inp.matrix[r][c]" ~ op ~ "matrix[r][c]");
            }
        }
    }

    Matrix!(mt, rows, T.cols) opBinary(string op : "*", T)(T inp) const if(isCompatibleMatrix!T && (T.rows == cols)) {
        Matrix!(mt, rows, T.cols) ret;

        foreach(r; TupleRange!(0, rows)) {
            foreach(c; TupleRange!(0, T.cols)) {
                ret.matrix[r][c] = 0;

                foreach(c2; TupleRange!(0, cols)) {
                    ret.matrix[r][c] += matrix[r][c2] * inp.matrix[c2][c];
                }
            }
        }

        return ret;
    }

    Vector!(mt, rows) opBinary(string op : "*", T : Vector!(mt, cols))(T inp) const {
        Vector!(mt, rows) ret;
        ret.clear(0);

        foreach(c; TupleRange!(0, cols)) {
            foreach(r; TupleRange!(0, rows)) {
                ret.vector[r] += matrix[r][c] * inp.vector[c];
            }
        }

        return ret;
    }

    Matrix opBinary(string op : "*")(mt inp) const {
        Matrix ret;
        mms(inp, ret);
        return ret;
    }

    Matrix opBinaryRight(string op : "*")(mt inp) const {
        return this.opBinary!(op)(inp);
    }

    Matrix opBinary(string op)(Matrix inp) const if((op == "+") || (op == "-")) {
        Matrix ret;
        masm!(op)(inp, ret);
        return ret;
    }

    void opOpAssign(string op : "*")(mt inp) {
        mms(inp, this);
    }

    void opOpAssign(string op)(Matrix inp) if((op == "+") || (op == "-")) {
        masm!(op)(inp, this);
    }

    void opOpAssign(string op)(Matrix inp) if(op == "*") {
        this = this * inp;
    }

    unittest {
        mat2 m2 = mat2(1.0f, 2.0f, 3.0f, 4.0f);
        vec2 v2 = vec2(2.0f, 2.0f);
        assert((m2*2).matrix == [[2.0f, 4.0f], [6.0f, 8.0f]]);
        assert((2*m2).matrix == (m2*2).matrix);
        m2 *= 2;
        assert(m2.matrix == [[2.0f, 4.0f], [6.0f, 8.0f]]);
        assert((m2*v2).vector == [12.0f, 28.0f]);
        assert((v2*m2).vector == [16.0f, 24.0f]);
        assert((m2*m2).matrix == [[28.0f, 40.0f], [60.0f, 88.0f]]);
        assert((m2-m2).matrix == [[0.0f, 0.0f], [0.0f, 0.0f]]);
        assert((m2+m2).matrix == [[4.0f, 8.0f], [12.0f, 16.0f]]);
        m2 += m2;
        assert(m2.matrix == [[4.0f, 8.0f], [12.0f, 16.0f]]);
        m2 -= m2;
        assert(m2.matrix == [[0.0f, 0.0f], [0.0f, 0.0f]]);

        mat3 m3 = mat3(1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f, 9.0f);
        vec3 v3 = vec3(2.0f, 2.0f, 2.0f);
        assert((m3*2).matrix == [[2.0f, 4.0f, 6.0f], [8.0f, 10.0f, 12.0f], [14.0f, 16.0f, 18.0f]]);
        assert((2*m3).matrix == (m3*2).matrix);
        m3 *= 2;
        assert(m3.matrix == [[2.0f, 4.0f, 6.0f], [8.0f, 10.0f, 12.0f], [14.0f, 16.0f, 18.0f]]);
        assert((m3*v3).vector == [24.0f, 60.0f, 96.0f]);
        assert((v3*m3).vector == [48.0f, 60.0f, 72.0f]);
        assert((m3*m3).matrix == [[120.0f, 144.0f, 168.0f], [264.0f, 324.0f, 384.0f], [408.0f, 504.0f, 600.0f]]);
        assert((m3-m3).matrix == [[0.0f, 0.0f, 0.0f], [0.0f, 0.0f, 0.0f], [0.0f, 0.0f, 0.0f]]);
        assert((m3+m3).matrix == [[4.0f, 8.0f, 12.0f], [16.0f, 20.0f, 24.0f], [28.0f, 32.0f, 36.0f]]);
        m3 += m3;
        assert(m3.matrix == [[4.0f, 8.0f, 12.0f], [16.0f, 20.0f, 24.0f], [28.0f, 32.0f, 36.0f]]);
        m3 -= m3;
        assert(m3.matrix == [[0.0f, 0.0f, 0.0f], [0.0f, 0.0f, 0.0f], [0.0f, 0.0f, 0.0f]]);

        // test opOpAssign for matrix multiplication
        auto m4 = mat4.translation(0,1,2);
        m4 *= mat4.translation(0,-1,2);
        assert(m4 == mat4.translation(0,0,4));

        //TODO: tests for mat4, mat34
    }

    unittest {
        assert(mat2(1.0f, 2.0f, 1.0f, 1.0f) == mat2(1.0f, 2.0f, 1.0f, 1.0f));
        assert(mat2(1.0f, 2.0f, 1.0f, 1.0f) != mat2(1.0f, 1.0f, 1.0f, 1.0f));

        assert(mat3(1.0f) == mat3(1.0f));
        assert(mat3(1.0f) != mat3(2.0f));

        assert(mat4(1.0f) == mat4(1.0f));
        assert(mat4(1.0f) != mat4(2.0f));

        assert(!mat4(float.nan).isFinite);
        if(mat4(1.0f).isFinite) { }
        else { assert(false); }
    }
}

/// Pre-defined matrix types, the first number represents the number of rows
/// and the second the number of columns, if there's just one it's a nxn matrix.
/// All of these matrices are floating-point matrices.
alias mat2 = Matrix!(float, 2, 2);
alias mat3 = Matrix!(float, 3, 3);
alias mat34 = Matrix!(float, 3, 4);
alias mat4 = Matrix!(float, 4, 4);

private unittest {
    Matrix!(float,  1, 1) A = 1;
    Matrix!(double, 1, 1) B = 1;
    Matrix!(real,   1, 1) C = 1;
    Matrix!(int,    1, 1) D = 1;
    Matrix!(float,  5, 1) E = 1;
    Matrix!(double, 5, 1) F = 1;
    Matrix!(real,   5, 1) G = 1;
    Matrix!(int,    5, 1) H = 1;
    Matrix!(float,  1, 5) I = 1;
    Matrix!(double, 1, 5) J = 1;
    Matrix!(real,   1, 5) K = 1;
    Matrix!(int,    1, 5) L = 1;
}

/// Base template for all quaternion-types.
/// Params:
///  type = all values get stored as this type
struct Quaternion(type) {
    alias qt = type; /// Holds the internal type of the quaternion.

    qt[4] quaternion; /// Holds the w, x, y and z coordinates.

    /// Returns a pointer to the quaternion in memory, it starts with the w coordinate.
    auto ptr() const { return quaternion.ptr; }

    /// Returns the current vector formatted as string, useful for printing the quaternion.
    string toString() const {
        return format("%s", quaternion);
    }

    /// Gets a hash of this item
    size_t toHash() const { return typeid(this).getHash(&this); }

@safe pure nothrow:
    qt get_(char coord)() const {
        return quaternion[coordToIndex!coord];
    }
    void set_(char coord)(qt value) {
        quaternion[coordToIndex!coord] = value;
    }

    alias w = get_!'w'; /// static properties to access the values.
    alias w = set_!'w';
    alias x = get_!'x'; /// ditto
    alias x = set_!'x';
    alias y = get_!'y'; /// ditto
    alias y = set_!'y';
    alias z = get_!'z'; /// ditto
    alias z = set_!'z';

    /// Constructs the quaternion.
    /// Takes a 4-dimensional vector, where vector.x = the quaternions w coordinate,
    /// or a w coordinate of type $(I qt) and a 3-dimensional vector representing the imaginary part,
    /// or 4 values of type $(I qt).
    this(qt w_, qt x_, qt y_, qt z_) {
        w = w_;
        x = x_;
        y = y_;
        z = z_;
    }

    /// ditto
    this(qt w_, Vector!(qt, 3) vec) {
        w = w_;
        quaternion[1..4] = vec.vector[];
    }

    /// ditto
    this(Vector!(qt, 4) vec) {
        quaternion[] = vec.vector[];
    }

    /// Returns true if all values are not nan and finite, otherwise false.
    bool isFinite() const {
        foreach(q; quaternion) {
            if(isNaN(q) || isInfinity(q)) {
                return false;
            }
        }
        return true;
    }
    deprecated("Use isFinite instead of ok") alias ok = isFinite;

    unittest {
        quat q1 = quat(0.0f, 0.0f, 0.0f, 1.0f);
        assert(q1.quaternion == [0.0f, 0.0f, 0.0f, 1.0f]);
        assert(q1.quaternion == quat(0.0f, 0.0f, 0.0f, 1.0f).quaternion);
        assert(q1.quaternion == quat(0.0f, vec3(0.0f, 0.0f, 1.0f)).quaternion);
        assert(q1.quaternion == quat(vec4(0.0f, 0.0f, 0.0f, 1.0f)).quaternion);

        assert(q1.isFinite);
        q1.x = float.infinity;
        assert(!q1.isFinite);
        q1.x = float.nan;
        assert(!q1.isFinite);
        q1.x = 0.0f;
        assert(q1.isFinite);
    }

    template coordToIndex(char c) {
        static if(c == 'w') {
            enum coordToIndex = 0;
        } else static if(c == 'x') {
            enum coordToIndex = 1;
        } else static if(c == 'y') {
            enum coordToIndex = 2;
        } else static if(c == 'z') {
            enum coordToIndex = 3;
        } else {
            static assert(false, "accepted coordinates are x, y, z and w not " ~ c ~ ".");
        }
    }

    /// Returns the squared magnitude of the quaternion.
    real magnitudeSquared() const {
        return to!real(w^^2 + x^^2 + y^^2 + z^^2);
    }

    /// Returns the magnitude of the quaternion.
    real magnitude() const {
        return sqrt(magnitudeSquared);
    }

    /// Returns an identity quaternion (w=1, x=0, y=0, z=0).
    static Quaternion identity() {
        return Quaternion(1, 0, 0, 0);
    }

    /// Makes the current quaternion an identity quaternion.
    void makeIdentity() {
        w = 1;
        x = 0;
        y = 0;
        z = 0;
    }

    /// Inverts the quaternion.
    void invert() {
        x = -x;
        y = -y;
        z = -z;
    }
    alias conjugate = invert; /// ditto

    /// Returns an inverted copy of the current quaternion.
    Quaternion inverse() const {
        return Quaternion(w, -x, -y, -z);
    }
    alias conjugated = inverse; /// ditto

    unittest {
        quat q1 = quat(1.0f, 1.0f, 1.0f, 1.0f);

        assert(q1.magnitude == 2.0f);
        assert(q1.magnitudeSquared == 4.0f);
        assert(q1.magnitude == quat(0.0f, 0.0f, 2.0f, 0.0f).magnitude);

        quat q2 = quat.identity;
        assert(q2.quaternion == [1.0f, 0.0f, 0.0f, 0.0f]);
        assert(q2.x == 0.0f);
        assert(q2.y == 0.0f);
        assert(q2.z == 0.0f);
        assert(q2.w == 1.0f);

        assert(q1.inverse.quaternion == [1.0f, -1.0f, -1.0f, -1.0f]);
        q1.invert();
        assert(q1.quaternion == [1.0f, -1.0f, -1.0f, -1.0f]);

        q1.makeIdentity();
        assert(q1.quaternion == q2.quaternion);

    }

    /// Creates a quaternion from a 3x3 matrix.
    /// Params:
    ///  matrix = 3x3 matrix (rotation)
    /// Returns: A quaternion representing the rotation (3x3 matrix)
    static Quaternion fromMatrix(Matrix!(qt, 3, 3) matrix) {
        Quaternion ret;

        auto mat = matrix.matrix;
        qt trace = mat[0][0] + mat[1][1] + mat[2][2];

        if(trace > 0) {
            real s = 0.5 / sqrt(trace + 1.0f);

            ret.w = to!qt(0.25 / s);
            ret.x = to!qt((mat[2][1] - mat[1][2]) * s);
            ret.y = to!qt((mat[0][2] - mat[2][0]) * s);
            ret.z = to!qt((mat[1][0] - mat[0][1]) * s);
        } else if((mat[0][0] > mat[1][1]) && (mat[0][0] > mat[2][2])) {
            real s = 2.0 * sqrt(1.0 + mat[0][0] - mat[1][1] - mat[2][2]);

            ret.w = to!qt((mat[2][1] - mat[1][2]) / s);
            ret.x = to!qt(0.25f * s);
            ret.y = to!qt((mat[0][1] + mat[1][0]) / s);
            ret.z = to!qt((mat[0][2] + mat[2][0]) / s);
        } else if(mat[1][1] > mat[2][2]) {
            real s = 2.0 * sqrt(1 + mat[1][1] - mat[0][0] - mat[2][2]);

            ret.w = to!qt((mat[0][2] - mat[2][0]) / s);
            ret.x = to!qt((mat[0][1] + mat[1][0]) / s);
            ret.y = to!qt(0.25f * s);
            ret.z = to!qt((mat[1][2] + mat[2][1]) / s);
        } else {
            real s = 2.0 * sqrt(1 + mat[2][2] - mat[0][0] - mat[1][1]);

            ret.w = to!qt((mat[1][0] - mat[0][1]) / s);
            ret.x = to!qt((mat[0][2] + mat[2][0]) / s);
            ret.y = to!qt((mat[1][2] + mat[2][1]) / s);
            ret.z = to!qt(0.25f * s);
        }

        return ret;
    }

    /// Returns the quaternion as matrix.
    /// Params:
    ///  rows = number of rows of the resulting matrix (min 3)
    ///  cols = number of columns of the resulting matrix (min 3)
    Matrix!(qt, rows, cols) toMatrix(int rows, int cols)() const if((rows >= 3) && (cols >= 3)) {
        static if((rows == 3) && (cols == 3)) {
            Matrix!(qt, rows, cols) ret;
        } else {
            Matrix!(qt, rows, cols) ret = Matrix!(qt, rows, cols).identity;
        }

        qt xx = x^^2;
        qt xy = x * y;
        qt xz = x * z;
        qt xw = x * w;
        qt yy = y^^2;
        qt yz = y * z;
        qt yw = y * w;
        qt zz = z^^2;
        qt zw = z * w;

        ret.matrix[0][0] = 1 - 2 * (yy + zz);
        ret.matrix[0][1] = 2 * (xy - zw);
        ret.matrix[0][2] = 2 * (xz + yw);

        ret.matrix[1][0] = 2 * (xy + zw);
        ret.matrix[1][1] = 1 - 2 * (xx + zz);
        ret.matrix[1][2] = 2 * (yz - xw);

        ret.matrix[2][0] = 2 * (xz - yw);
        ret.matrix[2][1] = 2 * (yz + xw);
        ret.matrix[2][2] = 1 - 2 * (xx + yy);

        return ret;
    }

    unittest {
        quat q1 = quat(4.0f, 1.0f, 2.0f, 3.0f);

        assert(q1.toMatrix!(3, 3).matrix == [[-25.0f, -20.0f, 22.0f], [28.0f, -19.0f, 4.0f], [-10.0f, 20.0f, -9.0f]]);
        assert(q1.toMatrix!(4, 4).matrix == [[-25.0f, -20.0f, 22.0f, 0.0f],
                                              [28.0f, -19.0f, 4.0f, 0.0f],
                                              [-10.0f, 20.0f, -9.0f, 0.0f],
                                              [0.0f, 0.0f, 0.0f, 1.0f]]);
        assert(quat.identity.toMatrix!(3, 3).matrix == Matrix!(qt, 3, 3).identity.matrix);
        assert(q1.quaternion == quat.fromMatrix(q1.toMatrix!(3, 3)).quaternion);

        assert(quat(1.0f, 0.0f, 0.0f, 0.0f).quaternion == quat.fromMatrix(mat3.identity).quaternion);

        quat q2 = quat.fromMatrix(mat3(1.0f, 3.0f, 2.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f));
        assert(q2.x == 0.0f);
        assert(almostEqual(q2.y, 0.7071067f));
        assert(almostEqual(q2.z, -1.060660));
        assert(almostEqual(q2.w, 0.7071067f));
    }
    
    /// Returns the quaternion as a vec3 (axis / angle representation).
    vec3 toAxisAngle() {
        vec3 ret;
        quat this_normalized = this.normalized();
        real angle = 2 * acos(this_normalized.w);
        real denominator = sqrt(1.0 - (this_normalized.w)^^2);

        if (almostEqual(denominator, 0)) { // Avoid divide by 0
            ret.x = 1;
            ret.y = ret.z = 0;
        }
        else {
            ret.x = x / denominator;
            ret.y = y / denominator;
            ret.z = z / denominator;
        }

        return ret * angle;
    }
    unittest {
        // See https://www.energid.com/resources/orientation-calculator

        void testPair(quat q, vec3 v) {
            vec3 q2v = q.toAxisAngle();
            assert(almostEqual(q2v.x, q2v.x) // epsilon = 0.000001f
                   && almostEqual(q2v.y, q2v.y)
                   && almostEqual(q2v.z, q2v.z),
                   "to_axisAngle does not yield a correct vector.");
        }

        quat q1 = quat(0.5, 0.5, 0.5, 0.5);
        vec3 v1 = vec3(1.2091996);
        testPair(q1, v1);

        q1 = quat(0.1825742, 0.3651484, 0.5477226, 0.7302967);
        v1 = vec3(1.030380653, 1.54557098, 2.060761024);
        testPair(q1, v1);

        q1 = quat(0, 0, 0, 1);
        v1 = vec3(0, 0, 0);
        testPair(q1, v1);

        q1 = quat(-0.1961161, 0.4358136, 0.8716273, 0.1089534);
        v1 = vec3(-1.220800604, -2.441601488, -0.305200151);
        testPair(q1, v1);
    }

    /// Normalizes the current quaternion.
    void normalize() {
        qt m = to!qt(magnitude);

        if(m != 0) {
            w = w / m;
            x = x / m;
            y = y / m;
            z = z / m;
        }
    }

    /// Returns a normalized copy of the current quaternion.
    Quaternion normalized() const {
        Quaternion ret;
        qt m = to!qt(magnitude);

        if(m != 0) {
            ret.w = w / m;
            ret.x = x / m;
            ret.y = y / m;
            ret.z = z / m;
        } else {
            ret = Quaternion(w, x, y, z);
        }

        return ret;
    }

    unittest {
        quat q1 = quat(1.0f, 2.0f, 3.0f, 4.0f);
        quat q2 = quat(1.0f, 2.0f, 3.0f, 4.0f);

        q1.normalize();
        assert(q1.quaternion == q2.normalized.quaternion);
        //assert(q1.quaternion == q1.normalized.quaternion);
        assert(almostEqual(q1.magnitude, 1.0));
    }

    /// Returns the yaw.
    real yaw() const {
        return atan2(to!real(2.0*(w*z + x*y)), to!real(1.0 - 2.0*(y*y + z*z)));
    }

    /// Returns the pitch.
    real pitch() const {
        return asin(to!real(2.0*(w*y - z*x)));
    }

    /// Returns the roll.
    real roll() const {
        return atan2(to!real(2.0*(w*x + y*z)), to!real(1.0 - 2.0*(x*x + y*y)));
    }

    unittest {
        quat q1 = quat.identity;
        assert(q1.pitch == 0.0f);
        assert(q1.yaw == 0.0f);
        assert(q1.roll == 0.0f);

        quat q2 = quat(1.0f, 1.0f, 1.0f, 1.0f);
        assert(almostEqual(q2.yaw, q2.roll));
        //assert(almostEqual(q2.yaw, 1.570796f));
        assert(q2.pitch == 0.0f);

        quat q3 = quat(0.1f, 1.9f, 2.1f, 1.3f);
        //assert(almostEqual(q3.yaw, 2.4382f));
        assert(isNaN(q3.pitch));
        //assert(almostEqual(q3.roll, 1.67719f));
    }

    /// Returns a quaternion with applied rotation around the x-axis.
    static Quaternion xRotation(real alpha) {
        Quaternion ret;

        alpha /= 2;
        ret.w = to!qt(cos(alpha));
        ret.x = to!qt(sin(alpha));
        ret.y = 0;
        ret.z = 0;

        return ret;
    }

    /// Returns a quaternion with applied rotation around the y-axis.
    static Quaternion yRotation(real alpha) {
        Quaternion ret;

        alpha /= 2;
        ret.w = to!qt(cos(alpha));
        ret.x = 0;
        ret.y = to!qt(sin(alpha));
        ret.z = 0;

        return ret;
    }

    /// Returns a quaternion with applied rotation around the z-axis.
    static Quaternion zRotation(real alpha) {
        Quaternion ret;

        alpha /= 2;
        ret.w = to!qt(cos(alpha));
        ret.x = 0;
        ret.y = 0;
        ret.z = to!qt(sin(alpha));

        return ret;
    }

    /// Returns a quaternion with applied rotation around an axis.
    static Quaternion axisRotation(real alpha, Vector!(qt, 3) axis) {
        if(alpha == 0) {
            return Quaternion.identity;
        }
        Quaternion ret;

        alpha /= 2;
        qt sinaqt = to!qt(sin(alpha));

        ret.w = to!qt(cos(alpha));
        ret.x = axis.x * sinaqt;
        ret.y = axis.y * sinaqt;
        ret.z = axis.z * sinaqt;

        return ret;
    }

    /// Creates a quaternion from an euler rotation.
    static Quaternion eulerRotation(real roll, real pitch, real yaw) {
        Quaternion ret;

        auto cr = cos(roll / 2.0);
        auto cp = cos(pitch / 2.0);
        auto cy = cos(yaw / 2.0);
        auto sr = sin(roll / 2.0);
        auto sp = sin(pitch / 2.0);
        auto sy = sin(yaw / 2.0);

        ret.w = cr * cp * cy + sr * sp * sy;
        ret.x = sr * cp * cy - cr * sp * sy;
        ret.y = cr * sp * cy + sr * cp * sy;
        ret.z = cr * cp * sy - sr * sp * cy;

        return ret;
    }

    @("quat eulerRotation")
    unittest {
        enum startPitch = 0.1;
        enum startYaw = -0.2;
        enum startRoll = 0.6;
        
        auto q = quat.eulerRotation(startRoll,startPitch,startYaw);
        
        assert(almostEqual(q.pitch,startPitch));
        assert(almostEqual(q.yaw,startYaw));
        assert(almostEqual(q.roll,startRoll));
    }

    /**
        Makes a quaternion with the specified forward and upwards directions in a Unity compatible manner.
        Implementation from http://answers.unity3d.com/questions/467614/what-is-the-source-code-of-quaternionlookrotation.html
    */
    static Quaternion lookRotation(vec3 forward, vec3 up) {
        forward = forward.normalized;
        vec3 right = cross(up, forward).normalized;
        up = cross(forward, right);

        qt m00 = right.x;
        qt m01 = right.y;
        qt m02 = right.z;
        qt m10 = up.x;
        qt m11 = up.y;
        qt m12 = up.z;
        qt m20 = forward.x;
        qt m21 = forward.y;
        qt m22 = forward.z;

        qt num8 = (m00 + m11) + m22;
        Quaternion!qt quat = Quaternion!qt.identity;
        if (num8 > 0) {
            qt num = sqrt(num8 + 1.0);
            quat.w = num * 0.5f;
            num = 0.5f/num;
			quat.x = (m12 - m21) * num;
			quat.y = (m20 - m02) * num;
			quat.z = (m01 - m10) * num;
            return quat;
        }
        if ((m00 >= m11) && (m00 >= m22)) {
            qt num7 = sqrt(((1f + m00) - m11) - m22);
			qt num4 = 0.5f / num7;
			quat.x = 0.5f * num7;
			quat.y = (m01 + m10) * num4;
			quat.z = (m02 + m20) * num4;
			quat.w = (m12 - m21) * num4;
			return quat;
        }
		if (m11 > m22)
		{
			qt num6 = sqrt(((1f + m11) - m00) - m22);
			qt num3 = 0.5f / num6;
			quat.x = (m10 + m01) * num3;
			quat.y = 0.5f * num6;
			quat.z = (m21 + m12) * num3;
			quat.w = (m20 - m02) * num3;
			return quat;
		}
		qt num5 = sqrt(((1f + m22) - m00) - m11);
		qt num2 = 0.5f / num5;
		quat.x = (m20 + m02) * num2;
		quat.y = (m21 + m12) * num2;
		quat.z = 0.5f * num5;
		quat.w = (m01 - m10) * num2;
		return quat;
    }

    @("quat lookRotation")
    unittest {
        // TODO: add unittest
    }

    /// Rotates the current quaternion around the x-axis and returns $(I this).
    Quaternion rotateX(real alpha) {
        this = xRotation(alpha) * this;
        return this;
    }

    /// Rotates the current quaternion around the y-axis and returns $(I this).
    Quaternion rotateY(real alpha) {
        this = yRotation(alpha) * this;
        return this;
    }

    /// Rotates the current quaternion around the z-axis and returns $(I this).
    Quaternion rotateZ(real alpha) {
        this = zRotation(alpha) * this;
        return this;
    }

    /// Rotates the current quaternion around an axis and returns $(I this).
    Quaternion rotateAxis(real alpha, Vector!(qt, 3) axis) {
        this = axisRotation(alpha, axis) * this;
        return this;
    }

    /// Applies an euler rotation to the current quaternion and returns $(I this).
    Quaternion rotateEuler(real heading, real attitude, real bank) {
        this = eulerRotation(heading, attitude, bank) * this;
        return this;
    }

    @("quat x/y/z/axis rotation")
    unittest {
        assert(quat.xRotation(PI).quaternion[1..4] == [1.0f, 0.0f, 0.0f]);
        assert(quat.yRotation(PI).quaternion[1..4] == [0.0f, 1.0f, 0.0f]);
        assert(quat.zRotation(PI).quaternion[1..4] == [0.0f, 0.0f, 1.0f]);
        assert((quat.xRotation(PI).w == quat.yRotation(PI).w) && (quat.yRotation(PI).w == quat.zRotation(PI).w));
        //assert(quat.rotateX(PI).w == to!(quat.qt)(cos(PI)));
        assert(quat.xRotation(PI).quaternion == quat.identity.rotateX(PI).quaternion);
        assert(quat.yRotation(PI).quaternion == quat.identity.rotateY(PI).quaternion);
        assert(quat.zRotation(PI).quaternion == quat.identity.rotateZ(PI).quaternion);

        assert(quat.axisRotation(PI, vec3(1.0f, 1.0f, 1.0f)).quaternion[1..4] == [1.0f, 1.0f, 1.0f]);
        assert(quat.axisRotation(PI, vec3(1.0f, 1.0f, 1.0f)).w == quat.xRotation(PI).w);
        assert(quat.axisRotation(PI, vec3(1.0f, 1.0f, 1.0f)).quaternion ==
               quat.identity.rotateAxis(PI, vec3(1.0f, 1.0f, 1.0f)).quaternion);

        quat q1 = quat.eulerRotation(PI, PI, PI);
        assert((q1.x > -1e-16) && (q1.x < 1e-16));
        assert((q1.y > -1e-16) && (q1.y < 1e-16));
        assert((q1.z > -1e-16) && (q1.z < 1e-16));
        //assert(q1.w == -1.0f);
        assert(quat.eulerRotation(PI, PI, PI).quaternion == quat.identity.rotateEuler(PI, PI, PI).quaternion);
    }

    Quaternion opBinary(string op : "*")(Quaternion inp) const {
        Quaternion ret;

        ret.w = -x * inp.x - y * inp.y - z * inp.z + w * inp.w;
        ret.x = x * inp.w + y * inp.z - z * inp.y + w * inp.x;
        ret.y = -x * inp.z + y * inp.w + z * inp.x + w * inp.y;
        ret.z = x * inp.y - y * inp.x + z * inp.w + w * inp.z;

        return ret;
    }

    auto opBinaryRight(string op, T)(T inp) const if(!isQuaternion!T) {
        return this.opBinary!(op)(inp);
    }

    Quaternion opBinary(string op)(Quaternion inp) const  if((op == "+") || (op == "-")) {
        Quaternion ret;

        mixin("ret.w = w" ~ op ~ "inp.w;");
        mixin("ret.x = x" ~ op ~ "inp.x;");
        mixin("ret.y = y" ~ op ~ "inp.y;");
        mixin("ret.z = z" ~ op ~ "inp.z;");

        return ret;
    }

    Vector!(qt, 3) opBinary(string op : "*")(Vector!(qt, 3) inp) const {
        Vector!(qt, 3) ret;

        qt ww = w^^2;
        qt w2 = w * 2;
        qt wx2 = w2 * x;
        qt wy2 = w2 * y;
        qt wz2 = w2 * z;
        qt xx = x^^2;
        qt x2 = x * 2;
        qt xy2 = x2 * y;
        qt xz2 = x2 * z;
        qt yy = y^^2;
        qt yz2 = 2 * y * z;
        qt zz = z * z;

        ret.vector =  [ww * inp.x + wy2 * inp.z - wz2 * inp.y + xx * inp.x +
                       xy2 * inp.y + xz2 * inp.z - zz * inp.x - yy * inp.x,
                       xy2 * inp.x + yy * inp.y + yz2 * inp.z + wz2 * inp.x -
                       zz * inp.y + ww * inp.y - wx2 * inp.z - xx * inp.y,
                       xz2 * inp.x + yz2 * inp.y + zz * inp.z - wy2 * inp.x -
                       yy * inp.z + wx2 * inp.y - xx * inp.z + ww * inp.z];

       return ret;
    }

    Quaternion opBinary(string op : "*")(qt inp) const {
        return Quaternion(w*inp, x*inp, y*inp, z*inp);
    }

    void opOpAssign(string op : "*")(Quaternion inp) {
        qt w2 = -x * inp.x - y * inp.y - z * inp.z + w * inp.w;
        qt x2 = x * inp.w + y * inp.z - z * inp.y + w * inp.x;
        qt y2 = -x * inp.z + y * inp.w + z * inp.x + w * inp.y;
        qt z2 = x * inp.y - y * inp.x + z * inp.w + w * inp.z;
        w = w2; x = x2; y = y2; z = z2;
    }

    void opOpAssign(string op)(Quaternion inp) if((op == "+") || (op == "-")) {
        mixin("w = w" ~ op ~ "inp.w;");
        mixin("x = x" ~ op ~ "inp.x;");
        mixin("y = y" ~ op ~ "inp.y;");
        mixin("z = z" ~ op ~ "inp.z;");
    }

    void opOpAssign(string op : "*")(qt inp) {
        quaternion[0] *= inp;
        quaternion[1] *= inp;
        quaternion[2] *= inp;
        quaternion[3] *= inp;
    }

    unittest {
        quat q1 = quat.identity;
        quat q2 = quat(3.0f, 0.0f, 1.0f, 2.0f);
        quat q3 = quat(3.4f, 0.1f, 1.2f, 2.3f);

        assert((q1 * q1).quaternion == q1.quaternion);
        assert((q1 * q2).quaternion == q2.quaternion);
        assert((q2 * q1).quaternion == q2.quaternion);
        quat q4 = q3 * q2;
        assert((q2 * q3).quaternion != q4.quaternion);
        q3 *= q2;
        assert(q4.quaternion == q3.quaternion);
        assert(almostEqual(q4.x, 0.4f));
        assert(almostEqual(q4.y, 6.8f));
        assert(almostEqual(q4.z, 13.8f));
        assert(almostEqual(q4.w, 4.4f));

        quat q5 = quat(1.0f, 2.0f, 3.0f, 4.0f);
        quat q6 = quat(3.0f, 1.0f, 6.0f, 2.0f);

        assert((q5 - q6).quaternion == [-2.0f, 1.0f, -3.0f, 2.0f]);
        assert((q5 + q6).quaternion == [4.0f, 3.0f, 9.0f, 6.0f]);
        assert((q6 - q5).quaternion == [2.0f, -1.0f, 3.0f, -2.0f]);
        assert((q6 + q5).quaternion == [4.0f, 3.0f, 9.0f, 6.0f]);
        q5 += q6;
        assert(q5.quaternion == [4.0f, 3.0f, 9.0f, 6.0f]);
        q6 -= q6;
        assert(q6.quaternion == [0.0f, 0.0f, 0.0f, 0.0f]);

        quat q7 = quat(2.0f, 2.0f, 2.0f, 2.0f);
        assert((q7 * 2).quaternion == [4.0f, 4.0f, 4.0f, 4.0f]);
        assert((2 * q7).quaternion == (q7 * 2).quaternion);
        q7 *= 2;
        assert(q7.quaternion == [4.0f, 4.0f, 4.0f, 4.0f]);

        vec3 v1 = vec3(1.0f, 2.0f, 3.0f);
        assert((q1 * v1).vector == v1.vector);
        assert((v1 * q1).vector == (q1 * v1).vector);
        assert((q2 * v1).vector == [-2.0f, 36.0f, 38.0f]);
    }

    int opCmp(ref const Quaternion qua) const {
        foreach(i, a; quaternion) {
            if(a < qua.quaternion[i]) {
                return -1;
            } else if(a > qua.quaternion[i]) {
                return 1;
            }
        }

        // Quaternions are the same
        return 0;
    }

    bool opEquals(const Quaternion qu) const {
        return quaternion == qu.quaternion;
    }

    unittest {
        assert(quat(1.0f, 2.0f, 3.0f, 4.0f) == quat(1.0f, 2.0f, 3.0f, 4.0f));
        assert(quat(1.0f, 2.0f, 3.0f, 4.0f) != quat(1.0f, 2.0f, 3.0f, 3.0f));

        assert(!quat(float.nan, float.nan, float.nan, float.nan).isFinite);
        if(quat(1.0f, 1.0f, 1.0f, 1.0f).isFinite) { }
        else { assert(false); }
    }
}

/// Pre-defined quaternion of type float.
alias quat = Quaternion!(float);

struct Rect(type) {
    alias rt = type; /// Holds the internal type of the rect.
    alias RectT = Rect!type;

    union {
        struct {
            rt x;
            rt y;
            rt width;
            rt height;
        }
        rt[4] elements; /// Holds the x, y, width and height
    }

    /// Returns a pointer to the quaternion in memory, it starts with the w coordinate.
    auto ptr() const { return elements.ptr; }
@safe pure nothrow:

    /**
        returns identity rect
    */
    static auto identity() {
        return RectT(0, 0, 0, 0);
    }

    /**
        Left coordinate of the rect
    */
    rt left() const {
        return x;
    }

    /**
        Right coordinate of the rect
    */
    rt right() const {
        return x+width;
    }

    /**
        Top coordinate of the rect
    */
    rt top() const {
        return y;
    }

    /**
        Bottom coordinate of the rect
    */
    rt bottom() const {
        return y+height;
    }

    /**
        Gets the center of a rect
    */
    vec2 center() const {
        return vec2(this.x + (this.width/2), this.y + (this.height/2));
    }

    /**
        Gets whether this rect intersects another rect
    */
    bool intersects(rtype)(rtype other) const if (isRect!rtype) {
        return !(other.left >= this.right || other.right <= this.left || other.top >= this.bottom || other.bottom <= this.top);
    }

    /**
        Gets whether this rect intersects a vector
    */
    bool intersects(vtype)(vtype other) const if (isVector!vtype) {
        return !(other.x >= this.right || other.x <= this.left || other.y >= this.bottom || other.y <= this.top);
    }

    /**
        Displaces the rect by the specified amount
    */
    void displace(vtype)(vtype other) const if (isVector!vtype && vtype.dimension == 2) {
        this.x += other.x;
        this.y += other.y;
    }

    /**
        Gets a rect that has been displaced by the specified amount
    */
    RectT displaced(vtype)(vtype other) const if (isVector!vtype && vtype.dimension == 2) {
        return RectT(this.x+other.x, this.y+other.y, this.width, this.height);
    }

    /**
        Expands the rect by the specified amount
    */
    void expand(vtype)(vtype other) const if (isVector!vtype && vtype.dimension == 2) {
        this.x -= other.x;
        this.y -= other.y;
        this.width += other.x*2;
        this.height += other.y*2;
    }

    /**
        Gets a rect that has been expanded by the specified amount
    */
    RectT expanded(vtype)(vtype other) const if (isVector!vtype && vtype.dimension == 2) {
        return RectT(this.x-other.x, this.y-other.y, this.width+(other.x*2), this.height+(other.y*2));
    }

    /**
        Gets the UV coordinates of each corner and returns them as a vec4 of the rect's type
    */
    Vector!(rt, 4) uvs() const {
        return Vector!(rt, 4)(this.left, this.top, this.right, this.bottom);
    }
}

alias rect = Rect!(float);
alias rectd = Rect!(double);

@("rect intersects")
unittest {
    rect a = rect(0, 0, 32, 32);
    rectd b = rectd(16, 16, 32, 32);
    rect c = rect(0, 32, 32, 32);
    vec2 p = vec2(8, 8);

    assert(a.intersects(b));
    assert(b.intersects(c));
    assert(!a.intersects(c));
    
    assert(a.intersects(p));
    assert(!b.intersects(p));
    assert(!c.intersects(p));
}

@("rect uvs")
unittest {
    assert(rect(16, 16, 32, 64).uvs == vec4(16, 16, 16+32, 16+64));
}

@("rect center")
unittest {
    rect a = rect(0, 0, 32, 32);
    assert(a.center == vec2(16, 16));
}