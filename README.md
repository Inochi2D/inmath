InMath
====

InMath, forked from gl3n, is a mathematics library for D providing math functions useful for multimedia software.

InMath provides all the math you need to work with OpenGL and Metal.

* Linear algebra
    * Vectors
    * Matrices
        * Metal-compatible orthographic and projection matrix functions! (`perspective01`, `orthographic01`)
    * Quaternions
    * Rectangles
    * AABBs
    * Planes,
    * Frustrums
* Interpolation
    * Linear interpolation (`lerp`)
    * Spherical linear interpolation (`slerp`)
    * Hermite interpolation
    * Catmull rom interpolation
    * Cubic interpolation
* Dampening
    * Basic distance based dampening (`dampen`)
    * Dampened Harmonic Oscillator (`smoothDamp`)
* Noise
    * OpenSimplex2
* Colors
    * HSV<->RGB conversion
    * Hexstring to RGB(A) conversion
    * RGBA32<->RGBAF conversion
* No unexpected side effects
    * No unexpected casts
    * Vector multiplication is _not_ dot product.
    * InMath uses camelCase for the naming of things
    * InMath simplifies parts of gl3n by removing excessive use of aliases.
    * InMath is usable without the GC!

Installation
============

You can use inmath in your project via the dub package database.

Examples
========

```D
vec4 v4 = vec4(1.0f, vec3(2.0f, 3.0f, 4.0f));
vec4 v4_2 = vec4(1.0f, vec4(1.0f, 2.0f, 3.0f, 4.0f).xyz); // "dynamic" swizzling with opDispatch
vec4 v4_3 = v4_2.xxyz; // opDispatch returns a static array which you can pass directly to the ctor of a vector!

vec3 v3 = my_3dvec.rgb;
vec3 foo = v4.xyzzzwzyyxw.xyz // not useful but possible!

mat4 m4fv = mat4.translation(-0.5f, -0.54f, 0.42f).rotateX(PI).rotateZ(PI/2);
glUniformMatrix4fv(location, 1, GL_TRUE, m4fv.ptr); // yes they are row major!

alias Matrix!(double, 4, 4) mat4d;
mat4d projection;
glGetDoublev(GL_PROJECTION_MATRIX, projection.ptr);

mat3 inv_view = view.rotation;
mat3 inv_view = mat3(view);

mat4 m4 = mat4(vec4(1.0f, 2.0f, 3.0f, 4.0f), 5.0f, 6.0f, 7.0f, 8.0f, vec4(...) ...); 
```

```D
    void strafeLeft(float delta) { // A
        vec3 vcross = cross(up, forward).normalized;
        _position = _position + vcross.dot(delta);
    }

    void strafeRight(float delta) { // D
        vec3 vcross = cross(up, forward).normalized;
        _position = _position - vcross.dot(delta);
    }
```
