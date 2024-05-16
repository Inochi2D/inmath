/**
    inmath.hsv

    Authors: David Herberth, Inochi2D Project
    License: MIT
*/
module inmath.color;

private {
    import std.conv : to;
    
    import inmath.linalg : vec3, vec4;
    import inmath.math : min, max, floor;

    version(unittest) {
        import inmath.math : almostEqual;
    }
}

@nogc:

/// Converts a 3 dimensional color-vector from the RGB to the HSV colorspace.
/// The function assumes that each component is in the range [0, 1].
@safe pure nothrow vec3 rgb2hsv(vec3 inp) {
    vec3 ret = vec3(0.0f, 0.0f, 0.0f);
    
    float h_max = max(inp.r, inp.g, inp.b);
    float h_min = min(inp.r, inp.g, inp.b);
    float delta = h_max - h_min;

   
    // h
    if(delta == 0.0f) {
        ret.x = 0.0f;
    } else if(inp.r == h_max) {
        ret.x = (inp.g - inp.b) / delta; // h
    } else if(inp.g == h_max) {
        ret.x = 2 + (inp.b - inp.r) / delta; // h
    } else {
        ret.x = 4 + (inp.r - inp.g) / delta; // h
    }

    ret.x = ret.x * 60;
    if(ret.x < 0) {
        ret.x = ret.x + 360;
    }

    // s
    if(h_max == 0.0f) {
        ret.y = 0.0f;
    } else {
        ret.y = delta / h_max;
    }

    // v
    ret.z = h_max;

    return ret;
}

/// Converts a 4 dimensional color-vector from the RGB to the HSV colorspace.
/// The alpha value is not touched. This function also assumes that each component is in the range [0, 1].
@safe pure nothrow vec4 rgb2hsv(vec4 inp) {
    return vec4(rgb2hsv(vec3(inp.rgb)), inp.a);
}

unittest {
    assert(rgb2hsv(vec3(0.0f, 0.0f, 0.0f)) == vec3(0.0f, 0.0f, 0.0f));
    assert(rgb2hsv(vec3(1.0f, 1.0f, 1.0f)) == vec3(0.0f, 0.0f, 1.0f));

    vec3 hsv = rgb2hsv(vec3(100.0f/255.0f, 100.0f/255.0f, 100.0f/255.0f));    
    assert(hsv.x == 0.0f && hsv.y == 0.0f && almostEqual(hsv.z, 0.392157, 0.000001));
    
    assert(rgb2hsv(vec3(0.0f, 0.0f, 1.0f)) == vec3(240.0f, 1.0f, 1.0f));
}

/// Converts a 3 dimensional color-vector from the HSV to the RGB colorspace.
/// RGB colors will be in the range [0, 1].
/// This function is not marked es pure, since it depends on std.math.floor, which
/// is also not pure.
@safe nothrow vec3 hsv2rgb(vec3 inp) {
    if(inp.y == 0.0f) { // s
        return vec3(inp.zzz); // v
    } else {
        float var_h = inp.x * 6;
        float var_i = to!float(floor(var_h));
        float var_1 = inp.z * (1 - inp.y);
        float var_2 = inp.z * (1 - inp.y * (var_h - var_i));
        float var_3 = inp.z * (1 - inp.y * (1 - (var_h - var_i)));

        if(var_i == 0.0f)      return vec3(inp.z, var_3, var_1);
        else if(var_i == 1.0f) return vec3(var_2, inp.z, var_1);
        else if(var_i == 2.0f) return vec3(var_1, inp.z, var_3);
        else if(var_i == 3.0f) return vec3(var_1, var_2, inp.z);
        else if(var_i == 4.0f) return vec3(var_3, var_1, inp.z);
        else                   return vec3(inp.z, var_1, var_2);
    }
}

/// Converts a 4 dimensional color-vector from the HSV to the RGB colorspace.
/// The alpha value is not touched and the resulting RGB colors will be in the range [0, 1].
@safe nothrow vec4 hsv2rgb(vec4 inp) {
    return vec4(hsv2rgb(vec3(inp.xyz)), inp.w);
}

unittest {
    assert(hsv2rgb(vec3(0.0f, 0.0f, 0.0f)) == vec3(0.0f, 0.0f, 0.0f));
    assert(hsv2rgb(vec3(0.0f, 0.0f, 1.0f)) == vec3(1.0f, 1.0f, 1.0f));

    vec3 rgb = hsv2rgb(vec3(0.0f, 0.0f, 0.392157f));
    assert(rgb == vec3(0.392157f, 0.392157f, 0.392157f));

    assert(hsv2rgb(vec3(300.0f, 1.0f, 1.0f)) == vec3(1.0f, 0.0f, 1.0f));
}

/**
    Parses a hex string to a RGB color
*/
vec3 hex2rgb(const(char)[] input) {
    return hex2rgba(input).xyz;
}

/**
    Parses a hex string to a RGBA color
*/
vec4 hex2rgba(const(char)[] input) {
    ubyte[4] colors = hex2rgbau(input);

    return vec4(
        cast(float)colors[3]/255.0,
        cast(float)colors[2]/255.0,
        cast(float)colors[1]/255.0,
        cast(float)colors[0]/255.0,
    );
}

/**
    Parses a hex string to a RGBA color

    Returns the color as a static 8 bit RGB(A) array.
*/
ubyte[4] hex2rgbau(const(char)[] input) {
    import core.stdc.stdio : printf;
    const(char)[] r = input;

    // storage for color values
    union u { int rgba; ubyte[4] colors; }
    u color;    // color

    // Skip initial if any
    if (input[0] == '#') r = r[1..$];
    else if (input[0] == '0' && input[1] == 'x') r = r[2..$];

    // Iterate through colors
    foreach(i; 0..min(r.length, 8)) {
        ubyte c;
        
        if (r[i] >= '0' && r[i] <= '9') {

            c = cast(ubyte)(r[i]-'0');
        } else if (r[i] >= 'A' && r[i] <= 'F') {

            c = cast(ubyte)(r[i]-'A'+10);
        } else if (r[i] >= 'a' && r[i] <= 'f') {

            c = cast(ubyte)(r[i]-'a'+10);
        } else {
            break;
        }

        color.rgba = (color.rgba << 4) | c;
    }

    return color.colors;
}

@("Hex color")
unittest {
    string col = "#FF0000FF";
    assert(hex2rgb(col) == vec3(1, 0, 0));
    assert(hex2rgba(col) == vec4(1, 0, 0, 1));

    string col2 = "0xFF00FF00";
    assert(hex2rgb(col2) == vec3(1, 0, 1));
    assert(hex2rgba(col2) == vec4(1, 0, 1, 0));

    string col3 = "0x00FF00FF";
    assert(hex2rgb(col3) == vec3(0, 1, 0));
    assert(hex2rgba(col3) == vec4(0, 1, 0, 1));

    string col4 = "#0A0A0A0A";
    ubyte[4] col4u = [10, 10, 10, 10];
    assert(hex2rgbau(col4) == col4u);

    string col5 = "#0A0AXXXX";
    ubyte[4] col5u = [10, 10, 0, 0];
    assert(hex2rgbau(col5) == col5u);

    string col6 = "#FF";
    ubyte[4] col6u = [255, 0, 0, 0];
    assert(hex2rgbau(col6) == col6u);
}

/**
    Returns RGBA color from 0..1 from colors in range 0..255
*/
vec4 rgbau2rgba(ubyte[4] colors) {
    return vec4(
        cast(float)colors[0]/255.0,
        cast(float)colors[1]/255.0,
        cast(float)colors[2]/255.0,
        cast(float)colors[3]/255.0,
    );
}

/**
    Returns RGB color from 0..1 from colors in range 0..255
*/
vec3 rgbu2rgb(ubyte[3] colors) {
    return vec3(
        cast(float)colors[0]/255.0,
        cast(float)colors[1]/255.0,
        cast(float)colors[2]/255.0,
    );
}