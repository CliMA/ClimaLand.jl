#version 300 es
precision mediump int;
precision mediump float;
precision mediump sampler2D;
precision mediump sampler3D;
precision mediump isampler2D;
precision mediump isampler3D;
precision mediump usampler2D;
precision mediump usampler3D;

out vec4 fragment_color;

// Uniforms: 
uniform mat4 view;
mat4 get_view(){return view;}
uniform mat4 model;
mat4 get_model(){return model;}
uniform mat4 projection;
mat4 get_projection(){return projection;}

void main(){}
