#version 460 core

out vec4 color;

layout (std140, binding = 1) uniform info {
    float time;
    int width;
    int height;
};

uniform sampler2D img;

void main() {
   color = texture(img, gl_FragCoord.xy / vec2(width, height));
}
