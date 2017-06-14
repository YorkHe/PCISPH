#version 430 core

out vec4 color;

uniform vec3 lightDir;

float Ns = 250;
vec4 mat_specular = vec4(1);
vec4 light_specular = vec4(1);

void main() {
	vec3 N;
	N.xy = gl_PointCoord * 2.0 - vec2(1.0);
	float mag = dot(N.xy, N.xy);
	if (mag > 1.0) discard;
	N.z = sqrt(1.0 - mag);

	float diffuse = max(0.0, dot(lightDir, N));

	vec3 eye = vec3(1.0f, 0.25f, 1.5f);
	vec3 halfVector = normalize(eye + lightDir);

	float spec = max(pow(dot(N, halfVector), Ns), 0.0f);

	vec4 S = light_specular * mat_specular * spec;

	color = vec4(0.0f, 0.0f, 1.0f, 1.0f) * diffuse + S;
}