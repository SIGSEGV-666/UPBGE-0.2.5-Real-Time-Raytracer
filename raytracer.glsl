#define PI 3.1415926535897932384626433832795
#define EPSILON 0.000000000000001
uniform sampler2D bgl_RenderedTexture;
uniform sampler2D bgl_DepthTexture;
uniform sampler2D sphere_tex;
uniform sampler2D plane_tex;
uniform sampler2D sky_tex;
uniform float bgl_RenderedTextureWidth;
uniform float bgl_RenderedTextureHeight;
uniform vec3 cam_pos;
uniform vec3 cam_orn_0;
uniform vec3 cam_orn_1;
uniform vec3 cam_orn_2;
uniform float cam_near;
uniform float cam_far;
uniform float sim_time;
uniform int use_perspective;
mat3 cam_orn = mat3(1.0);
uniform float cam_fov_y;
float aspect_ratio = bgl_RenderedTextureWidth/bgl_RenderedTextureHeight;
vec3 X_AXIS = vec3(1.0, 0.0, 0.0);
vec3 Y_AXIS = vec3(0.0, 1.0, 0.0);
vec3 Z_AXIS = vec3(0.0, 0.0, 1.0);
int NUM_SPHERES = 300;
vec3 SPHERE_POSITIONS[150];
float SPHERE_RADIUSES[150];
struct isect_result {
    vec3 hit;
    vec3 normal;
    int cube_face;
};
vec3 BOX_FACE_COLORS[6] = {
    vec3(1.0, 0.0, 0.0),
    vec3(0.0, 1.0, 0.0),
    vec3(0.0, 0.0, 1.0),
    vec3(1.0, 1.0, 0.0),
    vec3(0.0, 1.0, 1.0),
    vec3(1.0, 1.0, 1.0)
    };
int rayPlaneIntersect(vec3 rayP, vec3 rayD, vec3 planeP, vec3 planeN, out vec3 hit, out vec3 normal)
{
    float d = dot(planeP, -planeN);
    float e = dot(rayD, -planeN);
    //if (d >= -EPSILON && d <= EPSILON){return 0;}
    float t = -(d + rayP.z * planeN.z + rayP.y * planeN.y + rayP.x * planeN.x) / (rayD.z * planeN.z + rayD.y * planeN.y + rayD.x * planeN.x);
    if (t < 0){return 0;}
    hit = rayP + t * rayD;
    normal = (e >= 0) ? (planeN) : (-planeN);
    return 1;
}
mat4 trsMatrix(vec3 t, mat3 r, vec3 s)
{
    mat4 translation = mat4(
        1.0, 0.0, 0.0, t.x,
        0.0, 1.0, 0.0, t.y,
        0.0, 0.0, 1.0, t.z,
        0.0, 0.0, 0.0, 1.0
        );
    mat4 orientation = mat4(r);
    mat4 scale = mat4(
        s.x, 0.0, 0.0, 0.0,
        0.0, s.y, 0.0, 0.0,
        0.0, 0.0, s.z, 0.0,
        0.0, 0.0, 0.0, 1.0
        );
    return translation*orientation*scale;
}
mat3 rotationMatrix(vec3 axis, float angle) //http://www.neilmendoza.com/glsl-rotation-about-an-arbitrary-axis/
{
    axis = normalize(axis);
    float s = sin(angle);
    float c = cos(angle);
    float oc = 1.0 - c;
    
    return mat3(oc * axis.x * axis.x + c,           oc * axis.x * axis.y - axis.z * s,  oc * axis.z * axis.x + axis.y * s,
                oc * axis.x * axis.y + axis.z * s,  oc * axis.y * axis.y + c,           oc * axis.y * axis.z - axis.x * s,
                oc * axis.z * axis.x - axis.y * s,  oc * axis.y * axis.z + axis.x * s,  oc * axis.z * axis.z + c);
}
vec3 getRayVector_angles(float zangle, float xangle)
{
    mat3 zrot = rotationMatrix(Y_AXIS, radians(zangle));
    mat3 xrot = rotationMatrix(X_AXIS, radians(xangle));
    return (xrot*(zrot*-Z_AXIS));
}
vec3 getRayVectorSphere()
{
    vec2 tc = gl_TexCoord[0].st;
    float half_fov = cam_fov_y*0.5;
    float zangle = mix(-half_fov, half_fov, tc.x);
    float xangle = mix(-half_fov/aspect_ratio, half_fov/aspect_ratio, 1.0-(tc.y));
    mat3 zrot = rotationMatrix(Y_AXIS, radians(zangle));
    mat3 xrot = rotationMatrix(X_AXIS, radians(xangle));
    return (xrot*(zrot*-Z_AXIS));
}
vec3 getRayVectorPerspective()
{
    vec2 tc = gl_TexCoord[0].st;
    float half_fov = cam_fov_y*0.5;
    vec3 plane_pos = vec3(0.0, 0.0, -cam_near);
    vec3 plane_vec = vec3(0.0, 0.0, 1.0);
    vec3 tlvec = getRayVector_angles(-half_fov, half_fov);
    vec3 brvec = getRayVector_angles(half_fov, -half_fov);
    vec3 tlpos, brpos;
    vec3 nullnorm;
    rayPlaneIntersect(vec3(0.0), tlvec, plane_pos, plane_vec, tlpos, nullnorm);
    rayPlaneIntersect(vec3(0.0), brvec, plane_pos, plane_vec, brpos, nullnorm);
    vec3 pixelpos = vec3(mix(tlpos.x, brpos.x, tc.x), mix(tlpos.y, brpos.y, tc.y), -cam_near);
    return normalize(pixelpos);
}
int longest_vec3_axis(vec3 v)
{
    bool first = true;
    float longest_value;
    float cvalue;
    int laxis = 0;
    for (int i = 0; i < 3; i++)
    {
        cvalue = abs(v[i]);
        if (cvalue >= longest_value || first)
        {
            longest_value = cvalue;
            laxis = i;
        }
        first = false;
    }
    return laxis;
}
int rayBoxIntersect(vec3 ro, vec3 rv, vec3 bmin, vec3 bmax, out vec3 hit, out vec3 normal, out int face)
{
    float tmin, tmax, tymin, tymax, tzmin, tzmax, t;
    vec3 invdir = vec3(1.0)/rv;
    vec3 center = (bmin+bmax)*0.5;
    int signs[3] = {0, 0, 0};
    vec3 bounds[2] = {bmin, bmax};
    for (int i = 0; i < 3; i++){signs[i] = (invdir[i] < 0) ? (1) : (0);}
    tmin = (bounds[signs[0]][0] - ro[0])*(invdir[0]);
    tmax = (bounds[1-signs[0]][0] - ro[0])*(invdir[0]);
    tymin = (bounds[signs[1]][1] - ro[1])*(invdir[1]);
    tymax = (bounds[1-signs[1]][1] - ro[1])*(invdir[1]);
    if (tmin > tymax || tymin > tmax){return 0;}
    if (tymin > tmin){tmin = tymin;}
    if (tymax < tmax){tmax = tymax;}
    tzmin = (bounds[signs[2]][2] - ro[2])*(invdir[2]);
    tzmax = (bounds[1-signs[2]][2] - ro[2])*(invdir[2]);
    if (tmin > tzmax || tzmin > tmax){return 0;}
    if (tzmin > tmin){tmin = tzmin;}
    if (tzmax < tmax){tmax = tzmax;}
    t = tmin;
    if (t < 0)
    {
        t = tmax;
        if (t < 0){return 0;}
    }
    hit = ro+(rv*t);
    vec3 hitvec = normalize(hit-center);
    vec3 bscale = (bmax-bmin)*0.5;
    vec3 svec = normalize(hitvec/bscale);
    int laxis = longest_vec3_axis(svec);
    normal = vec3(0.0);
    normal[laxis] = (svec[laxis] >= 0) ? (1.0) : (-1.0);
    face = (2*laxis)+int(svec[laxis] > 0);
    return 1;
}
int raySphereIntersect(vec3 ro, vec3 rv, vec3 center, float radius, out vec3 hit, out vec3 normal)
{
    float t0, t1, temp, t; // solutions for t if the ray intersects 
    float radius2 = pow(radius, 2);
    // geometric solution
    vec3 L = center - ro; 
    float tca = dot(L, rv); 
    // if (tca < 0) return false;
    float d2 = dot(L, L) - tca * tca; 
    if (d2 > radius2){return 0;}
    float thc = sqrt(radius2 - d2); 
    t0 = tca - thc; 
    t1 = tca + thc;
    
    if (t0 > t1)
    {
        temp = t0;
        t0 = t1;
        t1 = temp;
    }
 
    if (t0 < 0)
    { 
        t0 = t1; // if t0 is negative, let's use t1 instead 
        if (t0 < 0) return 0; // both t0 and t1 are negative 
    } 
 
    t = t0;
    hit = ro+(rv*t);
    normal = normalize(hit-center);
    return 1;
}
float dist2depth(float dist)
{
    return (dist-cam_near)/(cam_far-cam_near);
}
float sine(float v)
{
    return sin(radians(180)*v);
}
float cosine(float v)
{
    return cos(radians(180)*v);
}
float checker_channel(float v)
{
    if (v < 0.5){return 0.0;}
    else        {return 1.0;}
}
vec4 checker(vec3 v, vec3 gridscale)
{
    vec3 modv = mod(v, gridscale)/gridscale;
    return vec4(checker_channel(modv.x), checker_channel(modv.y), checker_channel(modv.z), 1.0);
}
vec2 sphere2uv(vec3 vec)
{
    vec3 v = normalize(vec);
    vec3 d = vec3(v.x, -v.z, v.y);
    return vec2(
        ((0.5+(atan(d.z, d.x)/(2*PI)))),
        0.5-(asin(d.y)/PI)
    );
}
float spot_light_calc(vec3 point, vec3 normal, vec3 lightpos, vec3 spotdir, float max_angle, float maxdistance, float smoothing)
{
    vec3 lightvec = normalize(point-lightpos);
    float lightdist = distance(point, lightpos);
    float angledot = dot(lightvec, spotdir);
    float angle = degrees(acos(angledot));
    if (angle < max_angle)
    {
        return max(dot(lightvec, -normal), 0.0)*pow((max(angledot, 0.0)), smoothing)*(max(1.0-min(lightdist/maxdistance, 1.0), 0.0));
    }
    else
    {
        return 0.0;
    }
}
bool test_box(vec3 ro, vec3 rv, vec3 box_min, vec3 box_max, float cdistance, out vec3 color, out float afterdistance)
{
    vec3 hit, normal;
    int face;
    float dist, lightdist;
    vec3 center = (box_min+box_max)*0.5;
    vec3 hitvec;
    vec3 lightvec;
    float lightfac, specular;
    vec3 lightpos = cam_pos;
    //isect_result res;
    if (rayBoxIntersect(ro, rv, box_min, box_max, hit, normal, face) == 1)
    {
        dist = dist2depth(distance(hit, ro));
        if (dist <= cdistance)
        {
            afterdistance = dist;
            hitvec = normalize(hit-center);
            //lightvec = normalize(hit-lightpos);
            //lightdist = distance(hit, lightpos);
            //lightfac = (max(dot(lightvec, -normal), 0.0)*pow(max(dot(lightvec, cam_orn*-Z_AXIS), 0.0), 4.0)*(1.0-min(lightdist/100.0, 1.0)));
            lightfac = spot_light_calc(hit, normal, lightpos, cam_orn*-Z_AXIS, 90.0, 100.0, 8.0);
            specular = pow(lightfac, 16.0);
            //color = (texture2D(sphere_tex, sphere2uv(hitvec)).rgb*lightfac)+(vec3(1.0)*specular);
            color = (texture2D(sphere_tex, sphere2uv(hitvec)).rgb);
            return true;
        }
    }
    return false;
}
bool test_boxen(vec3 ro, vec3 rv, float cdistance, out vec3 color, out float afterdistance) //I BOUGHT TWO BOXEN OF DONUTS!
{
    vec3 box_pos;
    vec3 vec;
    vec3 offset = vec3(0.0, 0.0, 1.0);
    float radius = 7.5+(sine(sim_time)*2.5);
    int num_boxes = 30;
    float angle_incr = 360.0/num_boxes;
    float cangle = 0.0;
    mat3 rotmat;
    vec3 bmin = vec3(-0.5, -0.5, -0.5);
    vec3 bmax = vec3(0.5, 0.5, 0.5);
    float dist = cdistance;
    float adist;
    bool gotit = false;
    vec3 color2;
    //isect_result res;
    for (int i = 0; i < num_boxes; i++)
    {
        cangle = (angle_incr*i)+(sim_time*50.0);
        rotmat = rotationMatrix(Z_AXIS, radians(cangle));
        vec = rotmat*Y_AXIS;
        box_pos = offset+(vec*radius);
        if (test_box(ro, rv, bmin+box_pos, bmax+box_pos, dist, color, adist))
        {
            dist = adist;
            gotit = true;
        }
    }
    if (test_box(ro, rv, vec3(-10.0, -10.0, 3.0), vec3(10.0, 10.0, 23.0), dist, color2, adist))
    {
        dist = adist;
        gotit = true;
        color = color2;
    }
    afterdistance = dist;
    return gotit;
}
void main()
{
    cam_orn[0] = cam_orn_0;
    cam_orn[1] = cam_orn_1;
    cam_orn[2] = cam_orn_2;
    int row;
    vec3 rv = (use_perspective == 1) ? (getRayVectorPerspective()) : (getRayVectorSphere());
    vec3 sphere_origin;
    float sphere_radius;
    float lightlevel, specular;
    vec3 hit, normal;
    float cdistance = texture2D(bgl_DepthTexture, gl_TexCoord[0].st).r;
    float dist;
    vec4 bgColor = texture2D(sky_tex, sphere2uv(cam_orn*rv));
    vec4 fragColor = bgColor;
    vec3 lightpos = cam_pos;
    float lightdist;
    float lightmaxdist = 200.0;
    vec3 lightvec;
    float cz = 0.0;
    int count = 0;
    vec2 uv;
    vec3 texvec;
    vec3 box_min = vec3(-10.0, -10.0, -10.0);
    vec3 box_max = vec3(10.0, 10.0, 10.0);
    int face;
    vec3 plane_co = vec3(0.0, 0.0, 0.0);
    vec3 plane_n = vec3(0.0, 0.0, 1.0);
    vec3 plane_refl;
    vec3 plane_color = vec3(0.5);
    //vec3 plane_refl_hit, plane_refl_n;
    vec3 box_outcolor;
    float afterdist;
    if (test_boxen(cam_pos, cam_orn*rv, cdistance, box_outcolor, afterdist))
    {
        cdistance = afterdist;
        fragColor = vec4(box_outcolor, 1.0);
    }
    float nulldist1 = 1.0;
    vec3 plane_outnormal;
    vec3 plane_outcolor;
    if (rayPlaneIntersect(cam_pos, cam_orn*rv, plane_co, plane_n, hit, normal) == 1)
    {
        dist = dist2depth(distance(hit, cam_pos));
        if (dist <= cdistance)
        {
            cdistance = dist;
            //fragColor = texture2D(sphere_tex, hit.xy/10);
            //fragColor = vec4(normal, 1.0);
            //fragColor = (texture2D(plane_tex, hit.xy/100.0)*lightlevel)+(vec4(vec3(1.0)*specular, 1.0));
            plane_outnormal = normalize(normal+vec3(sine(hit.y+(sim_time*2))*0.005, cosine(hit.x+sim_time)*0.01, 0.0));
            plane_refl = reflect(cam_orn*rv, plane_outnormal);
            //lightvec = normalize(hit-lightpos);
            //lightdist = distance(hit, lightpos);
            //lightlevel = dot(-plane_outnormal, lightvec)*(1.0-clamp(distance(lightpos, hit)/lightmaxdist, 0.0, 1.0));
            //specular = pow(lightlevel, 16);
            //lightlevel = (max(dot(lightvec, -plane_refl), 0.0)*max(dot(lightvec, cam_orn*-Z_AXIS), 0.0)*(1.0-min(lightdist/1000.0, 1.0)));
            lightlevel = spot_light_calc(hit, plane_refl, lightpos, cam_orn*-Z_AXIS, 45.0, 1000.0, 2.0);
            specular = pow(lightlevel, 4.0);
            if (test_boxen(hit, plane_refl, nulldist1, box_outcolor, nulldist1))
            {
                plane_outcolor = mix(plane_color*lightlevel, box_outcolor, 0.75);
            }
            else
            {
                plane_outcolor = mix(plane_color*lightlevel, texture2D(sky_tex, sphere2uv(plane_refl)).rgb, 0.75);
            }
            fragColor = vec4((plane_outcolor)+(vec3(1.0)*specular), 1.0);
        }
    }
    gl_FragColor = mix(fragColor, bgColor, pow(cdistance, 2.0)); //fade into background with distance.
}