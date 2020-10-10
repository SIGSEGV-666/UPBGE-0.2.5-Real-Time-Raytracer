import bge, mathutils, math, bgl
from bgl import *
rot90 = mathutils.Matrix.Rotation(math.radians(-90), 4, mathutils.Vector((1, 0, 0))).to_3x3()
MOVEMENT_SPEED = 0.2
def update_movement(cont):
    owner = cont.owner
    values = dict((k, float(cont.sensors[k.upper()].positive)) for k in "wasdqz")
    owner.applyMovement((0, 0, -MOVEMENT_SPEED*(values['w']-values['s'])), True)
    owner.applyMovement((MOVEMENT_SPEED*(values['d']-values['a']), 0.0, 0.0), True)
    owner.applyMovement((0, 0, MOVEMENT_SPEED*(values['q']-values['z'])), False)
def update_shader(cont):
    owner = cont.owner
    rtshader = owner['rtshader']
    rtshader.setUniform3f("cam_pos", *owner.worldPosition)
    #rtshader.setUniformMatrix3("cam_orn", owner.worldOrientation, True)
    worn = owner.worldOrientation.transposed()
    for i, vec in enumerate(worn):
        rtshader.setUniform3f("cam_orn_"+str(i), *vec)
    rtshader.setUniform1f("cam_near", owner.near)
    rtshader.setUniform1f("cam_far", owner['far'])
    rtshader.setUniform1f("sim_time", bge.logic.getClockTime())
    rtshader.setUniform1f("cam_fov_y", owner['fov'])
    rtshader.setTexture(0, owner['textures']['Texture'].bindCode, "sphere_tex")
    rtshader.setTexture(1, owner['textures']['grid'].bindCode, "plane_tex")
    rtshader.setTexture(2, owner['textures']['sky'].bindCode, "sky_tex")
    rtshader.setUniform1i("use_perspective", int(owner['perspective']))
def setTextureNearest(texture):
    stid = texture.bindCode
    glBindTexture(GL_TEXTURE_2D, stid)
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR)
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR)
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE)
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE)
    glBindTexture(GL_TEXTURE_2D, 0)
def init(cont):
    owner = cont.owner
    scene = owner.scene
    cube = scene.objects['Cube']
    cube_textures = dict((t.name, t) for t in cube.meshes[0].materials[0].textures if t is not None)
    #ct0 = cube_textures['Texture']
    owner['textures'] = cube_textures
    st = owner['textures']['Texture']
    setTextureNearest(st)
    setTextureNearest(owner['textures']['sky'])
    with open(bge.logic.expandPath("//raytracer.glsl"), "r") as shaderfile:
        _RT_src = shaderfile.read()
    print(_RT_src)
    owner['rtshader'] = scene.filterManager.addFilter(0, bge.logic.RAS_2DFILTER_CUSTOMFILTER, _RT_src)
    update_shader(cont)
def run(cont):
    update_movement(cont)
    update_shader(cont)
def main(cont):
    owner = cont.owner
    if "init" not in owner:
        owner['init'] = False
        init(cont)
        owner['init'] = True
    elif owner['init']:
        run(cont)
