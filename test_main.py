import main

def test_cube_is_rigid():
    assert main.check_model_rigidity("models/cube.stl")

def test_plane_is_rigid():
    assert main.check_model_rigidity("models/plane.stl")

def test_paperplane_is_flexible():
    assert not main.check_model_rigidity("models/paperplane.stl") 

# Sphere is very slow due to "big" matrix
#def test_sphere_is_rigid():
    #assert main.check_model_rigidity("models/sphere.stl")