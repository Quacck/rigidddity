import networkx as nx
from sympy import Matrix, pprint, floor, binomial
import sympy
from linkages import *
import openmesh as om
import numpy as np

def graph_to_matrix(G):
    M = Matrix()
    for edge in G.edges:
        row = []
        d = edge[0] - edge[1]

        for vertex in G.nodes:
            if vertex == edge[0]  : row.extend(d)
            elif vertex == edge[1]: row.extend(-d)
            else:                   row.extend([0 for _ in range(DIM)])
        M=Matrix([M,row])
    return M

def pin_face(mesh: om.PolyMesh, M: Matrix):
    face_to_be_pinned = next(mesh.faces())
    pins = [vertex.idx() for vertex in mesh.fv(face_to_be_pinned)]
    return set_pinning(pins, M)
    

def set_pinning(pins, M: Matrix):
    if type(pins) is int: pins = [pins]

    for pin in pins:
        for dim in range(DIM):
            row = np.zeros(M.cols, int)
            row[pin*DIM + dim] = 1
            M = Matrix([M, list(row)])

    return M

def close_to_zero(x): #unused
    if isinstance(x, sympy.core.numbers.Rational):
        return abs(x) < 1e-30
    else:
        return x.is_zero

# helper function to convert nullspace to a list of motions 
def get_motions(N):
    if type(N) is Matrix: N = N.nullspace()
    return [[format(float(val), '.15f') + "*v" + str(floor(i/DIM)+1) + ("xyzwrüdiger"[i%DIM]) for i,val in enumerate(vector) if val != 0] for vector in N]

#just a cute function to convert detected motions into a human readible string
def motions_to_string(motions):
    string = ""
    for v in motions:
        if len(v) !=0:
            for val in v:
                if val == v[0] and len(v)>1: word = " depends on "
                elif val != v[len(v)-1]    : word = ", and "
                elif len(v) == 1           : word = " is free\n"
                else                       : word ="\n"
                string += val + word
    return string

def model_to_graph(mesh: om.PolyMesh) -> nx.Graph:
    mesh.update_normals() #lol
    graph = nx.Graph()

    points = mesh.points()

    rng = np.random.default_rng()
    wiggled_points = [Point([int((coord + rng.random() * 1e-3)  * 100000) for coord in point]) for point in points]

    for edge in mesh.edge_vertex_indices():
        graph.add_edge(wiggled_points[edge[0]], wiggled_points[edge[1]])

    for face in mesh.faces():
        neighbouring_faces = mesh.ff(face)
        for neighbour in neighbouring_faces:
            if np.allclose(mesh.normal(face), mesh.normal(neighbour)):
                vertices_face = mesh.fv(face)
                vertices_neighbour = mesh.fv(neighbour)
                first = [vertex for vertex in vertices_face if vertex not in vertices_neighbour]
                second = [vertex for vertex in vertices_neighbour if vertex not in vertices_face]
                first_point = wiggled_points[first[0].idx()]
                second_point = wiggled_points[second[0].idx()]
                if not graph.has_edge(first_point, second_point):
                    graph.add_edge(first_point, second_point)

    return graph    

def check_rigidity(M, pinned: bool):
    rank = M.rank()
    if pinned:
        return rank == M.cols
    
    return rank == M.cols - (DIM+1) * DIM / 2

def model_to_matrix(meshname):
    mesh = om.read_polymesh(meshname)

    graph = model_to_graph(mesh)

    A = graph_to_matrix(graph)
    A = pin_face(mesh, A)
    return A

def check_model_rigidity(meshname):
    return check_rigidity(model_to_matrix(meshname), True)

def get_model_motion_string(meshname):
   return motions_to_string(get_motions(model_to_matrix(meshname)))

if __name__ == "__main__":
    A = model_to_matrix("models/paperplane.stl")
    pprint(A)
    print(A.shape)
    print ("the linkage is infinitesimally rigid!" if check_rigidity(A, True) else "the linkage is infinitesimally flexible")
    print(motions_to_string(get_motions(A)))
