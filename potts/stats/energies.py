
def Hamiltonian(cocycle, cubes, faceIndex, field):
    return -sum(
        1 if not sum(int(cocycle[faceIndex[face]]) for face in cube.faces) % field else 0
        for cube in cubes
    )
