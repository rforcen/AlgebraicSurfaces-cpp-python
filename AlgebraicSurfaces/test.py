from cpp.AlgebraicSurfaces import algebraic_surface, func_name

resolution=200

def test():
    coords, textures, normals = algebraic_surface(0, 5)
    print(coords, textures, normals)
    print(f'sizes: coords:{len(coords)}, textures:{len(textures)}, normals:{len(normals)} ')

def test_all():
    for nf in range(0, 27):
        coords, textures, normals = algebraic_surface(nf, resolution)
        print(f'func:{nf:2}:{func_name(nf):22} sizes: coords:{len(coords)}, textures:{len(textures)}, normals:{len(normals)}')
    print('ok!')

test_all()