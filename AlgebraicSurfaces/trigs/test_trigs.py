from AlgebraicSurfaces import algebraic_surface, func_name

resolution=16

def test():
    print(f'resolution:{resolution}, r^2={resolution**2}')
    indexes, coords, textures, normals = algebraic_surface(0, resolution)
    print(indexes, coords, textures, normals, sep='\n\n')
    print(f'resolution:{resolution}\nsizes: indexes:{len(indexes)}, coords:{len(coords)}, textures:{len(textures)}, normals:{len(normals)} ')

def test_all():
    for nf in range(0, 27):
        indexes, coords, textures, normals = algebraic_surface(nf, resolution)
        print(f'func:{nf:2}:{func_name(nf):22} sizes: coords:{len(coords)}, textures:{len(textures)}, normals:{len(normals)}')
    print('ok!')

test()
# test_all()