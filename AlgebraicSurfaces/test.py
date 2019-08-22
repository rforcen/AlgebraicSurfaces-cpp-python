from AlgebraicSurfaces import algebraic_surface, func_name, func_names

resolution = 16


def test():
    print(f'resolution:{resolution}, r^2={resolution**2}')
    indexes, coords, textures, normals, colors = algebraic_surface(0, resolution, 2)
    print(indexes, coords, textures, normals, colors, sep='\n\n')
    print(
        f'resolution:{resolution}\nsizes: indexes:{len(indexes)}, coords:{len(coords)}, textures:{len(textures)}, normals:{len(normals)}, colors:{len(colors)}')


def test_all():
    for nf in range(0, 27):
        indexes, coords, textures, normals, colors = algebraic_surface(nf, resolution)
        print(
            f'func:{nf:2}:{func_name(nf):22} sizes: coords:{len(coords)}, textures:{len(textures)}, normals:{len(normals)}')
    print('ok!')


def print_func_names():
    print('supported funcs:')
    for i, fn in enumerate(func_names()):
        print(f'{i:2}:{fn}')


test()
print_func_names()
# test_all()
