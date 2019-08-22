import sys

from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import (QApplication, QMainWindow)
from rendererGL import RendererGL
from AlgebraicSurfaces import algebraic_surface, func_name, func_names
from array import array


# triangle strip drawing version

class AlgebraicSurfaces_widget(RendererGL):
    indexes = None
    coords = None
    textures = None
    normals = None
    colors = None

    color_bronze = (200, 132, 102)

    scale = 0.7
    win = None

    def __init__(self, mesh, win):
        super(AlgebraicSurfaces_widget, self).__init__()

        self.win = win

        self.indexes, self.coords, self.textures, self.normals, self.colors = mesh  # mesh expand (coords, textures, normals)

        self.setFocusPolicy(Qt.StrongFocus)  # accepts key events

    def init(self, gl):
        def set_draw(gl):
            def np2array(type, np_vect):  # array = numpy, take the fast lane
                arr = array(type)
                arr.frombytes(np_vect.tobytes())
                return arr

            # define draw components (indexes, vertex & normals)
            gl.glEnableClientState(gl.GL_VERTEX_ARRAY)
            gl.glEnableClientState(gl.GL_NORMAL_ARRAY)
            gl.glEnableClientState(gl.GL_COLOR_ARRAY)

            self.index_array = np2array('I', self.indexes)  # 'I' -> unsigned int 32, GL_UNSIGNED_INT

            gl.glVertexPointer(3, gl.GL_FLOAT, 0, np2array('f', self.coords))  # set (vertex, normals) pointers
            gl.glNormalPointer(gl.GL_FLOAT, 0, np2array('f', self.normals))
            gl.glColorPointer(3, gl.GL_FLOAT, 0, np2array('f', self.colors))

        self.sceneInit(gl)
        gl.glCullFace(gl.GL_FRONT)
        set_draw(gl)

        gl.glEnable(gl.GL_RESCALE_NORMAL)
        gl.glColor3ubv(self.color_bronze)  # for all render

    def draw(self, gl):
        gl.glScalef(self.scale, self.scale, self.scale)
        gl.glDrawElements(gl.GL_TRIANGLE_STRIP, len(self.index_array), gl.GL_UNSIGNED_INT, self.index_array)


class Main(QMainWindow):
    def __init__(self, n_func, resolution, color_map, *args):
        super(Main, self).__init__(*args)

        self.setWindowTitle(f'Algebraic Surfaces: {func_name(n_func)}')
        self.setCentralWidget(AlgebraicSurfaces_widget(algebraic_surface(n_func, resolution, color_map), self))
        self.show()


if __name__ == '__main__':
    def print_func_names():
        print('supported funcs:')
        for i, fn in enumerate(func_names()):
            print(f'{i:2}:{fn}')


    # print_func_names()

    n_func, resolution, color_map = 16, 800, 13

    app = QApplication(sys.argv)
    Main(n_func, resolution, color_map)
    app.exec_()
