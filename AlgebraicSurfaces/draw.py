import sys

from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import (QApplication, QMainWindow)
from rendererGL import RendererGL
import random
from cpp.AlgebraicSurfaces import algebraic_surface, func_name, func_names


class AlgebraicSurfaces_widget(RendererGL):
    coords = None
    faces = None
    textures = None
    normals = None

    scale = 0.6
    win = None
    needs_compile = True
    gl_compiled_list = 1

    def __init__(self, mesh, win):
        super(AlgebraicSurfaces_widget, self).__init__()

        self.win = win

        self.coords, self.textures, self.normals = mesh  # mesh expand (coords, textures, normals)

        self.setFocusPolicy(Qt.StrongFocus)  # accepts key events

    def init(self, gl):
        self.sceneInit(gl)
        gl.glCullFace(gl.GL_FRONT)

    def draw(self, gl):

        def draw_lines(gl):
            gl.glLineWidth(2)
            gl.glColor4f(1, 1, 1, 0.1)

            gl.glBegin(gl.GL_LINES)
            for coord in self.coords:
                gl.glVertex3fv(coord)
            gl.glEnd()

        def draw_mesh(gl):
            gl.glColor3ub(200, 132, 102)
            gl.glEnable(gl.GL_NORMALIZE)
            gl.glBegin(gl.GL_QUADS)

            for ic, c in enumerate(self.coords):
                gl.glVertex3fv(c)
                gl.glNormal3fv(self.normals[ic])

            gl.glEnd()

        def compile(gl):
            if self.needs_compile:
                gl.glNewList(self.gl_compiled_list, gl.GL_COMPILE)

                draw_mesh(gl)
                # draw_lines(gl)

                gl.glEndList()
                self.needs_compile = False

        def draw_list(gl):
            compile(gl)
            gl.glCallList(self.gl_compiled_list)

        gl.glScalef(self.scale, self.scale, self.scale)

        draw_list(gl)


class Main(QMainWindow):
    def __init__(self, n_func, resolution, *args):
        super(Main, self).__init__(*args)

        self.setWindowTitle(f'Algebraic Surfaces: {func_name(n_func)}')
        self.setCentralWidget(AlgebraicSurfaces_widget(algebraic_surface(n_func, resolution), self))
        self.show()


if __name__ == '__main__':
    def print_funciotn_names():
        print('supported funcs:')
        for i, fn in enumerate(func_names()):
            print(f'{i:2}:{fn}')

    print_funciotn_names()

    n_func, resolution = 6, 400
    app = QApplication(sys.argv)
    Main(n_func, resolution)
    app.exec_()
