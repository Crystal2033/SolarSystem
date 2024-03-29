import math
import sys

import numpy as np
import pyqtgraph as pg
import pyqtgraph.opengl as gl
from pyqtgraph.Qt import QtCore, QtGui, QtWidgets


class Rotator:
    # angles -- [x, y, z]
    @staticmethod
    def get_rotate_matrixes(angles):
        x = [[1, 0, 0],
             [0, math.cos(angles[0]), math.sin(angles[0])],
             [0, -math.sin(angles[0]), math.cos(angles[0])]]

        y = [[math.cos(angles[1]), 0, math.sin(angles[1])],
             [0, 1, 0],
             [-math.sin(angles[1]), 0, math.cos(angles[1])]]

        z = [[math.cos(angles[2]), math.sin(angles[2]), 0],
             [-math.sin(angles[2]), math.cos(angles[2]), 0],
             [0, 0, 1]]
        return x, y, z


class MovementManager:

    def update(self):
        self.make_movement(sun, earth)
        self.make_movement(earth, moon)

    def make_movement(self, planet, satellite):
        planet_orbit = planet.get_orbit()

        # Находим время, удобное для связи с углом.
        satellite.rotate_timer = int((satellite.rotate_timer + satellite.get_delta_rot_angle()) % 360)
        angle = math.radians(satellite.rotate_timer)  # angle in radians

        planet_orbit.move_satellite(angle, planet)


class CelestialObject:
    def __init__(self):
        self.sphere_meshdata = None  # Data about sphere. Info about it.
        self.sphere_meshitem = None  # Sphere itself
        self.radius = None
        self.plots = []
        self.plots_pos = []
        self.pos_x = None
        self.pos_y = None
        self.pos_z = None
        self.color = None
        self.orbit = None
        self.delta_rotate_angle = 1
        self.rotate_timer = 0

        self.axle_rotate_angle_x = 0
        self.axle_rotate_angle_y = 0
        self.axle_rotate_angle_z = 0

    def create_sphere(self, radius, sphere_color, global_position, _window):
        # Метод создает небесную сферу и устанавливает долготы и широты для того, чтобы видеть вращение.
        self.radius = radius
        self.color = sphere_color
        self.sphere_meshdata = gl.MeshData.sphere(radius=radius, rows=100, cols=100)
        self.sphere_meshitem = gl.GLMeshItem(meshdata=self.sphere_meshdata, smooth=False, color=self.color)
        self.pos_x = global_position[0]
        self.pos_y = global_position[1]
        self.pos_z = global_position[2]
        self.add_long_and_lati(_window)
        self.move_sphere_to_pos(self.pos_x, self.pos_y, self.pos_z, [0, 0, 0])
        self.sphere_meshitem.setGLOptions('opaque')  # opaque
        _window.addItem(self.sphere_meshitem)

    def set_orbit_params(self, _window, orbit_a, orbit_b, orbit_tilt_angle=0, orbit_color=(0.015, 0.67, 0.82, 1.)):
        # Метод устанавливает параметры орбиты и создает орбиту. Также поворачивает орбиту на заданный угол
        self.orbit = Orbit()
        self.orbit.set_orbit_tilt(orbit_tilt_angle)
        position = [self.pos_x, self.pos_y, self.pos_z]
        self.orbit.create_orbit(_window, orbit_a, orbit_b, orbit_color, position)

    def set_delta_rot_angle(self, angle_delta):
        self.delta_rotate_angle = angle_delta

    def get_delta_rot_angle(self):
        return self.delta_rotate_angle

    def set_satellite(self, satelite):
        self.orbit.set_celestian(satelite, [self.pos_x, self.pos_y, self.pos_z])

    def add_long_and_lati(self, _window):
        # longitudes = []
        # latitudes = []

        r = self.radius + 0.1
        phi_rng = np.linspace(0., 360., 360, endpoint=True)
        theta_rng = np.linspace(10., 170., 5, endpoint=True)
        cad = 0

        # ШИРОТЫ
        for theta in theta_rng:
            angle = math.pi * (theta / 180)
            theta_sin = math.sin(angle)
            theta_cos = math.cos(angle)

            i = 0
            # пересоздать массив - иначе будет только последний график
            longitudes = np.ndarray(shape=(phi_rng.size, 3), dtype=np.float32)

            for phi in phi_rng:
                angle2 = (math.pi * phi) / 180
                x = r * math.cos(angle2) * theta_sin
                y = r * math.sin(angle2) * theta_sin
                z = r * theta_cos
                longitudes[i] = [x, y, z]
                i = i + 1

            cad = cad + int(255 / theta_rng.size)
            plt = gl.GLLinePlotItem(pos=longitudes, color=pg.glColor(250, cad, cad))
            plt.setGLOptions('opaque')
            self.plots.append(plt)
            self.plots_pos.append(plt.pos)
            _window.addItem(plt)

        phi_rng = np.linspace(0, 180, 2, endpoint=False)
        theta_rng = np.linspace(0, 360, 360, endpoint=True)
        cad = 0
        # ДОЛГОТЫ
        for phi in phi_rng:
            angle = (math.pi * phi) / 180
            phi_sin = math.sin(angle)
            phi_cos = math.cos(angle)

            i = 0
            # пересоздать массив - иначе будет только последний график
            latitudes = np.ndarray((theta_rng.size, 3), dtype=np.float32)

            for theta in theta_rng:
                angle = (math.pi * theta) / 180
                x = r * math.sin(angle) * phi_cos  # CHANGED  +++ self.position_x, self.position_y, self.position_z
                y = r * math.sin(angle) * phi_sin
                z = r * math.cos(angle)
                latitudes[i] = [x, y, z]
                i = i + 1

            cad = cad + (1. / phi_rng.size)
            plt = gl.GLLinePlotItem(pos=latitudes, color=(cad, 1., cad, 1.))
            plt.setGLOptions('opaque')
            self.plots.append(plt)
            self.plots_pos.append(plt.pos)
            _window.addItem(plt)

    def move_long_and_lati(self, pos_x, pos_y, pos_z, angles, parent_planet=None):
        """
        Функция просто двигает широты и долготы за планетой. (широты и долготы нужны, чтобы видеть вращение).
        :param parent_planet:
        :return:
        """
        parent_orbit = None
        if parent_planet:
            parent_orbit = parent_planet.get_orbit()
        self.pos_x = pos_x if not parent_planet else pos_x + parent_planet.pos_x
        self.pos_y = pos_y if not parent_planet else pos_y + parent_planet.pos_y
        self.pos_z = pos_z if not parent_planet else pos_z + parent_planet.pos_z
        rotator = Rotator()
        x, y, z = rotator.get_rotate_matrixes(angles)
        rot_mat = np.dot(np.dot(z, x), z)

        for plot, plot_pos in zip(self.plots, self.plots_pos):
            position = plot_pos.copy()
            # Вращение спутника вокруг оси
            position = np.dot(position, rot_mat)
            # смещение спутника в центр.
            position = np.add(position, [pos_x, pos_y, pos_z])
            # Наклон орбиты
            if parent_orbit:
                position = np.dot(position, parent_orbit.get_orbit_tilt_matrix())
                position = np.add(position, [parent_planet.pos_x, parent_planet.pos_y, parent_planet.pos_z])
            # Установка новых координат
            plot.setData(pos=position)

    def move_sphere_to_pos(self, gl_pos_x, gl_pos_y, gl_pos_z, angles, parent_planet=None):
        parent_orbit = None
        if parent_planet:
            parent_orbit = parent_planet.get_orbit()

        self.pos_x = gl_pos_x if not parent_planet else gl_pos_x + parent_planet.pos_x
        self.pos_y = gl_pos_y if not parent_planet else gl_pos_y + parent_planet.pos_y
        self.pos_z = gl_pos_z if not parent_planet else gl_pos_z + parent_planet.pos_z
        rotator = Rotator()
        x, y, z = rotator.get_rotate_matrixes(angles)
        rot_mat = np.dot(np.dot(z, x), z)
        verts = self.sphere_meshdata.vertexes().copy()
        md = gl.MeshData(vertexes=verts, faces=self.sphere_meshdata.faces(),
                         edges=self.sphere_meshdata.edges(),
                         vertexColors=self.sphere_meshdata.vertexColors(),
                         faceColors=self.sphere_meshdata.faceColors())

        verts = np.dot(verts, rot_mat)
        # Смещаем спутник (планету-спутник) в центр планеты-родителя
        verts = np.add(verts, [gl_pos_x, gl_pos_y, gl_pos_z])
        # Наклон орбиты
        if parent_orbit:
            verts = np.dot(verts, parent_orbit.get_orbit_tilt_matrix())
            verts = np.add(verts, [parent_planet.pos_x, parent_planet.pos_y, parent_planet.pos_z])

        md.setVertexes(verts)
        self.sphere_meshitem.setMeshData(meshdata=md)

        if self.orbit:
            self.orbit.move_orbit_to_pos(self.pos_x, self.pos_y, self.pos_z)

        self.move_long_and_lati(gl_pos_x, gl_pos_y, gl_pos_z, angles, parent_planet)

    def get_orbit(self):
        return self.orbit


class Orbit:
    def __init__(self):
        self.phi_rng = np.linspace(0, 360, 37)
        self.orbit_points = np.ndarray((self.phi_rng.size, 3), dtype=np.float32)

        self.orbit_a = None
        self.orbit_b = None
        self.linear_eccentricity = None
        self.orbit_tilt = None  # Euler rotate matrix
        self.angle = None

        self.pos_x = None
        self.pos_y = None
        self.pos_z = None
        self.orbit_obj = None
        self.celestian_obj = None

    def get_celestian(self):
        return self.celestian_obj

    def move_orbit_to_pos(self, gl_pos_x, gl_pos_y, gl_pos_z):
        tmp = np.add(self.orbit_points, [gl_pos_x - self.linear_eccentricity, gl_pos_y, gl_pos_z])
        self.orbit_obj.setData(pos=tmp)

    def set_celestian(self, new_satellite, parent_planet_pos):
        """
        Ставит спутник на орбиту (не принципиально, но для статичной картинки без таймера -- хорошо)
        Также орбите задает дочерний элемент (объект, движущийся по орбите).
        """
        self.celestian_obj = new_satellite
        if self.orbit_obj:
            new_satellite.move_sphere_to_pos(self.orbit_points[0][0] + parent_planet_pos[0],
                                             self.orbit_points[0][1] + parent_planet_pos[1],
                                             self.orbit_points[0][2] + parent_planet_pos[2],
                                             [0, 0, 0])

    def get_orbit_obj(self):
        return self.orbit_obj

    def move_satellite(self, angle_delta, parent_planet):
        """
        Двигает свой дочерний объект по своим координатам на малый заданный угол rotate_angle.
        :param angle_delta:  Угол поворота в таймере.
        :param parent_planet: родительская планета (по орбите которой движется спутник)
        """
        x_orbit = self.orbit_a * math.cos(angle_delta)
        y_orbit = self.orbit_b * math.sin(angle_delta)
        # TODO:  Думаю, что углы должны быть не тут
        axle_rotation_angle_per_t = 1

        if self.celestian_obj:
            self.celestian_obj.axle_rotate_angle_z = (self.celestian_obj.axle_rotate_angle_z +
                                                      math.radians(axle_rotation_angle_per_t)) % (2 * math.pi)
            angles = [self.celestian_obj.axle_rotate_angle_x,
                      self.celestian_obj.axle_rotate_angle_y,
                      self.celestian_obj.axle_rotate_angle_z]

            self.celestian_obj.move_sphere_to_pos(x_orbit - self.linear_eccentricity, y_orbit, 0, angles, parent_planet)

    def set_orbit_tilt(self, angle):
        """
        Функция составляет матрицу трех поворотов Эйлера для заданного угла.
        :param angle: Угол поворота орбиты.
        """
        orbit_tilt_angles = [math.radians(angle), 0, 0]
        self.angle = angle
        rotator = Rotator()
        x, y, z = rotator.get_rotate_matrixes(orbit_tilt_angles)

        # Поворот на угол Эйлера, хотя в нашей задаче проблемы схлопывания плоскоcтей нет.
        self.orbit_tilt = np.dot(np.dot(z, x), z)

    def create_orbit(self, _window, orbit_a, orbit_b, color, celestian_pos):
        """
        Функция инициализирует орбиту. Задает положение в пространстве, устанавливает инфомрацию об орбите,
        поворачивает точки орбиты в соответствии с заданным углом и уже имеющейся матрицы поворота на углы Эйлера
        orbit_tilt.
        :param celestian_pos: Позиция объекта, к которогу привязана орбита.
        """
        if orbit_a < orbit_b:
            print("The value of orbit_a has to be greater than orbit_b.")
            exit(-1)
        i = 0

        self.pos_x = celestian_pos[0]
        self.pos_y = celestian_pos[1]
        self.pos_z = celestian_pos[2]

        self.orbit_a = orbit_a
        self.orbit_b = orbit_b
        self.linear_eccentricity = math.sqrt(float(math.pow(self.orbit_a, 2)) - float(math.pow(self.orbit_b, 2)))

        # Отрисовка орбиты
        for phi in self.phi_rng:
            angl = math.radians(phi)
            # Переменные без учета смещения в пространстве.
            x = self.orbit_a * math.cos(angl)
            y = self.orbit_b * math.sin(angl)
            z = 0
            self.orbit_points[i] = [x, y, z]
            i = i + 1

        # Получим значения точек орбиты после поворота.
        self.orbit_points = np.dot(self.orbit_points, self.orbit_tilt)

        temp_points = np.add(self.orbit_points, [self.pos_x - self.linear_eccentricity,
                                                 self.pos_y,  self.pos_z])

        self.orbit_obj = gl.GLLinePlotItem(pos=temp_points, color=color)
        _window.addItem(self.orbit_obj)

    def get_orbit_tilt_matrix(self):
        return self.orbit_tilt


if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)

    # surface grid z axis
    grid_z = gl.GLGridItem()
    grid_z.setSize(500, 500, 0)
    grid_z.setSpacing(10, 10, 10)
    grid_z.translate(0, 0, -150)

    window = gl.GLViewWidget()
    window.opts['distance'] = 150
    window.show()
    window.addItem(grid_z)

    #  Оси XYZ
    size = QtGui.QVector3D(100, 100, 100)
    axis = gl.GLAxisItem(size, antialias=False)
    window.addItem(axis)

    earth = CelestialObject()
    earth_color = (0.098, 0.49, 0.345, 1.0)
    earth.create_sphere(6, earth_color, (-20, -20, -20), window)
    earth.set_orbit_params(window, 25, 22, 30)
    earth.set_delta_rot_angle(1)

    moon = CelestialObject()
    moon_color = (0.26, 0.26, 0.26, 1.0)
    moon.create_sphere(1, moon_color, (50, 50, 50), window)
    moon.set_delta_rot_angle(3)
    earth.set_satellite(moon)

    sun = CelestialObject()
    sun_color = (1., 1., 0., 1.0)
    sun.create_sphere(19, sun_color, (40, 40, 40), window)
    sun.set_orbit_params(window, 200, 170, 0)
    sun.set_satellite(earth)

    movement_manager = MovementManager()

    t = QtCore.QTimer()
    t.timeout.connect(movement_manager.update)
    t.start(30)
    sys.exit(app.exec())
