import math
import sys

import numpy as np
import pyqtgraph as pg
import pyqtgraph.opengl as gl
from pyqtgraph.Qt import QtCore, QtGui, QtWidgets


# CONSTANTS
time_global = 0
time_global_delta_angle = 1
axle_rotate_angle_z = 0
#  CONSTANTS


class ObjectsCollector:
    def __init__(self):
        self.list_of_plots = []
        self.list_of_plots_pos = []

    def insert_item_and_pos(self, plot, plot_pos):
        self.list_of_plots.append(plot)
        self.list_of_plots_pos.append(plot_pos)

    def get_plots(self):
        return self.list_of_plots

    def get_plots_poses(self):
        return self.list_of_plots_pos


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
    @staticmethod
    def make_movement(planet, satell_planet):
        global time_global
        global time_global_delta_angle
        global axle_rotate_angle_z
        # Находим данные об орбите небесного тела.
        orbit_a = planet.get_orbit().get_orbit_a()
        orbit_b = planet.get_orbit().get_orbit_b()
        # planet.get_orbit().set_celestian(satell_planet)
        # Находим время, удобное для связи с углом.
        time_global = int((time_global + time_global_delta_angle) % 360)
        angle1 = (time_global * math.pi) / 180  # angle1 in radians

        # Находим точки орбиты.
        x_orbit = orbit_a * math.cos(angle1) - planet.get_orbit().linear_eccentricity
        y_orbit = orbit_b * math.sin(angle1)

        # Задаем углы вращения небесных тел.
        axle_rotate_angle_x = 0
        axle_rotate_angle_y = 0
        axle_rotation_angle_per_t = 5
        axle_rotate_angle_z = (axle_rotate_angle_z + math.radians(axle_rotation_angle_per_t)) % (2 * math.pi)
        angles = [axle_rotate_angle_x, axle_rotate_angle_y, axle_rotate_angle_z]
        rotator = Rotator()
        x, y, z = rotator.get_rotate_matrixes(angles)

        rot_mat = np.dot(np.dot(z, x), z)
        # движение широт и долгот
        for plt1, plt1_pos in zip(satell_planet.collector.get_plots(), satell_planet.collector.get_plots_poses()):
            pos1 = plt1_pos.copy()
            # 2 - Вращение спутника вокруг оси
            pos1 = np.dot(pos1, rot_mat)

            # 4 - Движение спутника по орбите
            pos1 = np.add(pos1, [x_orbit, y_orbit, 0])

            # 5 - Наклон орбиты
            pos1 = np.dot(pos1, planet.get_orbit().get_orbit_tilt_matrix())
            # Установка новых координат

            # pos1 = np.add(pos1, [satell_planet.position_x, satell_planet.position_y, satell_planet.position_z])
            plt1.setData(pos=pos1)

        # сферы
        verts = satell_planet.sphere_meshdata.vertexes().copy()
        md = gl.MeshData(vertexes=verts, faces=satell_planet.sphere_meshdata.faces(),
                         edges=satell_planet.sphere_meshdata.edges(),
                         vertexColors=satell_planet.sphere_meshdata.vertexColors(),
                         faceColors=satell_planet.sphere_meshdata.faceColors())
        verts = np.dot(verts, rot_mat)
        # Движение Земли по орбите
        verts = np.add(verts, [x_orbit, y_orbit, 0])
        # Наклон орбиты
        verts = np.dot(verts, planet.get_orbit().get_orbit_tilt_matrix())
        md.setVertexes(verts)
        satell_planet.sphere_meshitem.setMeshData(meshdata=md)


class CelestialObject:
    def __init__(self):
        self.sphere_meshdata = None  # Data about sphere. Info about it.
        self.sphere_meshitem = None  # Sphere itself
        self.radius = None
        self.latitudes = []  # широта
        self.longitudes = []  # долгота
        self.position_x = None
        self.position_y = None
        self.position_z = None
        self.collector = ObjectsCollector()
        self.color = None
        self.orbit = None

    def create_sphere(self, radius, sphere_color, global_position, _window):
        # Метод создает небесную сферу и устанавливает долготы и широты для того, чтобы видеть вращение.
        self.radius = radius
        self.color = sphere_color
        self.sphere_meshdata = gl.MeshData.sphere(radius=radius, rows=100, cols=100)
        self.sphere_meshitem = gl.GLMeshItem(meshdata=self.sphere_meshdata, smooth=False, color=self.color)
        self.position_x = global_position[0]
        self.position_y = global_position[1]
        self.position_z = global_position[2]
        self.sphere_meshitem.translate(self.position_x, self.position_y, self.position_z)
        self.sphere_meshitem.setGLOptions('additive')  # opaque
        self.add_long_and_lati(_window)

    def set_orbit_params(self, _window, orbit_a, orbit_b, orbit_tilt_angle=0, orbit_color=(0.015, 0.67, 0.82, 1.)):
        # Метод устанавливает параметры орбиты и создает орбиту. Также поворачивает орбиту на заданный угол
        self.orbit = Orbit()
        self.orbit.set_orbit_tilt(orbit_tilt_angle)
        position = [self.position_x, self.position_y, self.position_z]
        self.orbit.create_orbit(_window, orbit_a, orbit_b, orbit_color, position)

    def get_collector(self):
        return self.collector

    def get_meshitem(self):
        return self.sphere_meshitem

    def get_pos_z(self):
        return self.position_z

    def get_pos_y(self):
        return self.position_y

    def get_pos_x(self):
        return self.position_x

    def set_pos(self, pos):
        self.position_x = pos[0]
        self.position_y = pos[1]
        self.position_z = pos[2]

    def set_satellite(self, satellite):
        # Данный метод позволяет небесному телу установить себе спутник. Причем центр спутника изначально
        # Причем центр спутника будет переведен в центр планеты. А потом уже все будет отрисовано в зависимости от
        # расположения орбиты. Перемещает на разницу между двумя объектами
        # (очевидно, для перемещения спутнкка в центр планеты).

        if self.get_orbit():
            self.get_orbit().set_celestian_on_orbit(satellite)

        shifted_pos = [self.position_x - satellite.get_pos_x(),
                       self.position_y - satellite.get_pos_y(),
                       self.position_z - satellite.get_pos_z()]
        self_pos = [self.position_x, self.position_y, self.position_z]
        satellite.set_pos(self_pos)

        plots = satellite.get_collector().get_plots()
        # Переносим широты и долготы в нужное место: центр планеты.
        for plot, plot_pos in zip(satellite.collector.get_plots(), satellite.collector.get_plots_poses()):
            position = plot_pos.copy()
            # смещение спутника в центр.
            position = np.add(position, [shifted_pos[0], shifted_pos[1], shifted_pos[2]])
            # Установка новых координат
            plot.setData(pos=position)

        # Переносим центр спутника в центр планеты, потом уже перерисуется спутник в нужном месте при движении.
        # satellite.get_meshitem().translate(shifted_pos[0], shifted_pos[1], shifted_pos[2])

        verts = satellite.sphere_meshdata.vertexes().copy()
        md = gl.MeshData(vertexes=verts, faces=satellite.sphere_meshdata.faces(),
                         edges=satellite.sphere_meshdata.edges(),
                         vertexColors=satellite.sphere_meshdata.vertexColors(),
                         faceColors=satellite.sphere_meshdata.faceColors())

        # Смещаем спутник (планету-спутник) в центр планеты-родителя
        verts = np.add(verts, [shifted_pos[0], shifted_pos[1], shifted_pos[2]])
        md.setVertexes(verts)
        # satellite.sphere_meshdata = md
        satellite.sphere_meshitem.setMeshData(meshdata=md)

        satellite_orbit = satellite.get_orbit()
        if satellite_orbit:
            satellite_from_satellite = satellite_orbit.get_satellite()
            if satellite_from_satellite:
                satellite.set_satellite(satellite_orbit.get_satellite())
            # satellite_orbit.get_satellite()

    def add_long_and_lati(self, _window):
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
            self.longitudes = np.ndarray(shape=(phi_rng.size, 3), dtype=np.float32)

            for phi in phi_rng:
                angle2 = (math.pi * phi) / 180
                x = r * math.cos(angle2) * theta_sin  ## CHANGED  +++ self.position_x, self.position_y, self.position_z
                y = r * math.sin(angle2) * theta_sin
                z = r * theta_cos
                self.longitudes[i] = [x, y, z]
                i = i + 1

            cad = cad + int(255 / theta_rng.size)
            plt = gl.GLLinePlotItem(pos=self.longitudes, color=pg.glColor(250, cad, cad))
            plt.translate(self.position_x, self.position_y, self.position_z)
            plt.setGLOptions('opaque')
            self.collector.insert_item_and_pos(plt, plt.pos)
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
            self.latitudes = np.ndarray((theta_rng.size, 3), dtype=np.float32)

            for theta in theta_rng:
                angle = (math.pi * theta) / 180
                x = r * math.sin(angle) * phi_cos  # CHANGED  +++ self.position_x, self.position_y, self.position_z
                y = r * math.sin(angle) * phi_sin
                z = r * math.cos(angle)
                self.latitudes[i] = [x, y, z]
                i = i + 1

            cad = cad + (1. / phi_rng.size)
            plt = gl.GLLinePlotItem(pos=self.latitudes, color=(cad, 1., cad, 1.))
            plt.translate(self.position_x, self.position_y, self.position_z)
            plt.setGLOptions('opaque')
            self.collector.insert_item_and_pos(plt, plt.pos)
            _window.addItem(plt)

    def get_orbit(self):
        return self.orbit


class Orbit:
    def __init__(self):
        self.phi_rng = np.linspace(0, 360, 37)
        self.orbit_a = None
        self.orbit_b = None
        self.linear_eccentricity = None
        self.orbit_tilt = None  # rotate matrix
        self.angle = None
        self.orbit_points = np.ndarray((self.phi_rng.size, 3), dtype=np.float32)
        self.center_pos_x = None
        self.center_pos_y = None
        self.center_pos_z = None
        self.orbit_obj = None
        self.celestian_obj = None

    def get_orbit_a(self):
        return self.orbit_a

    def get_orbit_b(self):
        return self.orbit_b

    def get_satellite(self):
        return self.celestian_obj

    def set_celestian_on_orbit(self, new_satellite):
        self.celestian_obj = new_satellite

    def get_orbit_obj(self):
        return self.orbit_obj

    def set_orbit_tilt(self, angle):
        orbit_tilt_angles = [math.radians(angle), 0, 0]
        self.angle = angle
        rotator = Rotator()
        x, y, z = rotator.get_rotate_matrixes(orbit_tilt_angles)

        # Поворот на угол Эйлера, хотя в нашей задаче проблемы схлопывания плоскотей нет.
        self.orbit_tilt = np.dot(np.dot(z, x), z)

    def create_orbit(self, _window, orbit_a, orbit_b, color, celestian_pos):
        # assert orbit_a < orbit_b, "The value of orbit_a has to be greater than orbit_b."
        if orbit_a < orbit_b:
            print("The value of orbit_a has to be greater than orbit_b.")
            exit(-1)
        i = 0
        # self.celestian_obj = celestian

        self.center_pos_x = celestian_pos[0]
        self.center_pos_y = celestian_pos[1]
        self.center_pos_z = celestian_pos[2]
        self.orbit_a = orbit_a
        self.orbit_b = orbit_b
        self.linear_eccentricity = math.sqrt(float(math.pow(self.orbit_a, 2)) - float(math.pow(self.orbit_b, 2)))

        for phi in self.phi_rng:
            angl = math.radians(phi)
            x = self.orbit_a * math.cos(angl)
            y = self.orbit_b * math.sin(angl)
            z = 0
            self.orbit_points[i] = [x, y, z]
            # print(phi, angl, x, y)
            i = i + 1
        self.orbit_points = np.dot(self.orbit_points, self.orbit_tilt)

        for point in self.orbit_points:
            point[0] = point[0] + self.center_pos_x - self.linear_eccentricity  # xshift На нас, не трогать X
            point[1] = point[1] + self.center_pos_y  # y
            point[2] = point[2] + self.center_pos_z  # z_shift Z

        self.orbit_obj = gl.GLLinePlotItem(pos=self.orbit_points, color=color)
        _window.addItem(self.orbit_obj)

    def get_orbit_tilt_matrix(self):
        return self.orbit_tilt


# Слот для timer
def update():
    movement_manager = MovementManager()
    movement_manager.make_movement(earth, moon)
    # movement_manager.make_movement(sun, earth)


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

    earth.create_sphere(6, earth_color, (100, 100, 100), window)
    earth.set_orbit_params(window, 30, 27, 0)
    window.addItem(earth.sphere_meshitem)

    moon = CelestialObject()
    moon_color = (0.26, 0.26, 0.26, 1.0)
    moon.create_sphere(1, moon_color, (20, 20, 20), window)

    earth.set_satellite(moon)
    window.addItem(moon.sphere_meshitem)

    # sun = CelestialObject()
    # sun_color = (1., 0.5, 0., 1.0)
    # sun.create_sphere(19, sun_color, (0, 0, 0), window)
    # sun.set_orbit_params(window, 80, 70, 45)
    #
    # sun.set_satellite(earth)
    # window.addItem(sun.sphere_meshitem)

    t = QtCore.QTimer()
    t.timeout.connect(update)
    t.start(50)
    sys.exit(app.exec())
