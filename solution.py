import meshio
import numpy as np
from dolfin import *
import matplotlib.pyplot as plt

# Загрузите вашу mesh-файл
mesh = meshio.read("house.msh")

# Конвертируем mesh в формат, понятный FEniCS
points = mesh.points[:, :2]  # Только x и y координаты
cells = mesh.cells_dict["triangle"]
cell_data = mesh.cell_data_dict["gmsh:physical"]["triangle"]

# Создаем mesh в FEniCS
mesh_fenics = Mesh()
editor = MeshEditor()
editor.open(mesh_fenics, "triangle", 2, 2)
editor.init_vertices(len(points))
editor.init_cells(len(cells))

for i, point in enumerate(points):
    editor.add_vertex(i, point)

for i, cell in enumerate(cells):
    editor.add_cell(i, cell)

editor.close()

# Определяем функцию пространства
V = FunctionSpace(mesh_fenics, "P", 1)

# Задаем начальное условие
T0 = Constant(253.0)
T = interpolate(T0, V)

# Параметры
rho = 2000
C = 880
k_furnace = 0.67
k_air = 0.024
k_brick = 0.72
Q_furnace = 30000
h = 15
Text = 273
dt = 0.1  # Временной шаг

# Определяем коэффициенты и источники для различных групп
k_values = {1: k_furnace, 4: k_air, 2: k_brick}  # Печь, воздух и кирпич
Q_values = {1: Q_furnace}  # Печь

class ThermalConductivity(UserExpression):
    def __init__(self, cell_data, k_values, **kwargs):
        self.cell_data = cell_data
        self.k_values = k_values
        super().__init__(**kwargs)

    def eval_cell(self, values, x, cell):
        tag = self.cell_data[cell.index]
        values[0] = self.k_values.get(tag, 0.72)  # Default is brick

    def value_shape(self):
        return ()

class HeatSource(UserExpression):
    def __init__(self, cell_data, Q_values, **kwargs):
        self.cell_data = cell_data
        self.Q_values = Q_values
        super().__init__(**kwargs)

    def eval_cell(self, values, x, cell):
        tag = self.cell_data[cell.index]
        values[0] = self.Q_values.get(tag, 0.0)  # Default is no heat source

    def value_shape(self):
        return ()

# Создаем экземпляры пользовательских классов
k_expr = ThermalConductivity(cell_data, k_values, degree=1)
Q_expr = HeatSource(cell_data, Q_values, degree=1)

# Определяем вариационную форму
u = TrialFunction(V)
v = TestFunction(V)
a = (rho * C / dt) * u * v * dx + k_expr * dot(grad(u), grad(v)) * dx + h * u * v * ds
L = (rho * C / dt) * T * v * dx + Q_expr * v * dx + h * Text * v * ds

# Решаем уравнение
T_new = Function(V)
solve(a == L, T_new)

# Построение результата
plot_t = plot(T_new)  # сохраняем результат в переменную
plt.colorbar(plot_t)  # передаем переменную в colorbar
plt.show()
