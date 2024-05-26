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
k = 0.67
Q = 30000
h = 15
Text = 273

# Определяем источники и коэффициенты теплопроводности
k_values = np.ones(len(cell_data)) * k
Q_values = np.zeros(len(cell_data))

# Присвоение физических параметров различным группам
for i, tag in enumerate(cell_data):
    if tag == 1:  # Печь
        k_values[i] = k
        Q_values[i] = Q
    elif tag == 4:  # Воздух
        k_values[i] = 0.024
    else:  # Кирпич
        k_values[i] = 0.72

# Определяем вариационную форму
u = TrialFunction(V)
v = TestFunction(V)
a = rho * C * u * v * dx + dt * k * dot(grad(u), grad(v)) * dx
L = (rho * C * T * v + dt * Q * v) * dx + dt * h * Text * v * ds - dt * h * u * v * ds

# Решаем уравнение
T = Function(V)
solve(a == L, T)

# Построение результата
plot(T)
plt.show()
