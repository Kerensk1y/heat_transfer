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
rho_brick = 2000
C_brick = 880
k_brick = 0.67

rho_snow = 200
C_snow = 0.15
k_snow = 2090

rho_wood = 450
C_wood = 0.15
k_wood = 2700

rho_fou = 2500
C_fou = 1.7
k_fou = 840

rho_air = 1.2
C_air = 1005
k_air = 0.026

Q_furnace = 30000
h = 15
Text = 273
dt = 3600  # Временной шаг
t_end = 48*3600  # Конечное время

# Определяем коэффициенты и источники для различных групп
k_values = {5:k_brick, 3: k_air, 1: k_wood, 2: k_wood, 4: k_snow, 6: k_fou, 7: k_brick, 8: k_brick}  # Печь, воздух и кирпич
Q_values = {5: Q_furnace}  # Печь
rho_values = {5:rho_brick, 3: rho_air, 1: rho_wood, 2: rho_wood, 4: rho_snow, 6: rho_fou, 7: rho_brick, 8: rho_brick}
C_values = {5:C_brick, 3: C_air, 4: C_snow, 1: C_wood, 2: C_wood, 8: C_brick, 7: C_brick, 6: C_fou}

class ThermalProperties(UserExpression):
    def init(self, cell_data, property_values, default_value=0, **kwargs):
        self.cell_data = cell_data
        self.property_values = property_values
        self.default_value = default_value
        super().init(**kwargs)

    def eval_cell(self, values, x, cell):
        tag = self.cell_data[cell.index]
        values[0] = self.property_values.get(tag, self.default_value)

    def value_shape(self):
        return ()

# Создаем экземпляры пользовательских классов
k_expr = ThermalProperties(cell_data, k_values, default_value=k_brick, degree=1)
Q_expr = ThermalProperties(cell_data, Q_values, default_value=0, degree=1)
rho_expr = ThermalProperties(cell_data, rho_values, default_value=rho_brick, degree=1)
C_expr = ThermalProperties(cell_data, C_values, default_value=C_brick, degree=1)

# Определяем вариационную форму
u = TrialFunction(V)
v = TestFunction(V)
a = (rho_expr * C_expr / dt) * u * v * dx + k_expr * dot(grad(u), grad(v)) * dx + h * u * v * ds
L = (rho_expr * C_expr / dt) * T * v * dx + Q_expr * v * dx + h * Text * v * ds

# Временной цикл
T_new = Function(V)
t = 0.0

# Диагностика начальных условий
print(f"Initial temperature: {T.vector().get_local().min()} to {T.vector().get_local().max()}")

while t < t_end:
    t += dt
    
    # Диагностика значений выражений через интерполяцию
    k_interpolated = interpolate(k_expr, V)
    Q_interpolated = interpolate(Q_expr, V)
    rho_interpolated = interpolate(rho_expr, V)
    C_interpolated = interpolate(C_expr, V)
    
    # Вывод значений для диагностики
    print(f"Time step: {t:.2f}")
    print(f"Thermal conductivity: {k_interpolated.vector().get_local().min()} to {k_interpolated.vector().get_local().max()}")
    print(f"Heat source: {Q_interpolated.vector().get_local().min()} to {Q_interpolated.vector().get_local().max()}")
    print(f"Density: {rho_interpolated.vector().get_local().min()} to {rho_interpolated.vector().get_local().max()}")
    print(f"Heat capacity: {C_interpolated.vector().get_local().min()} to {C_interpolated.vector().get_local().max()}")
    
    solve(a == L, T_new)
    T.assign(T_new)  # Обновляем температуру
    # Проверка на наличие нечисловых значений
    if np.isnan(T.vector().get_local()).any():
        raise ValueError("Результат содержит некорректные значения (NaN)")
    
    plot_t = plot(T_new)  # сохраняем результат в переменную
    plt.colorbar(plot_t)  # передаем переменную в colorbar
    plt.title(f'Temperature at t={t:.2f}')
    plt.draw()
    plt.pause(0.1)  # пауза для обновления графика
    plt.clf()  # очистка графика для следующего шага

plt.show()  # отображаем конечный результат