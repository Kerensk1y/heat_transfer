import meshio
import matplotlib.pyplot as plt
import numpy as np

# Загрузите вашу mesh-файл
mesh = meshio.read("house.msh")

# Инициализация цветов для каждой физической группы
physical_groups = {
    1: 'red',       # Furnace
    2: 'blue',      # Walls
    3: 'green',     # Wooden ceiling
    4: 'cyan',      # Air
    5: 'magenta',   # Snow
    6: 'yellow',    # Foundation
    7: 'orange',    # Ground
    8: 'purple',    # Roof
    9: 'black'      # Outer Lines
}

# Получить метки физических групп для треугольников и линий
triangle_data = mesh.cell_data_dict.get("gmsh:physical", {}).get("triangle", None)
line_data = mesh.cell_data_dict.get("gmsh:physical", {}).get("line", None)

# Инициализация графика
fig, ax = plt.subplots()
ax.set_aspect('equal')

# Рисуем поверхности (треугольники)
if triangle_data is not None:
    for group_id in np.unique(triangle_data):
        color = physical_groups.get(group_id, 'grey')
        indices = np.where(triangle_data == group_id)[0]
        for idx in indices:
            cell = mesh.cells_dict["triangle"][idx]
            coords = mesh.points[cell][:, :2]  # Только x и y координаты
            polygon = plt.Polygon(coords, closed=True, fill=True, edgecolor='k', facecolor=color, alpha=0.5)
            ax.add_patch(polygon)

# Рисуем линии
if line_data is not None:
    for group_id in np.unique(line_data):
        color = physical_groups.get(group_id, 'grey')
        indices = np.where(line_data == group_id)[0]
        for idx in indices:
            cell = mesh.cells_dict["line"][idx]
            coords = mesh.points[cell][:, :2]  # Только x и y координаты
            ax.plot(coords[:, 0], coords[:, 1], color=color, linewidth=2)

# Установить метки и легенду
handles = [plt.Line2D([0], [0], color=color, lw=4) for color in physical_groups.values()]
labels = list(physical_groups.keys())
ax.legend(handles, labels)
ax.set_title('2D Layout of the House')
ax.set_xlabel('X-coordinate')
ax.set_ylabel('Y-coordinate')

# Показать график
plt.show()
