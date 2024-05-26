import gmsh
import json
import os

def coord_to_json(name, zero_point, width, height):
    current_dir = os.getcwd()
    print(current_dir)
    A = (zero_point[0], zero_point[1], 0)
    B = (zero_point[0]+width, zero_point[1], 0)
    C = (zero_point[0]+width, zero_point[1]+height, 0)
    D = (zero_point[0], zero_point[1]+height, 0)
    coordinates = (A, B, C, D)
    print(coordinates)
    try:
        with open('points.json', 'r+') as f:
            data = json.load(f)
    except FileNotFoundError:
        with open('points.json', 'w+') as f:
            data = {}
    # Добавить ключ и значение в словарь
    data[name] = coordinates
    # Записать обновленный словарь в файл
    with open('points.json', 'w') as f:
        json.dump(data, f)

def json_to_points():
    walls = {}
    try:
        with open('points.json', 'r') as f:
            data = json.load(f)
        # Открытие и загрузка данных из JSON файла
        filename = 'points.json'
        with open(filename, 'r') as file:
            data = json.load(file)
            print(data)
        #Итерация через координаты и добавление точек
        for name, points in data.items():
            point_ids = []
            for coords in points:
                x, y, z = coords
                point_id = gmsh.model.geo.addPoint(x, y, z)
                point_ids.append(point_id)
            walls[name] = point_ids
        return(walls)
    except FileNotFoundError:
        pass

def house_drawing():
    gmsh.initialize()
    gmsh.model.add("house")
    # Размеры дома
    width = 5.0
    height = 2.0
    wall_thickness = 0.2
    # Размеры печки
    furn_width = 0.8
    furn_height = 1
    # Толщина потолка в комнате
    ceeling_thikness = 0.5
    roof_height = 2.4
    roof_thikness = 0.2
    snow_thikness = 0.2

    coord_to_json('l_wall', [0, 0], wall_thickness, height)
    coord_to_json('r_wall', [wall_thickness + width, 0], wall_thickness, height)
    coord_to_json('furnace', [wall_thickness + width/2-furn_width/2, 0], furn_width, furn_height)
    coord_to_json('ceeling', [-0.2, height], width+0.2*2 + wall_thickness*2, ceeling_thikness)
    walls = json_to_points()
    print(walls)

    # Линии стен
    wall_l1 = gmsh.model.geo.addLine(2, 1)
    wall_l2 = gmsh.model.geo.addLine(1, 4)
    wall_l3 = gmsh.model.geo.addLine(4, 3)
    wall_l4 = gmsh.model.geo.addLine(3, 2)
    wall_r1 = gmsh.model.geo.addLine(6, 5)
    wall_r2 = gmsh.model.geo.addLine(5, 8)
    wall_r3 = gmsh.model.geo.addLine(8, 7)
    wall_r4 = gmsh.model.geo.addLine(7, 6)
    #Линии Печи
    furn_1 = gmsh.model.geo.addLine(10, 9)
    furn_2 = gmsh.model.geo.addLine(9, 12)
    furn_3 = gmsh.model.geo.addLine(12, 11)
    furn_4 = gmsh.model.geo.addLine(11, 10)
    # Земля
    Gnd_A = gmsh.model.geo.addPoint(-0.7, -2, 0)
    Gnd_B = gmsh.model.geo.addPoint(2*wall_thickness + width+0.7, -2, 0)
    Gnd_C = gmsh.model.geo.addPoint(2*wall_thickness + width+0.7, 0, 0)
    Gnd_D = gmsh.model.geo.addPoint(2*wall_thickness + width+0.2, 0, 0)
    Gnd_E = gmsh.model.geo.addPoint(2*wall_thickness + width+0.2, -0.5, 0)
    Gnd_F = gmsh.model.geo.addPoint(-0.2, -0.5, 0)
    Gnd_G = gmsh.model.geo.addPoint(-0.2, 0, 0)
    Gnd_J = gmsh.model.geo.addPoint(-0.7, 0, 0)
    # Фундамент
    Fou_A = Gnd_F
    Fou_B = Gnd_E
    Fou_C = Gnd_D
    # Gnd_D -> 6
    # wall_r1
    # 12 -> 5
    # furn_1
    Fou_D = Gnd_G

    # Создаем линии для фундамента
    fou_l1 = gmsh.model.geo.addLine(Fou_A, Fou_B)
    fou_l2 = gmsh.model.geo.addLine(Fou_B, Fou_C)
    fou_l3 = gmsh.model.geo.addLine(Fou_C, 6)
    fou_l4 = wall_r1 # 6->5
    fou_l5 = gmsh.model.geo.addLine(5, 10)
    fou_l6 = furn_1 # 10->9
    fou_l7 = gmsh.model.geo.addLine(9, 2)
    fou_l8 = wall_l1 # 2->1
    fou_l9 = gmsh.model.geo.addLine(1, Fou_D)
    fou_l10 = gmsh.model.geo.addLine(Fou_D, Fou_A)
    fou_loop = gmsh.model.geo.addCurveLoop([fou_l1, fou_l2, fou_l3, fou_l4, fou_l5, fou_l6, fou_l7, fou_l8, fou_l9, fou_l10])
    gmsh.model.geo.synchronize()
    fou_surface = gmsh.model.geo.addPlaneSurface([fou_loop])

    # Создаем линии для земли
    gnd_l1 = gmsh.model.geo.addLine(Gnd_A, Gnd_B)
    gnd_l2 = gmsh.model.geo.addLine(Gnd_B, Gnd_C)
    gnd_l3 = gmsh.model.geo.addLine(Gnd_C, Gnd_D)
    gnd_l4 = gmsh.model.geo.addLine(Gnd_D, Gnd_E)
    gnd_l5 = gmsh.model.geo.addLine(Gnd_E, Gnd_F)
    gnd_l6 = gmsh.model.geo.addLine(Gnd_F, Gnd_G)
    gnd_l7 = gmsh.model.geo.addLine(Gnd_G, Gnd_J)
    gnd_l8 = gmsh.model.geo.addLine(Gnd_J, Gnd_A)
    gnd_loop = gmsh.model.geo.addCurveLoop([gnd_l1, gnd_l2, gnd_l3, gnd_l4, gnd_l5, gnd_l6, gnd_l7, gnd_l8])
    # Синхронизируем геометрию
    gmsh.model.geo.synchronize()
    gnd_surface = gmsh.model.geo.addPlaneSurface([gnd_loop])
    gmsh.model.geo.synchronize()

    # Поверхность левой стены
    wall_loop = gmsh.model.geo.addCurveLoop([wall_l1, wall_l2, wall_l3, wall_l4])
    gmsh.model.geo.synchronize()
    wall_surface = gmsh.model.geo.addPlaneSurface([wall_loop])
    gmsh.model.geo.synchronize()

    # Поверхность правой стены
    wallr_loop = gmsh.model.geo.addCurveLoop([wall_r1, wall_r2, wall_r3, wall_r4])
    gmsh.model.geo.synchronize()
    wallr_surface = gmsh.model.geo.addPlaneSurface([wallr_loop])
    gmsh.model.geo.synchronize()

    # Поверхность печи
    furn_loop = gmsh.model.geo.addCurveLoop([furn_1, furn_2, furn_3, furn_4])
    gmsh.model.geo.synchronize()
    furn_surface = gmsh.model.geo.addPlaneSurface([furn_loop])
    gmsh.model.geo.synchronize()

    # Крыша внешний слой
    Roof_A = 16
    Roof_B = 15
    Roof_C = gmsh.model.geo.addPoint((width+0.2*2 + wall_thickness*2)/2, height+ceeling_thikness+roof_height, 0)

    # Крыша внутреннний слой
    Roof_D = gmsh.model.geo.addPoint(-0.2+roof_thikness, height+ceeling_thikness, 0)
    Roof_E = gmsh.model.geo.addPoint(width+0.2 + wall_thickness*2-roof_thikness, height+ceeling_thikness, 0)
    Roof_F = gmsh.model.geo.addPoint((width+0.2*2 + wall_thickness*2)/2, height+ceeling_thikness+roof_height-roof_thikness, 0)
    roof_1 = gmsh.model.geo.addLine(Roof_D, Roof_A)
    roof_2 = gmsh.model.geo.addLine(Roof_A, Roof_C)
    roof_3 = gmsh.model.geo.addLine(Roof_C, Roof_B)
    roof_4 = gmsh.model.geo.addLine(Roof_B, Roof_E)
    roof_5 = gmsh.model.geo.addLine(Roof_E, Roof_F)
    roof_6 = gmsh.model.geo.addLine(Roof_F, Roof_D)
    roof_loop = gmsh.model.geo.addCurveLoop([roof_1, roof_2, roof_3, roof_4, roof_5, roof_6])
    roof_surface = gmsh.model.geo.addPlaneSurface([roof_loop])
    gmsh.model.geo.synchronize()

    # Крыша внутреннний слой
    loft_1 = gmsh.model.geo.addLine(Roof_D, Roof_E)
    loft_2 = roof_5
    loft_3 = roof_6
    loft_loop = gmsh.model.geo.addCurveLoop([loft_1, loft_2, loft_3])
    loft_surface = gmsh.model.geo.addPlaneSurface([loft_loop])
    gmsh.model.geo.synchronize()

    # Cнег
    Snow_A = gmsh.model.geo.addPoint(-0.2, snow_thikness+height+ceeling_thikness, 0)
    Snow_B = gmsh.model.geo.addPoint(width+0.2 + wall_thickness*2, snow_thikness+height+ceeling_thikness, 0)
    Snow_C = gmsh.model.geo.addPoint((width+0.2*2 + wall_thickness*2)/2, snow_thikness+height+ceeling_thikness+roof_height, 0)
    snow_1 = gmsh.model.geo.addLine(16, Snow_A)
    snow_2 = gmsh.model.geo.addLine(Snow_A, Snow_C)
    snow_3 = gmsh.model.geo.addLine(Snow_C, Snow_B)
    snow_4 = gmsh.model.geo.addLine(Snow_B, 15)
    snow_5 = gmsh.model.geo.addLine(15, Roof_C)
    snow_6 = gmsh.model.geo.addLine(Roof_C, 16)
    snow_loop = gmsh.model.geo.addCurveLoop([snow_1, snow_2, snow_3, snow_4, snow_5, snow_6])
    snow_surface = gmsh.model.geo.addPlaneSurface([snow_loop])
    gmsh.model.geo.synchronize()

    # Потолок
    # Линии потолка
    ceel_1 = gmsh.model.geo.addLine(13, 4)
    ceel_2 = wall_l3 # 4 -> 3
    ceel_3 = gmsh.model.geo.addLine(3, 8)
    ceel_4 = wall_r3 # 8->7
    ceel_5 = gmsh.model.geo.addLine(7, 14)
    ceel_6 = gmsh.model.geo.addLine(14, 15)
    ceel_7 = roof_4 # 15 -> E
    ceel_8 = gmsh.model.geo.addLine(Roof_E, Roof_D)
    ceel_9 = roof_1 # -> 16
    ceel_10 = gmsh.model.geo.addLine(16, 13)
    ceel_loop = gmsh.model.geo.addCurveLoop([ceel_1, ceel_2, ceel_3, ceel_4, ceel_5, ceel_6, ceel_7, ceel_8, ceel_9, ceel_10])
    gmsh.model.geo.synchronize()
    ceel_surface = gmsh.model.geo.addPlaneSurface([ceel_loop])
    gmsh.model.geo.synchronize()

    room_1 = fou_l5
    room_2 = gmsh.model.geo.addLine(10, 11)
    room_3 = gmsh.model.geo.addLine(11, 12)
    room_4 = gmsh.model.geo.addLine(12, 9)
    room_5 = fou_l7
    room_6 = gmsh.model.geo.addLine(2, 3)
    room_7 = ceel_3
    room_8 = gmsh.model.geo.addLine(8, 5)
    room_loop = gmsh.model.geo.addCurveLoop([room_1, room_2, room_3, room_4, room_5, room_6, room_7, room_8])
    gmsh.model.geo.synchronize()
    room_surface = gmsh.model.geo.addPlaneSurface([room_loop])
    gmsh.model.geo.synchronize()

    gmsh.model.addPhysicalGroup(2, [wall_surface, wallr_surface], name="Walls")
    gmsh.model.addPhysicalGroup(2, [ceel_surface], name="Wooden ceeling")
    gmsh.model.addPhysicalGroup(2, [loft_surface, room_surface], name="Air")
    gmsh.model.addPhysicalGroup(2, [snow_surface], name="Snow")
    gmsh.model.addPhysicalGroup(2, [furn_surface], name="Furnace")
    gmsh.model.addPhysicalGroup(2, [fou_surface], name="Foundament")
    gmsh.model.addPhysicalGroup(2, [gnd_surface], name="Ground")
    gmsh.model.addPhysicalGroup(2, [roof_surface], name="Roof")
    gmsh.model.addPhysicalGroup(1, [gnd_l7, fou_l9, wall_l2, ceel_1, ceel_10, snow_1, snow_2, snow_3, snow_4, ceel_5, ceel_6, wall_r4, fou_l3, gnd_l3], name="Outer Lines")
    

    # Генерация сетки второго порядка (2D)
    gmsh.model.mesh.generate(2)

    # Визуализация геометрии и сетки
    # gmsh.fltk.run()

    # Сохранение геометрии в файл .msh
    gmsh.write("house.msh")

    # Закрыть Gmsh
    gmsh.finalize()

house_drawing()