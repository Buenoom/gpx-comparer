import os
import gpxpy
import folium
from shapely.geometry import LineString, Point
from math import radians, sin, cos, sqrt, atan2
import pandas as pd


# Функция для расчета расстояния между двумя точками (формула гаверсинуса)
def haversine_distance(coord1, coord2):
    """
    Вычисляет расстояние между двумя точками на сфере (формула Haversine).
    :param coord1: tuple (latitude, longitude) первой точки.
    :param coord2: tuple (latitude, longitude) второй точки.
    :return: Расстояние в метрах.
    """
    lat1, lon1 = coord1
    lat2, lon2 = coord2
    R = 6371000  # Радиус Земли в метрах

    phi1, phi2 = radians(lat1), radians(lat2)
    delta_phi = radians(lat2 - lat1)
    delta_lambda = radians(lon2 - lon1)

    a = sin(delta_phi / 2) ** 2 + cos(phi1) * cos(phi2) * sin(delta_lambda / 2) ** 2
    c = 2 * atan2(sqrt(a), sqrt(1 - a))

    return R * c


# Функция для загрузки GPX-файлов
def load_gpx_files(file_paths):
    """
    Загружает координаты из GPX-файлов.
    :param file_paths: Список путей к GPX-файлам.
    :return: Список маршрутов, где каждый маршрут — список кортежей (lat, lon, elevation, time).
    """
    tracks = []
    for file_path in file_paths:
        try:
            with open(file_path, 'r') as gpx_file:
                gpx = gpxpy.parse(gpx_file)
                track_points = [
                    (point.latitude, point.longitude, point.elevation, point.time)
                    for track in gpx.tracks
                    for segment in track.segments
                    for point in segment.points
                ]
                if not track_points:
                    raise ValueError(f"Файл {file_path} не содержит точек.")
                tracks.append(track_points)
        except Exception as e:
            print(f"Ошибка при загрузке файла {file_path}: {e}")
    return tracks


# Функция для анализа маршрута
def analyze_track(track):
    """
    Анализирует маршрут и возвращает его характеристики.
    :param track: Маршрут как список кортежей (lat, lon, elevation, time).
    :return: Словарь с характеристиками маршрута.
    """
    total_distance = 0
    elevations = [point[2] for point in track if point[2] is not None]
    times = [point[3] for point in track if point[3] is not None]

    for i in range(1, len(track)):
        total_distance += haversine_distance(track[i - 1][:2], track[i][:2])

    total_time = (times[-1] - times[0]).total_seconds() if times else 0
    min_elevation = min(elevations) if elevations else None
    max_elevation = max(elevations) if elevations else None
    ascent = sum(max(0, track[i][2] - track[i - 1][2]) for i in range(1, len(track)) if track[i][2] and track[i - 1][2])
    descent = sum(max(0, track[i - 1][2] - track[i][2]) for i in range(1, len(track)) if track[i][2] and track[i - 1][2])

    return {
        "total_distance": total_distance,
        "total_time": total_time,
        "min_elevation": min_elevation,
        "max_elevation": max_elevation,
        "ascent": ascent,
        "descent": descent
    }


# Функция для создания сводной таблицы
def create_summary_table(tracks):
    """
    Создает сводную таблицу с характеристиками маршрутов.
    :param tracks: Список маршрутов.
    :return: DataFrame с характеристиками маршрутов.
    """
    data = []
    for i, track in enumerate(tracks):
        stats = analyze_track(track)
        data.append({
            "Track": f"Track {i + 1}",
            "Total Distance (m)": round(stats["total_distance"], 2),
            "Total Time (s)": round(stats["total_time"], 2),
            "Min Elevation (m)": stats["min_elevation"],
            "Max Elevation (m)": stats["max_elevation"],
            "Ascent (m)": round(stats["ascent"], 2),
            "Descent (m)": round(stats["descent"], 2)
        })

    return pd.DataFrame(data)


# Функция для визуализации маршрутов на карте
def visualize_tracks(tracks, output_map_path="tracks_map.html"):
    """
    Создает карту Folium с визуализацией маршрутов.
    :param tracks: Список маршрутов, где каждый маршрут — список кортежей (lat, lon, elevation, time).
    :param output_map_path: Путь для сохранения HTML-карты.
    """
    if not tracks:
        print("Нет данных для визуализации.")
        return

    # Создаем карту Folium. Берем первую точку первого маршрута для центрирования карты
    first_track = tracks[0]
    m = folium.Map(location=[first_track[0][0], first_track[0][1]], zoom_start=13)

    colors = ['red', 'blue', 'green', 'purple', 'orange']
    for i, track in enumerate(tracks):
        # Оставляем только широту и долготу для PolyLine
        simplified_track = [(point[0], point[1]) for point in track]
        color = colors[i % len(colors)]
        folium.PolyLine(
            locations=simplified_track,
            color=color,
            weight=5,
            opacity=0.8,
            tooltip=f"Маршрут {i + 1}"
        ).add_to(m)

    # Сохраняем карту в HTML файл
    m.save(output_map_path)
    print(f"Карта сохранена в {output_map_path}")


# Функция для расчета пересечений
def calculate_intersections(base_track, other_track):
    """
    Рассчитывает точки пересечения двух маршрутов.
    :param base_track: Базовый маршрут (LineString).
    :param other_track: Вторичный маршрут (LineString).
    :return: Список точек пересечения [(lat, lon)].
    """
    intersections = base_track.intersection(other_track)
    if intersections.is_empty:
        return []
    return list(intersections.coords)


# Функция для расчета среднего отклонения
def calculate_average_deviation(base_track, other_track):
    """
    Рассчитывает среднее отклонение вторичного маршрута от базового.
    :param base_track: Базовый маршрут (LineString).
    :param other_track: Вторичный маршрут (LineString).
    :return: Среднее отклонение в метрах.
    """
    total_deviation = sum(
        Point(point).distance(base_track) * 1000 for point in other_track.coords
    )
    return total_deviation / len(other_track.coords) if len(other_track.coords) > 0 else 0


# Функция для генерации отчета
def generate_report(differences, output_report_path="report.json"):
    """
    Генерирует отчет о различиях между маршрутами.
    :param differences: Список словарей с информацией о различиях.
    :param output_report_path: Путь для сохранения отчета.
    """
    import json
    report_data = {
        "tracks": [
            {
                "track_index": diff["track_index"],
                "intersections": diff["intersections"],
                "avg_deviation": round(diff["avg_deviation"], 2)
            }
            for diff in differences
        ]
    }

    with open(output_report_path, 'w') as report_file:
        json.dump(report_data, report_file, indent=4)
    print(f"Отчет сохранен в {output_report_path}")


# Основная функция
def main():
    print("Введите пути к GPX-файлам по одному. Для завершения ввода введите 'done'.")

    # Список для хранения путей к файлам
    gpx_files = []

    while True:
        file_path = input("Введите путь к GPX-файлу (или 'done' для завершения): ").strip()

        if file_path.lower() == 'done':
            break

        if not os.path.isfile(file_path):
            print(f"Файл '{file_path}' не найден. Попробуйте снова.")
            continue

        gpx_files.append(file_path)

    if not gpx_files:
        print("Не было добавлено ни одного файла. Завершение программы.")
        return

    print(f"Загружено {len(gpx_files)} файлов: {', '.join(gpx_files)}")

    # Загружаем треки
    tracks = load_gpx_files(gpx_files)

    if not tracks:
        print("GPX-файлы не содержат данных или произошла ошибка при загрузке.")
        return

    # Визуализируем маршруты
    visualize_tracks(tracks)

    # Создаем сводную таблицу
    summary_table = create_summary_table(tracks)
    print("\nСводная таблица характеристик маршрутов:")
    print(summary_table)

    # Сохраняем сводную таблицу в CSV
    summary_table.to_csv("summary_table.csv", index=False)
    print("Сводная таблица сохранена в summary_table.csv")

    # Рассчитываем различия
    base_track = LineString([(point[0], point[1]) for point in tracks[0]])  # Базовый маршрут
    differences = []
    for i, track in enumerate(tracks[1:], start=1):
        other_track = LineString([(point[0], point[1]) for point in track])
        intersections = calculate_intersections(base_track, other_track)
        avg_deviation = calculate_average_deviation(base_track, other_track)
        differences.append({
            "track_index": i,
            "intersections": intersections,
            "avg_deviation": avg_deviation
        })

    # Генерируем отчет
    generate_report(differences)


if __name__ == "__main__":
    main()