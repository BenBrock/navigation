import numpy as np

import csv
import math

class Coord:
    def __init__(self, lat, lon, code=None):
        self.lat = lat
        self.lon = lon
        self.code = code

    def __str__(self):
        if self.code is None:
            return "(%.7f, %.7f)" % (self.lat, self.lon)
        else:
            return "%s (%.7f, %.7f)" % (self.code, self.lat, self.lon)

    def __repr__(self):
        return 'Coord object: ' + str(self)

class Vor:
    def __init__(self, code, coord, freq, range=60, dme=True):
        self.code = code
        self.coord = coord
        self.freq = freq
        self.range = range
        self.dme = dme
        self.lat = coord.lat
        self.lon = coord.lon

    def __str__(self):
        return '%s (%.7f, %.7f) %s MHz' % (self.code, self.lat, self.lon, self.freq)

    def __repr__(self):
        return 'VOR Object: ' + str(self)

vors = [
        # SF Bay Area
        Vor('SFO', Coord(37.619499, -122.374001), 115.8),
        Vor('OAK', Coord(37.726, -122.223667), 116.8),
        Vor('SAU', Coord(37.855333, -122.522833), 116.2),
        Vor('CCR', Coord(38.045000, -122.045167), 117),
        Vor('OSI', Coord(37.392500, -122.281333), 113.9),
        Vor('SJC', Coord(37.374667, -121.944667), 114.1),
        Vor('PYE', Coord(38.079833, -122.867833), 113.7),
        Vor('SNS', Coord(36.663833, -121.603167), 117.3),
        # Canada, BC
        Vor('YVR', Coord(49.077333, -123.149000), 115.9),
        Vor('YPK', Coord(49.215833, -122.715000), 112.4),
        Vor('YYJ', Coord(48.727000, -123.484333), 113.7),
        Vor('NUW', Coord(48.355000, -122.661833), 113.8),
        Vor('CVV', Coord(48.244667, -122.724500), 117.2),
        Vor('PAE', Coord(47.919833, -122.277833), 110.6),
        Vor('SEA', Coord(47.435333, -122.309667), 116.8),
        Vor('HUH', Coord(48.945333, -122.579333), 113),
        Vor('YDC', Coord(49.381500, -120.373667), 113.9),
        Vor('YWL', Coord(52.237000, -122.169667), 113.6),
        Vor('IBL', Coord(49.955167, -125.330333), 109.1),
        Vor('YZT', Coord(50.684500, -127.366000), 112),
        # Melbourne Area
        Vor('ML', Coord(-37.660000, 144.841667), 114.1),
        Vor('LVT', Coord(-37.859500, 144.749167), 113.1),
        Vor('AV', Coord(-38.048833, 144.459000), 116.1),
        Vor('MNG', Coord(-36.886667, 145.195000), 113.2)
        ]

normal_waypoints = [
        Coord( -38.019667, 144.763333 , 'MPZSF'),
        Coord( -37.819167, 144.932167 , 'BLTB'),
        Coord( -37.720000, 144.963333 , 'PASCO'),
        Coord( -37.802500, 145.151333 , 'MENEC'),
        Coord( -37.718333, 145.466667 , 'BADGR'),
        Coord( -37.516667, 145.750000 , 'MVL'),
        Coord( -37.302500, 145.669500 , 'BOYSE'),
        Coord( -37.206333, 145.835000 , 'TEMIS'),
        Coord( -37.066667, 146.116667 , 'YMFD'),
        Coord( -37.825167, 144.943167 , 'BOLTY'),
        Coord( -37.791667, 144.770000 , 'STEVO')]

waypoints = vors + normal_waypoints

def get_waypoint(code):
    for waypoint in waypoints:
        if waypoint.code == code:
            return waypoint

def vor_rnav(vor: Vor, bearing, distance, dest, altitude = 0, mag_dev = 0, vor_mdev = None):
    if vor_mdev == None:
        vor_mdev = mag_dev
    bearing += vor_mdev

    v = get_waypoint(vor)
    print('Using %s' % (v))

    horizontal_distance = np.sqrt(distance**2 - (altitude / 6076.115)**2)
    current_pos = offset_coord(v, ((bearing + 180) % 360), horizontal_distance)

    d = get_waypoint(dest)

    print('Your pos estimated at %s' % (current_pos))

    magnetic_bearing = (np.round(get_bearing(current_pos, d) - mag_dev) + 360) % 360

    # return int(np.round(get_bearing(current_pos, d)))
    print('%3.0f° at %4.1fnm' % (magnetic_bearing, get_distance(current_pos, d)))

# Distance in nautical miles from `orig` to `dest`
def get_distance(orig: Coord, dest: Coord, radius_of_earth = 3440.0647948):
    import pyproj
    geodesic = pyproj.Geod(ellps='WGS84')
    fwd_azimuth,back_azimuth,distance = geodesic.inv(orig.lon, orig.lat, dest.lon, dest.lat)
    return distance / 1852
    # from geopy.distance import distance
    # return distance((orig.lat, orig.lon), (dest.lat, dest.lon)).nautical

# Bearing in degrees from `orig` to `dest`
def get_bearing(orig: Coord, dest: Coord):
    import pyproj
    geodesic = pyproj.Geod(ellps='WGS84')
    fwd_azimuth,back_azimuth,distance = geodesic.inv(orig.lon, orig.lat, dest.lon, dest.lat)
    return (fwd_azimuth + 360) % 360

def get_closest_vor(vors, coord: Coord):
    best_vor = None
    best_dist = np.inf
    for vor in vors:
        dist = get_distance(coord, vor.coord)
        if dist < best_dist:
            best_vor = vor
            best_dist = dist
    return best_vor

def get_active_vors(vors, coord: Coord):
    active_vors = []
    for vor in vors:
        dist = get_distance(coord, vor.coord)
        if dist < vor.range:
            active_vors.append(vor)
    return sorted(active_vors, key=lambda x: get_distance(coord, x.coord))

def get_offsets(vors, coord):
    return [get_offset(vor.coord, coord) for vor in vors]

# Offset from `orig` to `dest` (distance in nautical miles and bearing in degrees)
def get_offset(orig: Coord, dest: Coord):
    dist = get_distance(orig, dest)
    a = get_bearing(orig, dest)
    return (a, dist)

# Return coordinates of the point at `angle` degrees,
# `dist` nautical miles away from `orig`.
def offset_coord(orig: Coord, angle, dist):
    import pyproj
    geodesic = pyproj.Geod(ellps='WGS84')
    lon,lat,azimuth = geodesic.fwd(orig.lon, orig.lat, angle, dist*1852)
    return Coord(lat, lon)
    import geopy.distance

    lat2, lon2, d = geopy.distance.distance(nautical=dist).destination((orig.lat, orig.lon), angle)
    return Coord(lat2, lon2)

'''
with open('/home/xiii/Downloads/vor_list.csv') as csv_file:
    csv = csv.reader(csv_file, delimiter=',')

    vor_names = list(map(lambda x: x[0], csv))

for vor in vor_names:
    print(vor)
'''

def print_active_vors(waypoint: Coord):
    my_vors = get_active_vors(vors, waypoint)
    print("VORs in Range for waypoint %s:" % (waypoint,))
    for vor in my_vors:
        offset = get_offset(waypoint, vor.coord)
        print('  %5.1f° at %4.1fnm from %s (%s MHz)' % (offset[0], offset[1], vor.code, vor.freq))


def generate_approach(touchdown_coord, runway_degrees, runway_altitude, fix_altitude, increments = 3.5, glideslope_angle=3):
    total_drop = fix_altitude - runway_altitude
    # Total horizontal length of the glidepath, in ft
    glidepath_length = total_drop / np.tan(np.radians(glideslope_angle))

    # Bearing towards the glideslope fix (opposite runway bearing)
    # towards_fix = (runway_degrees + 180 + 80) % 360
    towards_fix = (runway_degrees + 180) % 360

    # Distance of the initial fix, in nm
    fix_distance = glidepath_length / 6076.12

    waypoints = []

    while fix_distance > 0:
        waypoint = offset_coord(touchdown_coord, towards_fix, fix_distance)
        waypoint.code = "TR%.1f" % (fix_distance,)
        waypoint_altitude = 6076.12 * fix_distance * np.tan(np.radians(glideslope_angle)) + runway_altitude
        print("%s, %.1fnm RW. Alt: %dft" % (waypoint, get_distance(waypoint, touchdown_coord), waypoint_altitude))
        waypoints.append(waypoint)
        fix_distance -= increments

whistler_touchdown = Coord(50.147275, -122.941241)
offset = get_offset(Coord(50.147275, -122.941241), Coord(50.153882, -122.927808))
print("%s" % (offset,))

print("Whistler touchdown")
print("%s" % (whistler_touchdown))

generate_approach(whistler_touchdown, 48, 2100, 6000)

def descent_rate_at_angle_speed(glideslope_angle, ground_speed):
    descent_speed = ground_speed*np.tan(np.radians(glideslope_angle))
    # Convert knots to feet/s
    print(descent_speed)
    return descent_speed * 6076.12 / 60

def test_random():
    lat1 = np.random.random()*180 - 90
    lon1 = np.random.random()*180

    lat2 = np.random.random()*180 - 90
    lon2 = np.random.random()*180
    coord = offset_coord(Coord(lat1, lon1), *get_offset(Coord(lat1, lon1), Coord(lat2, lon2)))
    print('%.4f, %.4f' % (lat2, lon2))
    print('%.4f, %.4f' % (coord.lat, coord.lon))
    print(get_distance(Coord(lat2, lon2), coord))
