import pymol
from pymol import cmd
from pymol.cgo import *
import numpy as np
import math as m
from typing import List, Dict, Union, Optional, Tuple
from pathlib import Path
from enum import Enum
from dataclasses import dataclass

def unit_vector(vector: np.array) -> np.array:
    """unit_vector generates a unit vector given any numpy array

    :return: normalised vector
    """ 
    return vector / np.linalg.norm(vector)

def bohr_to_angstrom(n):
    return n*0.529177

def closestNumber(n, m):
    q = n%m

    if n%m != 0:
        return int(n - q)
    else:
        return n
class CriticalPointType(Enum):
    BondCriticalPoint = "BCP"
    RingCriticalPoint = "RCP"
    CageCriticalPoint = "CCP"


class CriticalPoint:
    def __init__(self, ty_: CriticalPointType, x: float, y: float, z: float, rho: float, other_atoms: List[str]):
        self.type= ty_
        self.coordinates = np.array([x, y, z])
        self.rho = rho
        self.other_atoms = other_atoms

@dataclass
class Nuclei:
    atomic_number: float
    x: float
    y: float
    z: float

    @property
    def coordinates(self) -> np.ndarray:
        return np.array([self.x, self.y, self.z])

class GradientPath:
    def __init__(self, idx: int, npoints: int, fuck_knows: float) -> None:
        self.idx = idx
        self.npoints = npoints
        self.fuck_knows = fuck_knows  # TODO: rename
        self.points = np.empty((npoints, 4))
    
    @property
    def coords(self) -> np.ndarray:
        return self.points[:,:3]

    @property
    def rho(self) -> np.ndarray:
        return self.points[:,3]


class InterAtomicSurface(list):
    @property
    def points(self) -> np.ndarray:
        return np.vstack(path.coords for path in self)

    @property
    def rho(self) -> np.ndarray:
        return np.hstack([path.rho for path in self])


class IsoDensitySurface:
    def __init__(self, npoints: int):
        self.npoints = npoints
        self.points = np.empty((npoints, 3))
        self.order = np.empty((self.npoints), dtype=int)
        self.intersections = np.empty((self.npoints, 5))
        self.rho_scaling: Dict[float, np.ndarray] = {}

    def points_for_rho_value(self, rho_value: float) -> np.ndarray:
        scaling = self.rho_scaling[rho_value]
        select = (scaling != -1).flatten()
        return (self.points * scaling)[select]


class IasViz:
    def __init__(self, path: Union[str, Path]):
        self.path = Path(path)
        self.atom_name: str = ""
        self.nuclei: Dict[str, Nuclei] = {}
        self.critical_points: List[CriticalPoint] = []
        self.inter_atomic_surfaces: Dict[int, InterAtomicSurface] = {}
        self.iso_density_surface: Optional[IsoDensitySurface] = None
        self.parse()

    @property
    def ias_points(self) -> np.ndarray:
        return np.vstack(ias.points for ias in self.inter_atomic_surfaces.values())

    @property
    def ias_rho(self) -> np.ndarray:
        return np.hstack(ias.rho for ias in self.inter_atomic_surfaces.values())
    
    @property
    def index(self) -> int:
        return int("".join(c for c in self.atom_name if c.isdigit()))

    def iso_density_surface_for_rho(self, rho_value: float) -> np.ndarray:
        return self.iso_density_surface.points_for_rho_value(rho_value) + self.nuclei[self.atom_name].coordinates

    def get_color(self, color: Optional[str], color_list: List[Tuple[int, int, int]]):
        if color is None:
            return color_list[self.index - 1]
        elif '#' in color:
            return [int(color.lstrip('#')[i:i+2], 16) for i in (0, 2, 4)]
        else:
            from PIL import ImageColor
            return ImageColor.getrgb(color)
    
    def parse(self):
        with open(self.path, 'r') as f:
            for line in f:
                if line.startswith("<Atom>"):
                    self.atom_name = next(f).strip()
                    next(f) # </Atom>
                if line.startswith("<Nuclei of Molecule>"):
                    natoms = int(next(f))
                    for _ in range(natoms):
                        record = next(f).split()
                        self.nuclei[record[0]] = Nuclei(*map(float, record[1:]))
                    next(f) # </Nuclei of Molecule>
                if line.startswith("<Electron Density Critical Points in Atomic Surface>"):
                    ncps = int(next(f))
                    for i in range(ncps):
                        record = next(f).split()
                        other_atoms = [] if len(record) <= 5 else record[5:]
                        critical_point = CriticalPoint(CriticalPointType(record[0]), *map(float, record[1:5]), other_atoms)
                        self.critical_points.append(critical_point)
                        if critical_point.type is CriticalPointType.BondCriticalPoint:
                            self.inter_atomic_surfaces[i+1] = InterAtomicSurface()
                    next(f) # </Electron Density Critical Points in Atomic Surface>
                if line.startswith("<IAS Path>"):
                    record = next(f).split()
                    gradient_path = GradientPath(int(record[0]), int(record[1]), float(record[2]))
                    for i in range(gradient_path.npoints):
                        gradient_path.points[i,:] = [x for x in map(float, next(f).split())]
                    next(f) # </IAS Path>
                    self.inter_atomic_surfaces[gradient_path.idx].append(gradient_path)
                if line.startswith("<Intersections of Integration Rays with Atomic Surface>"):
                    record = next(f).split()
                    self.iso_density_surface = IsoDensitySurface(int(record[-1]))
                    for i in range(self.iso_density_surface.npoints):
                        record = next(f).split()
                        self.iso_density_surface.points[i,:] = np.array([x for x in map(float, record[3:6])])
                        self.iso_density_surface.order[i] = int(record[2])
                        self.iso_density_surface.intersections[i,:] = np.array([x for x in map(float, record[6:])])
                    next(f) # </Intersections of Integration Rays with Atomic Surface>
                if line.startswith("<Intersections of Integration Rays With IsoDensity Surfaces>"):
                    record = next(f).split()
                    nrho = int(record[0])
                    rho_values = list(map(float, record[1:1+nrho]))
                    for rho in rho_values:
                        self.iso_density_surface.rho_scaling[rho] = np.empty((self.iso_density_surface.npoints, 1))
                    for i in range(self.iso_density_surface.npoints):
                        record = map(float, next(f).split())
                        for rho_key, rho_value in zip(rho_values, record):
                            self.iso_density_surface.rho_scaling[rho_key][i] = rho_value
                    next(f) # </Intersections of Integration Rays With IsoDensity Surfaces>
       

def create_obj_points(coords, colors, define_normals=False):
    obj = [BEGIN, POINTS]
    div_by_3 = int(closestNumber(len(coords),3))
    coords = coords[:div_by_3]
    colors = colors[:div_by_3]
    for i in range(0, len(coords), 3):
        tri = np.array([coords[i], coords[i+1], coords[i+2]])
        v = tri[1] - tri[0]
        w = tri[2] - tri[0]
        normal = unit_vector(np.cross(v, w))

        obj.append(COLOR)
        obj.extend(colors[i])
        if define_normals:
            obj.append(NORMAL)
            obj.extend(normal)
        obj.append(VERTEX)
        obj.extend(coords[i])

        obj.append(COLOR)
        obj.extend(colors[i+1])
        if define_normals:
            obj.append(NORMAL)
            obj.extend(normal)
        obj.append(VERTEX)
        obj.extend(coords[i+1])

        obj.append(COLOR)
        obj.extend(colors[i+2])
        if define_normals:
            obj.append(NORMAL)
            obj.extend(normal)
        obj.append(VERTEX)
        obj.extend(coords[i+2])
    obj.append(END)
    return obj

def create_obj_wrapper(coords,colors,meshtype='POINTS',define_normals=False):
    '''
    Function that wraps around cgo create objects functions. 
    
    NOTE: ONLY WORKS FOR POINTS (DEFAULT) IN THIS VERSION.
    '''
    if meshtype == 'DOT_LINES':
        obj = create_obj_dot_lines(coords,colors)
    elif meshtype == 'TRIANGLES':
        obj = create_obj_triangles(coords,colors,define_normals)
    else:
        obj = create_obj_points(coords,colors, define_normals)
    return obj

def qtaim_visualise_iasviz(selection = "(all)", file = None, color=None, meshtype = 'POINTS', define_normals = False, rho = 1e-3, transparency = 0.0, *args, **kwargs):
    pymol.color_list = []
    cmd.iterate(selection, 'pymol.color_list.append(color)')
    pymol.color_list = [cmd.get_color_tuple(color) for color in pymol.color_list]

    iasviz = IasViz(file)
    point_color = np.array(iasviz.get_color(color, pymol.color_list))
    #Time to create the mesh object for the iasviz
    rho = float(rho)
    ias_viz_points = iasviz.ias_points
    ias_viz_points = ias_viz_points[iasviz.ias_rho > rho]
    iso_surface_points = iasviz.iso_density_surface_for_rho(rho)
    ias_coords = bohr_to_angstrom(ias_viz_points)
    iso_surface_coords = bohr_to_angstrom(iso_surface_points)
    ias_color = np.full(ias_coords.shape, point_color)
    iso_surface_color = np.full(iso_surface_coords.shape, point_color)
    ias_obj = create_obj_wrapper(ias_coords,ias_color, meshtype, define_normals)
    iso_surface_obj = create_obj_wrapper(iso_surface_coords,iso_surface_color, meshtype, define_normals)
    ias_name = f"{iasviz.atom_name}_ias_{selection}"
    iso_surface_name = f"{iasviz.atom_name}_iso_{selection}"
    cmd.load_cgo(ias_obj, ias_name)
    cmd.load_cgo(iso_surface_obj, iso_surface_name)
    grouped_name = f"{iasviz.atom_name}_iasviz_{selection}"
    cmd.group(grouped_name, f"{ias_name} {iso_surface_name}")
    cmd.group(f"{selection}_iasviz", grouped_name)
    v = cmd.get_view()
    cmd.set('cgo_transparency', transparency, grouped_name)
    cmd.set_view(v)

def qtaim_visualiser(selection = "(all)", file = None, *args, **kwargs):
    file = Path(file)

    if file.is_dir():
        file = [f for f in file.iterdir() if f.is_file() and f.suffix == ".iasviz"]
    else:
        file = [file]

    for f in file:
        qtaim_visualise_iasviz(selection, f, *args, **kwargs)
    
cmd.extend('qtaim_visualiser',qtaim_visualiser)

# tab-completion of arguments
cmd.auto_arg[0]["qtaim_visualiser"] = [cmd.object_sc, "selection=", ", "]