from vector import *
from matrix import *
from objects3d import *

class Rasterizer:
    def __init__(self, surf):
        self.surf = surf
        self.sw = surf.get_width()
        self.sh = surf.get_height()

    def render(self):
        pass