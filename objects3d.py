import vector as vector
from matrix import *
import pygame
import math
import random

class Light:
    """There was darkness, then there was the Light class, then it just needed ambient, diffuse, and specular values."""
    def __init__(self, pos, diff, spec):
        self.pos = pos
        self.diff = diff
        self.spec = spec

class Camera:
    """Camer class makes a camera position, along with some other variables to help with 3d rendering"""
    def __init__(self, pos, coi, up, fov, near, surf):
        self.pos = pos
        self.coi = coi
        up = up.normalized()
        self.fov = fov
        self.surf = surf
        self.aspect_ratio = self.surf.get_width() / self.surf.get_height()
        self.local_z = (coi - pos).normalized()
        self.local_x = vector.cross(up, self.local_z).normalized()
        self.local_y = vector.cross(self.local_z, self.local_x)
        self.VPH = 2 * near * math.tan((math.radians(fov/2)))
        self.VPW = self.VPH * self.aspect_ratio
        Q = pos + self.local_z * near
        R = Q - self.local_x * (self.VPW / 2)
        S = R + self.local_y * (self.VPH / 2)
        self.viewplane_origin = S

    def getpixelpos(self, other):
        pcnt_x = other.x / (self.surf.get_width() - 1)
        pcnt_y = other.y / (self.surf.get_height() - 1)
        P = self.viewplane_origin - self.local_y * (pcnt_y * self.VPH) + self.local_x * (pcnt_x * self.VPW)
        return P

    def rotate(self, radians):
        self.local_z *= math.sin(radians)

class Raytracer:
    """Ray trace class puts 3d shapes on a 2d viewplane."""
    def __init__(self, surf, camera, shapelist, lightlist, scene):
        self.camera = camera
        self.shapelist = shapelist
        self.surf = surf
        self.lightlist = lightlist
        self.scene = scene

    def renderOneLine(self, y):
        """Where most of the raytracer happens, the z in range determines the quality of anti aliasing"""
        for x in range(self.surf.get_width()):
            total = 0
            color = vector.Vector3(0, 0, 0)
            temp_color = vector.Vector3(0, 0, 0)
            for z in range(1):
                L = []
                dx = random.uniform(-1, 1)
                dy = random.uniform(-1, 1)
                pixelvector = vector.Vector2(x + dx, y + dy)
                otherpixel = vector.Vector2(x, y)
                e = math.exp(-((dx ** 2) / (2 * (0.333 ** 2)) + (dy ** 2) / (2 * (0.333 ** 2))))
                total += e
                direction = (self.camera.getpixelpos(pixelvector) - self.camera.pos).normalized()
                ray = Ray(color, self.camera.getpixelpos(pixelvector), direction, 1)
                b = self.raycast(ray)
                temp_color += b * e

            avg_color = temp_color / total
            color = avg_color.clamp()
            color = (color * 255).i
            self.surf.set_at(otherpixel.i, color)

    def raycast(self, ray):
        L = []
        maxi = 0
        weight = ray.weight
        temp_color = vector.Vector3(0, 0, 0)
        for s in range(len(self.shapelist)):
            a = self.shapelist[s].rayTest(ray)
            L += a
        if L == []:
            background = self.scene.pairwise_mult(vector.Vector3(0.5, 0.5, 0.5))
            temp_color += background
        else:
            L.sort()
            amb_c = self.scene.pairwise_mult(L[0].otherobject.material["amb"])
            temp_color += amb_c
            for l in self.lightlist:
                n = (l.pos - L[0].point).normalized()
                o = vector.dot(n, L[0].normal)
                if o < 0:
                    pass
                else:
                    r = 2 * o * L[0].normal - n
                    v = (self.camera.pos - L[0].point).normalized()
                    m = vector.dot(v, r)
                    if m > maxi:
                        maxi = m
                    shadow = False
                    c = Ray(L[0].otherobject.material, L[0].point + 0.0001 * L[0].normal, n, weight)
                    abba = ((L[0].point + 0.0001 * L[0].normal) - l.pos).magnitude()
                    for s in self.shapelist:
                        temp = s.rayTest(c)
                        for z in temp:
                            if z.distance <= abba:
                                shadow = True

                            elif z.distance >= abba and shadow is not True:
                                shadow = False

                    if shadow is False:
                        diff_c = o * l.diff.pairwise_mult(L[0].otherobject.material["diff"])
                        if m < 0:
                            spec_c = vector.Vector3(0, 0, 0)
                        else:
                            spec_c = m ** L[0].otherobject.material["shiny"] * \
                                     (l.spec.pairwise_mult(L[0].otherobject.material["spec"]))
                        temp_color += spec_c + diff_c

            if maxi > 0:
                weight -= maxi
                if weight < 0:
                    weight = 0
                z = 2 * (vector.dot(-ray.Dhat, L[0].normal) * L[0].normal - (-ray.Dhat))
                newray = Ray(L[0].otherobject.material, L[0].point + (z), z, 1)
                o = self.raycast(newray)
                temp_color = weight * temp_color + (1 - weight) * o
        return temp_color

class Baseobjects:
    """The class from which all further objects made will derive their color and origin point attributes from."""
    def __init__(self, material, origin):
        """Takes two Vector3 and raises exception if you attempt to pass it something else."""
        self.material = material
        self.origin = origin


class Raycollision:
    """This class will make new instances of itself, every time a ray intersects an object, this will help sort them"""
    def __init__(self, ray, otherobject, distance, point, normal):
        """This takes a ray, another object, and the distance and point at which they intersect to make a new
         Raycollision object"""
        self.ray = ray
        self.otherobject = otherobject
        self.distance = distance
        self.point = point
        self.normal = normal.normalized()

    def __lt__(self, other):
        """This method will sort out the collisions by their distance values."""
        if self.distance < other.distance:
            return True
        else:
            return False


class Ray(Baseobjects):
    """A ray consists of a origin point and a direction which goes off into infinity,
    (a little bigger then screen diagonal in actual application)"""
    def __init__(self, material, origin, direction, weight):
        """Derives origin and color attributes from baseobjects class, but also has direction for Ray class."""
        super().__init__(material, origin)
        self.Dhat = direction.normalized()
        self.weight = weight
        if self.Dhat[1] > 1:
            self.Dhat[1] = 1
        elif self.Dhat[1] < -1:
            self.Dhat[1] = -1

    def getPoint(self, distance):
        """Takes a positive distance and returns a point that far along the ray."""
        z = self.origin + distance * self.Dhat
        return z

    def pygameDraw(self, surf):
        """Draws the ray to screen, need to pass the surf in order for it to work."""
        p = self.getPoint(10000)
        pygame.draw.line(surf, (self.material["diff"] * 255).i, (self.origin.x, self.origin.y), (p.x, p.y))
        pygame.draw.circle(surf, (self.material["diff"] * 255).i, (int(self.origin.x), int(self.origin.y)), 5, 0)


class Plane(Baseobjects):
    """An infinite plane, which in Lab05 will be represented as a line."""
    def __init__(self, material, origin, norm, width, height):
        """Derives the base attributes from baseobject, but also takes a normal direction and point on the plane."""
        super().__init__(material, origin)
        self.norm = norm.normalized()
        self.distance = vector.dot(self.origin, self.norm)
        self.height = height
        self.width = width

    def pygameDraw(self, surf):
        """Draws the plane line to the surf."""
        if abs(self.norm[1]) >= abs(self.norm[0]):
            # Horizontal
            self.B = vector.Vector3(self.width, ((self.distance - self.width * self.norm.x) / self.norm.y), 0)
            self.startpos = vector.Vector3(0, (self.distance / self.norm.y), 0)
        else:
            # Vertical
            self.B = vector.Vector3((self.distance - self.height * self.norm.y) / self.norm.x, self.height, 0)
            self.startpos = vector.Vector3((self.distance / self.norm.x), 0, 0)

        pygame.draw.line(surf, (self.material["diff"] * 255).i, (int(self.startpos.x), int(self.startpos.y)),
                         (int(self.B.x), int(self.B.y)), 1)

    def rayTest(self, R):
        """Returns a list of sorted Raycollisions or none if its empty."""
        Q = vector.dot(R.Dhat, self.norm)
        if Q != 0:
            t = (self.distance - vector.dot(R.origin, self.norm)) / Q
            if t < 0:
                return []
            else:
                p = R.getPoint(t)
                return [Raycollision(R, self, t, p, self.norm)]
        else:
            return[]


class Sphere(Baseobjects):
    """Will draw a circle, but it also really be a sphere mathematically speaking.."""
    def __init__(self, material, origin, radius):
        """Draws from the baseobject's attributes, but also gets a radius for the pygame circle command"""
        super().__init__(material, origin)
        self.radius = radius
        self.pos = origin

    def pygameDraw(self, surf):
        """Draws the sphere to the surf"""
        pygame.draw.circle(surf, (self.material["diff"] * 255).i, (int(self.pos.x), int(self.pos.y)), self.radius, 1)

    def rayTest(self, R):
        """Returns a list of sorted Raycollisions or none if its empty."""
        Q = self.pos - R.origin
        Qpara = (vector.dot(Q, R.Dhat)) * R.Dhat
        Qperp = Q - Qpara
        paradist = vector.dot(Q, R.Dhat)
        perpdist = Qperp.magnitude()
        if perpdist > self.radius:
            # Case I
            return []

        elif perpdist < self.radius and vector.dot(Q, Q) > (self.radius ** 2) and vector.dot(R.Dhat, Q) > 0:
            # Case II
            offset = (self.radius ** 2 - perpdist ** 2) ** (1/2)
            t1 = paradist - offset
            P = (R.getPoint(t1) - self.pos).normalized()
            t2 = paradist + offset
            P2 = (R.getPoint(t2) - self.pos).normalized()
            return[Raycollision(R, self, t1, R.getPoint(t1), P), Raycollision(R, self, t2, R.getPoint(t2), P2)]

        elif perpdist == self.radius and vector.dot(Q, R.Dhat) > 0:
            # Case III
            offset = (self.radius ** 2 - perpdist ** 2) ** (1 / 2)
            t1 = paradist - offset
            P = (R.getPoint(t1) - self.pos).normalized()
            return [Raycollision(R, self, t1, R.getPoint(t1), P)]

        if vector.dot(Q, Q) < (self.radius ** 2):
            # Case IV
            offset = (self.radius ** 2 - perpdist ** 2) ** (1 / 2)
            t2 = paradist + offset
            P = (R.getPoint(t2) - self.pos).normalized()
            return [Raycollision(R, self, t2, R.getPoint(t2), P)]
        else:
            return []


class Triangle(Baseobjects):
    """Will draw a triangle, which will also have 3d to it from the ray detection point of view."""
    def __init__(self, material, origin, point1, point2, n1=0, n2=0, n3=0):
        """Draws from the baseobject's attributes, but also has two more points from which to draw a triangle."""
        super().__init__(material, origin)
        self.point1 = point1
        self.point2 = point2
        self.origin_norm = n1
        self.point1_norm = n2
        self.point2_norm = n3

    def pygameDraw(self, surf):
        """Draws the triangle to the surf"""
        pygame.draw.polygon(surf, (self.material["diff"] * 255).i, ((self.origin.x, self.origin.y), (self.point1.x, self.point1.y),
                                                         (self.point2.x, self.point2.y)), 1)

    def rayTest(self, R):
        """Returns a list of sorted Raycollisions or none if its empty."""
        v = self.point2 - self.origin
        w = self.point1 - self.origin
        l = self.point1 - self.point2
        norm = vector.cross(v, w).normalized()
        P = Plane(self.material, self.origin, norm, 1000, 1000)
        z = P.rayTest(R)
        if z == []:
            return []
        A = (vector.cross(v, w) / 2).magnitude()
        PBC = (vector.cross(z[0].point - self.point2, l) / 2 / A).magnitude()
        PBA = (vector.cross(z[0].point - self.origin, v) / 2 / A).magnitude()
        PAC = (vector.cross(z[0].point - self.origin, w) / 2 / A).magnitude()
        n1 = self.origin_norm * PBC
        n2 = self.point1_norm * PBA
        n3 = self.point2_norm * PAC
        new_norm = (n1 + n2 + n3)
        if PBC + PBA + PAC <= 1.0000001:
            return [Raycollision(R, self, z[0].distance, z[0].point, new_norm)]
        else:
            return []


class Box(Baseobjects):
    """Will draw a Box, which will also have 3d to it from the ray detection point of view."""
    def __init__(self, material, origin, hext):
        """Draws from the baseobject's attributes, but also has a vector representing the half extents."""
        super().__init__(material, origin)
        self.hext = hext
        self.Xhat = vector.Vector3(1, 0, 0)
        self.Yhat = vector.Vector3(0, 1, 0)
        self.Zhat = vector.Vector3(0, 0, 1)

    def Xrotation(self, angle):
        self.Xhat = self.Xhat.rotateX(angle)
        self.Yhat = self.Yhat.rotateX(angle)
        self.Zhat = self.Zhat.rotateX(angle)

    def Yrotation(self, angle):
        self.Xhat = self.Xhat.rotateY(angle)
        self.Yhat = self.Yhat.rotateY(angle)
        self.Zhat = self.Zhat.rotateY(angle)

    def Zrotation(self, angle):
        self.Xhat = self.Xhat.rotateZ(angle)
        self.Yhat = self.Yhat.rotateZ(angle)
        self.Zhat = self.Zhat.rotateZ(angle)

    def pygameDraw(self, surf):
        """Draws the box to the surf"""
        # First Box
        A_local = vector.Vector3(-self.hext.x, -self.hext.y, -self.hext.z)
        B_local = vector.Vector3(self.hext.x, -self.hext.y, -self.hext.z)
        C_local = vector.Vector3(self.hext.x, self.hext.y, -self.hext.z)
        D_local = vector.Vector3(-self.hext.x, self.hext.y, -self.hext.z)
        A_world = self.origin + A_local.x * self.Xhat + A_local.y * self.Yhat + A_local.z * self.Zhat
        B_world = self.origin + B_local.x * self.Xhat + B_local.y * self.Yhat + B_local.z * self.Zhat
        C_world = self.origin + C_local.x * self.Xhat + C_local.y * self.Yhat + C_local.z * self.Zhat
        D_world = self.origin + D_local.x * self.Xhat + D_local.y * self.Yhat + D_local.z * self.Zhat

        # Second Box
        A2_local = vector.Vector3(-self.hext.x, -self.hext.y, self.hext.z)
        B2_local = vector.Vector3(self.hext.x, -self.hext.y, self.hext.z)
        C2_local = vector.Vector3(self.hext.x, self.hext.y, self.hext.z)
        D2_local = vector.Vector3(-self.hext.x, self.hext.y, self.hext.z)
        A2_world = self.origin + A2_local.x * self.Xhat + A2_local.y * self.Yhat + A2_local.z * self.Zhat
        B2_world = self.origin + B2_local.x * self.Xhat + B2_local.y * self.Yhat + B2_local.z * self.Zhat
        C2_world = self.origin + C2_local.x * self.Xhat + C2_local.y * self.Yhat + C2_local.z * self.Zhat
        D2_world = self.origin + D2_local.x * self.Xhat + D2_local.y * self.Yhat + D2_local.z * self.Zhat

        # Drawing the boxes
        pygame.draw.line(surf, (self.material["diff"] * 255).i, (A_world.x, A_world.y), (A2_world.x, A2_world.y), 1)
        pygame.draw.line(surf, (self.material["diff"] * 255).i, (B_world.x, B_world.y), (B2_world.x, B2_world.y), 1)
        pygame.draw.line(surf, (self.material["diff"] * 255).i, (C_world.x, C_world.y), (C2_world.x, C2_world.y), 1)
        pygame.draw.line(surf, (self.material["diff"] * 255).i, (D_world.x, D_world.y), (D2_world.x, D2_world.y), 1)
        pygame.draw.polygon(surf,(self.material["diff"] * 255).i, ((A_world.x, A_world.y), (B_world.x, B_world.y),
                                                        (C_world.x, C_world.y), (D_world.x, D_world.y)), 1)
        pygame.draw.polygon(surf, (self.material["diff"] * 255).i, ((A2_world.x, A2_world.y), (B2_world.x, B2_world.y),
                                                         (C2_world.x, C2_world.y), (D2_world.x, D2_world.y)), 1)

    def rayTest(self, R):
        """Returns a list of sorted Raycollisions or none if its empty."""
        hits = []
        hat_rack = [self.Xhat, -self.Xhat, self.Yhat, -self.Yhat, self.Zhat, -self.Zhat]
        hext_test = [self.hext.x, self.hext.x, self.hext.y, self.hext.y, self.hext.z, self.hext.z]
        for i in range(6):
            norm = hat_rack[i]
            d = vector.dot((self.origin + (hext_test[i] * norm)), norm)
            if vector.dot(R.Dhat, norm) != 0:
                t = (d - vector.dot(R.origin, norm)) / (vector.dot(R.Dhat, norm))
                if t < 0:
                    continue
                Pw = R.getPoint(t)
                Q = Pw - self.origin
                Pl = vector.Vector3(vector.dot(Q, self.Xhat), vector.dot(Q, self.Yhat), vector.dot(Q, self.Zhat))
                if -self.hext.x - 0.0001 <= Pl.x and Pl.x <= self.hext.x + 0.0001:
                    if -self.hext.y - 0.0001 <= Pl.y and Pl.y <= self.hext.y + 0.0001:
                        if -self.hext.z - 0.0001 <= Pl.z and Pl.z <= self.hext.z + 0.0001:
                            hits += [Raycollision(R, self, t, R.getPoint(t), norm)]
        return hits


class Objmesh:
    """Objmesh is not limited to triangles like triangle mesh is."""
    def __init__(self, children, objfile=None):
        """Draws from baseobject's attributes, but also takes an objfile exported from blender."""
        mode = None
        self.transform = identity(4)
        self.facelist = []
        self.mdict = {}
        self.children = children
        print(self.children)
        self.vlist = []
        self.nlist = []
        self.plist = []
        d = {}
        fp = open(objfile + str(".obj"))
        mp = open(objfile + str(".mtl"))
        materials = None
        for line in mp:
            if line[0:6] == "newmtl":
                mode = line.split(" ")
                mode = mode[1]
                mode = mode.split("\n")[0]
                d[str(mode)] = {}

            if line[0:2] == "Ks":
                temp_spec = []
                deparsed = line.split("Ks ")
                deparsed.remove(deparsed[0])
                refined = deparsed[0].split("\n")
                refined.remove(refined[-1])
                almost = refined[0].split(" ")
                for i in range(len(almost)):
                    temp_spec += [float(almost[i])]
                x = vector.VectorN(*temp_spec)
                d[str(mode)]["spec"] = x

            elif line[0:2] == "Ka":
                temp_amb = []
                deparsed = line.split("Ka ")
                deparsed.remove(deparsed[0])
                refined = deparsed[0].split("\n")
                refined.remove(refined[-1])
                almost = refined[0].split(" ")
                for i in range(len(almost)):
                    temp_amb += [float(almost[i])]
                x = vector.VectorN(*temp_amb)
                d[str(mode)]["amb"] = x

            elif line[0:2] == "Kd":
                temp_diff = []
                deparsed = line.split("Kd ")
                deparsed.remove(deparsed[0])
                refined = deparsed[0].split("\n")
                refined.remove(refined[-1])
                almost = refined[0].split(" ")
                for i in range(len(almost)):
                    temp_diff += [float(almost[i])]
                x = vector.VectorN(*temp_diff)
                d[str(mode)]["diff"] = x

            elif line[0:2] == "Ns":
                temp_hard = 0
                deparsed = line.split("Ns ")
                deparsed.remove(deparsed[0])
                refined = deparsed[0].split("\n")
                refined.remove(refined[-1])
                almost = refined[0].split(" ")
                for i in range(len(almost)):
                    temp_hard += float(almost[i])
                x = temp_hard
                d[str(mode)]["hard"] = x
        self.mdict = d

        for line in fp:
            if line[0] == "v" and line[1] != "n":
                temp_v = []
                deparsed = line.split("v ")
                deparsed.remove(deparsed[0])
                refined = deparsed[0].split("\n")
                refined.remove(refined[-1])
                almost = refined[0].split(" ")
                for i in range(len(almost)):
                    temp_v += [float(almost[i])]
                temp_v += [float(1)]
                x = vector.VectorN(*temp_v)
                self.vlist.append(x)
            elif line[0] == "f":
                temp_p = []
                deparsed = line.split("f ")
                deparsed.remove(deparsed[0])
                refined = deparsed[0].split("\n")
                refined.remove(refined[-1])
                almost = refined[0].split(" ")
                for i in range(len(almost)):
                    temp_p += [int(almost[i]) - 1]
                self.plist.append([temp_p, materials])
                b = self.NormalFace(self.vlist, temp_p)
                c = self.getCenter(self.vlist, temp_p)
                self.facelist.append(Face(temp_p, materials,b , c, self.mdict))

            elif line[0:6] == "usemtl":
                z = line.split(" ")[1]
                z = z.strip()
                materials = str(z)



    def NormalFace(self, vects, point):
        vectlist = []
        for p in range(len(point)):
            index = point[p]
            vectlist += [vects[index]]
        origin = vectlist[0]
        p1 = vectlist[1]
        p2 = vectlist[-1]
        v = p1 - origin
        v = VectorN(v.x, v.y, v.z)
        w = p2 - origin
        w = VectorN(w.x, w.y, w.z)
        n = cross(v, w).normalized()
        n = VectorN(n.x, n.y, n.z, 0)
        return n

    def getCenter(self, vects, point):
        vectlist = []
        total = VectorN(0, 0, 0, 0)
        for p in range(len(point)):
            index = point[p]
            vectlist += [vects[index]]
        for z in range(len(vectlist)):
            total += vectlist[z]
        total /= len(vectlist)
        return total

    def pygameDraw(self, surf, camera, lights, M):
        C = self.transform * M
        v = []
        for c in self.children:
            v += c.pygameDraw(surf, camera, lights, C)

        trilist = []
        for f in range(len(self.facelist)):
            plist = []
            n = self.facelist[f].normal * C
            newcent = self.facelist[f].center * C
            d4 = VectorN(camera.local_z.x, camera.local_z.y, camera.local_z.z, 0)
            b = dot(n, d4)
            if b < 0:
                for p in range(len(self.facelist[f].plist)):
                    z = self.vlist[self.facelist[f].plist[p]] * C
                    first = self.facelist[f].material
                    final_color = self.facelist[f].lighting(lights, camera, C)
                trilist.append([final_color, self.facelist[f], newcent.z, self.vlist, C])
        trilist += v
        return trilist


class Face:
    def __init__(self, plist, material, normal, center, md):
        self.md = md
        self.plist = plist
        self.material = material
        self.normal = normal
        self.center = center
        self.scene = Vector3(0.04, 0.04, 0.04)

    def __lt__(self, other):
        if other.center.z > self.center.z:
            return True
        else:
            return False

    def lighting(self, lights, camera, transform):
        temp_color = Vector3(0, 0, 0)
        supern = (self.normal * transform).normalized()
        for l in range(len(lights)):
            center = self.center * transform
            amb = self.scene.pairwise_mult(self.md[self.material]["amb"])
            temp_color += amb
            n = (lights[l].pos - center).normalized()
            o = vector.dot(n, supern)
            if o < 0:
                pass
            else:

                r = 2 * o * supern - n
                i = VectorN(camera.pos.x, camera.pos.y, camera.pos.z, 1)
                v = (i - center).normalized()
                m = vector.dot(v, r)
                diff_c = o * lights[l].diff.pairwise_mult(self.md[self.material]["diff"])
                if m < 0:
                    spec_c = vector.Vector3(0, 0, 0)
                else:
                    spec_c = m ** self.md[self.material]["hard"] * \
                            (lights[l].spec.pairwise_mult(self.md[self.material]["spec"]))
                temp_color += spec_c + diff_c
        return temp_color.clamp()

class TriangleMesh(Baseobjects):
    """TriangleMesh class is for drawing objfiles exported from blender."""
    def __init__(self, material, origin, objfile):
        """Draws from the baseobject's attributes, but also takes an objfile exported from blender."""
        super().__init__(material, origin)
        self.Xhat = vector.Vector3(1, 0, 0)
        self.Yhat = vector.Vector3(0, 1, 0)
        self.Zhat = vector.Vector3(0, 0, 1)
        self.vlist = []
        self.plist = []
        self.nlist = []
        self.max_x = 0
        self.max_y = 0
        self.max_z = 0
        fp = open(objfile)
        for line in fp:
            line.strip()
            if line[0] == "v" and line[1] == "n":
                deparsed = line.split("vn ")
                deparsed.remove(deparsed[0])
                for i in range(1):
                    refined = deparsed[i].split("\n")
                    refined.remove(refined[-1])
                almost = refined[0].split(" ")
                temp_n = []
                for i in range(len(almost)):
                    temp_n += [float(almost[i])]
                x = vector.Vector3(temp_n[0], temp_n[1], temp_n[2])
                self.nlist.append(x)

            if line[0] == "v" and line[1] != "n":
                parsed = line.split(" ")
                temp_v = []
                for i in parsed:
                    if i != "v":
                        temp_v.append(float(i))
                x = vector.Vector3(temp_v[0], temp_v[1], temp_v[2])
                self.vlist.append(x)
            if line[0] == "f":
                unparsed = line.split("f ")
                for i in unparsed:
                    if i == "":
                        unparsed.remove(i)
                    refined = unparsed[0].split(" ")
                    refined[-1] = refined[-1].rstrip("\n")
                    rl = []
                    for i in refined:
                        almost = i.split("//")
                        rl += almost
                    temp_p = []
                    for r in rl:
                        temp_p.append(int(r) - 1)
                    self.plist.append(temp_p)
        for i in self.vlist:
            if abs(i.x) > self.max_x:
                self.max_x = abs(i.x)
            if abs(i.y) > self.max_y:
                self.max_y = abs(i.y)
            if abs(i.z) > self.max_z:
                self.max_z = abs(i.z)
        hext = vector.Vector3(self.max_x, self.max_y, self.max_z)
        self.boundbox = Box(self.material, self.origin, hext)
        self.tlist = []
        for i in range(len(self.plist)):
            a = self.plist[i][0]
            b = self.plist[i][2]
            c = self.plist[i][4]
            a1 = self.plist[i][1]
            b2 = self.plist[i][3]
            c2 = self.plist[i][5]
            n1 = self.nlist[a1]
            n2 = self.nlist[b2]
            n3 = self.nlist[c2]
            v1 = self.vlist[a]
            v2 = self.vlist[b]
            v3 = self.vlist[c]
            p1 = self.origin + v1.x * self.Xhat + v1.y * self.Yhat + v1.z * self.Zhat
            p2 = self.origin + v2.x * self.Xhat + v2.y * self.Yhat + v2.z * self.Zhat
            p3 = self.origin + v3.x * self.Xhat + v3.y * self.Yhat + v3.z * self.Zhat
            triangle = Triangle(self.material, p1, p2, p3, n1, n2, n3)
            self.tlist.append(triangle)

    def pygameDraw(self, surf):
        """Draws the triangle mesh to the surf."""
        self.boundbox.pygameDraw(surf)

        for i in self.tlist:
            i.pygameDraw(surf)

    def Xrotation(self, angle):
        self.Xhat = self.Xhat.rotateX(angle)
        self.Yhat = self.Yhat.rotateX(angle)
        self.Zhat = self.Zhat.rotateX(angle)
        self.boundbox.Xrotation(angle)
        self.tlist = []
        for i in range(len(self.plist)):
            a = self.plist[i][0]
            b = self.plist[i][2]
            c = self.plist[i][4]
            a1 = self.plist[i][1]
            b2 = self.plist[i][3]
            c2 = self.plist[i][5]
            n1 = self.nlist[a1]
            n2 = self.nlist[b2]
            n3 = self.nlist[c2]
            v1 = self.vlist[a]
            v2 = self.vlist[b]
            v3 = self.vlist[c]
            p1 = self.origin + v1.x * self.Xhat + v1.y * self.Yhat + v1.z * self.Zhat
            p2 = self.origin + v2.x * self.Xhat + v2.y * self.Yhat + v2.z * self.Zhat
            p3 = self.origin + v3.x * self.Xhat + v3.y * self.Yhat + v3.z * self.Zhat
            triangle = Triangle(self.material, p1, p2, p3, n1, n2, n3)
            self.tlist.append(triangle)

    def Yrotation(self, angle):
        self.Xhat = self.Xhat.rotateY(angle)
        self.Yhat = self.Yhat.rotateY(angle)
        self.Zhat = self.Zhat.rotateY(angle)
        self.boundbox.Yrotation(angle)
        self.tlist = []
        for i in range(len(self.plist)):
            a = self.plist[i][0]
            b = self.plist[i][2]
            c = self.plist[i][4]
            a1 = self.plist[i][1]
            b2 = self.plist[i][3]
            c2 = self.plist[i][5]
            n1 = self.nlist[a1]
            n2 = self.nlist[b2]
            n3 = self.nlist[c2]
            v1 = self.vlist[a]
            v2 = self.vlist[b]
            v3 = self.vlist[c]
            p1 = self.origin + v1.x * self.Xhat + v1.y * self.Yhat + v1.z * self.Zhat
            p2 = self.origin + v2.x * self.Xhat + v2.y * self.Yhat + v2.z * self.Zhat
            p3 = self.origin + v3.x * self.Xhat + v3.y * self.Yhat + v3.z * self.Zhat
            triangle = Triangle(self.material, p1, p2, p3, n1, n2, n3)
            self.tlist.append(triangle)

    def Zrotation(self, angle):
        self.Xhat = self.Xhat.rotateZ(angle)
        self.Yhat = self.Yhat.rotateZ(angle)
        self.Zhat = self.Zhat.rotateZ(angle)
        self.boundbox.Zrotation(angle)
        self.tlist = []
        for i in range(len(self.plist)):
            a = self.plist[i][0]
            b = self.plist[i][2]
            c = self.plist[i][4]
            a1 = self.plist[i][1]
            b2 = self.plist[i][3]
            c2 = self.plist[i][5]
            n1 = self.nlist[a1]
            n2 = self.nlist[b2]
            n3 = self.nlist[c2]
            v1 = self.vlist[a]
            v2 = self.vlist[b]
            v3 = self.vlist[c]
            p1 = self.origin + v1.x * self.Xhat + v1.y * self.Yhat + v1.z * self.Zhat
            p2 = self.origin + v2.x * self.Xhat + v2.y * self.Yhat + v2.z * self.Zhat
            p3 = self.origin + v3.x * self.Xhat + v3.y * self.Yhat + v3.z * self.Zhat
            triangle = Triangle(self.material, p1, p2, p3, n1, n2, n3)
            self.tlist.append(triangle)

    def rayTest(self, R):
        """Returns a list of sorted Raycollisions or none if its empty."""
        total = []
        if self.boundbox.rayTest(R) != []:
            for t in self.tlist:
                if t.rayTest(R) != []:
                    total += t.rayTest(R)
            return total
        else:
            return[]