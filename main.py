from objects3d import *
import pygame
from vector import *
import operator
pygame.init()
screen = pygame.display.set_mode((800, 600))
obj1 = "obj//sun"
obj2 = "obj//saturn"
obj3 = "obj//ship"
littleship = Objmesh([], obj3)
bigship = Objmesh([], obj3)
saturn = Objmesh([littleship], obj2)
sun = Objmesh([saturn, bigship], obj1)
done = False
clock = pygame.time.Clock()
coi = VectorN(400, 300, 0)
pos = VectorN(400, 300,500)
up = VectorN(0, 1, 0)
c = Camera(pos, coi, up, 90, 0.1, screen)
theta = 0
zeta = 0
beta = 0
yeta = 0
pos = VectorN(0, 300, 50, 1)
diff = Vector3(1, 1, 1)
spec = Vector3(1, 1, 1)
L1 = Light(pos, diff, spec)
lights = [L1]
geta = 0
while not done:
    evt = pygame.event.poll()
    dt = clock.tick() / 1000
    theta += 10 * dt
    yeta += 10 * dt
    zeta += 100 * dt

    keys = pygame.key.get_pressed()
    sun.transform = Scale(50, 50, 50, False) * RotateX(math.radians(90), False) * RotateZ(math.radians(180), False) * RotateY(math.radians(geta), False) * Translate(400, 300, 0, False)
    saturn.transform = Scale(0.5, 0.5, 0.5, False) * RotateZ(math.radians(yeta), False) * Translate(4, 0, 0, False) * RotateY(math.radians(theta), False)
    littleship.transform = Scale(0.2, 0.2, 0.2, False) * RotateX(math.radians(-90), False) * Translate(2, 0, 0, False) * RotateY(math.radians(zeta), False)
    bigship.transform = Scale(0.3, 0.3, 0.3, False) * RotateX(math.radians(90), False) * Translate(-5, 0, 0, False) * RotateY(math.radians(theta), False)
    if evt.type == pygame.QUIT:
        done = True
    elif evt.type == pygame.KEYDOWN:
        if evt.key == pygame.K_ESCAPE:
            done = True
    if keys[pygame.K_d]:
        geta += 10 * dt
    elif keys[pygame.K_a]:
        geta -= 10 * dt
    screen.fill((0, 0, 0))
    a = sun.pygameDraw(screen, c, lights, identity(4))
    superbiglist = a
    superbiglist.sort(key=operator.itemgetter(2))
    for b in range(len(superbiglist)):
        vlist = []
        temp = []
        for p in range(len(superbiglist[b][1].plist)):
            vlist.append(superbiglist[b][3][superbiglist[b][1].plist[p]] * superbiglist[b][4])
        for v in vlist:
            z = Vector2(v.x, v.y)
            temp.append(z)
        pygame.draw.polygon(screen, (superbiglist[b][0] * 255), temp, 0)
    pygame.display.flip()
pygame.quit()