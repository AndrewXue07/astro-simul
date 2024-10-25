# random walk until it reaches the edge, where it becomes immobile
# second particles starts at center, does random walk until it sticks to either edge or another particle (if pos == s. pos or t.pos, etc., break loop)
# third loop, etc.
# list storing values; loop through position list every time the position is updated. 

from vpython import *
import numpy as np
import random

def up(pos, posList):
    for position in posList:
        if pos.y + 1 == position:
            # posList.append(pos)
            return -1

    return pos.y + 1

def down(pos, posList):
    for position in posList:
        if pos.y - 1 == position:
            # posList.append(pos)
            return -1

    return pos.y - 1

def left(pos, posList):
    for position in posList:
        if pos.x - 1 == position:
            # posList.append(pos)
            return -1
    return pos.x - 1

def right(pos, posList):
    for position in posList:
        if pos.x + 1 == position:
            # posList.append(pos)
            return -1
    return pos.x + 1


def main():
    posList = []
    canvas(center = vector(50, 50, 0), xmin = 0, xmax = 101, ymin = 0, ymax = 101)
    n = 1000000

    while True: 
        s = vector(50, 50, 0) # sphere(pos = vector(0, 0, 0))
        if s in posList:
            break

        for i in range(n):
            # rate(10)

            randomDir = random.random()
            if randomDir < 0.25: # move up
                proposedDir = up(s, posList)
                if proposedDir == -1: 
                    a = sphere(pos = vector(s.x, s.y - 1, s.z))
                    posList.append(a.pos)
                    break
                elif proposedDir == 101: 
                    a = sphere(pos = vector(s.x, s.y, s.z))
                    posList.append(s)
                    break
                elif proposedDir < 101: # within y bounds
                    s.y += 1
                    
            
            elif randomDir >= 0.25 and randomDir < 0.5: # move down
                proposedDir = down(s, posList)
                if proposedDir == -1: 
                    a = sphere(pos = vector(s.x, s.y + 1, s.z))
                    posList.append(a.pos)
                    break
                elif proposedDir == 0: 
                    a = sphere(pos = vector(s.x, s.y, s.z))
                    posList.append(s)
                    break
                elif proposedDir > 0: # within y bounds
                    s.y -= 1

            elif randomDir > 0.5 and randomDir <= 0.75: # move left
                proposedDir = left(s, posList)
                if proposedDir == -1:
                    a = sphere(pos = vector(s.x + 1, s.y, s.z)) 
                    posList.append(a.pos)
                    break
                elif proposedDir == 0: 
                    a = sphere(pos = vector(s.x, s.y, s.z))
                    posList.append(s)
                    break
                elif proposedDir > 0: # within x bounds
                    s.x -= 1

            elif randomDir > 0.75 and randomDir <= 1: # move right
                proposedDir = right(s, posList)
                if proposedDir == -1: 
                    a = sphere(pos = vector(s.x - 1, s.y, s.z))
                    posList.append(a.pos)
                    break
                elif proposedDir == 101: 
                    a = sphere(pos = vector(s.x, s.y, s.z))
                    posList.append(s)
                    break
                elif proposedDir < 101: # within x bounds
                    s.x += 1



if __name__ == "__main__":
    main()