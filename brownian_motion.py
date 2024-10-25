from vpython import *
import numpy as np
import random

def up(pos):
    return pos.y + 1

def down(pos):
    return pos.y - 1

def left(pos):
    return pos.x - 1

def right(pos):
    return pos.x + 1


def main():
    canvas(center = vector(50, 50, 0), xmin = 0, xmax = 101, ymin = 0, ymax = 101)
    n = 1000000
    s = sphere(pos = vector(50, 50, 0))

    for i in range(n):
        rate(10)

        randomDir = random.random()
        if randomDir < 0.25: # move up
            proposedDir = up(s.pos)
            if proposedDir < 101: # within bounds
                s.pos.y += 1
        
        elif randomDir >= 0.25 and randomDir < 0.5: # move down
            proposedDir = down(s.pos)
            if proposedDir > 0: # within bounds
                s.pos.y -= 1

        elif randomDir > 0.5 and randomDir <= 0.75: # move left
            proposedDir = left(s.pos)
            if proposedDir > 0: # within bounds
                s.pos.x -= 1

        elif randomDir > 0.75 and randomDir <= 1: # move left
            proposedDir = right(s.pos)
            if proposedDir < 101: # within bounds
                s.pos.x += 1





# particle starts at middle of square grid
# random direction (up, down, left, right)
# if new (proposed) direction outside of bounds vertically or horizontally, don't go forth and instead choose a new direction until valid
# repeat until n = 1,000,000

if __name__ == "__main__":
    main()