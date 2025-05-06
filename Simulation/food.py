"""
Baer Hodge
Jan 7, 2025
"""


class Food:
    """
    Object acting as food for Organism objects. They diffuse through the space
    and respawn when they go out of bounds or if they age past the shelfLife.

    Attributes
    ----------
    spawnPoint : int
        Location in which Food starts and to which Food respawns.
    xCoord : int
        Location in 1-Dimensional space of the simulation.
    age : int
        Number of timesteps since Food has respawned.

    """
    def __init__(self, spawnPoint):
        self.spawnPoint = spawnPoint
        self.xCoord = spawnPoint
        self.age = 0
