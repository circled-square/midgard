# Midgard

## Steps

- Generate biomes: 
  - Desert
  - Forest
  - Mountain(s)
  - Water
  - Plains
- Generate altitude (with biome bias)
- Generate tiles based on biome and altitude
  - water biome: ShallowWater/DeepWater 
  - mountain biome: Hill/Mountain/Snow/Lava
  - desert biome: Sand/Lava
  - forest biome: Grass
  - plains biome: Grass
- Generate road network
- Generate Teleports
- Generate Contents, based on biome
  - water biome: Water, Fish
  - mountain biome: Rock, Tree
  - desert biome: Fire
  - forest biome: Tree
  - plains biome: Rock
- Generate random Contents
  - Garbage, Coin, Bin, Crate, Bank, Market
  - Not on: DeepWater, ShallowWater


## TO-DO

- Fill water as content
- Generate roads
- Generate random Contents
  - Garbage, Coin, Bin, Crate, Bank, Market
  - Not on: DeepWater, ShallowWater
