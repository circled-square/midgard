mod fall_off_noise_function;

use noise::NoiseFn;
use fall_off_noise_function::FallOffNoiseFunction;
use robotics_lib::world::worldgenerator::Generator;
use robotics_lib::world::tile::Tile;
use robotics_lib::world::environmental_conditions::EnvironmentalConditions;
use robotics_lib::world::environmental_conditions::WeatherType;
use robotics_lib::world::tile::TileType;
use robotics_lib::world::tile::Content;


pub struct WorldGenerator {
    world_size : usize, seed : u32
}

impl WorldGenerator {
    pub fn new(world_size : usize, seed : u32) -> Self {
        Self { world_size, seed }
    }
}

impl Generator for WorldGenerator {
    fn gen(&mut self) -> (Vec<Vec<Tile>>, (usize, usize), EnvironmentalConditions, f32) {
        // parameters
        let scale = 1./20.; // scaling of the noise texture
        let deep_water_level = -0.9;
        let shallow_water_level = -0.65;
        let sand_level = -0.50;
        let hill_level = 0.15;
        let mountain_level = 0.55;
        let snow_level = 0.85;

        let noise_function =
            FallOffNoiseFunction::new(
                noise::ScalePoint::new(
                    noise::Perlin::new(self.seed)
                ).set_scale(scale),
                (self.world_size/2) as f64, (self.world_size/2) as f64
            );

        let mut world = vec![];
        for i in 0..self.world_size {
            let mut row = vec![];
            for j in 0..self.world_size {
                let height = noise_function.get([i as f64, j as f64]);
                let tile = match height {
                    h if h < deep_water_level => Tile {
                        tile_type: TileType::DeepWater,
                        content: Content::None,
                        elevation: 0,
                    },
                    h if h < shallow_water_level => Tile {
                        tile_type: TileType::ShallowWater,
                        content: Content::None,
                        elevation: 0,
                    },
                    h if h < sand_level => Tile {
                        tile_type: TileType::Sand,
                        content: Content::None,
                        elevation: 0,
                    },
                     h if h < hill_level => Tile {
                        tile_type: TileType::Grass,
                        content: Content::None,
                        elevation: 0,
                    },
                    h if h < mountain_level => Tile {
                        tile_type: TileType::Hill,
                        content : Content::None,
                        elevation: 0,
                    },
                    h if h < snow_level => Tile {
                        tile_type: TileType::Mountain,
                        content : Content::None,
                        elevation: 0,
                    },
                    _ => Tile {
                        tile_type: TileType::Snow,
                        content : Content::None,
                        elevation: 0,
                    },
                };
                row.push(tile);
            }
            world.push(row);
        }

        let pos = (self.world_size/2,self.world_size/2); // should check if it is water or not
        let env_conditions = EnvironmentalConditions::new(&[WeatherType::Sunny], 30, 8);
        let max_score = 100.;

        (world, pos, env_conditions, max_score)
    }
}
