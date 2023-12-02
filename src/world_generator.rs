use robotics_lib::world::{worldgenerator::Generator, tile::Tile, environmental_conditions::{EnvironmentalConditions, WeatherType}};

pub struct WorldGenerator {
    seed: u32,
    world_size: usize,
}

impl WorldGenerator {
    pub fn new(seed: u32, world_size: usize) -> Self {
        Self { seed, world_size }
    }

    fn generate_weather() -> EnvironmentalConditions {
        EnvironmentalConditions::new(
            &[WeatherType::Sunny, WeatherType::Rainy],
            1,
            0
        ).unwrap()
    }

    fn generate_altitude(world: &mut Vec<Vec<Tile>>) {

    }

    fn generate_biomes(world: &mut Vec<Vec<Tile>>) {

    }

    fn generate_spawnpoint(&self) -> (usize, usize) {
        (0, 0)
    }
}

impl Generator for WorldGenerator {
    fn gen(&mut self) -> (Vec<Vec<Tile>>, (usize, usize), EnvironmentalConditions, f32) {        
        let mut world = vec![];
        
        WorldGenerator::generate_altitude(&mut world);
        WorldGenerator::generate_biomes(&mut world);

        let weather = WorldGenerator::generate_weather();
        let spawnpoint = self.generate_spawnpoint();
        let score = 100.0;

        (world, spawnpoint, weather, score)
    }
}
