use robotics_lib::world::worldgenerator::Generator;
use robotics_lib::world::tile::Tile;
use robotics_lib::world::environmental_conditions::EnvironmentalConditions;

pub struct WorldGenerator {
    seed : u32,
    world_size : usize,
}

impl Generator for WorldGenerator {
    fn gen(&mut self) -> (Vec<Vec<Tile>>, (usize, usize), EnvironmentalConditions, f32) {
        todo!()
    }
}