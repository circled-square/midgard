mod world_generator;

use robotics_lib::world::worldgenerator::Generator;
use world_generator::WorldGenerator;

fn main() {
    let mut world_generator = WorldGenerator::new(1, 800);
    let (world, (spawn_x, spawn_y), weather, max_score) = world_generator.gen();
    world_generator.visualize();
}