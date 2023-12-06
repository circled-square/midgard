mod world_generator;

use robotics_lib::world::worldgenerator::Generator;
use world_generator::WorldGenerator;

fn main() {
    let mut world_generator = WorldGenerator::new(rand::random::<u32>(), 100);
    let (world, (_spawn_x, _spawn_y), _weather, _max_score) = world_generator.gen();
    WorldGenerator::visualize(world, 800);
}