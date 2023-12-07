mod world_generator;
mod world_visualizer;

use robotics_lib::world::worldgenerator::Generator;
use world_generator::WorldGenerator;
use world_visualizer::WorldVisualizer;


fn main() {
    let mut world_generator = WorldGenerator::new(1, 100);
    let (world, (_spawn_x, _spawn_y), _weather, _max_score) = world_generator.gen();
    WorldVisualizer::visualize(world, 600);
}