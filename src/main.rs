mod world_generator;
mod world_visualizer;

use robotics_lib::world::worldgenerator::Generator;
use world_generator::WorldGenerator;
use world_visualizer::WorldVisualizer;


fn main() {
    let mut world_generator = WorldGenerator::new(rand::random(), 300);
    let (world, (_spawn_x, _spawn_y), _weather, _max_score, _score_table) = world_generator.gen();
    WorldVisualizer::visualize(world, 1200, 2);
}