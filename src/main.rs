pub mod world_generator;
pub mod world_visualizer;

use crate::world_generator::WorldGeneratorParameters;
use robotics_lib::world::worldgenerator::Generator;
use world_generator::WorldGenerator;
use world_visualizer::WorldVisualizer;

fn main() {
    let world_generator_parameters = WorldGeneratorParameters {
        world_size: 300,
        ..Default::default()
    };
    let mut world_generator = WorldGenerator::new(world_generator_parameters);
    let (world, (_spawn_x, _spawn_y), _weather, _max_score, _score_table) = world_generator.gen();
    WorldVisualizer::visualize(world, 600, 2);
}
