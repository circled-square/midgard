use midgard::params::WorldGeneratorParameters;
use robotics_lib::world::world_generator::Generator;
use midgard::WorldGenerator;
use midgard::WorldVisualizer;

fn main() {
    WorldVisualizer::visualize_realtime(| | {
        let params = WorldGeneratorParameters {
            world_size: 500,
            amount_of_streets: Some(0.7),
            ..Default::default()
        };
        let mut world_generator = WorldGenerator::new(params);
        let (world, _spawn_point, _weather, _max_score, _score_table) = world_generator.gen();
        world
    }, 1000)
}
