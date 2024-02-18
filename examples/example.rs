use midgard::params::WorldGeneratorParameters;
use robotics_lib::world::world_generator::Generator;
use midgard::WorldGenerator;
use midgard::WorldVisualizer;

fn main() {
    let params = WorldGeneratorParameters {
        seed: 15, // fixed seed
        world_size: 200, // smaller world
        amount_of_rivers: Some(1.2), // more rivers
        amount_of_streets: None, // disable streets
        ..Default::default() // the rest of the parameters keep their default value
    };
    let mut world_generator = WorldGenerator::new(params);
    let (world, _spawn_point, _weather, _max_score, _score_table) = world_generator.gen();
    WorldVisualizer::visualize(world, 800);
}
