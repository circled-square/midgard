use midgard::world_generator::WorldGeneratorParameters;
use robotics_lib::world::world_generator::Generator;
use midgard::world_generator::WorldGenerator;
use midgard::world_visualizer::WorldVisualizer;
use std::{thread, time};

// NOT WORKING
fn main() {
    let params = WorldGeneratorParameters {
        world_size: 300,
        ..Default::default()
    };

    let ten_seconds = time::Duration::from_secs(10);

    loop {
        let mut world_generator = WorldGenerator::new(params.clone());
        let (world, (_spawn_x, _spawn_y), _weather, _max_score, _score_table) = world_generator.gen();
        WorldVisualizer::visualize(world, 800); // This closes the program :(
        thread::sleep(ten_seconds);
    }
}
