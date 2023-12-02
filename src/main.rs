mod world_generator;

use bevy::prelude::*;
use bevy_pixels::prelude::*;
use robotics_lib::world::worldgenerator::Generator;
use world_generator::WorldGenerator;

fn main() {
    App::new()
        .add_plugins((DefaultPlugins, PixelsPlugin::default()))
        .add_systems(Update, bevy::window::close_on_esc)
        // Add systems that draw to the buffer to `Draw` schedule
        // to ensure they are rendered in the current frame.
        .add_systems(Draw, draw)
        .run();
}

/// Draw solid background to window buffer.
fn draw(mut wrapper_query: Query<&mut PixelsWrapper>) {
    //Bevy pixels stuff
    let Ok(mut wrapper) = wrapper_query.get_single_mut() else { return };
    let frame = wrapper.pixels.frame_mut();

    //World generation stuff
    let mut world_generator = WorldGenerator::new(1, 16);
    let (world, (spawn_x, spawn_y), weather, max_score) = world_generator.gen();

    frame.copy_from_slice(&[0x48, 0xb2, 0xe8, 0xff].repeat(frame.len() / 4));
}