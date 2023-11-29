mod generator;

extern crate core;

use bevy::prelude::*;
use bevy::window::WindowResolution;
use bevy_pixels::prelude::*;

use robotics_lib::world::tile::Tile;
use robotics_lib::world::tile::TileType;
use robotics_lib::world::worldgenerator::Generator;
use generator::WorldGenerator;

// a bevy resource is a single object accessible by bevy systems, in this case the resource is the world and the system is the function which draws it onto the screen
#[derive(Resource)]
struct WorldMatrixResource {
    world_matrix: Vec<Vec<Tile>>
}

const WINDOW_SIZE : usize = 800;

fn main() {
    let mut world_generator = WorldGenerator::new(200, 0);
    let (world, _pos, _env_conditions, _max_score) = world_generator.gen();

    let window_plugin = WindowPlugin {
        primary_window: Some(Window {
            title: "MIDGARD".into(),
            resolution: WindowResolution::new(WINDOW_SIZE as f32,WINDOW_SIZE as f32),
            resizable: false,
            ..default()
        }),
        ..default()
    };

    App::new()
        .add_plugins((DefaultPlugins.set(window_plugin), PixelsPlugin::default()))
        .add_systems(Update, bevy::window::close_on_esc)
        .add_systems(Draw, draw)
        .insert_resource(WorldMatrixResource{ world_matrix: world })
        .run();
}



fn draw(mut wrapper_query: Query<&mut PixelsWrapper>, world_matrix: Res<WorldMatrixResource>) {
    let world_size = world_matrix.world_matrix.len();

    let Ok(mut wrapper) = wrapper_query.get_single_mut() else { return };
    wrapper.pixels.resize_buffer(world_size as u32, world_size as u32).unwrap();
    let frame = wrapper.pixels.frame_mut();
    assert_eq!(frame.len(), world_size * world_size * 4);

    let to_world_coords = |i| (i%world_size, i/world_size);
    for i in 0..frame.len()/4 {
        let (x,y) = to_world_coords(i);

        let color = match (*world_matrix).world_matrix[x][y] {
            Tile{ tile_type: TileType::Grass, ..} => { [0, 128, 0] }
            Tile{ tile_type: TileType::ShallowWater, ..} => { [0, 0, 128] }
            _ => todo!(),
        };
        frame[i*4..i*4+3].copy_from_slice(&color);
        frame[i*4 + 3] = 0xff; // alpha channel
    }
}