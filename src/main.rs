mod generator;

extern crate core;

use std::collections::{HashMap, HashSet};
use bevy::prelude::*;
use bevy::window::WindowResolution;
use bevy_pixels::prelude::*;
use noise::{Constant, Multiply, NoiseFn, Perlin, ScalePoint};

use robotics_lib::world::tile::Tile;
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

    let scale = 1. / 50.;
    let scaled_perlin = |seed, scale| -> Multiply<f64, ScalePoint<Perlin>, Constant, 2> {
        Multiply::new(
            ScalePoint::new(
                Perlin::new(seed)
            ).set_scale(scale),
            Constant::new(10.)
        )
    };
    let noise_function =
        noise::Displace::new(
            noise::ScalePoint::new(
                noise::Worley::new(2),
            ).set_scale(scale),
            scaled_perlin(0, scale * 2.),
            scaled_perlin(1, scale * 2.),
            // noise::Constant::new(0.), noise::Constant::new(0.),
            noise::Constant::new(0.), noise::Constant::new(0.),
        );

    for i in 0..frame.len()/4 {
        let (x, y) = to_world_coords(i);
        let noise = noise_function.get([x as f64, y as f64]);

        let color : [u8; 3]= match noise {
            n if n < -0.75 => [0, 0, 128],// sea
            n if n < -0.35 => [128, 128, 0], // desert
            n if n < 0.55 => [0, 128, 0], // plains
            n if n < 0.8 => [0, 64, 0], // forest
            _ => [64, 64, 64], // mountain
        };
        /*
        let color = (noise*127.+127.) as u8;
        let color = [color, color, color];
        */
        /*
        let color = match (*world_matrix).world_matrix[x][y] {
            Tile{ tile_type: TileType::ShallowWater, ..} => { [0, 0, 128] }
            Tile{ tile_type: TileType::DeepWater, ..} => { [0, 0, 64] }
            Tile{ tile_type: TileType::Grass, ..} => { [0, 128, 0] }
            Tile{ tile_type: TileType::Sand, ..} => { [128, 128, 0] }
            Tile{ tile_type: TileType::Hill, ..} => { [92, 92, 0] }
            Tile{ tile_type: TileType::Mountain, ..} => { [64, 64, 64] }
            Tile{ tile_type: TileType::Snow, ..} => { [164, 164, 164] }
            _ => todo!(),
        };
         */
        frame[i*4..i*4+3].copy_from_slice(&color);

        frame[i*4 + 3] = 0xff; // alpha channel
    }
}