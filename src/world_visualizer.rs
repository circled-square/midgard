use bevy::app::{App, Update};
use bevy::DefaultPlugins;
use bevy::prelude::*;
use bevy::window::WindowResolution;
use bevy_pixels::{PixelsPlugin, PixelsWrapper};
use bevy_pixels::prelude::*;
use robotics_lib::world::tile::{Content, Tile, TileType};

#[derive(Resource)]
struct WorldMatrixResource {
    matrix: Vec<Vec<Tile>>
}

#[derive(Resource)]
struct PixelScalingResource {
    pixel_scaling: usize
}

pub struct WorldVisualizer {}

impl WorldVisualizer {
    fn draw_window(mut wrapper_query: Query<&mut PixelsWrapper>, world: Res<WorldMatrixResource>, pixel_scaling: Res<PixelScalingResource>) {
        let pixel_scaling = pixel_scaling.pixel_scaling;

        //Bevy pixels stuff
        let Ok(mut wrapper) = wrapper_query.get_single_mut() else { return };
        let frame_row_len = world.matrix.len() * pixel_scaling;
        wrapper.pixels.resize_buffer(frame_row_len as u32, frame_row_len as u32).unwrap();
        let frame = wrapper.pixels.frame_mut();

        assert_eq!(frame_row_len * frame_row_len * 4, frame.len());

        for i in 0..world.matrix.len() {
            for j in 0..world.matrix.len() {
                for x in 0..pixel_scaling {
                    for y in 0..pixel_scaling {
                        let color = if x == pixel_scaling/2 && y == pixel_scaling/2 {
                            Self::color_tile_content(&world.matrix[i][j]).unwrap_or(Self::color_tile(&world.matrix[i][j]))
                        } else {
                            Self::color_tile(&world.matrix[i][j])
                        };
                        let pixel_coords = (i*pixel_scaling + x, j*pixel_scaling + y);
                        let pixel_index = pixel_coords.0 + frame_row_len * pixel_coords.1;
                        frame[pixel_index * 4..pixel_index * 4 + 4].copy_from_slice(&color);
                    }
                }
            }
        }
    }
    pub fn visualize(world: Vec<Vec<Tile>>, resolution: usize, pixel_scaling: usize) {
        let mut resolution = WindowResolution::new(resolution as f32, resolution as f32);
        resolution.set_scale_factor_override(Some(1.0));

        let window_plugin = WindowPlugin {
            primary_window: Some(Window {
                title: "MIDGARD".into(),
                resolution,
                resizable: false,
                ..default()
            }),
            ..default()
        };


        App::new()
            .add_plugins((DefaultPlugins.set(window_plugin), PixelsPlugin::default()))
            .add_systems(Update, bevy::window::close_on_esc)
            .add_systems(Draw, Self::draw_window)
            .insert_resource(WorldMatrixResource { matrix: world })
            .insert_resource(PixelScalingResource { pixel_scaling })
            .run();
    }
    fn color_tile_content(tile: &Tile) -> Option<[u8; 4]> {
        return match tile.content{
            Content::Fish(_) => Some([255, 153, 102, 255]),
            Content::Tree(_) => Some([0, 51, 0, 255]),
            Content::Rock(_) => Some([128, 128, 128, 255]),
            Content::Fire => Some([255, 102, 0, 255]),
            _ => None
        }
    }
    fn color_tile(tile : &Tile) -> [u8; 4] {
        return match tile.tile_type {
            TileType::DeepWater => [0, 0, 127, 255],
            TileType::ShallowWater => [0, 0, 255, 255],
            TileType::Grass => [0, 255, 0, 255],
            TileType::Sand => [255, 255, 0, 255],
            TileType::Lava => [255, 0, 0, 255],
            TileType::Hill => [0, 127, 0, 255],
            TileType::Mountain => [153, 102, 51, 255],
            TileType::Snow => [255, 255, 255, 255],
            _ => [0, 0, 0, 255]
        }
    }
}
