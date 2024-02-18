use std::ops::Deref;
use bevy::app::{App, Update};
use bevy::prelude::*;
use bevy::window::WindowResolution;
use bevy::DefaultPlugins;
use bevy_pixels::prelude::*;
use bevy_pixels::{PixelsPlugin, PixelsWrapper};
use robotics_lib::world::tile::{Content, Tile, TileType};
use embed_doc_image::embed_doc_image;

#[derive(Resource)]
struct WorldMatrixResource {
    matrix: Vec<Vec<Tile>>,
}

#[derive(Resource)]
struct PixelScalingResource {
    pixel_scaling: usize,
}

#[derive(Resource)]
struct RegenWorldFuncResource {
    func : fn() -> Vec<Vec<Tile>>,
}

#[derive(Resource)]
struct TimeCounterResource {
    t : f64,
    last_t : f64,
}

#[embed_doc_image("world_render", "misc/world_render.png")]
/// Provides the `visualize` method to render the world
/// 
/// # Render example
/// 
/// ![World render image][world_render]
/// 
/// <style>
/// .legend {
///     font-family: Arial, sans-serif;
///     padding: 10px;
///     border: 1px solid #ccc;
///     border-radius: 5px;
/// }
/// 
/// .legend-item {
///     display: flex;
///     align-items: center;
///     margin-bottom: 5px;
/// }
/// 
/// .color-box {
///     border-style: solid;
///     border-width: 1px;
///     border-color: white;
///     box-shadow: 2px 2px 2px 0px #b3b3b3;
///     width: 20px;
///     height: 20px;
///     margin-right: 10px;
/// }
/// </style>
/// 
/// # Legend
/// 
/// ## Tiles
/// 
/// <div class="legend">
///     <div class="legend-item">
///         <div class="color-box" style="background-color: #000066;"></div>
///         Deep Water
///     </div>
///     <div class="legend-item">
///         <div class="color-box" style="background-color: #3333ff;"></div>
///         Shallow Water
///     </div>
///     <div class="legend-item">
///         <div class="color-box" style="background-color: #009933;"></div>
///         Grass
///     </div>
///     <div class="legend-item">
///         <div class="color-box" style="background-color: #ffcc00;"></div>
///         Sand
///     </div>
///     <div class="legend-item">
///         <div class="color-box" style="background-color: #cf1020;"></div>
///         Lava
///     </div>
///     <div class="legend-item">
///         <div class="color-box" style="background-color: #669900;"></div>
///         Hill
///     </div>
///     <div class="legend-item">
///         <div class="color-box" style="background-color: #996633;"></div>
///         Mountain
///     </div>
///     <div class="legend-item">
///         <div class="color-box" style="background-color: #ccffff;"></div>
///         Snow
///     </div>
///     <div class="legend-item">
///         <div class="color-box" style="background-color: #66ffff;"></div>
///         Teleport
///     </div>
///     <div class="legend-item">
///         <div class="color-box" style="background-color: #404040;"></div>
///         Street
///     </div>
/// </div>
/// 
/// ## Contents
/// 
/// <div class="legend">
///     <div class="legend-item">
///         <div class="color-box" style="background-color: #ff9933;"></div>
///         Fish
///     </div>
///     <div class="legend-item">
///         <div class="color-box" style="background-color: #336600;"></div>
///         Tree
///     </div>
///     <div class="legend-item">
///         <div class="color-box" style="background-color: #999966;"></div>
///         Rock
///     </div>
///     <div class="legend-item">
///         <div class="color-box" style="background-color: #000000;"></div>
///         Bush (missing visualization)
///     </div>
///     <div class="legend-item">
///         <div class="color-box" style="background-color: #ff6600;"></div>
///         Fire
///     </div>
///     <div class="legend-item">
///         <div class="color-box" style="background-color: #333300;"></div>
///         Garbage
///     </div>
///     <div class="legend-item">
///         <div class="color-box" style="background-color: #cc9900ff;"></div>
///         Coin
///     </div>
///     <div class="legend-item">
///         <div class="color-box" style="background-color: #ff6600;"></div>
///         Bin
///     </div>
///     <div class="legend-item">
///         <div class="color-box" style="background-color: #666633;"></div>
///         Crate
///     </div>
///     <div class="legend-item">
///         <div class="color-box" style="background-color: #cc00ff;"></div>
///         Market
///     </div>
///     <div class="legend-item">
///         <div class="color-box" style="background-color: #000000;"></div>
///         Bank (missing visualization)
///     </div>
///     <div class="legend-item">
///         <div class="color-box" style="background-color: #000000;"></div>
///         Building (missing visualization)
///     </div>
///     <div class="legend-item">
///         <div class="color-box" style="background-color: #000000;"></div>
///         Scarecrow (missing visualization)
///     </div>
///     <div class="legend-item">
///         <div class="color-box" style="background-color: #000000;"></div>
///         Jolly Block (missing visualization)
///     </div>
/// </div>
/// 
/// # Examples
/// 
/// ```
/// use robotics_lib::world::world_generator::Generator;
/// use midgard::WorldGenerator;
/// use midgard::params::WorldGeneratorParameters;
/// use midgard::WorldVisualizer;
///
/// # fn main() {
/// let mut world_generator = WorldGenerator::new(WorldGeneratorParameters::default());
/// // Generate the world
/// let (world, _spawn_point, _weather, _max_score, _score_table) = world_generator.gen();
///
/// // Use 'WorldVisualizer::visualize' to render the world at the specified resolution
/// WorldVisualizer::visualize(world, 600);
/// # }
/// ```
pub struct WorldVisualizer {}

impl WorldVisualizer {
    fn draw_window(mut wrapper_query: Query<&mut PixelsWrapper>, world: Res<WorldMatrixResource>, pixel_scaling: Res<PixelScalingResource>) {
        let pixel_scaling = pixel_scaling.pixel_scaling;

        //Bevy pixels stuff
        let Ok(mut wrapper) = wrapper_query.get_single_mut() else {
            return;
        };
        let frame_row_len = world.matrix.len() * pixel_scaling;
        wrapper
            .pixels
            .resize_buffer(frame_row_len as u32, frame_row_len as u32)
            .unwrap();
        let frame = wrapper.pixels.frame_mut();

        assert_eq!(frame_row_len * frame_row_len * 4, frame.len());

        for i in 0..world.matrix.len() {
            for j in 0..world.matrix.len() {
                for x in 0..pixel_scaling {
                    for y in 0..pixel_scaling {
                        let color = if x <= pixel_scaling / 2 && y <= pixel_scaling / 2 {
                            Self::color_tile_content(&world.matrix[i][j])
                                .unwrap_or(Self::color_tile(&world.matrix[i][j]))
                        } else {
                            Self::color_tile(&world.matrix[i][j])
                        };
                        let pixel_coords = (i * pixel_scaling + x, j * pixel_scaling + y);
                        let pixel_index = pixel_coords.0 + frame_row_len * pixel_coords.1;
                        frame[pixel_index * 4..pixel_index * 4 + 4].copy_from_slice(&color);
                    }
                }
            }
        }
    }

    /// This method creates a window and renders the world in it.
    ///
    /// # Arguments
    /// - `world` - The world you want to render
    /// - `resolution` - The resolution of the output window
    ///
    /// # Examples
    /// ```
    /// # use robotics_lib::world::world_generator::Generator;
    /// # use midgard::{ *, params::*};
    /// # let mut world_generator = WorldGenerator::new(Default::default());
    /// let (world, _spawn_point, _weather, _max_score, _score_table) = world_generator.gen();
    /// WorldVisualizer::visualize(world, 600);
    /// ```
    /// 
    /// # Notes
    /// 
    /// The function panics if the resolution is lower than the world size.
    pub fn visualize(world: Vec<Vec<Tile>>, resolution: usize) {
        assert!(resolution >= world.len(), "WorldVisualizer::visualize must be called with resolution >= world_size ({resolution} < {})", world.len());

        let pixel_scaling = resolution / world.len();

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
        let color: u32 = match tile.content {
            // Content::Water(_) => Some()
            Content::Fish(_) => 0xff9933ff,
            Content::Tree(_) => 0x336600ff,
            Content::Rock(_) => 0x999966ff,
            Content::Fire => 0xff6600ff,
            Content::Garbage(_) => 0x333300ff,
            Content::Coin(_) => 0xcc9900ff,
            Content::Bin(_) => 0xff6600ff,
            Content::Crate(_) => 0x666633ff,
            Content::Market(_) => 0xcc00ffff,
            _ => return None,
        };
        return Some((color as u32).to_be_bytes());
    }
    fn color_tile(tile: &Tile) -> [u8; 4] {
        let color: u32 = match tile.tile_type {
            TileType::DeepWater => 0x000066ff,
            TileType::ShallowWater => 0x3333ffff,
            TileType::Grass => 0x009933ff,
            TileType::Sand => 0xffcc00ff,
            TileType::Lava => 0xcf1020ff,
            TileType::Hill => 0x669900ff,
            TileType::Mountain => 0x996633ff,
            TileType::Snow => 0xccffffff,
            TileType::Teleport(_) => 0x66ffffff,
            TileType::Street => 0x404040ff,
            _ => 0x000000ff,
        };
        return (color as u32).to_be_bytes();
    }


    pub fn visualize_realtime(generator_function: fn() -> Vec<Vec<Tile>>, resolution: usize) {
        let world = generator_function();
        assert!(resolution >= world.len(), "WorldVisualizer::visualize_realtime must be called with resolution >= world_size ({resolution} < {})", world.len());

        let pixel_scaling = resolution / world.len();

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
            .add_systems(Update, Self::regen_world)
            .insert_resource(WorldMatrixResource { matrix: world })
            .insert_resource(PixelScalingResource { pixel_scaling })
            .insert_resource(RegenWorldFuncResource { func: generator_function })
            .insert_resource(TimeCounterResource { t: 0.0, last_t: 0.0 })
            .run();
    }
    fn regen_world(mut world: ResMut<WorldMatrixResource>, gen_func: Res<RegenWorldFuncResource>, time: Res<Time>, mut time_counter_resource: ResMut<TimeCounterResource>) {
        time_counter_resource.t += time.delta_seconds_f64();
        if time_counter_resource.t > time_counter_resource.last_t + 1.0 {
            world.matrix = (gen_func.deref().func)();
            time_counter_resource.last_t = time_counter_resource.t;
        }
    }
}
