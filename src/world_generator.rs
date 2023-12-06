mod multi_octave_noise;

use bevy::{prelude::*, window::WindowResolution};
use bevy_pixels::prelude::*;
use noise::*;
use multi_octave_noise::Multi;
use robotics_lib::world::{worldgenerator::Generator, tile::Tile, tile::{TileType, Content}, environmental_conditions::{EnvironmentalConditions, WeatherType}};

pub struct WorldGenerator {
    seed: u32,
    world_size: usize,
}

#[derive(Resource)]
struct WorldMatrixResource {
    matrix: Vec<Vec<Tile>>
}
impl WorldGenerator {
    pub fn new(seed: u32, world_size: usize) -> Self {
        Self { seed, world_size }
    }

    fn generate_weather(&self) -> EnvironmentalConditions {
        EnvironmentalConditions::new(
            &[WeatherType::Sunny, WeatherType::Rainy],
            1,
            0
        ).unwrap()
    }

    fn generate_altitude(&self) -> Vec<Vec<f64>> {
        let noise_function = Multi::new(Perlin::new(self.seed), 7, 1.0/90.0);

        let noise_function = Multiply::new(Constant::new(1.5), noise_function);

        // clamp values in [-1,1] range
        // let noise_function = Clamp::new(noise_function).set_bounds(-1.0, 1.0);

        let mut elevation_map = vec![vec![0.0; self.world_size]; self.world_size];

        let mut min = 10000.;
        let mut max = -10000.;
        for x in 0..self.world_size {
            for y in 0..self.world_size {
                let noise = noise_function.get([x as f64, y as f64]);
                if noise > max { max = noise; }
                if noise < min { min = noise; }
                elevation_map[x][y] = noise;
            }
        }
        println!("altitude range: = [{min}, {max}]");

        return elevation_map;
    }

    fn generate_biomes(&self, elevation_map: Vec<Vec<f64>>) -> Vec<Vec<Tile>> {
        let deep_water_tile = Tile { tile_type: TileType::DeepWater, content: Content::None, elevation: 0};
        let shallow_water_tile = Tile { tile_type: TileType::ShallowWater, content: Content::None, elevation: 0};
        let grass_tile = Tile { tile_type: TileType::Grass, content: Content::None, elevation: 0};
        let sand_tile = Tile { tile_type: TileType::Sand, content: Content::None, elevation: 0};
        let lava_tile = Tile { tile_type: TileType::Lava, content: Content::None, elevation: 0};
        let hill_tile = Tile { tile_type: TileType::Hill, content: Content::None, elevation: 0};
        let mountain_tile = Tile { tile_type: TileType::Mountain, content: Content::None, elevation: 0};
        let snow_tile = Tile { tile_type: TileType::Snow, content: Content::None, elevation: 0};

        let mut world = vec![vec![deep_water_tile.clone(); self.world_size]; self.world_size];
        let mut temperature_map: Vec<Vec<f64>> = vec![vec![0.0; self.world_size]; self.world_size];

        let scaling = 1.0/35.0;
        let perlin = Perlin::new(self.seed + 42);
        let amplitudes = [1.0, 1.0/2.0, 1.0/4.0, 1.0/8.0, 1.0/16.0, 1.0/32.0];
        let frequencies = [1.0, 2.0, 4.0, 8.0, 16.0, 32.0].map(|x| x * scaling);
        for x in 0..self.world_size {
            for y in 0..self.world_size {
                temperature_map[x][y] =
                    (
                        amplitudes[0] * perlin.get([frequencies[0] * x as f64, frequencies[0] * y as f64])
                        + amplitudes[1] * perlin.get([frequencies[1] * x as f64, frequencies[1] * y as f64])
                        + amplitudes[2] * perlin.get([frequencies[2] * x as f64, frequencies[2] * y as f64])
                        + amplitudes[3] * perlin.get([frequencies[3] * x as f64, frequencies[3] * y as f64])
                        + amplitudes[4] * perlin.get([frequencies[4] * x as f64, frequencies[4] * y as f64])
                        + amplitudes[5] * perlin.get([frequencies[5] * x as f64, frequencies[5] * y as f64])
                    ) / amplitudes.iter().sum::<f64>()
            }
        }

        let get_biome = |x: usize, y: usize| -> Tile {
            //println!("{}", temperature_map[x][y]);

            match temperature_map[x][y] {
                //plains
                v if v > -0.5 && v <= 0.0 => grass_tile.clone(),
                //desert
                v if v >  0.0 && v < 0.5 => sand_tile.clone(),
                //forest
                _ => lava_tile.clone()
            }
        };

        for x in 0..self.world_size {
            for y in 0..self.world_size {
                let tile_type = match elevation_map[x][y] {
                    h if h < -0.65 => TileType::DeepWater,
                    h if h < -0.40 => TileType::ShallowWater,
                    h if h < -0.30 => TileType::Sand,
                    h if h <  0.35 => TileType::Grass,
                    h if h <  0.55 => TileType::Hill,
                    h if h <  0.85 => TileType::Mountain,
                    _ => TileType::Snow,
                };
                world[x][y] = Tile { tile_type, content: Content::None, elevation: 0 };
            }
        }

        return world;
    }

    fn generate_spawnpoint(&self) -> (usize, usize) {
        (0, 0)
    }

    fn color_tile(tile: &Tile) -> [u8; 4] {
        return match tile.tile_type {
            TileType::DeepWater => [0, 0, 127, 255],
            TileType::ShallowWater => [0, 0, 255, 255],
            TileType::Grass => [0, 255, 0, 255],
            TileType::Sand => [255, 255, 0, 255],
            TileType::Lava => [255, 0, 0, 255],
            TileType::Hill => [0, 127, 0, 255],
            TileType::Mountain => [153, 102, 51, 255],
            TileType::Snow => [255, 255, 255, 255],
            TileType::Sand => [192, 192, 0, 255],
            _ => [0, 0, 0, 255]
        }
    }

    fn draw_window(mut wrapper_query: Query<&mut PixelsWrapper>, world: Res<WorldMatrixResource>) {
        //Bevy pixels stuff
        let Ok(mut wrapper) = wrapper_query.get_single_mut() else { return };
        wrapper.pixels.resize_buffer(world.matrix.len() as u32, world.matrix.len() as u32).unwrap();
        let frame = wrapper.pixels.frame_mut();

        assert_eq!(world.matrix.len() * world.matrix.len() * 4, frame.len());

        for i in 0..(frame.len() / 4) {
            let color = Self::color_tile(&world.matrix[i % world.matrix.len()][i / world.matrix.len()]);
            frame[i*4..i*4+4].copy_from_slice(&color);
        }
    }

    pub fn visualize(world: Vec<Vec<Tile>>, resolution : usize) {
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
            .insert_resource(WorldMatrixResource{ matrix: world })
            .run();
    }
}

impl Generator for WorldGenerator {
    fn gen(&mut self) -> (Vec<Vec<Tile>>, (usize, usize), EnvironmentalConditions, f32) {    
        let altitude_map = self.generate_altitude();
        let world = self.generate_biomes(altitude_map);

        let weather = self.generate_weather();
        let spawnpoint = self.generate_spawnpoint();
        let score = 100.0;

        (world, spawnpoint, weather, score)
    }
}
