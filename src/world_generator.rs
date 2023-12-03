use bevy::{prelude::*, window::WindowResolution};
use bevy_pixels::prelude::*;
use robotics_lib::world::{worldgenerator::Generator, tile::Tile, tile::{TileType, Content}, environmental_conditions::{EnvironmentalConditions, WeatherType}};
use noise::{NoiseFn, Perlin};

pub struct WorldGenerator {
    seed: u32,
    world_size: usize,
    world: Vec<Vec<Tile>>,
}

#[derive(Resource)]
struct WorldMatrixResource {
    matrix: Vec<Vec<Tile>>
}

//forse non dovrei inizializzarlo. perdita di tempo e poco significativo.
impl WorldGenerator {
    pub fn new(seed: u32, world_size: usize) -> Self {
        let default_tile = Tile { 
            tile_type: TileType::DeepWater, 
            content: Content::None, 
            elevation: 0 
        };
        let mut world = vec![vec![default_tile; world_size]; world_size];

        Self { seed, world_size, world }
    }

    fn generate_weather() -> EnvironmentalConditions {
        EnvironmentalConditions::new(
            &[WeatherType::Sunny, WeatherType::Rainy],
            1,
            0
        ).unwrap()
    }

    /**
     * Se usiamo la elevation serve trovare una soluzione per convertire i valori
     * da f64 -1,+1 a usize. si può fare un mapping manuale enorme ma non so quanto sia una buona idea.
     * possiamo fare arrotondamento ma è praticamente inutile con valori così piccoli.
     * il casting non ci aiuta.
     * 
     * altrimenti ce ne sbattiamo e teniamo la elevation come vogliamo noi per poi in caso mapparla
     * ad elevation delle tiles come si era detto. volevo evitare di utilizzare altre strutture per nulla ma sembra l'unica.
     */
    fn generate_altitude(&self) -> Vec<Vec<f64>> {
        let perlin = Perlin::new(1);
        let mut elevation_map = vec![vec![0.0; self.world_size]; self.world_size];
        
        for x in 0..self.world_size {
            for y in 0..self.world_size {
                elevation_map[x][y] = perlin.get([x as f64, y as f64]);
            }
        }

        return elevation_map;
    }

    fn generate_biomes(&mut self, elevation_map: Vec<Vec<f64>>) {
        let deep_water_tile = Tile { tile_type: TileType::DeepWater, content: Content::None, elevation: 0};
        let shallow_water_tile = Tile { tile_type: TileType::ShallowWater, content: Content::None, elevation: 0};
        let grass_tile = Tile { tile_type: TileType::Grass, content: Content::None, elevation: 0};
        let hill_tile = Tile { tile_type: TileType::Hill, content: Content::None, elevation: 0};
        let mountain_tile = Tile { tile_type: TileType::Mountain, content: Content::None, elevation: 0};
        let snow_tile = Tile { tile_type: TileType::Snow, content: Content::None, elevation: 0};
        
        for x in 0..self.world_size {
            for y in 0..self.world_size {
                self.world[x][y] = match elevation_map[x][y] {
                    h if h < -0.75 => deep_water_tile.clone(),
                    h if h < -0.50 => shallow_water_tile.clone(),
                    h if h < -0.25 => grass_tile.clone(),
                    h if h >  0.25 => hill_tile.clone(),
                    h if h >  0.50 => mountain_tile.clone(),
                    h if h >  0.75 => snow_tile.clone(),
                    _ => deep_water_tile.clone(),
                }
            }
        }
    }

    fn generate_spawnpoint(&self) -> (usize, usize) {
        (0, 0)
    }

    fn color_tile(tile: &Tile) -> [u8; 4] {
        return match tile.tile_type {
            TileType::DeepWater => [0, 0, 255, 127],
            TileType::ShallowWater => [0, 0, 255, 255],
            TileType::Grass => [0, 255, 0, 255],
            TileType::Hill => [0, 127, 0, 255],
            TileType::Mountain => [153, 102, 51, 255],
            TileType::Snow => [255, 255, 255, 255],
            _ => [0, 0, 255, 255]
        }
    }

    fn draw_window(mut wrapper_query: Query<&mut PixelsWrapper>, world: Res<WorldMatrixResource>) {
        //Bevy pixels stuff
        let Ok(mut wrapper) = wrapper_query.get_single_mut() else { return };
        let frame = wrapper.pixels.frame_mut();

        //non ho capito al massimo come funzioni la cosa dei pixel
        for i in 0..(frame.len() / 4) {
            frame[i..i+4].copy_from_slice(&Self::color_tile(&world.matrix[i % 800][i / 800]));
        }   
    }

    pub fn visualize(&self) {
        let window_plugin = WindowPlugin {
            primary_window: Some(Window {
                title: "MIDGARD".into(),
                resolution: WindowResolution::new(self.world_size as f32, self.world_size as f32),
                resizable: false,
                ..default()
            }),
            ..default()
        };
        
        App::new()
            .add_plugins((DefaultPlugins.set(window_plugin), PixelsPlugin::default()))
            .add_systems(Update, bevy::window::close_on_esc)
            .add_systems(Draw, Self::draw_window)
            .insert_resource(WorldMatrixResource{ matrix: self.world.clone() }) //si qua non va bene clonare bisognerebbe passare reference
            .run();
    }
}

impl Generator for WorldGenerator {
    fn gen(&mut self) -> (Vec<Vec<Tile>>, (usize, usize), EnvironmentalConditions, f32) {    
        let altitude_map = self.generate_altitude();
        self.generate_biomes(altitude_map);

        let weather = WorldGenerator::generate_weather();
        let spawnpoint = self.generate_spawnpoint();
        let score = 100.0;

        (self.world.to_owned(), spawnpoint, weather, score)
    }
}
