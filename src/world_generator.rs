mod multi_octave_noise;
mod isize_index_matrix;
mod performance_profiler;

use std::collections::{HashMap, HashSet};
use std::time::SystemTime;
use noise::*;
use multi_octave_noise::Multi;
use robotics_lib::world::{worldgenerator::Generator, tile::Tile, tile::{TileType, Content}, environmental_conditions::{EnvironmentalConditions, WeatherType}};
use fast_poisson::Poisson2D;
use rand::{Rng, SeedableRng};
use rand::rngs::StdRng;
use isize_index_matrix::*;
use crate::world_generator::performance_profiler::PerformanceProfiler;

pub struct WorldGenerator {
    seed: u32,
    world_size: usize,
}

#[derive(Eq, Hash, PartialEq, Debug, Clone)]
enum Biomes {
    Deepwater,
    ShallowWater,
    Beach,
    Desert,
    Plain,
    Forest,
    Hill,
    Mountain,
    SnowyMountain
}

impl WorldGenerator {
    pub fn new(seed: u32, world_size: usize) -> Self {
        Self { seed, world_size }
    }

    fn generate_weather(&self) -> EnvironmentalConditions {
        let weather_types = vec![WeatherType::Sunny, WeatherType::Foggy, WeatherType::Rainy, WeatherType::TrentinoSnow, WeatherType::TropicalMonsoon];
        let mut weather_forecast: Vec<WeatherType> = vec![];

        let mut rng = StdRng::seed_from_u64(self.seed as u64);

        for _i in 0..rng.gen_range(1..10) {
            weather_forecast.push(weather_types[rng.gen_range(0..weather_types.len())]);
        }

        print!("Weather forecast: ");
        for weather_type in &weather_forecast {
            print!("{weather_type:?} ");
        }
        println!();

        return EnvironmentalConditions::new(
            &weather_forecast,
            1,
            0,
        ).unwrap();
    }

    fn generate_altitude(&self, octaves: u8) -> Vec<Vec<f64>> {
        let mut noise_function = Multi::new(Perlin::new(self.seed), octaves, 1.0 / 180.0);
        noise_function.set_ampl_decay(0.5);

        let noise_function = Multiply::new(Constant::new(1.5), noise_function);

        // clamp values in [-1,1] range
        // let noise_function = Clamp::new(noise_function).set_bounds(-1.0, 1.0);

        let mut elevation_map = vec![vec![0.0; self.world_size]; self.world_size];

        for x in 0..self.world_size {
            for y in 0..self.world_size {
                elevation_map[x][y] = noise_function.get([x as f64, y as f64]);
            }
        }

        return elevation_map;
    }

    fn generate_temperature_map(&self) -> Vec<Vec<f64>> {
        let mut temperature_map = vec![vec![0.0; self.world_size]; self.world_size];

        let noise_function = Multiply::new(Constant::new(1.5), Multi::new(Perlin::new(self.seed + 42), 7, 1.0 / 100.0));

        for x in 0..self.world_size {
            for y in 0..self.world_size {
                temperature_map[x][y] = noise_function.get([x as f64, y as f64]);
            }
        }

        return temperature_map;
    }

    fn generate_biomes(&self, world: &mut Vec<Vec<Tile>>, elevation_map: &Vec<Vec<f64>>) -> HashMap<Biomes, Vec<(usize, usize)>> {
        let mut biomes_map: HashMap<Biomes, Vec<(usize, usize)>> = HashMap::new();
        let temperature_map = self.generate_temperature_map();

        for x in 0..self.world_size {
            for y in 0..self.world_size {
                match elevation_map[x][y] {
                    h if h < -0.65 => biomes_map.entry(Biomes::Deepwater).or_insert(vec![(x, y)]).push((x, y)),
                    h if h < -0.40 => biomes_map.entry(Biomes::ShallowWater).or_insert(vec![(x, y)]).push((x, y)),
                    h if h < -0.30 => biomes_map.entry(Biomes::Beach).or_insert(vec![(x, y)]).push((x, y)),
                    h if h < 0.35 => {
                        match temperature_map[x][y] {
                            t if t < -0.3 => biomes_map.entry(Biomes::Forest).or_insert(vec![(x, y)]).push((x, y)),
                            t if t > 0.2 => biomes_map.entry(Biomes::Desert).or_insert(vec![(x, y)]).push((x, y)),
                            _ => biomes_map.entry(Biomes::Plain).or_insert(vec![(x, y)]).push((x, y)),
                        }
                    }
                    h if h < 0.55 => biomes_map.entry(Biomes::Hill).or_insert(vec![(x, y)]).push((x, y)),
                    h if h < 0.85 => biomes_map.entry(Biomes::Mountain).or_insert(vec![(x, y)]).push((x, y)),
                    _ => biomes_map.entry(Biomes::SnowyMountain).or_insert(vec![(x, y)]).push((x, y)),
                };
            }
        }

        let mut fill_tiles_with = |coords: &Vec<(usize, usize)>, tile_type: TileType| {
            for (x, y) in coords {
                world[*x][*y] = Tile { tile_type, content: Content::None, elevation: 0 };
            }
        };

        for (biome, coords) in biomes_map.iter() {
            match biome {
                Biomes::Deepwater => fill_tiles_with(coords, TileType::DeepWater),
                Biomes::ShallowWater => fill_tiles_with(coords, TileType::ShallowWater),
                Biomes::Beach | Biomes::Desert => fill_tiles_with(coords, TileType::Sand),
                Biomes::Plain | Biomes::Forest => fill_tiles_with(coords, TileType::Grass),
                Biomes::Hill => fill_tiles_with(coords, TileType::Hill),
                Biomes::Mountain => fill_tiles_with(coords, TileType::Mountain),
                Biomes::SnowyMountain => fill_tiles_with(coords, TileType::Snow),
            }
        }

        return biomes_map;
    }

    fn generate_hellfire(&self, world: &mut Vec<Vec<Tile>>, biomes_map: &mut HashMap<Biomes, Vec<(usize, usize)>>) {
        if biomes_map.get(&Biomes::Desert).is_none() { return; }

        let mut lava_lakes_map = vec![vec![0.0; self.world_size]; self.world_size];
        let mut fire_zones_map = vec![vec![0.0; self.world_size]; self.world_size];

        let lava_noise_function = Multiply::new(Constant::new(1.5), Multi::new(Perlin::new(self.seed + 7), 7, 1.0 / 30.0));
        let fire_noise_function = Multiply::new(Constant::new(1.5), Multi::new(Perlin::new(self.seed + 77), 7, 1.0 / 30.0));
        for x in 0..self.world_size {
            for y in 0..self.world_size {
                lava_lakes_map[x][y] = lava_noise_function.get([x as f64, y as f64]);
                fire_zones_map[x][y] = fire_noise_function.get([x as f64, y as f64]);
            }
        }

        for (x, y) in biomes_map.get(&Biomes::Desert).unwrap().iter() {
            if world[*x][*y].tile_type != TileType::ShallowWater {
                if world[*x][*y].content != Content::Fire && lava_lakes_map[*x][*y] < -0.6 {
                    world[*x][*y].tile_type = TileType::Lava;
                }
                if world[*x][*y].tile_type != TileType::Lava && fire_zones_map[*x][*y] < -0.5 {
                    world[*x][*y].content = Content::Fire;
                }
            }
        }
    }

    fn generate_rivers(&self, world: &mut Vec<Vec<Tile>>, elevation : &Vec<Vec<f64>>, rivers_amount : f64) {
        let number_of_rivers = (world.len() as f64 * world.len() as f64 * rivers_amount / 1000000.0) as usize;
        let inertia_factor = 0.005;

        let rng_seed = {
            let mut rng_seed = [0u8; 32];
            rng_seed[0..4].copy_from_slice(&self.seed.to_le_bytes());
            rng_seed
        };

        let mut rng = StdRng::from_seed(rng_seed.into());

        let world_pos_distribution = rand::distributions::Uniform::new(0, world.len() as isize);
        let float_distribution = rand::distributions::Uniform::new(0.0, 1.0);
        let directions = [(-1, 0), (1, 0), (0, -1), (0, 1)];

        for _ in 0..number_of_rivers {
            let river_tiles = loop {

                let mut number_of_loops_looking_for_start_coord = 0; //safeguard against infinite looping while looking for a valid random start coordinates in worlds that lack it
                let start_coords = loop {
                    let coords = (rng.sample(world_pos_distribution), rng.sample(world_pos_distribution));

                    if *elevation.at(coords) > 0.35 {
                        break coords;
                    }

                    if number_of_loops_looking_for_start_coord > self.world_size * self.world_size * 2 {
                        if *elevation.at(coords) > 0.15 {
                            break coords;
                        }

                        if number_of_loops_looking_for_start_coord > self.world_size * self.world_size * 10 {
                            return;
                        }
                    }
                    number_of_loops_looking_for_start_coord += 1;
                };

                let mut river_tiles_stack = vec![start_coords];

                let mut avoid_tiles = HashSet::new();
                let mut river_tiles_set = HashSet::new();
                let mut river_inertia = (0.0, 0.0);

                loop {
                    if river_tiles_stack.is_empty() {
                        break;
                    }
                    let coords = *river_tiles_stack.last().unwrap();
                    if [0, world.len() as isize-1].contains(&coords.0) || [0, world.len() as isize-1].contains(&coords.1) {
                        match rng.sample(float_distribution) {
                            n if n < 0.6 => {
                                //discard this river
                                river_tiles_stack.clear();
                                break;
                            }
                            n if n < 0.7 => {
                                //flow off the edge of the world
                                break;
                            }
                            n if n < 0.9 => {
                                //backtrack to the previous tile
                                river_tiles_stack.pop().unwrap();
                                continue;
                            }
                            _ => {}
                        }
                    }

                    if [TileType::ShallowWater, TileType::DeepWater].contains(&world.at(coords).tile_type) {
                        break;
                    }

                    let candidates: Vec<_> =
                        directions.iter()
                            .filter_map(|(x, y)| {
                                let c = (coords.0 + x, coords.1 + y);

                                if avoid_tiles.contains(&c) {
                                    return None;
                                }
                                if !(0..world.len() as isize).contains(&c.0) || !(0..world.len() as isize).contains(&c.1) {
                                    return None;
                                }
                                //discard tiles which are adjacent to a river tile (other than the last)
                                for d in directions {
                                    let c_near = (c.0 + d.0, c.1 + d.1);
                                    if c_near == coords { continue; }

                                    if river_tiles_set.contains(&c_near) {
                                        return None;
                                    }
                                }
                                return Some(c);
                            }).collect();

                    if candidates.len() == 0 {
                        avoid_tiles.insert(coords);
                        let tile = river_tiles_stack.pop().unwrap();
                        river_tiles_set.remove(&tile);
                        continue;
                    }
                    let mut min_index = 0;

                    for i in 1..candidates.len() {
                        let c = candidates[i];
                        let m = candidates[min_index];

                        let direction_c = ((c.0 - coords.0) as f64, (c.1 - coords.1) as f64);
                        let inertia_c = (direction_c.0 * river_inertia.0 + direction_c.1 * river_inertia.1) * inertia_factor;

                        let direction_m = ((m.0 - coords.0) as f64, (m.1 - coords.1) as f64);
                        let inertia_m = (direction_m.0 * river_inertia.0 + direction_m.1 * river_inertia.1) * inertia_factor;

                        if elevation.at(c) - inertia_c < elevation.at(m) - inertia_m { min_index = i; }
                    }

                    let direction = (candidates[min_index].0 - coords.0, candidates[min_index].1 - coords.1);

                    river_inertia.0 = (river_inertia.0 * 1.0 + direction.0 as f64).clamp(-4.0, 4.0);
                    river_inertia.1 = (river_inertia.1 * 1.0 + direction.1 as f64).clamp(-4.0, 4.0);

                    avoid_tiles.insert(candidates[min_index]);
                    river_tiles_stack.push(candidates[min_index]);
                    river_tiles_set.insert(candidates[min_index]);
                }

                if river_tiles_stack.is_empty() {
                    continue;
                } else {
                    break river_tiles_stack;
                }
            };

            // for each river tile fill it and its neighbors with shallow water
            for coords in river_tiles {
                for (x,y) in directions.iter() {
                    if let Some(tile) = world.at_mut_checked((coords.0 + x, coords.1 + y)) {
                        tile.tile_type = TileType::ShallowWater;
                    }
                }
                world.at_mut(coords).tile_type = TileType::ShallowWater;
            }
        }
    }

    fn generate_teleports(&self, world: &mut Vec<Vec<Tile>>) {
        let mut teleport_coordinates = Poisson2D::new();
        let random_bias = StdRng::seed_from_u64(self.seed as u64).gen_range(0..100);
        teleport_coordinates.set_seed(self.seed as u64 + random_bias);
        teleport_coordinates.set_dimensions([self.world_size as f64, self.world_size as f64], 100.0);
        let tiles_types_to_avoid = vec![TileType::ShallowWater, TileType::DeepWater, TileType::Lava];

        for coordinate in teleport_coordinates.iter() {
            let (x, y) = (coordinate[0] as usize, coordinate[1] as usize);
            if tiles_types_to_avoid.iter().all(|t| *t != world[x][y].tile_type) && world[x][y].content != Content::Fire {
                world[x][y] = Tile { tile_type: TileType::Teleport(false), content: Content::None, elevation: 0 };
            }
        }
    }

    fn generate_content(&self, world: &mut Vec<Vec<Tile>>, biomes_map: &mut HashMap<Biomes, Vec<(usize, usize)>>, coordinates: &Vec<[f64; 2]>, allowed_biomes: &Vec<Biomes>, tiles_types_to_avoid: Option<&Vec<TileType>>, content: &Content) {
        for coordinate in coordinates.iter() {
            let (x, y) = (coordinate[0] as usize, coordinate[1] as usize);
            for biome in allowed_biomes {
                if biomes_map.get(biome).is_some() {
                    let coords = biomes_map.get(biome).unwrap();
                    if coords.contains(&(x, y)) {
                        match tiles_types_to_avoid {
                            Some(tiles) => {
                                if tiles.iter().all(|t| *t != world[x][y].tile_type) {
                                    world[x][y].content = content.clone();
                                }
                            },
                            None => world[x][y].content = content.clone()
                        }
                    }
                }
            }
        }
    }

    fn generate_contents(&self, world: &mut Vec<Vec<Tile>>, biomes_map: &mut HashMap<Biomes, Vec<(usize, usize)>>) {
        //Water
        if biomes_map.get(&Biomes::Deepwater).is_some() {
            for (x, y) in biomes_map.get(&Biomes::Deepwater).unwrap() {
                world[*x][*y].content = Content::Water(2);
            }
        }
        if biomes_map.get(&Biomes::Deepwater).is_some() {
            for (x, y) in biomes_map.get(&Biomes::ShallowWater).unwrap() {
                world[*x][*y].content = Content::Water(1);
            }
        }

        let mut poisson = Poisson2D::new();
        poisson.set_seed(self.seed as u64);
        poisson.set_samples(15);

        let tiles_types_to_avoid = vec![TileType::ShallowWater, TileType::DeepWater, TileType::Lava];
        let allowed_biomes = vec![Biomes::Beach, Biomes::Desert, Biomes::Plain, Biomes::Forest, Biomes::Hill];
        let configurations = vec![
            (3, vec![Biomes::Forest], Some(&tiles_types_to_avoid), Content::Tree(1)),
            (4, vec![Biomes::Hill], Some(&tiles_types_to_avoid), Content::Tree(1)),
            (5, vec![Biomes::Mountain], Some(&tiles_types_to_avoid), Content::Tree(1)),
            (5, vec![Biomes::Plain], Some(&tiles_types_to_avoid), Content::Rock(1)),
            (4, vec![Biomes::Hill], Some(&tiles_types_to_avoid), Content::Rock(1)),
            (3, vec![Biomes::Mountain], Some(&tiles_types_to_avoid), Content::Rock(1)),
            (3, vec![Biomes::SnowyMountain], Some(&tiles_types_to_avoid), Content::Rock(1)),
            (5, vec![Biomes::ShallowWater], None, Content::Fish(1)),
            (4, vec![Biomes::Deepwater], None, Content::Fish(1)),
            (10, allowed_biomes.clone(), Some(&tiles_types_to_avoid), Content::Garbage(1)),
            (20, allowed_biomes.clone(), Some(&tiles_types_to_avoid), Content::Coin(1)),
            (20, allowed_biomes.clone(), Some(&tiles_types_to_avoid), Content::Bin(0..5)),
            (40, allowed_biomes.clone(), Some(&tiles_types_to_avoid), Content::Crate(0..5)),
            (50, allowed_biomes.clone(), Some(&tiles_types_to_avoid), Content::Market(1)),
        ];

        let mut coords: HashMap<u64, Vec<[f64; 2]>> = HashMap::new();
        for (radius, test_allowed_biomes, tiles_types_to_avoid, content) in configurations {
            if coords.get(&radius).is_none() {
                poisson.set_dimensions([self.world_size as f64, self.world_size as f64], radius as f64);
                coords.insert(radius, poisson.generate());
            }

            self.generate_content(world, biomes_map, coords.get(&radius).unwrap(), &test_allowed_biomes, tiles_types_to_avoid, &content);
        }
    }

    fn generate_spawnpoint(&self, biomes_map: &mut HashMap<Biomes, Vec<(usize, usize)>>) -> (usize, usize) {
        let mut spawnpint_coords = Poisson2D::new();
        let random_bias = StdRng::seed_from_u64(self.seed as u64).gen_range(0..100);
        spawnpint_coords.set_seed(self.seed as u64 + random_bias);
        spawnpint_coords.set_dimensions([self.world_size as f64, self.world_size as f64], 50.0);

        let allowed_biomes = vec![Biomes::Plain, Biomes::Beach, Biomes::Forest];
        let mut spawnpoint = (0, 0);

        'outer: for coordinate in spawnpint_coords.iter() {
            let (x, y) = (coordinate[0] as usize, coordinate[1] as usize);
            for allowed_biome in allowed_biomes.iter() {
                if biomes_map.get(allowed_biome).is_some() {
                    let biome_coords = biomes_map.get(allowed_biome).unwrap();
                    if biome_coords.contains(&(x, y)) {
                        spawnpoint = (x, y);
                        break 'outer;
                    }
                }
            }
        }

        println!("Spawnpoint {spawnpoint:?}");
        return spawnpoint;
    }
}

impl Generator for WorldGenerator {
    fn gen(&mut self) -> (Vec<Vec<Tile>>, (usize, usize), EnvironmentalConditions, f32, Option<HashMap<Content, f32>>) {
        let mut profiler = PerformanceProfiler::new(SystemTime::now());

        println!("World seed: {}", self.seed);

        let mut world = vec![vec![Tile { tile_type: TileType::DeepWater, content: Content::None, elevation: 0 }; self.world_size]; self.world_size];
        profiler.print_elapsed_time_in_ms("water world generation time");

        let altitude_map = self.generate_altitude(5);
        profiler.print_elapsed_time_in_ms("altitude generation time");

        let mut biomes_map = self.generate_biomes(&mut world, &altitude_map);
        profiler.print_elapsed_time_in_ms("biomes generation time");

        self.generate_rivers(&mut world, &self.generate_altitude(7), 30.0);
        profiler.print_elapsed_time_in_ms("rivers generation time");

        self.generate_hellfire(&mut world, &mut biomes_map);
        profiler.print_elapsed_time_in_ms("hellfire generation time");

        self.generate_teleports(&mut world);
        profiler.print_elapsed_time_in_ms("teleports generation time");

        self.generate_contents(&mut world, &mut biomes_map);
        profiler.print_elapsed_time_in_ms("Total contents generation time");

        let weather_forecast = self.generate_weather();
        profiler.print_elapsed_time_in_ms("weather generation time");

        let spawnpoint = self.generate_spawnpoint(&mut biomes_map);
        profiler.print_elapsed_time_in_ms("spawnpoint generation time");

        let score = 1000.0;

        profiler.print_total_elapsed_time_in_ms("Total generation time");

        //TO-DO score table
        (world, spawnpoint, weather_forecast, score, None)
    }
}
