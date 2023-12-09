mod multi_octave_noise;
mod isize_index_matrix;
mod performance_telemetry;
mod vector_math;

use std::collections::{HashMap, HashSet};
use std::time::SystemTime;
use noise::*;
use multi_octave_noise::Multi;
use robotics_lib::world::{worldgenerator::Generator, tile::Tile, tile::{TileType, Content}, environmental_conditions::{EnvironmentalConditions, WeatherType}};
use fast_poisson::Poisson2D;
use rand::{distributions, Rng, SeedableRng};
use rand::prelude::*;
use isize_index_matrix::*;
use vector_math::*;
use performance_telemetry::PerformanceTelemetry;

pub struct WorldGenerator {
    seed: u32,
    world_size: usize,
}

#[derive(Clone, PartialEq)]
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
        EnvironmentalConditions::new(
            &[WeatherType::Sunny, WeatherType::Rainy],
            1,
            0
        ).unwrap()
    }

    fn generate_elevation(&self, octaves: u8) -> Vec<Vec<f64>> {
        let mut noise_function = Multi::new(Perlin::new(self.seed), octaves, 1.0 / 180.0);
        noise_function.set_ampl_decay(0.5);

        let noise_function = Multiply::new(Constant::new(1.5), noise_function);

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

    fn generate_biomes(&self, world: &mut Vec<Vec<Tile>>, elevation_map: &Vec<Vec<f64>>) -> Vec<Vec<Biomes>> {
        let mut biomes_map: Vec<Vec<Biomes>> = vec![vec![Biomes::Deepwater; self.world_size]; self.world_size];
        let temperature_map = self.generate_temperature_map();

        let get_biome = |temperature: f64| -> Biomes {
            return match temperature {
                t if t < -0.3 => Biomes::Forest,
                t if t > 0.2 => Biomes::Desert,
                _ => Biomes::Plain
            }
        };

        for x in 0..self.world_size {
            for y in 0..self.world_size {
                biomes_map[x][y] = match elevation_map[x][y] {
                    h if h < -0.65 => Biomes::Deepwater,
                    h if h < -0.40 => Biomes::ShallowWater,
                    h if h < -0.30 => Biomes::Beach,
                    h if h < 0.35 => get_biome(temperature_map[x][y]),
                    h if h < 0.55 => Biomes::Hill,
                    h if h < 0.75 => Biomes::Mountain,
                    _ => Biomes::SnowyMountain,
                };
            }
        }

        for x in 0..self.world_size {
            for y in 0..self.world_size {
                let tile_type = match biomes_map[x][y] {
                    Biomes::Deepwater => TileType::DeepWater,
                    Biomes::ShallowWater => TileType::ShallowWater,
                    Biomes::Beach | Biomes::Desert => TileType::Sand,
                    Biomes::Plain | Biomes::Forest => TileType::Grass,
                    Biomes::Hill => TileType::Hill,
                    Biomes::Mountain => TileType::Mountain,
                    Biomes::SnowyMountain => TileType::Snow,
                };
                world[x][y] = Tile { tile_type, content: Content::None, elevation: 0 };
            }
        }

        return biomes_map;
    }

    fn generate_hellfire(&self, world: &mut Vec<Vec<Tile>>, biomes_map: &Vec<Vec<Biomes>>) {
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

        for x in 0..self.world_size {
            for y in 0..self.world_size {
                if biomes_map[x][y] == Biomes::Desert && world[x][y].tile_type != TileType::ShallowWater {
                    if world[x][y].content != Content::Fire && lava_lakes_map[x][y] < -0.6 {
                        world[x][y].tile_type = TileType::Lava;
                    }
                    if world[x][y].tile_type != TileType::Lava && fire_zones_map[x][y] < -0.5 {
                        world[x][y].content = Content::Fire;
                    }
                }
            }
        }
    }

    fn generate_content(&self, world: &mut Vec<Vec<Tile>>, biomes_map: &Vec<Vec<Biomes>>, coordinates: &mut Poisson2D, radius: f64, biome: Biomes, to_avoid: Option<&Vec<TileType>>, content: Content) {
        coordinates.set_dimensions([self.world_size as f64, self.world_size as f64], radius);

        for coordinate in coordinates.iter() {
            let (x, y) = (coordinate[0] as usize, coordinate[1] as usize);
            match to_avoid {
                Some(tiles) => {
                    if biomes_map[x][y] == biome && tiles.iter().all(|b| *b != world[x][y].tile_type) {
                        world[x][y].content = content.clone();
                    }
                },
                None => {
                    if biomes_map[x][y] == biome {
                        world[x][y].content = content.clone();
                    }
                }
            }
        }
    }

#[allow(dead_code)]
    fn generate_trees(&self, world: &mut Vec<Vec<Tile>>, biomes_map: &Vec<Vec<Biomes>>) {
        let mut trees_coordinates = Poisson2D::new();
        trees_coordinates.set_seed(self.seed as u64);
        let to_avoid = vec![TileType::ShallowWater, TileType::DeepWater];

        self.generate_content(world, biomes_map, &mut trees_coordinates, 2.5, Biomes::Forest, Some(&to_avoid), Content::Tree(1));
        self.generate_content(world, biomes_map, &mut trees_coordinates, 3.5, Biomes::Hill, Some(&to_avoid), Content::Tree(1));
        self.generate_content(world, biomes_map, &mut trees_coordinates, 4.5, Biomes::Mountain, Some(&to_avoid), Content::Tree(1));
    }

#[allow(dead_code)]
    fn generate_rocks(&self, world: &mut Vec<Vec<Tile>>, biomes_map: &Vec<Vec<Biomes>>) {
        let mut rocks_coordinates = Poisson2D::new();
        rocks_coordinates.set_seed(self.seed as u64);
        let to_avoid = vec![TileType::ShallowWater, TileType::DeepWater];

        self.generate_content(world, biomes_map, &mut rocks_coordinates, 5.0, Biomes::Plain, Some(&to_avoid), Content::Rock(1));
        self.generate_content(world, biomes_map, &mut rocks_coordinates, 4.0, Biomes::Hill, Some(&to_avoid), Content::Rock(1));
        self.generate_content(world, biomes_map, &mut rocks_coordinates, 3.0, Biomes::Mountain, Some(&to_avoid), Content::Rock(1));
        self.generate_content(world, biomes_map, &mut rocks_coordinates, 2.0, Biomes::SnowyMountain, Some(&to_avoid), Content::Rock(1));
    }
#[allow(dead_code)]
    fn generate_fishes(&self, world: &mut Vec<Vec<Tile>>, biomes_map: &Vec<Vec<Biomes>>, ) {
        let mut fishes_coordinates = Poisson2D::new();
        fishes_coordinates.set_seed(self.seed as u64);

        self.generate_content(world, biomes_map, &mut fishes_coordinates, 4.0, Biomes::ShallowWater, None, Content::Fish(1));
        self.generate_content(world, biomes_map, &mut fishes_coordinates, 3.0, Biomes::Deepwater, None, Content::Fish(1));
    }

    fn generate_rivers(&self, world: &mut Vec<Vec<Tile>>, elevation : &Vec<Vec<f64>>, amount : f64, inertia_factor : f64, inertia_decay : f64, max_inertia : f64) {
        let number_of_rivers = (world.len() as f64 * world.len() as f64 * amount / 1000000.0) as usize;
        let random_jitter = 0.12;

        let mut rng = StdRng::seed_from_u64(self.seed as u64 + 1);//seed incremented by one because otherwise we would end up trying to place roads on all the same places where we placed rivers

        let world_pos_distribution = distributions::Uniform::new(1, world.len() as isize-1);
        let float_distribution = distributions::Uniform::new(0.0, 1.0);
        let directions = [(-1, 0), (1, 0), (0, -1), (0, 1)];

        let poisson = Poisson2D::new()
            .with_seed((self.seed+1) as u64)
            .with_dimensions([world.len() as f64, world.len() as f64], world.len() as f64 / (number_of_rivers as f64));
        let mut poisson_coords_iterator = poisson.into_iter();

        for _ in 0..number_of_rivers {
            let river_tiles = loop {
                let mut number_of_loops_looking_for_start_coord = 0; //safeguard against infinite looping while looking for a valid random start coordinates in worlds that lack it
                let start_coords = loop {
                    let coords = poisson_coords_iterator.next()
                        .map(|p| (p[0] as isize, p[1] as isize))
                        .unwrap_or_else(|| (rng.sample(world_pos_distribution), rng.sample(world_pos_distribution)));

                    if *elevation.at(coords) > 0.35 {
                        let mut should_discard = false;
                        'outer: for x in -3..3 {
                            for y in -3..3 {
                                if let Some(Tile {tile_type: TileType::ShallowWater,..}) = world.at_checked((coords.0+x, coords.1+y)) {
                                    should_discard = true;
                                    break 'outer;
                                }
                            }
                        }

                        if !should_discard {
                            break coords;
                        }
                    }

                    if number_of_loops_looking_for_start_coord > self.world_size * self.world_size * 10 {
                        return;
                    }
                    number_of_loops_looking_for_start_coord += 1;
                };

                let mut river_tiles_stack = vec![start_coords];
                let mut avoid_tiles = HashSet::from([start_coords]);
                let mut river_tiles_set = HashSet::from([start_coords]);
                let mut inertia = (0.0, 0.0);
                let mut movement_debt = (0.0, 0.0);

                loop {
                    if river_tiles_stack.is_empty() {
                        break;
                    }
                    let coords = *river_tiles_stack.last().unwrap();

                    // backtrack to prev tile if the current is close to lava or close to a tile in avoid tiles
                    {
                        let mut should_backtrack = false;
                        'outer:
                        for x in -3..=3 {
                            for y in -3..=3 {
                                let c = (coords.0 + x, coords.1 + y);
                                let last_tiles_size = river_tiles_stack.len().min(7);
                                let last_tiles = &river_tiles_stack[river_tiles_stack.len() - last_tiles_size..];
                                if river_tiles_set.contains(&c) && !last_tiles.contains(&c) {
                                    should_backtrack = true;
                                    break 'outer;
                                }
                                if let Some(Tile { tile_type: TileType::Lava, .. }) = world.at_checked(c) {
                                    should_backtrack = true;
                                    break 'outer;
                                }
                            }
                        }
                        if should_backtrack {
                            let tile = river_tiles_stack.pop().unwrap();//backtrack to the previous tile
                            river_tiles_set.remove(&tile);
                            continue;
                        }
                    }

                    if [0, world.len() as isize-1].contains(&coords.0) || [0, world.len() as isize-1].contains(&coords.1) {
                        match rng.sample(float_distribution) {
                            n if n < 0.5 => {
                                break;//flow off the edge of the world
                            }
                            _ => {
                                let tile = river_tiles_stack.pop().unwrap();//backtrack to the previous tile
                                river_tiles_set.remove(&tile);
                                continue;
                            }
                        }
                    }

                    if [TileType::ShallowWater, TileType::DeepWater].contains(&world.at(coords).tile_type) {
                        break; // flow into another river or into a lake
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
                                return Some(c);
                            }).collect();

                    if candidates.len() == 0 {
                        avoid_tiles.insert(coords);
                        let tile = river_tiles_stack.pop().unwrap();
                        river_tiles_set.remove(&tile);
                        continue;
                    }


                    let gradient_of_elevation = get_gradient(elevation, coords).unwrap();
                    let direction =
                        vec_normalize(
                            vec_sum(
                                vec_sum(
                                    vec_mul_by_scalar(gradient_of_elevation, -1.0),
                                    vec_mul_by_scalar(inertia, inertia_factor)
                                ),
                                vec_mul_by_scalar(((rng.sample(float_distribution) - 0.5) * 2.0, (rng.sample(float_distribution) - 0.5) * 2.0), random_jitter)
                            )
                        );
                    inertia =
                        vec_clamp(
                            vec_sum(
                                vec_mul_by_scalar(inertia, inertia_decay),
                                direction
                            ),
                            max_inertia
                        );

                    movement_debt = vec_clamp(movement_debt, 1.0);
                    movement_debt.0 += direction.0;
                    movement_debt.1 += direction.1;
                    let move_along_x_axis = movement_debt.0.abs() > movement_debt.1.abs();

                    let x_candidate_dir = if direction.0 > 0.0 { 1 } else { -1 };
                    let x_candidate = (coords.0 + x_candidate_dir, coords.1);
                    let minus_x_candidate = (coords.0 - x_candidate_dir, coords.1);
                    let y_candidate_dir = if direction.1 > 0.0 { 1 } else { -1 };
                    let y_candidate = (coords.0, coords.1 + y_candidate_dir);
                    let minus_y_candidate = (coords.0, coords.1 - y_candidate_dir);

                    let target_dir = if move_along_x_axis {
                        if candidates.contains(&x_candidate) { (x_candidate_dir, 0) }
                        else if candidates.contains(&y_candidate) { (0, y_candidate_dir) }
                        else if candidates.contains(&minus_y_candidate) { (0, -y_candidate_dir) }
                        else { (-x_candidate_dir, 0) }
                    }
                    else {
                        if candidates.contains(&y_candidate) { (0, y_candidate_dir) }
                        else if candidates.contains(&x_candidate) { (x_candidate_dir, 0) }
                        else if candidates.contains(&minus_x_candidate) { (-x_candidate_dir, 0) }
                        else { (0, -y_candidate_dir) }
                    };

                    movement_debt.0 -= target_dir.0 as f64;
                    movement_debt.1 -= target_dir.1 as f64;

                    let target_tile = (coords.0 + target_dir.0, coords.1 + target_dir.1);

                    avoid_tiles.insert(target_tile);
                    river_tiles_stack.push(target_tile);
                    river_tiles_set.insert(target_tile);
                }

                if river_tiles_stack.is_empty() {
                    continue;
                } else {
                    break river_tiles_stack;
                }
            };

            // for each river tile fill it and its neighbors with shallow water
            {
                for c in river_tiles {
                    for (x, y) in directions.iter() {
                        if let Some(tile) = world.at_mut_checked((c.0 + x, c.1 + y)) {
                            tile.tile_type = TileType::ShallowWater;
                        }
                    }
                    world.at_mut(c).tile_type = TileType::ShallowWater;
                }
            }
        }
    }
    fn generate_street(&self, world: &mut Vec<Vec<Tile>>, elevation: &Vec<Vec<f64>>,
                       inertia_factor: f64, inertia_decay: f64, max_inertia: f64, inertia: &mut (f64, f64),
                       poi1: (isize, isize), poi2: (isize, isize))
    {
        let mut street_tiles_stack = vec![poi1];
        let mut street_tiles_set = HashSet::from([poi1]);
        let mut avoid_tiles = HashSet::from([poi1]);
        let mut movement_debt = (0.0, 0.0);

        loop {
            if street_tiles_stack.is_empty() {
                break;
            }

            let coords = *street_tiles_stack.last().unwrap();


            let dist_to_poi2 = ((poi2.0 - coords.0) as f64, (poi2.1 - coords.1) as f64);
            if vec_module(dist_to_poi2) < 1.2 {
                break;
            }

            let gradient_of_elevation = get_gradient(elevation, coords).unwrap_or((0.0, 0.0));

            let perpendiculars_to_gradient_of_elevation = [
                vec_normalize((-gradient_of_elevation.1, gradient_of_elevation.0)),
                vec_normalize((gradient_of_elevation.1, -gradient_of_elevation.0)),
            ];

            let desired_direction = vec_normalize(dist_to_poi2);

            let perpendicular_to_gradient_of_elevation_concordant_to_desired_direction = {
                let desired_direction_dot_perpendiculars_to_gradient_of_elevation =
                    perpendiculars_to_gradient_of_elevation.map(|p| vec_dot(p, desired_direction));
                let mut max = 0;
                for i in 1..desired_direction_dot_perpendiculars_to_gradient_of_elevation.len() {
                    if desired_direction_dot_perpendiculars_to_gradient_of_elevation[i] > desired_direction_dot_perpendiculars_to_gradient_of_elevation[max] {
                        max = i;
                    }
                }
                perpendiculars_to_gradient_of_elevation[max]
            };

            let proximity_threshold = (((poi1.0 - poi2.0).pow(2) + (poi1.1 - poi2.1).pow(2)) as f64).sqrt() / 2.0;
            let proximity_multiplier = f64::max(1.0, proximity_threshold / vec_module(dist_to_poi2));

            let direction = vec_sum(
                vec_sum(
                    vec_mul_by_scalar(desired_direction, proximity_multiplier),
                    vec_mul_by_scalar(*inertia, inertia_factor)
                ),
                vec_mul_by_scalar(perpendicular_to_gradient_of_elevation_concordant_to_desired_direction, 0.5)
            );
            let direction = vec_normalize(direction);

            movement_debt = vec_clamp(movement_debt, 5.0);
            movement_debt.0 += direction.0;
            movement_debt.1 += direction.1;
            let move_along_x_axis = movement_debt.0.abs() > movement_debt.1.abs();


            let x_candidate_dir = if direction.0 > 0.0 { 1 } else { -1 };
            let y_candidate_dir = if direction.1 > 0.0 { 1 } else { -1 };

            let target_dir = if move_along_x_axis {
                (x_candidate_dir, 0)
            } else {
                (0, y_candidate_dir)
            };

            movement_debt.0 -= target_dir.0 as f64;
            movement_debt.1 -= target_dir.1 as f64;

            *inertia = (inertia.0 + direction.0, inertia.1 + direction.1);
            *inertia = vec_mul_by_scalar(*inertia, inertia_decay);
            *inertia = vec_clamp(*inertia, max_inertia);

            let target_tile = (coords.0 + target_dir.0, coords.1 + target_dir.1);

            avoid_tiles.insert(target_tile);
            street_tiles_stack.push(target_tile);
            street_tiles_set.insert(target_tile);
        }

        for pos in street_tiles_stack {
            if let Some(tile) = world.at_mut_checked(pos) {
                if ![TileType::Lava, TileType::DeepWater, TileType::ShallowWater].contains(&tile.tile_type) {
                    tile.tile_type = TileType::Street;
                }
            }
        }
    }

    fn generate_streets(&self, world: &mut Vec<Vec<Tile>>, elevation : &Vec<Vec<f64>>, amount : f64, inertia_factor : f64, inertia_decay : f64, max_inertia : f64) {
        //this function generates points of interest in the map and tries to create roads connecting them
        let poi_distance= 50.0 / amount;
        let mut rng = StdRng::seed_from_u64(self.seed as u64 + 1);//seed incremented by one because otherwise we would end up trying to place roads on all the same places where we placed rivers

        let poisson = Poisson2D::new()
            .with_seed((self.seed+1) as u64)
            .with_dimensions([world.len() as f64 * 1.5, world.len() as f64 * 1.5], poi_distance);
        let poi_vec : Vec<_> =
            poisson.into_iter()
                .map(|p| (p[0] - world.len() as f64 * 0.25, p[1] - world.len() as f64 * 0.25))
                .map(|p| (p.0 as isize, p.1 as isize))
                .collect();

        let random_poi_distribution = distributions::Uniform::new(0, poi_vec.len());

        let mut street_map = HashSet::new();
        let number_of_streets = ((world.len() * world.len()) as f64/ (poi_distance * poi_distance)) as usize * 3 / 2;
        println!("number_of_streets {number_of_streets}");
        let mut number_of_empty_streets_produced = 0; // stop after number_of_streets empty streets have been produced
        while street_map.len() < number_of_streets * 2 {
            // generate a street
            let start_poi = rng.sample(random_poi_distribution);
            let number_of_pois_in_street = world.len() / poi_distance as usize * 2;
            let mut street = vec![start_poi];
            let mut avoid_set = HashSet::from([start_poi]);

            for _ in 0..number_of_pois_in_street {
                let last_poi_index = *street.last().unwrap();
                let last_poi = poi_vec[last_poi_index];

                let mut loop_iterations = 0;
                let next_poi_index = loop {
                    let next_index = rng.sample(random_poi_distribution);
                    if avoid_set.contains(&next_index) || next_index == last_poi_index {
                        continue;
                    }
                    if street_map.contains(&(last_poi_index, next_index)) {
                        continue;
                    }

                    let next = poi_vec[next_index];

                    let dist_vec = ((last_poi.0 - next.0) as f64, (last_poi.1 - next.1) as f64);
                    let dist = vec_module(dist_vec);

                    if dist < poi_distance * 1.7 {
                        break Some(next_index);
                    }

                    if loop_iterations as f64 > amount * amount * 10.0 {
                        //no valid next_poi could be found
                        break None;
                    }
                    loop_iterations += 1;
                };

                if let Some(next_poi_index) = next_poi_index {
                    street_map.insert((last_poi_index, next_poi_index));
                    street_map.insert((next_poi_index, last_poi_index));
                    avoid_set.insert(next_poi_index);
                    street.push(next_poi_index);
                } else {
                    break;
                }
            }

            let mut inertia = (0.0, 0.0);
            for i in 0..street.len()-1 {
                let poi1 = poi_vec[street[i]];
                let poi2 = poi_vec[street[i+1]];
                self.generate_street(world, elevation, inertia_factor, inertia_decay, max_inertia, &mut inertia, poi1, poi2);
                inertia = vec_mul_by_scalar(inertia, 5.0);
            }
            if street.len() == 1 {
                number_of_empty_streets_produced += 1;
            }
            if number_of_empty_streets_produced > number_of_streets {
                break;
            }
        }
    }
    fn generate_spawnpoint(&self) -> (usize, usize) {
        (0, 0)
    }
}

impl Generator for WorldGenerator {
    fn gen(&mut self) -> (Vec<Vec<Tile>>, (usize, usize), EnvironmentalConditions, f32, Option<HashMap<Content, f32>>) {
        let mut telemetry = PerformanceTelemetry::new(SystemTime::now());

        println!("World seed: {}", self.seed);

        let mut world = vec![vec![Tile { tile_type: TileType::DeepWater, content: Content::None, elevation: 0 }; self.world_size]; self.world_size];
        telemetry.print_elapsed_time_in_ms("water world generation time");

        let elevation_map = self.generate_elevation(6);
        telemetry.print_elapsed_time_in_ms("elevation generation time");

        let biomes_map = self.generate_biomes(&mut world, &elevation_map);
        telemetry.print_elapsed_time_in_ms("biomes generation time");

        self.generate_hellfire(&mut world, &biomes_map);
        telemetry.print_elapsed_time_in_ms("hellfire generation time");

        self.generate_streets(&mut world, &elevation_map, 0.5, 2.0, 0.8, 2.5);
        telemetry.print_elapsed_time_in_ms("streets generation time");

        self.generate_rivers(&mut world, &elevation_map, 80.0, 0.05, 0.9, 4.0);
        telemetry.print_elapsed_time_in_ms("rivers generation time");

        // self.generate_trees(&mut world, &biomes_map);
        telemetry.print_elapsed_time_in_ms("trees generation time");

        // self.generate_rocks(&mut world, &biomes_map);
        telemetry.print_elapsed_time_in_ms("rocks generation time");

        // self.generate_fishes(&mut world, &biomes_map);
        telemetry.print_elapsed_time_in_ms("fishes generation time");

        let weather = self.generate_weather();
        telemetry.print_elapsed_time_in_ms("weather generation time");
        let spawnpoint = self.generate_spawnpoint();
        telemetry.print_elapsed_time_in_ms("spawnpoint generation time");
        let score = 100.0;

        telemetry.print_total_elapsed_time_in_ms("Total generation time");

        (world, spawnpoint, weather, score, None)
    }
}
