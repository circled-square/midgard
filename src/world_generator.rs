mod multi_octave_noise;
mod isize_index_matrix;
mod performance_profiler;
mod vector_math;

use std::collections::{HashMap, HashSet};
use std::time::SystemTime;
use noise::*;
use multi_octave_noise::Multi;
use robotics_lib::world::{worldgenerator::Generator, tile::Tile, tile::{TileType, Content}, environmental_conditions::{EnvironmentalConditions, WeatherType}};
use fast_poisson::Poisson2D;
use rand::rngs::StdRng;
use rand::{distributions, Rng, SeedableRng};
use isize_index_matrix::*;
use vector_math::*;
use performance_profiler::PerformanceProfiler;

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
                if world[*x][*y].tile_type != TileType::Street && world[*x][*y].tile_type != TileType::Lava && fire_zones_map[*x][*y] < -0.5 {
                    world[*x][*y].content = Content::Fire;
                }
            }
        }
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
    fn generate_teleports(&self, world: &mut Vec<Vec<Tile>>) {
        let mut teleport_coordinates = Poisson2D::new();
        let random_bias = StdRng::seed_from_u64(self.seed as u64).gen_range(0..100);
        teleport_coordinates.set_seed(self.seed as u64 + random_bias);
        teleport_coordinates.set_dimensions([self.world_size as f64, self.world_size as f64], 100.0);
        let tile_types_to_avoid = vec![TileType::ShallowWater, TileType::DeepWater, TileType::Lava];

        for coordinate in teleport_coordinates.iter() {
            let (x, y) = (coordinate[0] as usize, coordinate[1] as usize);
            if tile_types_to_avoid.iter().all(|t| *t != world[x][y].tile_type) && world[x][y].content != Content::Fire {
                world[x][y] = Tile { tile_type: TileType::Teleport(false), content: Content::None, elevation: 0 };
            }
        }
    }

    fn generate_content(&self, world: &mut Vec<Vec<Tile>>, biomes_map: &mut HashMap<Biomes, Vec<(usize, usize)>>, coordinates: &Vec<[f64; 2]>, allowed_biomes: &Vec<Biomes>, tile_types_to_avoid: Option<&Vec<TileType>>, content: &Content) {
        for coordinate in coordinates.iter() {
            let (x, y) = (coordinate[0] as usize, coordinate[1] as usize);
            for biome in allowed_biomes {
                if biomes_map.get(biome).is_some() {
                    let coords = biomes_map.get(biome).unwrap();
                    if coords.contains(&(x, y)) {
                        match tile_types_to_avoid {
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

        let tile_types_to_avoid = vec![TileType::ShallowWater, TileType::DeepWater, TileType::Lava];
        let tile_types_to_avoid_for_trees = vec![TileType::ShallowWater, TileType::DeepWater, TileType::Lava, TileType::Street];
        let allowed_biomes = vec![Biomes::Beach, Biomes::Desert, Biomes::Plain, Biomes::Forest, Biomes::Hill];
        let configurations = vec![
            (3, vec![Biomes::Forest], Some(&tile_types_to_avoid_for_trees), Content::Tree(1)),
            (4, vec![Biomes::Hill], Some(&tile_types_to_avoid_for_trees), Content::Tree(1)),
            (5, vec![Biomes::Mountain], Some(&tile_types_to_avoid_for_trees), Content::Tree(1)),
            (5, vec![Biomes::Plain], Some(&tile_types_to_avoid), Content::Rock(1)),
            (4, vec![Biomes::Hill], Some(&tile_types_to_avoid), Content::Rock(1)),
            (3, vec![Biomes::Mountain], Some(&tile_types_to_avoid), Content::Rock(1)),
            (3, vec![Biomes::SnowyMountain], Some(&tile_types_to_avoid), Content::Rock(1)),
            (5, vec![Biomes::ShallowWater], None, Content::Fish(1)),
            (4, vec![Biomes::Deepwater], None, Content::Fish(1)),
            (10, allowed_biomes.clone(), Some(&tile_types_to_avoid), Content::Garbage(1)),
            (20, allowed_biomes.clone(), Some(&tile_types_to_avoid), Content::Coin(1)),
            (20, allowed_biomes.clone(), Some(&tile_types_to_avoid), Content::Bin(0..5)),
            (40, allowed_biomes.clone(), Some(&tile_types_to_avoid), Content::Crate(0..5)),
            (50, allowed_biomes.clone(), Some(&tile_types_to_avoid), Content::Market(1)),
        ];

        let mut coords: HashMap<u64, Vec<[f64; 2]>> = HashMap::new();
        for (radius, test_allowed_biomes, tile_types_to_avoid, content) in configurations {
            if coords.get(&radius).is_none() {
                poisson.set_dimensions([self.world_size as f64, self.world_size as f64], radius as f64);
                coords.insert(radius, poisson.generate());
            }

            self.generate_content(world, biomes_map, coords.get(&radius).unwrap(), &test_allowed_biomes, tile_types_to_avoid, &content);
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

        let elevation_map = self.generate_elevation(6);
        profiler.print_elapsed_time_in_ms("elevation generation time");

        let mut biomes_map = self.generate_biomes(&mut world, &elevation_map);
        profiler.print_elapsed_time_in_ms("biomes generation time");

        self.generate_rivers(&mut world, &elevation_map, 80.0,  0.05, 0.9, 4.0);
        profiler.print_elapsed_time_in_ms("rivers generation time");

        self.generate_streets(&mut world, &elevation_map, 0.5, 2.0, 0.8, 2.5);
        profiler.print_elapsed_time_in_ms("streets generation time");

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
