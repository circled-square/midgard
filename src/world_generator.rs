mod isize_index_matrix;
mod multi_octave_noise;
pub mod params;
mod performance_profiler;
mod vector_math;

use fast_poisson::Poisson2D;
use isize_index_matrix::*;
use multi_octave_noise::Multi;
use noise::*;
use num_traits::pow::Pow;
use performance_profiler::PerformanceProfiler;
use rand::rngs::StdRng;
use rand::{distributions, Rng, SeedableRng};
use robotics_lib::world::{environmental_conditions::{EnvironmentalConditions, WeatherType}, tile::Tile, tile::{Content, TileType}, world_generator::Generator};
use std::collections::{HashMap, HashSet};
use std::hash::{Hash, Hasher};
use std::collections::hash_map::DefaultHasher;
use std::time::SystemTime;
use vector_math::*;

use params::*;
use rayon::prelude::*;

macro_rules! call_with_seed {
    ($this:ident . $fn_name:ident ($($arg:expr),*) ) => {{
        let mut hasher = DefaultHasher::new();
        let s = stringify! { $fn_name };
        s.hash(&mut hasher);
        $this.$fn_name($this.params.seed.wrapping_add(hasher.finish()), $($arg),*)
    }}
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
    SnowyMountain,
}

const WORLD_SCALE_MULTIPLIER : f64 = 180.0;

pub struct WorldGenerator {
    params: WorldGeneratorParameters,
    poisson: Poisson2D
}
impl WorldGenerator {
    pub fn new(params: WorldGeneratorParameters) -> Self {
        let mut poisson = Poisson2D::new();
        poisson.set_samples(15);
        Self { params, poisson }
    }

    fn generate_environmental_conditions(&self, seed: u64) -> EnvironmentalConditions {
        let weather_forecast =
            if self.params.always_sunny {
                vec![WeatherType::Sunny]
            } else {
                let weather_types = vec![
                    WeatherType::Sunny,
                    WeatherType::Foggy,
                    WeatherType::Rainy,
                    WeatherType::TrentinoSnow,
                    WeatherType::TropicalMonsoon,
                ];
                let mut weather_forecast = vec![];

                let mut rng = StdRng::seed_from_u64(seed);

                assert_ne!(self.params.weather_forecast_length, 0, "weather_forecast_length must be > 0");

                for _ in 0..self.params.weather_forecast_length {
                    weather_forecast.push(weather_types[rng.gen_range(0..weather_types.len())]);
                }
                weather_forecast
            };

        return EnvironmentalConditions::new(&weather_forecast, self.params.time_progression_minutes, self.params.starting_hour).unwrap();
    }

    fn generate_elevation(&self, seed: u64) -> Vec<Vec<f64>> {
        let octaves = 6;
        let noise_function = Multi::new(Perlin::new(seed as u32), octaves, 1.0 / (self.params.world_scale * WORLD_SCALE_MULTIPLIER));
        let world_size = self.params.world_size;

        let noise_function = Multiply::new(Constant::new(1.5), noise_function);

        let mut elevation_map = vec![vec![0.0; world_size]; world_size];

        elevation_map
            .par_iter_mut()
            .enumerate()
            .for_each(|(x, col)|{
                for y in 0..world_size {
                    col[y] = noise_function.get([x as f64, y as f64]);
                }
            });

        return elevation_map;
    }

    fn set_elevation_on_tiles(&self, world: &mut Vec<Vec<Tile>>, elevation_map: &Vec<Vec<f64>>) {
        world
            .par_iter_mut()
            .enumerate()
            .for_each(|(x, col)| {
                for y in 0..self.params.world_size {
                    let elevation_normalized = (elevation_map[x][y].clamp(-1.0, 1.0) + 1.0) / 2.0;
                    col[y].elevation = (elevation_normalized * self.params.elevation_multiplier.unwrap()) as usize;
                }
            });
    }

    fn generate_temperature_map(&self, seed: u64) -> Vec<Vec<f64>> {
        let world_size = self.params.world_size;
        let mut temperature_map = vec![vec![0.0; world_size]; world_size];

        let noise_function = Multiply::new(
            Constant::new(1.5),
            Multi::new(Perlin::new(seed as u32), 7, 1.0 / (self.params.world_scale * 0.56 * WORLD_SCALE_MULTIPLIER)),
        );

        temperature_map
            .par_iter_mut()
            .enumerate()
            .for_each(|(x, col)| {
                for y in 0..world_size {
                    col[y] = noise_function.get([x as f64, y as f64]);
                }
            });

        return temperature_map;
    }

    fn generate_biomes(&self, seed: u64, world: &mut Vec<Vec<Tile>>, elevation_map: &Vec<Vec<f64>>) -> HashMap<Biomes, HashSet<(usize, usize)>> {
        let mut biomes_map: HashMap<Biomes, HashSet<(usize, usize)>> = HashMap::new();
        let temperature_map = self.generate_temperature_map(seed);

        for x in 0..self.params.world_size {
            for y in 0..self.params.world_size {
                let biome = match elevation_map[x][y]{
                    h if h < -0.65 => Biomes::Deepwater,
                    h if h < -0.40 => Biomes::ShallowWater,
                    h if h < -0.30 => Biomes::Beach,
                    h if h < 0.35 => match temperature_map[x][y] {
                        t if t < -0.3 => Biomes::Forest,
                        t if t > 0.2 => Biomes::Desert,
                        _ => Biomes::Plain,
                    },
                    h if h < 0.55 => Biomes::Hill,
                    h if h < 0.85 => Biomes::Mountain,
                    _ => Biomes::SnowyMountain,
                };
                biomes_map
                    .entry(biome)
                    .or_insert(HashSet::new())
                    .insert((x, y));
            }
        }

        let mut fill_tiles_with = |coords: &HashSet<(usize, usize)>, tile_type: TileType| {
            for (x, y) in coords {
                world[*x][*y] = Tile{ tile_type, content: Content::None, elevation: 0 };
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

    fn generate_hellfire(&self, seed: u64, world: &mut Vec<Vec<Tile>>, biomes_map: &mut HashMap<Biomes, HashSet<(usize, usize)>>) {
        // Lava
        if biomes_map.get(&Biomes::Desert).is_none() {
            return;
        }

        let lava_noise_function = Multiply::new(
            Constant::new(1.5),
            Multi::new(Perlin::new(seed as u32), 7, 1.0 / (self.params.world_scale * 0.17 * WORLD_SCALE_MULTIPLIER)),
        );
        let lava_noise = |x : usize, y : usize| lava_noise_function.get([x as f64, y as f64]);

        for (x, y) in biomes_map.get(&Biomes::Desert).unwrap() {
            if world[*x][*y].tile_type != TileType::ShallowWater && world[*x][*y].content != Content::Fire && lava_noise(*x,*y) < -0.6 {
                world[*x][*y].tile_type = TileType::Lava;
            }
        }

        // Fire
        if biomes_map.get(&Biomes::Plain).is_none() {
            return;
        }

        let fire_noise_function = Multiply::new(
            Constant::new(1.5),
            Multi::new(Perlin::new(seed as u32 + 1), 7, 1.0 / (self.params.world_scale * 0.17 * WORLD_SCALE_MULTIPLIER)),
        );
        let fire_noise = |x : usize, y : usize| fire_noise_function.get([x as f64, y as f64]);

        for (x, y) in biomes_map.get(&Biomes::Plain).unwrap() {
            if world[*x][*y].tile_type.properties().can_hold(&Content::Fire) && fire_noise(*x,*y) < -0.5 {
                world[*x][*y].content = Content::Fire;
            }
        }
    }

    fn generate_rivers(&mut self, seed: u64, world: &mut Vec<Vec<Tile>>, elevation: &Vec<Vec<f64>>) {
        if self.params.amount_of_rivers.is_none() {
            return;
        }

        let world_size = self.params.world_size;
        let amount = self.params.amount_of_rivers.unwrap() * 9.0;
        let number_of_rivers = (world.len() as f64 * world.len() as f64 * amount * amount / 1000000.0) as usize;
        let random_jitter = 0.12;
        let inertia_factor = 0.05;
        let inertia_decay = 0.9;
        let max_inertia = 4.0;

        let mut rng = StdRng::seed_from_u64(seed);

        let world_pos_distribution = distributions::Uniform::new(1, world.len() as isize - 1);
        let float_distribution = distributions::Uniform::new(0.0, 1.0);
        let directions = [(-1, 0), (1, 0), (0, -1), (0, 1)];

        self.poisson.set_seed(seed);
        self.poisson.set_dimensions([1.0, 1.0], 1.0 / (number_of_rivers as f64).sqrt());

        let mut poisson_coords_iterator = self.poisson.iter();

        for _ in 0..number_of_rivers {
            //get a vector containing a river's tiles
            let river_tiles = loop {
                let mut number_of_loops_looking_for_start_coord = 0; //safeguard against infinite looping while looking for a valid random start coordinates in worlds that lack it

                //sample the poisson distribution for a start coordinate, and discard it if it is not
                //Hill, Mountain or Snow or if it is close (within 3 tiles) to a ShallowWater tile
                let start_coords = loop {
                    let coords = poisson_coords_iterator
                        .next()
                        .map(|p| ((p[0] * world_size as f64) as isize, (p[1] * world_size as f64) as isize))
                        .unwrap_or_else(|| { (rng.sample(world_pos_distribution), rng.sample(world_pos_distribution)) });

                    if [TileType::Hill, TileType::Mountain, TileType::Snow].contains(&world.at(coords).tile_type) {
                        let mut should_discard = false;
                        'outer: for x in -3..3 {
                            for y in -3..3 {
                                if let Some(Tile { tile_type: TileType::ShallowWater,..}) = world.at_checked((coords.0 + x, coords.1 + y)) {
                                    should_discard = true;
                                    break 'outer;
                                }
                            }
                        }

                        if !should_discard {
                            break coords;
                        }
                    }

                    if number_of_loops_looking_for_start_coord > world_size * world_size * 10 {
                        return;
                    }
                    number_of_loops_looking_for_start_coord += 1;
                };

                let mut river_tiles_stack = vec![start_coords];
                let mut avoid_tiles = HashSet::from([start_coords]);
                let mut river_tiles_set = HashSet::from([start_coords]); // always contains the same as river_tiles_stack
                let mut inertia = (0.0, 0.0);
                let mut movement_debt = (0.0, 0.0);

                // build river_tiles_stack with tiles based on the inertia of the river, backtracking
                // when there is no valid direction that hasn't been visited or when
                loop {
                    if river_tiles_stack.is_empty() {
                        break;
                    }
                    let coords = *river_tiles_stack.last().unwrap();

                    // backtrack to a prev tile if the current is close to lava or close to a tile in river_tiles_set (but has not been added recently)
                    {
                        let mut should_backtrack = false;
                        let recent_tiles_size = river_tiles_stack.len().min(7);
                        let recent_tiles = &river_tiles_stack[river_tiles_stack.len() - recent_tiles_size..];
                        'outer: for x in -3..=3 {
                            for y in -3..=3 {
                                let c = (coords.0 + x, coords.1 + y);

                                if river_tiles_set.contains(&c) && !recent_tiles.contains(&c) {
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
                            let tile = river_tiles_stack.pop().unwrap(); //backtrack to the previous tile
                            river_tiles_set.remove(&tile);
                            continue;
                        }
                    }

                    // when the river is close to the edge of the world randomly decide if it should
                    // flow off of it or if it should backtrack
                    if [0, world.len() as isize - 1].contains(&coords.0) || [0, world.len() as isize - 1].contains(&coords.1) {
                        match rng.sample(float_distribution) {
                            n if n < 0.5 => {
                                break; //flow off the edge of the world
                            }
                            _ => {
                                let tile = river_tiles_stack.pop().unwrap(); //backtrack to the previous tile
                                river_tiles_set.remove(&tile);
                                continue;
                            }
                        }
                    }

                    // stop if river flowed into another river or into a lake
                    if [TileType::ShallowWater, TileType::DeepWater].contains(&world.at(coords).tile_type) {
                        break;
                    }

                    // candidates contains the possible future tiles for the river to flow in
                    let candidates: Vec<_> = directions
                        .iter()
                        .filter_map(|(x, y)| {
                            let c = (coords.0 + x, coords.1 + y);

                            if avoid_tiles.contains(&c) {
                                return None;
                            }
                            if !(0..world.len() as isize).contains(&c.0) || !(0..world.len() as isize).contains(&c.1) {
                                return None;
                            }
                            return Some(c);
                        })
                        .collect();

                    // backtrack if candidates is empty (nowhere to flow to)
                    if candidates.len() == 0 {
                        let tile = river_tiles_stack.pop().unwrap();
                        river_tiles_set.remove(&tile);
                        continue;
                    }


                    //calculate a theoretical direction of movement based on a weighed sum of:
                    // - gradient of elevation,
                    // - inertia
                    // - random jitter
                    let gradient_of_elevation = get_gradient(elevation, coords).unwrap();
                    let direction = vec_normalize(vec_sum(
                        vec_sum(
                            vec_mul_by_scalar(gradient_of_elevation, -1.0),
                            vec_mul_by_scalar(inertia, inertia_factor),
                        ),
                        vec_mul_by_scalar(
                            (
                                (rng.sample(float_distribution) - 0.5) * 2.0,
                                (rng.sample(float_distribution) - 0.5) * 2.0,
                            ),
                            random_jitter,
                        ),
                    ));
                    //update inertia based on movement
                    inertia = vec_clamp(
                        vec_sum(vec_mul_by_scalar(inertia, inertia_decay), direction),
                        max_inertia,
                    );

                    //update movement_debt
                    movement_debt = vec_clamp(movement_debt, 1.0);
                    movement_debt = vec_sum(movement_debt, direction);

                    // actually compute movement
                    let move_along_x_axis = movement_debt.0.abs() > movement_debt.1.abs();

                    let x_candidate_dir = if direction.0 > 0.0 { 1 } else { -1 };
                    let x_candidate = (coords.0 + x_candidate_dir, coords.1);
                    let minus_x_candidate = (coords.0 - x_candidate_dir, coords.1);
                    let y_candidate_dir = if direction.1 > 0.0 { 1 } else { -1 };
                    let y_candidate = (coords.0, coords.1 + y_candidate_dir);
                    let minus_y_candidate = (coords.0, coords.1 - y_candidate_dir);

                    let target_dir = if move_along_x_axis {
                        if candidates.contains(&x_candidate) {
                            (x_candidate_dir, 0)
                        } else if candidates.contains(&y_candidate) {
                            (0, y_candidate_dir)
                        } else if candidates.contains(&minus_y_candidate) {
                            (0, -y_candidate_dir)
                        } else {
                            (-x_candidate_dir, 0)
                        }
                    }
                    else {
                        if candidates.contains(&y_candidate) {
                            (0, y_candidate_dir)
                        } else if candidates.contains(&x_candidate) {
                            (x_candidate_dir, 0)
                        } else if candidates.contains(&minus_x_candidate) {
                            (-x_candidate_dir, 0)
                        } else {
                            (0, -y_candidate_dir)
                        }
                    };

                    // update movement_debt removing from the debt the movement that actually happened
                    movement_debt = vec_subtract(movement_debt, (target_dir.0 as f64, target_dir.1 as f64));

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

    fn generate_street(&self, world: &mut Vec<Vec<Tile>>, elevation: &Vec<Vec<f64>>, poi1: (isize, isize), poi2: (isize, isize), inertia: &mut (f64, f64)) {
        let inertia_factor = 2.0;
        let inertia_decay = 0.8;
        let max_inertia = 2.5;

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

            // use dot product to find which perpendicular to the gradient is more concordant with the desired direction
            let perpendicular_to_gradient_of_elevation_concordant_to_desired_direction = {
                let desired_direction_dot_perpendiculars_to_gradient_of_elevation =
                    perpendiculars_to_gradient_of_elevation.map(|p| vec_dot(p, desired_direction));
                let mut max = 0;
                for i in 1..desired_direction_dot_perpendiculars_to_gradient_of_elevation.len() {
                    if desired_direction_dot_perpendiculars_to_gradient_of_elevation[i]
                        > desired_direction_dot_perpendiculars_to_gradient_of_elevation[max]
                    {
                        max = i;
                    }
                }
                perpendiculars_to_gradient_of_elevation[max]
            };

            // when the dist to poi2 is smaller than dist(poi1, poi2)/2
            // we will start giving more weight to desired_direction and less to inertia.
            let proximity_threshold =
                (((poi1.0 - poi2.0).pow(2) + (poi1.1 - poi2.1).pow(2)) as f64).sqrt() / 2.0;
            let proximity_multiplier =
                f64::max(1.0, proximity_threshold / vec_module(dist_to_poi2));

            // movement direction will be a weighed sum of desired_direction, inertia and
            // perpendicular_to_gradient_of_elevation_concordant_to_desired_direction
            let direction = vec_sum(
                vec_sum(
                    vec_mul_by_scalar(desired_direction, proximity_multiplier),
                    vec_mul_by_scalar(*inertia, inertia_factor),
                ),
                vec_mul_by_scalar(
                    perpendicular_to_gradient_of_elevation_concordant_to_desired_direction,
                    0.5,
                ),
            );
            let direction = vec_normalize(direction);

            //update movement_debt
            movement_debt = vec_clamp(movement_debt, 5.0);
            movement_debt = vec_sum(movement_debt, direction);

            //actually compute movement
            let move_along_x_axis = movement_debt.0.abs() > movement_debt.1.abs();

            let x_candidate_dir = if direction.0 > 0.0 { 1 } else { -1 };
            let y_candidate_dir = if direction.1 > 0.0 { 1 } else { -1 };

            let target_dir = if move_along_x_axis {
                (x_candidate_dir, 0)
            } else {
                (0, y_candidate_dir)
            };

            //update movement_debt subtracting the movement
            movement_debt = vec_subtract(movement_debt, (target_dir.0 as f64, target_dir.1 as f64));

            *inertia = vec_sum(*inertia, direction);
            *inertia = vec_mul_by_scalar(*inertia, inertia_decay);
            *inertia = vec_clamp(*inertia, max_inertia);

            let target_tile = (coords.0 + target_dir.0, coords.1 + target_dir.1);

            avoid_tiles.insert(target_tile);
            street_tiles_stack.push(target_tile);
            street_tiles_set.insert(target_tile);
        }

        //fill tiles in street_tiles_stack with Street (if they aren't water or lava)
        for pos in street_tiles_stack {
            if let Some(tile) = world.at_mut_checked(pos) {
                if ![TileType::Lava, TileType::DeepWater, TileType::ShallowWater]
                    .contains(&tile.tile_type)
                {
                    tile.tile_type = TileType::Street;
                }
            }
        }
    }
    fn generate_streets(&mut self, seed: u64, world: &mut Vec<Vec<Tile>>, elevation: &Vec<Vec<f64>>) {
        if self.params.amount_of_streets.is_none() {
            return;
        }

        let amount = self.params.amount_of_streets.unwrap() * 0.3;
        //this function generates points of interest in the map and tries to create roads connecting them
        let poi_distance = 50.0 / amount;
        let mut rng = StdRng::seed_from_u64(seed);

        //generate random Points Of Interest on the map to connect with streets
        self.poisson.set_seed(seed);
        self.poisson.set_dimensions([1.0, 1.0],poi_distance / (1.5 * world.len() as f64));
        let poi_vec: Vec<_> = self.poisson
            .iter()
            .map(|p| { p.map(|c| ((c * 1.5 - 0.25) * world.len() as f64) as isize) })
            .map(|[x, y]| (x, y))
            .collect();
        if poi_vec.is_empty() {
            return;
        }
        let random_poi_distribution = distributions::Uniform::new(0, poi_vec.len());

        // for any 2 POIs A, B that are connected by a street street_map will contain (A,B) and (B,A)
        let mut street_map = HashSet::new();
        let number_of_streets = (world.len() as f64 / poi_distance).pow(2.0) as usize * 2;
        let mut number_of_empty_streets_produced = 0; // stop after number_of_streets empty streets have been produced (empty street = street generation failure)

        while street_map.len() < number_of_streets * 2 {
            // generate a street (a vec of POIs)
            let start_poi = rng.sample(random_poi_distribution);
            let number_of_pois_in_street = world.len() / poi_distance as usize;
            let mut street = vec![start_poi];
            let mut avoid_set = HashSet::from([start_poi]);

            for _ in 0..number_of_pois_in_street {
                let last_poi_index = *street.last().unwrap();
                let last_poi = poi_vec[last_poi_index];

                let mut loop_iterations = 0;
                let next_poi_index = loop {
                    loop_iterations += 1;
                    if loop_iterations as f64 > amount * amount * 100.0 {
                        break None;//no valid next_poi could be found
                    }
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

                    if dist < poi_distance * 2.0 {
                        break Some(next_index);
                    }
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

            //actually build the street and populate the tiles
            let mut inertia = (0.0, 0.0);
            for i in 0..street.len() - 1 {
                let poi1 = poi_vec[street[i]];
                let poi2 = poi_vec[street[i + 1]];
                self.generate_street(world, elevation, poi1, poi2, &mut inertia);
                inertia = vec_mul_by_scalar(inertia, 5.0); // exaggerate inertia when passing through poi
            }

            if street.len() == 1 {
                number_of_empty_streets_produced += 1;
            }
            if number_of_empty_streets_produced > number_of_streets {
                break;
            }
        }
    }

    fn generate_teleports(&mut self, seed: u64, world: &mut Vec<Vec<Tile>>) {
        if self.params.amount_of_teleports.is_none() {
            return;
        }

        let amount = self.params.amount_of_teleports.unwrap();
        let world_size = self.params.world_size;

        self.poisson.set_seed(seed);
        self.poisson.set_dimensions([world_size as f64, world_size as f64], 100.0 / amount);
        let tile_types_to_avoid = vec![TileType::ShallowWater, TileType::DeepWater, TileType::Lava];

        for coordinate in self.poisson.iter() {
            let (x, y) = (coordinate[0] as usize, coordinate[1] as usize);
            if world[x][y].content != Content::Fire
                && tile_types_to_avoid.iter().all(|t| *t != world[x][y].tile_type)
            {
                world[x][y] = Tile {
                    tile_type: TileType::Teleport(false),
                    content: Content::None,
                    elevation: 0,
                };
            }
        }
    }

    fn generate_content(&self, world: &mut Vec<Vec<Tile>>, biomes_map: &mut HashMap<Biomes, HashSet<(usize, usize)>>, coords: &Vec<[f64; 2]>, allowed_biomes: &Vec<Biomes>, content: &Content) {
        for coord in coords.iter() {
            let (x, y) = (coord[0] as usize, coord[1] as usize);

            for biome in allowed_biomes {
                let Some(biome_coords) = biomes_map.get(biome) else {
                    continue;
                };

                if biome_coords.contains(&(x, y)) {
                    if world[x][y].tile_type.properties().can_hold(&content) {
                        world[x][y].content = content.clone();
                    }
                }
            }
        }
    }

    fn generate_contents(&mut self, seed: u64, world: &mut Vec<Vec<Tile>>, biomes_map: &mut HashMap<Biomes, HashSet<(usize, usize)>>) {
        let world_size = self.params.world_size;
        //Water
        if biomes_map.get(&Biomes::Deepwater).is_some() {
            for (x, y) in biomes_map.get(&Biomes::Deepwater).unwrap() {
                world[*x][*y].content = Content::Water(2);
            }
        }
        if biomes_map.get(&Biomes::ShallowWater).is_some() {
            for (x, y) in biomes_map.get(&Biomes::ShallowWater).unwrap() {
                world[*x][*y].content = Content::Water(1);
            }
        }

        self.poisson.set_seed(seed);

        let allowed_biomes = vec![Biomes::Beach, Biomes::Desert, Biomes::Plain, Biomes::Forest, Biomes::Hill];
        let radii = &self.params.contents_radii;
        let configurations = vec![
            (radii.trees_in_forest, vec![Biomes::Forest], Content::Tree(1)),
            (radii.trees_in_hill, vec![Biomes::Hill], Content::Tree(1)),
            (radii.trees_in_mountain, vec![Biomes::Mountain], Content::Tree(1)),

            (radii.rocks_in_plains, vec![Biomes::Plain], Content::Rock(1)),
            (radii.rocks_in_hill, vec![Biomes::Hill], Content::Rock(1)),
            (radii.rocks_in_mountain, vec![Biomes::Mountain], Content::Rock(1)),
            (radii.rocks_in_mountain, vec![Biomes::SnowyMountain], Content::Rock(1)),

            (radii.bushes_in_plains, vec![Biomes::Plain], Content::Bush(1)),

            (radii.fish_in_shallow_water, vec![Biomes::ShallowWater], Content::Fish(1)),
            (radii.fish_in_deep_water, vec![Biomes::Deepwater], Content::Fish(1)),

            (radii.garbage, allowed_biomes.clone(), Content::Garbage(1)),
            (radii.coins, allowed_biomes.clone(), Content::Coin(1)),
            (radii.garbage_bins, allowed_biomes.clone(), Content::Bin(0..5)),
            (radii.crates, allowed_biomes.clone(), Content::Crate(0..5)),
            (radii.markets, allowed_biomes.clone(), Content::Market(1)),
            (radii.banks, allowed_biomes.clone(), Content::Bank(0..5)),
            (radii.buildings, allowed_biomes.clone(), Content::Building),
            (radii.scarecrows, allowed_biomes.clone(), Content::Scarecrow),
            (radii.jolly_blocks, allowed_biomes.clone(), Content::JollyBlock(1)),
        ];

        let mut coords: HashMap<u64, Vec<[f64; 2]>> = HashMap::new();
        for (radius, allowed_biomes, content) in configurations {
            if coords.get(&radius).is_none() {
                self.poisson.set_dimensions([world_size as f64, world_size as f64], radius as f64);
                coords.insert(radius, self.poisson.generate());
            }

            self.generate_content(world, biomes_map, coords.get(&radius).unwrap(), &allowed_biomes, &content);
        }
    }

    fn generate_spawnpoint(&mut self, seed: u64, biomes_map: &mut HashMap<Biomes, HashSet<(usize, usize)>>) -> (usize, usize) {
        self.poisson.set_seed(seed);
        self.poisson.set_dimensions([self.params.world_size as f64, self.params.world_size as f64], self.params.world_size as f64 / 10.0);

        let allowed_biomes = vec![Biomes::Plain, Biomes::Beach, Biomes::Forest];
        let mut spawnpoint = (0, 0);

        'outer: for coordinate in self.poisson.iter() {
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

        return spawnpoint;
    }
}


impl Generator for WorldGenerator {
    fn gen(&mut self) -> (Vec<Vec<Tile>>, (usize, usize), EnvironmentalConditions, f32, Option<HashMap<Content, f32>>) {
        assert!(self.params.world_size > 2, "world size must be 3 or more");

        let mut profiler = PerformanceProfiler::new(SystemTime::now());

        println!("World seed: {}, size {}", self.params.seed, self.params.world_size);

        let deep_water_tile = Tile {
            tile_type: TileType::DeepWater,
            content: Content::None,
            elevation: 0,
        };
        let mut world = vec![vec![deep_water_tile; self.params.world_size]; self.params.world_size];
        profiler.print_elapsed_time_in_ms("water world generation time");

        let elevation_map = call_with_seed!(self.generate_elevation());
        profiler.print_elapsed_time_in_ms("elevation generation time");

        let mut biomes_map = call_with_seed!(self.generate_biomes(&mut world, &elevation_map));
        profiler.print_elapsed_time_in_ms("biomes generation time");

        call_with_seed!(self.generate_rivers(&mut world, &elevation_map));
        profiler.print_elapsed_time_in_ms("rivers generation time");

        call_with_seed!(self.generate_streets(&mut world, &elevation_map));
        profiler.print_elapsed_time_in_ms("streets generation time");

        call_with_seed!(self.generate_hellfire(&mut world, &mut biomes_map));
        profiler.print_elapsed_time_in_ms("hellfire generation time");

        call_with_seed!(self.generate_teleports(&mut world));
        profiler.print_elapsed_time_in_ms("teleports generation time");

        call_with_seed!(self.generate_contents(&mut world, &mut biomes_map));
        profiler.print_elapsed_time_in_ms("Total contents generation time");

        let environmental_conditions = call_with_seed!(self.generate_environmental_conditions());
        profiler.print_elapsed_time_in_ms("weather generation time");

        let spawn_point = call_with_seed!(self.generate_spawnpoint(&mut biomes_map));
        println!("Spawn point {spawn_point:?}");
        profiler.print_elapsed_time_in_ms("spawn point generation time");

        if self.params.elevation_multiplier.is_some() {
            self.set_elevation_on_tiles(&mut world, &elevation_map);
            profiler.print_elapsed_time_in_ms("elevation setting time");
        }

        profiler.print_total_elapsed_time_in_ms("Total generation time");

        (world, spawn_point, environmental_conditions, self.params.max_score, self.params.score_table.clone())
    }
}